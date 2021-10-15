import pysam
from sys import path as sys_path
if 'source' in sys_path[0]:
    from re import search as re_search
    p_to_add=sys_path[0][:re_search('Binning_refiner',sys_path[0]).span()[1]]
    if not sys_path[0].endswith('Binning_refiner') :  sys_path[0]=p_to_add

from math import log10,log2,log,ceil


#for pre_processor
from itertools import product
import logging
from time import time,sleep
import numpy as np
import pandas as pd
from scipy.spatial import distance
from psutil import virtual_memory
#for contig to candidate contig distance
from tqdm import tqdm as progress_bar
from sklearn.metrics.pairwise import cosine_similarity
from multiprocessing import Manager

class Pre_Processor():
    def __init__(self,assembly_fasta, gff_annotation,bam_file, depth_tsv,binned_contigs,annotations_tsv,taxonomy_tsv, min_cont_size,min_annot_cont_size,kmer_sizes):
        #input files
        self.assembly_fasta=assembly_fasta
        self.gff_annotation=gff_annotation
        self.bam_file=bam_file
        self.depth_tsv=depth_tsv
        self.binned_contigs=binned_contigs
        self.annotations_tsv=annotations_tsv
        self.taxonomy_tsv=taxonomy_tsv

        #parameters
        self.min_cont_size = min_cont_size
        self.min_annot_cont_size = min_annot_cont_size
        self.kmer_sizes = kmer_sizes
        # cityblock for manhattan, cosine for cosine. we use cosine so it is scaled from 0 to 1
        self.distance_function=distance.cosine
        self.invalid_score=-1
        self.rounding_digits=4
        #to convert into sqlite db:
        self.contigs_to_cds={} #contigs and respective cds
        #cursor for sqlite connection, if ram is enough cursor remains none
        self.cursor=None
        self.bins={} #bins and respective contigs; probably fine in RAM
        self.unbinned=set() # just all unbinned contigs; probably fine in RAM
        self.contigs_taxonomy={}
        self.avg_vector_bin={} # bin and respective max distance; probably fine in RAM
        self.idf={} # probably fine in ram
        self.bin_details={} # bins with two sets, one for metabolites and other for annotations

        #saving instances to memory so we dont load them twice
        self.annotations_reactions={}
        self.annotations={}


    def clear_annotations_dicts(self):
        del self.annotations
        del self.annotations_reactions


    ####Loading assembly information

    def extract_features_gff(self):
        '''
        returns low complexity features e.g., rRNA, crispr
        '''
        low_comp_feature_dict = {}
        contigs_to_cds={}
        # Get all info for each line
        with open(self.gff_annotation, 'r') as file:
            for line in file:
                line = line.strip(' \t\n').split('\t')
                contig, feature_type, feature_start, feature_end, attributes = line[0], line[2], line[3], line[4], line[8].split(';')
                if contig not in contigs_to_cds:contigs_to_cds[contig]=set()
                # Check if feature has no start and/or stop position or length = 0
                if not feature_start or not feature_end or int(feature_start)-int(feature_end) == 0:
                    continue
                # Only add line if feature was annotated as essential gene
                for attribute in attributes:
                    attribute_name, attribute_id = attribute.split('=')[0], attribute.split('=')[1]
                    if feature_type=='CDS' and attribute_name=='ID':
                        contigs_to_cds[contig].add(attribute_id)
                    elif attribute_name == 'product' and 'ribosomal RNA' in attribute_id:
                        if contig not in low_comp_feature_dict: low_comp_feature_dict[contig]=[]
                        low_comp_feature_dict[contig].append((int(feature_start), int(feature_end)))
                    elif attribute_name == 'rpt_family':
                        if contig not in low_comp_feature_dict: low_comp_feature_dict[contig] = []
                        low_comp_feature_dict[contig].append((int(feature_start), int(feature_end)))
        self.contigs_to_cds=contigs_to_cds
        return low_comp_feature_dict

    def read_fasta_generator(self,fasta_path):
        query=None
        seq=[]
        with open(fasta_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    if query:
                        yield query,''.join(seq).upper()
                        seq=[]
                    query=line.replace('>','').strip()
                else:
                    seq.append(line.strip())
            if query:
                yield query, ''.join(seq).upper()

    def is_valid_annot_contig(self,contig_id,seq_str):
        contig_seqs = self.contigs_to_cds.get(contig_id,[])
        for seq_id in contig_seqs:
            if seq_id in self.annotations and len(seq_str)>= self.min_annot_cont_size:
                return True
        return False

    def clean_annotations(self,contig_seqs):
        '''
        will remove annotations for sequences that are not valid (below min_annot_contig_size)
        '''
        valid_seqs=set()
        for contig_id in contig_seqs:
            for seq_id in self.contigs_to_cds.get(contig_id,[]):
                valid_seqs.add(seq_id)
        to_remove=set()
        for seq_id in self.annotations:
            if seq_id not in valid_seqs:
                to_remove.add(seq_id)
        for seq_id in to_remove:
            self.annotations.pop(seq_id)


    def load_assembly_len(self):
        fasta_reader = self.read_fasta_generator(self.assembly_fasta)
        self.contigs_len = {}
        for contig_id, seq_str in fasta_reader:
            self.contigs_len[contig_id]=len(seq_str)

    def load_assembly(self,paired_contigs):
        fasta_reader = self.read_fasta_generator(self.assembly_fasta)
        contig_seqs = {}
        for contig_id, seq_str in fasta_reader:
            if contig_id in paired_contigs:
                contig_seqs[contig_id] = seq_str
            elif len(seq_str) >= self.min_cont_size:
                contig_seqs[contig_id] = seq_str
            elif self.is_valid_annot_contig(contig_id,seq_str):
                contig_seqs[contig_id] = seq_str
        self.clean_annotations(contig_seqs)
        return contig_seqs

    def mask_rep_features(self, contig_rep_feature_dict, contig_seqs):
        for contig_id in contig_seqs:
            if contig_id in contig_rep_feature_dict:
                regions = contig_rep_feature_dict[contig_id]
                for region in regions:
                    end_region=region[1]
                    if region[0] - 1 <= 0:
                        start_region = 0
                        start2feature = 0
                    else:
                        start_region = region[0] - 1
                        start2feature = region[0] - 2
                    contig_seq = contig_seqs[contig_id]
                    masked_contig_seq = contig_seq[:start2feature] \
                                 + contig_seq[start_region:end_region - 1].lower() \
                                 + contig_seq[end_region:]
                    contig_seqs[contig_id]=masked_contig_seq

    ####Kmer frequency counting

    def rec_comp(self,kmer):
        base_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        rc_kmer = "".join(base_pairs.get(base, base) for base in reversed(kmer))
        return rc_kmer

    def kmer_is_canonical(self,kmer):
        rc_kmer = self.rec_comp(kmer)
        if kmer <= rc_kmer:
            canonical_kmer = kmer
            is_canonical = True
        else:
            canonical_kmer = rc_kmer
            is_canonical = False
        return canonical_kmer, is_canonical

    def frequencyCount(self,string, substr):
       count = 0
       pos = 0
       while True:
           pos = string.find(substr, pos)
           if pos > -1:
               count = count + 1
               pos += 1
           else:
               break
       return count

    def kmer_counter(self,kmer_list_can, sequence):
        kmer_count = [(self.frequencyCount(sequence, kl_km) + self.frequencyCount(sequence, self.rec_comp(kl_km))) / len(sequence)
                      if kl_km != self.rec_comp(kl_km) else self.frequencyCount(sequence, kl_km) / len(sequence) for kl_km in kmer_list_can]
        return kmer_count

    def get_contig_kmer_freq(self, k_mers, assembly_chunk):
        chunk_list=[]
        for contig in assembly_chunk:
            contig_seq=assembly_chunk[contig]
            kmer_count_list = self.kmer_counter(k_mers, contig_seq)
            contig_kfreq = contig,kmer_count_list
            chunk_list.append(contig_kfreq)
        return chunk_list


    def chunk_generator(self,to_chunk, chunk_size):
        for i in range(0, len(to_chunk), chunk_size):
            yield to_chunk[i:i + chunk_size]

    def chunk_generator_load_balanced(self,list_ordered, chunk_size, time_limit):
        n_chunks = ceil(len(list_ordered) / chunk_size)
        res = []
        direction_chunks = {}
        for i in range(n_chunks):
            res.append([])
            direction_chunks[i] = True
        chunk_index = 0
        start = time()
        while list_ordered:
            if time_limit:
                if time() - start > time_limit: raise TimeoutError
            if direction_chunks[chunk_index]:
                chunk_val = list_ordered.pop(0)
                direction_chunks[chunk_index] = False
            else:
                chunk_val = list_ordered.pop(-1)
                direction_chunks[chunk_index] = True
            res[chunk_index].append(chunk_val)
            if chunk_index == n_chunks - 1:
                chunk_index = 0
            else:
                chunk_index += 1
        return res

    def get_contig_kmer_freq_worker_function(self,queue, master_pid):
        res=[]
        while True:
            record = queue.pop(0)
            if record is None: return
            arg1,arg2=record
            self.mp_results.extend(self.get_contig_kmer_freq(arg1,arg2))


    def launch_get_contig_kmer_freq(self,kmer_list_can,contig_chunks):
        for chunk in contig_chunks:
            self.queue.append([kmer_list_can,chunk])
        self.processes_handler(self.get_contig_kmer_freq_worker_function)
        while self.mp_results:
            yield self.mp_results.pop(0)

    def split_contigs_into_chunks(self,contig_seqs,load_balancing_limit=1000000,time_limit=120):
        chunk_size=ceil(len(contig_seqs)/self.worker_count)
        if len(contig_seqs) < load_balancing_limit:
            contig_seqs_keys_len = {i: len(contig_seqs[i]) for i in contig_seqs}
            list_ordered = sorted(contig_seqs_keys_len, key=contig_seqs_keys_len.__getitem__)
            contig_chunks = self.chunk_generator_load_balanced(list_ordered, chunk_size,time_limit=time_limit)
        else:
            contig_seqs_keys = list(contig_seqs.keys())
            contig_chunks = self.chunk_generator(contig_seqs_keys, chunk_size)
        res=[]
        for chunk in contig_chunks:
            yield {contig:contig_seqs[contig] for contig in chunk}


    def get_contig_kmer_matrix(self,contig_chunks):
        all_kmers=[]
        all_contigs=[]
        contig_kmer_freq_list = []

        for ksize in self.kmer_sizes:
            d = {i: ['A', 'T', 'G', 'C'] for i in range(ksize)}
            kmer_list = [''.join(combination) for combination in product(*[d[k] for k in sorted(d.keys())])]
            kmer_list_can = sorted({i for i in kmer_list if self.kmer_is_canonical(i)[1]})
            #order of kmer frequence per contig is the same as the all_kmers
            all_kmers.extend(kmer_list_can)
        contig_kmer_freq=self.launch_get_contig_kmer_freq(all_kmers,contig_chunks)

        for contig,kmer_freq in contig_kmer_freq:
            if sum(kmer_freq)>0:
                all_contigs.append(contig)
                contig_kmer_freq_list.append(kmer_freq)
        return all_contigs,contig_kmer_freq_list


    def load_depth_dict(self):
        depth_dict = {}
        with open(self.depth_tsv) as file:
            for line in file:
                line = line.strip('\n')
                seq_id,depth_value=line.split('\t')
                depth_dict[seq_id]=float(depth_value)
        return depth_dict

    def multiplicative_replacement(self, mat, delta=None):
        # from sciki-bio http://scikit-bio.org/
        mat = self.closure(mat)
        z_mat = (mat == 0)
        num_feats = mat.shape[-1]
        tot = z_mat.sum(axis=-1, keepdims=True)
        if delta is None:
            delta = (1. / num_feats) ** 2
        zcnts = 1 - tot * delta
        mat = np.where(z_mat, delta, zcnts * mat)
        return mat.squeeze()

    def closure(self, mat):
        # from sciki-bio http://scikit-bio.org/
        mat = np.atleast_2d(mat)
        if np.any(mat < 0):
            raise ValueError("Cannot have negative proportions")
        if mat.ndim > 2:
            raise ValueError("Input matrix can only have two dimensions or less")
        mat = mat / mat.sum(axis=1, keepdims=True)
        return mat.squeeze()

    def clr(self, mat):
        # from sciki-bio http://scikit-bio.org/
        mat = self.closure(mat)
        lmat = np.log(mat)
        gm = lmat.mean(axis=-1, keepdims=True)
        return (lmat - gm).squeeze()

    def assembly_matrix(self,paired_contigs):
        #loading assembly
        start=time()
        low_comp_feature_dict = self.extract_features_gff()
        self.parse_mantis()
        contig_seqs = self.load_assembly(paired_contigs)
        self.mask_rep_features(low_comp_feature_dict, contig_seqs)
        print('Assembly information loading took', time()-start,'seconds')
        #chunking
        contig_chunks=self.split_contigs_into_chunks(contig_seqs)
        #kmer frequency
        start=time()
        contigs_headers,kfreq_array = self.get_contig_kmer_matrix(contig_chunks)
        print('Kmer counting took',time()-start,'seconds')
        start=time()
        contig_kfreqs = np.array(kfreq_array)
        contig_kfreqs = self.multiplicative_replacement(contig_kfreqs)
        # Clr transform   -   justified with binny paper
        contig_kfreqs_scaled = self.clr(contig_kfreqs)
        print('Array transformation took',time()-start,'seconds')

        #depth
        start=time()
        depth_dict = self.load_depth_dict()
        #just for testing
        for i in contigs_headers:
            if i not in depth_dict:
                print('MISSING DEPTH',i)
                depth_dict[i]=1.0

        contig_seqs_depth = np.array([depth_dict[contig] for contig in contigs_headers])
        print('Depth loading took',time()-start,'seconds')

        #preparing array for distance analysis
        #maybe we need to properly scale depth and kmer counting?
        start=time()
        contig_depth_kfreqs_scaled = np.array([np.append(depths, kmers) for depths, kmers in zip(contig_seqs_depth, contig_kfreqs_scaled)])
        contig_depth_kfreq_dict = {contig_header: d_kfreq_array for contig_header, d_kfreq_array in zip(contigs_headers, contig_depth_kfreqs_scaled)}
        self.store_kmer_depth_vectors(contig_depth_kfreq_dict)
        print('Preparing and storing contig vectors took',time()-start,'seconds')

        return set(contigs_headers)

    def load_contigs_taxonomy(self):
        self.contigs_taxonomy={}
        with open(self.taxonomy_tsv) as file:
            line=file.readline()
            #contig	superkingdom	phylum	class	order	family	genus	species
            for line in file:
                line=line.strip('\n')
                line=line.split('\t')
                if len(line)>2:
                    for i in range(1,len(line)):
                        if line[i]=='unknown': line[i]=None
                        elif line[i].startswith('uc_'):line[i]=None
                    contig,superkingdom,phylum,class_str,order,family,genus,species=line
                    self.contigs_taxonomy[contig]=superkingdom,phylum,class_str,order,family,genus,species







    def calculate_contig_distance_vector(self,contigs_set,target_contig):
        contigs_distance=[]
        already_added=set()
        to_fetch=contigs_set.union({target_contig})
        contigs_vectors = self.fetch_kmer_depth_vectors(to_fetch)
        contigs_vectors_keys=list(contigs_vectors.keys())
        contigs_vectors_keys.remove(target_contig)
        contigs_array = np.array([contigs_vectors[i] for i in contigs_vectors])
        candidate_vector = contigs_vectors[target_contig].reshape(1, -1)
        similarity = cosine_similarity(contigs_array,candidate_vector)
        for c1 in range(len(contigs_vectors_keys)):
            contig1=contigs_vectors_keys[c1]
            contig_key=self.get_contig_key(contig1,target_contig)
            if contig_key not in already_added:
                taxa_distance = self.invalid_score
                if contig1 in self.contigs_taxonomy and target_contig in self.contigs_taxonomy:
                    taxa_distance = self.calculate_taxa_distance(self.contigs_taxonomy[contig1],self.contigs_taxonomy[target_contig])
                kmer_depth_distance=similarity[c1][0]
                contigs_distance.append([contig_key, kmer_depth_distance, taxa_distance])
                already_added.add(contig_key)
        return contigs_distance

    def calculate_contig_distance_matrix(self,contigs_set):
        #inverse distance, 0 for maximum distance, 1 for minimum distance
        to_fetch=set()
        contigs_distance=[]
        already_added=set()
        contigs_vectors=self.fetch_kmer_depth_vectors(contigs_set)
        contigs_vectors_keys=list(contigs_vectors.keys())
        contigs_array=np.array([contigs_vectors[i] for i in contigs_vectors])
        similarity=cosine_similarity(contigs_array)
        for c1 in range(len(contigs_vectors_keys)):
            contig1=contigs_vectors_keys[c1]
            for c2 in range(len(contigs_vectors_keys)):
                contig2 = contigs_vectors_keys[c2]
                if contig1!=contig2:
                    contig_key=self.get_contig_key(contig1,contig2)
                    if contig_key not in already_added:
                        taxa_distance = self.invalid_score
                        if contig1 in self.contigs_taxonomy and contig2 in self.contigs_taxonomy:
                            taxa_distance = self.calculate_taxa_distance(self.contigs_taxonomy[contig1],self.contigs_taxonomy[contig2])
                        kmer_depth_distance=similarity[c1][c2]
                        contigs_distance.append([contig_key, kmer_depth_distance, taxa_distance])
                        already_added.add(contig_key)
        return contigs_distance


    def calculate_contig_distance_worker_function(self,queue, master_pid):
        res=[]
        #since we need to access the sqlite db we open a connection to the db per worker
        self.restart_sqlite_cursor()
        while True:
            record = queue.pop(0)
            if record is None:
                self.stop_sqlite_cursor()
                return
            if len(record)==1:
                contigs_set=record[0]
                self.mp_results.append(self.calculate_contig_distance_matrix(contigs_set))
            elif len(record)==2:
                contigs_set,target_contig=record
                self.mp_results.append(self.calculate_contig_distance_vector(contigs_set,target_contig))



    def yield_contigs(self, contigs,target_contig=None):
        contigs=list(contigs)
        step=self.numpy_array_size
        if target_contig:
            for i in range(0, len(contigs), step):
                yield [set(contigs[i:i + step]),target_contig]
        else:
            for i in range(0, len(contigs), step):
                yield [set(contigs[i:i + step])]


    def launch_calculate_distance(self,contigs_keys,target_contig=None):
        clean_contigs_keys=set()
        for contig_key in contigs_keys:
            contig1,contig2=self.get_original_contigs(contig_key)
            clean_contigs_keys.add(contig1)
            clean_contigs_keys.add(contig2)
        if target_contig:
            if target_contig in clean_contigs_keys:
                clean_contigs_keys.remove(target_contig)
        for chunk in self.yield_contigs(clean_contigs_keys,target_contig):
            self.queue.append(chunk)
        self.processes_handler(self.calculate_contig_distance_worker_function)
        while self.mp_results:
            #contigs distance cant be a generator because of multiprocessing. maybe it would work if we would apply multiprocessing before?
            contigs_distance=self.mp_results.pop(0)
            contigs_keys_list, kmer_depth_distance_list, taxa_distance_list = [], [], []
            for contig_key, kmer_depth_distance, taxa_distance in contigs_distance:
                contigs_keys_list.append(contig_key)
                kmer_depth_distance_list.append(kmer_depth_distance)
                taxa_distance_list.append(taxa_distance)
            self.store_distance(contigs_keys_list,kmer_depth_distance_list,taxa_distance_list)





    ####paired reads

    def get_paired_reads(self):
        read_connected_contigs = {}
        paired_contigs=set()
        #to avoid missing index file message
        pysam.set_verbosity(0)
        samfile = pysam.AlignmentFile(self.bam_file, "rb")

        '''
        try:
            print(f'Checking if index is available for {self.bam_file}')
            fetch_generator=samfile.fetch()
            print(f'Index is available')
        except:
            print(f'Creating index')
            pysam.index(self.bam_file)
            print(f'Index created')
            samfile = pysam.AlignmentFile(self.bam_file, "rb")
            fetch_generator=samfile.fetch()
        '''
        for read in samfile.fetch(until_eof=True):
            if read.is_paired and read.is_read1 and not read.is_unmapped and not read.is_secondary and not read.mate_is_unmapped:
                read1_contig = read.reference_name
                read2_contig = read.next_reference_name
                if read1_contig != read2_contig:

                    try:            read_connected_contigs[read1_contig].add(read2_contig)
                    except:         read_connected_contigs[read1_contig] = {read2_contig}

                    try:            read_connected_contigs[read2_contig].add(read1_contig)
                    except:         read_connected_contigs[read2_contig] = {read1_contig}
                    paired_contigs.add(read1_contig)
                    paired_contigs.add(read2_contig)

        return read_connected_contigs,paired_contigs


    ####get_bins_and_unbinned

    def get_seqs_generator(self,fasta_path):
        with open(fasta_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    query=line.replace('>','').strip()
                    yield query

    def map_bins(self,bins):
        #some bins have very long names, so we map them to ints to reduce memory footprint
        c=1
        self.bins_map={}
        self.bins={}
        for bin_id in bins:
            self.bins_map[c]=bin_id
            self.bins[c]=bins[bin_id]
            c+=1

    def add_paired_reads_contigs(self,bins, unbinned, read_connected_contigs):
        newly_binned = set()
        before_len_bins = 0
        after_len_bins = sum([len(bins[i]) for i in bins])
        while before_len_bins != after_len_bins:
            before_len_bins = sum([len(bins[i]) for i in bins])
            for unbinned_contig in unbinned:
                if unbinned_contig not in newly_binned:
                    if unbinned_contig in read_connected_contigs:
                        pair_contig = read_connected_contigs[unbinned_contig]
                        pair_contig.add(unbinned_contig)
                        current_candidates = {}
                        for bin_id in bins:
                            bin_contigs = bins[bin_id]
                            # checking the amount of contig & paired contigs in the bin
                            intersection_unbinned = len(pair_contig.intersection(bin_contigs))
                            current_candidates[bin_id] = intersection_unbinned

                        best_candidates = []
                        max_value = max(current_candidates.values())
                        if max_value:
                            for c in current_candidates:
                                if current_candidates[c] == max_value: best_candidates.append(c)
                            if len(best_candidates) > 1:
                                print(f'There are ambiguous contigs pairing (n={len(best_candidates)}) for the contig {unbinned_contig} so we are ignoring this assignment')
                            else:
                                bins[best_candidates[0]].add(unbinned_contig)
                                newly_binned.add(unbinned_contig)
            after_len_bins = sum([len(bins[i]) for i in bins])

        return newly_binned

    def get_bins_and_unbinned(self,valid_contigs,read_connected_contigs):
        '''
        :param valid_contigs:  contigs above the minimum contig length threshold or minimum marker contig length threshold (if these are marker genes)

        '''

        contigs_generator=self.get_seqs_generator(self.assembly_fasta)
        unbinned=set()
        bins={}
        binned=set()
        with open(self.binned_contigs) as file:
            for line in file:
                line=line.strip('\n')
                if line and not line.startswith('@'):
                    line=line.split('\t')
                    contig_id,bin_id=line[0],line[1]
                    if bin_id not in bins: bins[bin_id]=set()
                    bins[bin_id].add(contig_id)
                    binned.add(contig_id)
        for contig_id in contigs_generator:
            if contig_id not in binned:
                unbinned.add(contig_id)
        print(f'There is a total of {len(bins)} bins in this assembly, with a total of {len(binned)} contigs already binned')
        print(f'There is a total of {len(unbinned)} contigs to bin')
        unbinned=unbinned.intersection(valid_contigs)
        print(f'There is a total of {len(unbinned)} VALID contigs to bin')
        for bin_id in bins:
            bins[bin_id]=bins[bin_id].intersection(valid_contigs)
        newly_binned=self.add_paired_reads_contigs(bins,unbinned,read_connected_contigs)
        unbinned=unbinned-newly_binned
        self.map_bins(bins)
        self.unbinned=unbinned


    ####distance calculation
    def calculate_taxa_distance(self,contig1_lineage,contig2_lineage):
        '''
        species to genus = 0.5
        genus to family = 1
        family to order = 1.5
        order to class = 2
        class to phylum = 2.5
        phylum to superkingdom = 3
        maximum distance = 10.5
        '''
        max_distance=10.5
        distance=0
        jump=3
        for i in range(len(contig1_lineage)):
            if contig1_lineage[i] and contig2_lineage[i]:
                if contig1_lineage[i]!=contig2_lineage[i]:
                    distance+=jump
            if not contig1_lineage[i] or not contig2_lineage[i]:
                distance+=jump
            jump-=0.5
        return 1*distance/max_distance

    def get_contig_distance(self,contig1,contig2):
        avg_vector1, taxa_vector1 = contig1
        avg_vector2, taxa_vector2 = contig2
        taxa_distance = None
        if taxa_vector1 and taxa_vector2:
            taxa_distance = self.calculate_taxa_distance(taxa_vector1,taxa_vector2)
        kmer_depth_distance = cosine_similarity(np.array([avg_vector1,avg_vector2]))[0][1]
        return kmer_depth_distance, taxa_distance

    def get_contigs_taxa_lineage(self,contigs):
        taxonomy={}
        contigs_taxa=None
        map_taxa={}
        for contig in contigs:
            if contig in self.contigs_taxonomy:
                contig_taxonomy=str(self.contigs_taxonomy[contig])
                if contig_taxonomy not in map_taxa: map_taxa[contig_taxonomy]=self.contigs_taxonomy[contig]
                if contig_taxonomy not in taxonomy: taxonomy[contig_taxonomy]=0
                taxonomy[contig_taxonomy]+=self.contigs_len[contig]
        if taxonomy:
            contigs_taxa = max(taxonomy, key=taxonomy.get)
            contigs_taxa=map_taxa[contigs_taxa]
        return contigs_taxa



    def get_average_vector(self,contigs):
        vectors_array=self.fetch_kmer_depth_vectors_np(contigs)
        avg_vector=vectors_array.mean(axis=0)
        contigs_taxa_lineage=self.get_contigs_taxa_lineage(contigs)
        return avg_vector,contigs_taxa_lineage

    def get_avg_vector_bin(self):
        print('Calculating average vector per bin')
        #just to allow correct prints for tqdm
        sleep(0.1)
        for bin_id in progress_bar(self.bins):
            average_vector,contigs_taxa_lineage = self.get_average_vector(self.bins[bin_id])
            self.avg_vector_bin[bin_id]=average_vector,contigs_taxa_lineage




    ####parse_mantis

    def parse_mantis(self):
        res={}
        with open(self.annotations_tsv) as file:
            line=file.readline()
            if 'Links' in line: line=file.readline()
            while line:
                line=line.strip('\n')
                seq=line.split('\t')[0]
                annotations=''.join(line.split('|')[-1])
                annotations=[i for i in annotations.split('\t') if i]
                for db_id_annot in annotations:
                    db_id=db_id_annot.split(':')[0]
                    if db_id in self.protein_ids or db_id in self.reaction_ids:
                        annot=''.join(db_id_annot.split(':')[1:])
                        if db_id=='enzyme_ec' and annot.endswith('-'): annot=''
                        if annot:
                            if seq not in res: res[seq] = {}
                            if db_id not in res[seq]: res[seq][db_id]=set()
                            res[seq][db_id].add(annot)
                line=file.readline()
        self.annotations=res

    def run_preprocessing(self):
        start_ram=virtual_memory().available
        start=time()
        read_connected_contigs, paired_contigs = self.get_paired_reads()
        valid_contigs = self.assembly_matrix(paired_contigs)
        self.load_assembly_len()
        self.load_contigs_taxonomy()
        self.get_bins_and_unbinned(valid_contigs,read_connected_contigs)
        self.get_avg_vector_bin()
        end_ram=virtual_memory().available
        ram_cost=(start_ram-end_ram)/ (1024.0 ** 3)
        print(f'Pre-processing took {round(time()-start,self.rounding_digits)} seconds, and occupied {round(ram_cost,self.rounding_digits)} GBs of RAM.')

