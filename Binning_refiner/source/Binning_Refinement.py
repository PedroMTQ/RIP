import os
import re
import pysam
from sys import path as sys_path
if 'source' in sys_path[0]:
    from re import search as re_search
    p_to_add=sys_path[0][:re_search('Binning_refiner',sys_path[0]).span()[1]]
    if not sys_path[0].endswith('Binning_refiner') :  sys_path[0]=p_to_add

from source.Neo4j_Connector import *
from source.Pre_Processor import *
from source.SQLITE_Connector import *
from source.util import check_environment_cores,MAIN_FOLDER,KEGG_MODULES
from math import log10,log2,log


#for pre_processor
import logging
from time import time,sleep
import numpy as np
import pandas as pd
from scipy.spatial import distance
from scipy.special import softmax
from psutil import virtual_memory
#for contig to candidate contig distance
from tqdm import tqdm as progress_bar
import pickle
from math import ceil





class Binning_Refinement(Neo4j_Connector,Pre_Processor,SQLITE_Connector):

    def __init__(self,assembly_fasta, gff_annotation,bam_file, depth_tsv,binned_contigs, annotations_tsv,taxonomy_tsv,output_folder,min_cont_size=1000,min_annot_cont_size=300,db_name=None):
        Neo4j_Connector.__init__(self,db_name=db_name)
        Pre_Processor.__init__(self,assembly_fasta=assembly_fasta,gff_annotation=gff_annotation,bam_file=bam_file,
                               depth_tsv=depth_tsv,binned_contigs=binned_contigs,annotations_tsv=annotations_tsv,taxonomy_tsv=taxonomy_tsv,
                               min_cont_size=min_cont_size,min_annot_cont_size=min_annot_cont_size,kmer_sizes = [2, 3, 4])
        self.ready_mp()
        self.connect_to_databases()
        #from kofam, uniprot_ec and rhea
        self.protein_ids=['enzyme_ec','kegg_ko']
        #from rhea
        self.reaction_ids=['biocyc_reaction','kegg_reaction']
        self.output_folder=output_folder
        self.output_file=f'{self.output_folder}output.tsv'
        self.kegg_modules_tree = load_metrics(KEGG_MODULES)
        SQLITE_Connector.__init__(self)


        self.penalty_round=0.05
        self.threshold_mcs=0.0
        self.threshold_distance=0.0
        self.threshold_candidate=0.0

    def ready_mp(self):
        self.manager = Manager()
        self.queue = self.manager.list()
        self.mp_results = self.manager.list()
        self.worker_count = check_environment_cores()




    ####calculate_idf

    def find_sign(self,string_to_search):
        # first search should be for one of these signs, which means directionality is set
        compiled_sign = re.compile(
            '( = | <=> | <= | => | → | ← | ↔ | ⇄ | <-> | <--> | <- | -> | [=-]?&gt; | &lt;=&gt; | &lt;[=-]? )')
        sign = re.search(compiled_sign, string_to_search)
        if sign:
            return sign
        else:
            # the second search should be for unknown directionality
            # this way it avoids recognizing unknown compounds "?" as the direction sign of a reaciton
            interrogation_sign = re.compile('\?')
            sign = re.search(interrogation_sign, string_to_search)
            if sign: return sign
        return None

    def uniform_sign(self,sign):
        both = re.compile('( = | <=> | ↔ | ⇄ | <-> | <--> | &lt;=&gt; )')
        right = re.compile('( => | → | -> | [=-]?&gt; )')
        left = re.compile('( <= | ← | <- | &lt;[=-]? )')
        if re.search(both, sign):  return '<=>'
        if re.search(right, sign):  return '=>'
        if re.search(left, sign):  return '<='

    def calculate_idf(self):
        res={}
        total_reactions= self.run_command_neo4j('MATCH (r:Reaction)--(c:Compound) return count(r) as totalR')[0]['totalR']
        total_connections= self.run_command_neo4j('MATCH (c:Compound)-[rel:HAS_CPD]-(r) RETURN id(c)')
        for i in total_connections:
            id_cpd=i['id(c)']
            if id_cpd not in res: res[id_cpd]=0
            res[id_cpd]+=1
        idf={}
        for id_cpd in res:
            #idf[id_cpd]=total_reactions/res[id_cpd]
            idf[id_cpd]=log(total_reactions/res[id_cpd],10)
            #idf[id_cpd]=round(log(total_reactions/res[id_cpd] ** (log(total_reactions/res[id_cpd],1.5) / log(100,1.5)),1.5),self.rounding_digits)
        self.idf=idf

    ####get_reactions_unassigned

    def find_matches_neo4j_protein(self,id_type,annot_id):
        neo4j_command=f'MATCH (p:Protein)--(i:Identifier) WHERE i.ID="{annot_id}" AND i.Database="{id_type}" WITH p' \
                      f' WITH p MATCH (p)--(r:Reaction) WITH r' \
                      f' MATCH (r)-[rc]-(c:Compound) WITH r,c,rc' \
                      f' MATCH (r)--(rs:Reaction_str) RETURN id(r),r.Reversible,rs.Reaction_str,id(c),rc.side'

        res=self.run_command_neo4j(neo4j_command)
        reaction_dict={}
        reactions_list=[]
        for reaction in res:
            reaction_id=reaction['id(r)']
            if reaction_id not in reaction_dict:
                reaction_dict[reaction_id]={'reversibility':reaction['r.Reversible'],'reaction_str':reaction['rs.Reaction_str'],'compounds':{'left':set(),'right':set()}}
            compound_id=reaction['id(c)']
            side=reaction['rc.side']
            reaction_dict[reaction_id]['compounds'][side].add(compound_id)

            reaction_str=reaction['rs.Reaction_str']
            sign = self.find_sign(reaction_str)
            if sign:
                sign = sign.group()
                sign = self.uniform_sign(sign)
            else: sign='<=>'
            reaction_dict[reaction_id]['sign']=sign

        for reaction_id in reaction_dict:
            is_reversible=reaction_dict[reaction_id]['reversibility']
            left_reaction=reaction_dict[reaction_id]['compounds']['left']
            right_reaction=reaction_dict[reaction_id]['compounds']['right']
            if left_reaction and right_reaction:
                reactions_list.append([reaction_id,is_reversible,left_reaction,right_reaction,sign])
        return reactions_list

    def find_matches_neo4j_reaction(self,id_type,annot_id):
        neo4j_command=f'MATCH (r:Reaction)--(i:Identifier) WHERE i.ID="{annot_id}" AND i.Database="{id_type}" WITH r' \
                      f' MATCH (r)-[rc]-(c:Compound) WITH r,c,rc' \
                      f' MATCH (r)--(rs:Reaction_str) RETURN id(r),r.Reversible,rs.Reaction_str,id(c),rc.side'
        res=self.run_command_neo4j(neo4j_command)
        reaction_dict={}
        reactions_list=[]
        for reaction in res:
            reaction_id=reaction['id(r)']
            if reaction_id not in reaction_dict:
                reaction_dict[reaction_id]={'reversibility':reaction['r.Reversible'],'reaction_str':reaction['rs.Reaction_str'],'compounds':{'left':set(),'right':set()}}
            compound_id=reaction['id(c)']
            side=reaction['rc.side']
            reaction_dict[reaction_id]['compounds'][side].add(compound_id)

            reaction_str=reaction['rs.Reaction_str']
            sign = self.find_sign(reaction_str)
            if sign:
                sign = sign.group()
                sign = self.uniform_sign(sign)
            else: sign='<=>'
            reaction_dict[reaction_id]['sign']=sign

        for reaction_id in reaction_dict:
            is_reversible=reaction_dict[reaction_id]['reversibility']
            left_reaction=reaction_dict[reaction_id]['compounds']['left']
            right_reaction=reaction_dict[reaction_id]['compounds']['right']
            if left_reaction and right_reaction:
                reactions_list.append([reaction_id,is_reversible,left_reaction,right_reaction,sign])
        return reactions_list



    def get_reactions_list(self,id_type,annot_id,instance_type):
        if instance_type=='Reaction':
            if id_type == 'biocyc_reaction':                id_type = 'biocyc'
            elif id_type == 'kegg_reaction':                id_type = 'kegg'

        if instance_type not in self.annotations_reactions: self.annotations_reactions[instance_type]={}
        if id_type not in self.annotations_reactions[instance_type]:
            self.annotations_reactions[instance_type][id_type]={}
        if annot_id in self.annotations_reactions[instance_type][id_type]:
            reactions_list = self.annotations_reactions[instance_type][id_type][annot_id]
            return reactions_list
        reactions_list=[]
        #maybe we can accelerate this with better syntax
        if instance_type=='Protein':
            reactions_list=self.find_matches_neo4j_protein(id_type,annot_id)
        elif instance_type=='Reaction':
            if id_type == 'biocyc_reaction':                id_type = 'biocyc'
            elif id_type == 'kegg_reaction':                id_type = 'kegg'
            reactions_list=self.find_matches_neo4j_reaction(id_type,annot_id)
        return reactions_list


    def get_reactions_unassigned(self):

        res={}
        unbinned_reactions=0
        print('Getting unbinned contigs reactions')
        sleep(0.1)
        for contig_id in progress_bar(self.unbinned):
            contig_seqs=self.contigs_to_cds.get(contig_id)
            if contig_seqs:
                if contig_id not in res:res[contig_id]={}
                for seq_id in contig_seqs:
                    if seq_id not in res[contig_id]: res[contig_id][seq_id]={}
                    seq_annotations=self.annotations.get(seq_id)
                    if seq_annotations:
                        for id_type in seq_annotations:
                            for annot_id in seq_annotations[id_type]:
                                current_id= f'{id_type}{self.key_split}{annot_id}'
                                if id_type in self.protein_ids: instance_type='Protein'
                                elif id_type in self.reaction_ids: instance_type='Reaction'
                                reactions = self.get_reactions_list(id_type, annot_id,instance_type)
                                res[contig_id][seq_id][current_id] = reactions
                                unbinned_reactions+=len(reactions)
            #export to sqlite here
        return res,unbinned_reactions



    ####get_candidates

    def calculate_bin_to_candidate_contig_distance(self,bin_id,candidate_contig_id):
        kmer_depth_distance, taxa_distance=self.get_contig_distance(self.avg_vector_bin[bin_id],self.get_average_vector({candidate_contig_id}))
        if kmer_depth_distance and taxa_distance:
            return kmer_depth_distance,taxa_distance
        elif kmer_depth_distance and not taxa_distance:
            return kmer_depth_distance,self.invalid_score
        else:
            print('Error here',kmer_depth_distance,taxa_distance,bin_id,candidate_contig_id)
            return self.invalid_score,self.invalid_score

    ####get_bin_details

    def merge_bin_annotations(self,bin_id):
        res={}
        for contig_id in self.bins[bin_id]:
            contig_seqs = self.contigs_to_cds.get(contig_id)
            if contig_seqs:
                for seq_id in contig_seqs:
                    seq_annotations = self.annotations.get(seq_id)
                    if seq_annotations:
                        for id_type in seq_annotations:
                            if id_type not in res: res[id_type]=set()
                            for annot_id in seq_annotations[id_type]:
                                res[id_type].add(annot_id)
        return res

    def get_metabolites_bin(self,all_annotations):
        bin_substrates=set()
        bin_products=set()
        set_annotations=set()
        for id_type in all_annotations:
            if id_type in self.protein_ids:                instance_type = 'Protein'
            elif id_type in self.reaction_ids:             instance_type = 'Reaction'

            for annot_id in all_annotations[id_type]:
                set_annotations.add(f'{id_type}{self.key_split}{annot_id}')
                reactions_list = self.get_reactions_list(id_type, annot_id, instance_type)
                for reaction_id,is_reversible,left_reaction,right_reaction,sign in reactions_list:
                    if sign=='<=>':
                        for cpd_id in left_reaction: bin_substrates.add(cpd_id)
                        for cpd_id in right_reaction: bin_substrates.add(cpd_id)
                        for cpd_id in left_reaction: bin_products.add(cpd_id)
                        for cpd_id in right_reaction: bin_products.add(cpd_id)
                    elif sign=='=>':
                        for cpd_id in left_reaction: bin_substrates.add(cpd_id)
                        for cpd_id in right_reaction: bin_products.add(cpd_id)

                    elif sign=='<=':
                        for cpd_id in right_reaction: bin_substrates.add(cpd_id)
                        for cpd_id in left_reaction: bin_products.add(cpd_id)
        return bin_substrates,bin_products,set_annotations

    def get_bin_details(self):
        self.bin_details={}
        print('Getting bins metabolic information')
        sleep(0.1)
        for bin_id in progress_bar(self.bins):
            start = time()
            #mcs= metabolite connectivity score
            #https://academic.oup.com/bioinformatics/article/32/6/867/1744247
            bin_annotations=self.merge_bin_annotations(bin_id)
            bin_substrates,bin_products, set_annotations = self.get_metabolites_bin(bin_annotations)
            self.bin_details[bin_id]={'bin_substrates':bin_substrates,'bin_products':bin_products,'set_annotations':set_annotations}




    ####assign_contigs


    def get_candidate_score(self,mcs,kegg_module_score,kmer_depth_dist,taxa_score):
        '''
        :param mcs: 0 to 1, 1 is best
        :param kmer_depth_dist: inverse kmer and depth distance from 0 to 1, 1 being very similar kmer and depth , 0 being very different
        :param taxa_score: inverse taxa score from 0 to 1, 1 being the same species, 0 different superkingdoms
        '''
        if taxa_score!=self.invalid_score:
            score= mcs * kmer_depth_dist * taxa_score
        else:
            score = mcs * kmer_depth_dist
        return score

    def clean_candidates_scores(self,mcs_score,kegg_modules_score,contig_id):

        res={}
        file= open(f'{test_file}_scores.txt','a+')
        for bin_id in mcs_score:
            if bin_id in mcs_score:

                #we should probably already have a threshold here to reduce runtime? maybe get the top mcs candidates instead of checking all??
                bin_mcs_score=mcs_score[bin_id]
                if bin_id in kegg_modules_score:
                    bin_kegg_module_score=kegg_modules_score[bin_id]
                else: bin_kegg_module_score=self.invalid_score

                if bin_mcs_score>=self.threshold_mcs:
                    bin_dist_kmer_depth, bin_dist_taxa = self.calculate_bin_to_candidate_contig_distance(bin_id, contig_id)

                    print(contig_id,self.bins_map[bin_id],round(bin_mcs_score, 3), round(bin_kegg_module_score, 3), round(bin_dist_kmer_depth, 3), round(bin_dist_taxa, 3),file=file, flush=True)



                    if bin_dist_kmer_depth >= self.threshold_distance:
                        candidate_score=self.get_candidate_score(bin_mcs_score,bin_kegg_module_score, bin_dist_kmer_depth, bin_dist_taxa)
                        #we need some threshold here

                        #print(self.bins_map[bin_id], contig_id, round(bin_mcs_score, self.rounding_digits),round(bin_dist_kmer_depth, self.rounding_digits),candidate_score)

                        if candidate_score:
                            res[bin_id]=candidate_score
        return res

    def get_mcs_bins_score(self,reaction_candidates):
        res={}
        for contig_id in reaction_candidates:
            sum_bin_score = {}
            n_reactions=0
            for seq_id in reaction_candidates[contig_id]:
                for id_type in reaction_candidates[contig_id][seq_id]:
                    for annot_id in reaction_candidates[contig_id][seq_id][id_type]:
                        for reaction_id in reaction_candidates[contig_id][seq_id][id_type][annot_id]:
                            candidates = reaction_candidates[contig_id][seq_id][id_type][annot_id][reaction_id]
                            n_reactions+=1
                            for bin_id,bin_score in candidates:
                                if bin_id not in sum_bin_score: sum_bin_score[bin_id]=0
                                sum_bin_score[bin_id]+=bin_score
            for bin_id in sum_bin_score:
                sum_bin_score[bin_id]/=n_reactions
            res[contig_id]=sum_bin_score
        return res

    def assign_contigs_to_best_bins(self,reaction_candidates,kegg_module_candidates):
        assigned_contigs={}
        all_contigs=set(reaction_candidates.keys()).union(kegg_module_candidates.keys())
        reaction_candidates=self.get_mcs_bins_score(reaction_candidates)

        #need to paralelyze
        print('Assigning contigs to best candidate bins')
        sleep(0.1)



        for contig_id in all_contigs:
            if contig_id in kegg_module_candidates:
                kegg_modules_score=self.calculate_softmax_scores(kegg_module_candidates[contig_id])
            else:
                kegg_modules_score={}

            if contig_id in reaction_candidates:
                mcs_score=self.calculate_softmax_scores(reaction_candidates[contig_id])
            else:
                mcs_score={}



            bins_score=self.clean_candidates_scores(mcs_score,kegg_modules_score,contig_id)

            #just in case we have ties, this is a list of best candidates
            if bins_score:
                best_candidates=[]
                max_value=max(bins_score.values())
                if max_value:
                    for c in bins_score:
                        if bins_score[c]==max_value: best_candidates.append(c)
                    if len(best_candidates)>1:
                        print(f'There are ambiguous assignments (n={len(best_candidates)}) for the contig {contig_id} so we are ignoring this assignment')
                    else:
                        assigned_contigs[contig_id] = best_candidates

        return assigned_contigs



    def calculate_scaled_idf(self,reaction_cpds):
        idfs_vals=[self.idf[cpd] for cpd in reaction_cpds]
        idfs_vals=softmax(np.array(idfs_vals))
        res={}
        for i in range(len(reaction_cpds)):
            res[reaction_cpds[i]] = idfs_vals[i]
        return res

    def calculate_softmax_scores(self,scores_dict):
        if not scores_dict: return {}
        scores_dict_keys=list(scores_dict.keys())
        dict_vals=[scores_dict[bin_id] for bin_id in scores_dict_keys]
        dict_vals=softmax(np.array(dict_vals))
        res={}
        for i in range(len(scores_dict_keys)):
            res[scores_dict_keys[i]] = dict_vals[i]
        return res



    def calculate_mcs(self,left_reaction,right_reaction,bin_substrates,bin_products):
        mcs_sub, mcs_prod = 0, 0
        list_cpds=list(left_reaction.union(right_reaction))
        reaction_idfs=self.calculate_scaled_idf(list_cpds)
        #for transport reactions
        cpds_in_both=left_reaction.intersection(right_reaction)
        for cpd in cpds_in_both:
            reaction_idfs[cpd]/=2
        #right_reaction_idfs=self.calculate_scaled_idf(right_reaction)
        for cpd in left_reaction:
            if cpd in bin_substrates:
                #mcs_sub += 1
                #print(cpd.get_most_common_synonym(),self.calculate_compound_weight(cpd))
                #mcs_sub += self.calculate_compound_weight(cpd)
                mcs_sub += reaction_idfs[cpd]
            #else:
            #    print(cpd.get_most_common_synonym(),self.calculate_compound_weight(cpd))
            #    mcs_sub -= self.calculate_compound_weight(cpd)/len(left_reaction)

        for cpd in right_reaction:
            if cpd in bin_products:
                #mcs_prod += 1
                #print(cpd.get_most_common_synonym(),self.calculate_compound_weight(cpd))
                mcs_prod += reaction_idfs[cpd]
            #else:
            #    print(cpd.get_most_common_synonym(),self.calculate_compound_weight(cpd))
            #    mcs_prod -= self.calculate_compound_weight(cpd)/len(right_reaction)
        #mcs_sub /= len(left_reaction)
        #mcs_prod /= len(right_reaction)
        mcs = mcs_sub + mcs_prod
        #if mcs<1/(len(left_reaction)+len(right_reaction)): mcs=0
        return mcs

    def get_valid_contigs_per_candidate(self,set_annotations,unassigned_reactions):
        valid_contigs=set()
        for contig_id in unassigned_reactions:
            valid_contig=True
            for seq_id in unassigned_reactions[contig_id]:
                for id_type_annot_id in unassigned_reactions[contig_id][seq_id]:
                    #we only add contigs which could be missing in the bins
                    if  id_type_annot_id in set_annotations:
                        valid_contig=False
                    #if id_type_annot_id.startswith('enzyme_ec') and id_type_annot_id in set_annotations:
                    #    valid_contig=False
            if valid_contig: valid_contigs.add(contig_id)
        return valid_contigs


    #kegg modules

    def get_best_sample_module_path(self,sample_kos,all_paths):
        best_score=0
        best_path=set()
        for current_path in all_paths:
            available_kos=current_path.intersection(sample_kos)
            current_score=len(available_kos)/len(current_path)
            if current_score>best_score:
                best_score=current_score
                best_path=current_path
        available_kos = best_path.intersection(sample_kos)
        missing_kos = best_path.difference(available_kos)
        return best_score,best_path,available_kos,missing_kos

    def get_kegg_module_score(self,contig_kos,bin_kos):
        res = {}
        bin_score=[]
        temp_kos= bin_kos.union(contig_kos)
        for main_module_name in self.kegg_modules_tree:
            for sub_module_name in self.kegg_modules_tree[main_module_name]:
                for module_key in self.kegg_modules_tree[main_module_name][sub_module_name]:
                    module_name,module_paths = self.kegg_modules_tree[main_module_name][sub_module_name][module_key]
                    module_score,module_path,available_kos,missing_kos = self.get_best_sample_module_path(bin_kos,module_paths)
                    res[module_key]=module_score,module_path,available_kos,missing_kos
        for module_key in res:
            module_score,module_path,available_kos,missing_kos=res[module_key]
            if module_score:
                score_before=len(available_kos)/len(module_path)
                intersection_bin_contig=missing_kos.intersection(contig_kos)
                if intersection_bin_contig:
                    score_after=(len(available_kos)+len(intersection_bin_contig))/len(module_path)
                    bin_score.append(score_after-score_before) # after 80% before 50% = improvement 30%
                else:
                    bin_score.append(0)
        avg_improvement=sum(bin_score)/len(bin_score)
        return avg_improvement


    def add_candidates_kegg_module(self,bin_id,res,bin_annotations,valid_contigs,unassigned_reactions,refine_round):
        sample_kos=set()
        for id_type_annot_id in bin_annotations:
            id_type, annot_id = id_type_annot_id.split(self.key_split)
            if id_type=='kegg_ko':
                sample_kos.add(annot_id)
        for contig_id in valid_contigs:
            if contig_id not in res: res[contig_id] = {}
            contig_kos=set()
            contig_kegg_module_score=[0,0]
            for seq_id in unassigned_reactions[contig_id]:
                for id_type_annot_id in unassigned_reactions[contig_id][seq_id]:
                    id_type,annot_id=id_type_annot_id.split(self.key_split)
                    if id_type == 'kegg_ko':
                        contig_kos.add(annot_id)
            score_kegg_module =self.get_kegg_module_score(contig_kos,sample_kos)
            if score_kegg_module:
                res[contig_id][bin_id]=score_kegg_module



    def add_candidates_mcs(self,bin_id,res,bin_annotations,valid_contigs,unassigned_reactions,bin_substrates,bin_products,refine_round):
        for contig_id in valid_contigs:
            for seq_id in unassigned_reactions[contig_id]:
                for id_type_annot_id in unassigned_reactions[contig_id][seq_id]:
                    if id_type_annot_id not in bin_annotations:
                        id_type,annot_id=id_type_annot_id.split(self.key_split)
                        reactions=unassigned_reactions[contig_id][seq_id][id_type_annot_id]
                        for reaction_id,is_reversible,left_reaction,right_reaction,sign in reactions:
                            mcs=self.calculate_mcs(left_reaction,right_reaction,bin_substrates,bin_products)
                            if is_reversible:
                                reversible_mcs=self.calculate_mcs(right_reaction,left_reaction,bin_substrates,bin_products)
                                if reversible_mcs>mcs: mcs=reversible_mcs
                            mcs=mcs#-mcs*self.penalty_round*refine_round
                            if mcs>self.threshold_mcs:
                                if contig_id not in res: res[contig_id]={}
                                if seq_id not in res[contig_id]: res[contig_id][seq_id]={}
                                if id_type not in res[contig_id][seq_id]: res[contig_id][seq_id][id_type]={}
                                if annot_id not in res[contig_id][seq_id][id_type]: res[contig_id][seq_id][id_type][annot_id]={}
                                if reaction_id not in res[contig_id][seq_id][id_type][annot_id]: res[contig_id][seq_id][id_type][annot_id][reaction_id]=[]
                                mcs=round(mcs,self.rounding_digits)
                                res[contig_id][seq_id][id_type][annot_id][reaction_id].append([bin_id,mcs])


    def add_candidates(self,reaction_candidates,kegg_module_candidates,bin_id,bin_substrates,bin_products, set_annotations,unassigned_reactions,refine_round):
        #the unassigned contigs we need to do one at a time, this way we can reuse the get_metabolites_bin function
        valid_contigs=self.get_valid_contigs_per_candidate(set_annotations,unassigned_reactions)
        self.add_candidates_kegg_module(bin_id,kegg_module_candidates,set_annotations,valid_contigs,unassigned_reactions,refine_round)
        self.add_candidates_mcs(bin_id,reaction_candidates,set_annotations,valid_contigs,unassigned_reactions,bin_substrates,bin_products,refine_round)





    def assign_contigs(self,unassigned_reactions,refine_round):
        reaction_candidates={}
        kegg_module_candidates={}
        print('Getting candidate bins for contigs')
        sleep(0.1)
        for bin_id in progress_bar(self.bins):
            bin_details=self.bin_details[bin_id]
            bin_substrates=bin_details['bin_substrates']
            bin_products=bin_details['bin_products']
            set_annotations=bin_details['set_annotations']
            self.add_candidates(reaction_candidates,kegg_module_candidates,bin_id,bin_substrates,bin_products, set_annotations,unassigned_reactions,refine_round)
        assigned_contigs=self.assign_contigs_to_best_bins(reaction_candidates,kegg_module_candidates)
        return assigned_contigs

    ####add_to_bin_details

    def add_to_bin_details(self,assigned_contigs,unassigned_reactions):
        print('Adding assigned contigs to bins')
        sleep(0.1)
        for contig_id in progress_bar(assigned_contigs):
            contig_info = unassigned_reactions[contig_id]
            contig_substrates=[]
            contig_products=[]
            contig_annotations=set()
            for seq_id in unassigned_reactions[contig_id]:
                for annotation in unassigned_reactions[contig_id][seq_id]:
                    contig_annotations.add(annotation)
                    for reaction_id,is_reversible,left_reaction,right_reaction,sign in unassigned_reactions[contig_id][seq_id][annotation]:
                        contig_substrates.extend(left_reaction)
                        contig_products.extend(right_reaction)
                        if is_reversible:
                            contig_substrates.extend(right_reaction)
                            contig_products.extend(left_reaction)
            for bin_id in assigned_contigs[contig_id]:
                self.bin_details[bin_id]['bin_substrates'].update(contig_substrates)
                self.bin_details[bin_id]['bin_products'].update(contig_products)
                self.bin_details[bin_id]['set_annotations'].update(contig_annotations)

    ####export_refined_bins

    def bins_to_contigs(self,bins):
        res={}
        for bin_id in bins:
            for contig_id in bins[bin_id]:
                res[contig_id]=[bin_id]
        return res


    def reformat_refined_bins(self,refined_bins):
        res={}
        #res = self.bins_to_contigs(self.bins)
        for refine_round in refined_bins:
            for contig_id in refined_bins[refine_round]:
                res[contig_id]=[]
                for bin_id in refined_bins[refine_round][contig_id]:
                    res[contig_id].append(bin_id)
        return res

    def export_refined_bins(self,refined_bins):
        print('Exporting refined bins')
        sleep(0.1)
        contigs_binned=0
        with open(self.output_file,'w+') as file:
            first_line=['Contig','Bin']
            file.write('\t'.join(first_line)+'\n')
            for contig_id in progress_bar(refined_bins):
                line=[]
                contigs_binned+=1
                for bin_id in refined_bins[contig_id]:
                    mapped_bin_id=self.bins_map[bin_id]
                    line.append(mapped_bin_id)
                if line:
                    line.insert(0,contig_id)

                    line=[str(i) for i in line]
                    file.write('\t'.join(line) + '\n')
        try:    print(f'We managed to bin {contigs_binned} contigs out of {len(self.unbinned)} ({round(100*contigs_binned/len(self.unbinned),self.rounding_digits)}%) previously unbinned contigs.')
        except: pass


    def update_avg_vector_bin(self,assigned_contigs):
        print('Updating average vectors per bin')
        sleep(0.1)
        bins_with_all_contigs={}
        for assigned_contig in progress_bar(assigned_contigs):
            assigned_to_bins=assigned_contigs[assigned_contig]
            for bin_id in assigned_to_bins:
                if bin_id not in bins_with_all_contigs: bins_with_all_contigs[bin_id]=self.bins[bin_id]
                bins_with_all_contigs[bin_id].add(assigned_contig)
        for bin_id in bins_with_all_contigs:
            average_vector,contigs_taxa = self.get_average_vector(bins_with_all_contigs[bin_id])
            self.avg_vector_bin[bin_id]=average_vector,contigs_taxa




    def refine_bins(self):
        start_time = time()
        start=time()
        self.run_preprocessing()
        res = {}

        if self.unbinned:
            self.calculate_idf()

            start=time()
            start_ram=virtual_memory().available
            unassigned_reactions,unbinned_reactions=self.get_reactions_unassigned()
            end_ram = virtual_memory().available
            ram_cost = (start_ram - end_ram) / (1024.0 ** 3)
            print(f'Unbinned reactions fetching took {round(time()-start,self.rounding_digits)} seconds, and occupied {round(ram_cost,self.rounding_digits)} GBs of RAM.')

            if unbinned_reactions:
                start_ram = virtual_memory().available
                start=time()
                self.get_bin_details()
                end_ram = virtual_memory().available
                ram_cost = (start_ram - end_ram) / (1024.0 ** 3)
                print(f'Bin metabolic information fetching took {round(time() - start, self.rounding_digits)} seconds, and occupied {round(ram_cost, self.rounding_digits)} GBs of RAM.')

                #we dont need the annotations anymore so we just clear them
                self.clear_annotations_dicts()

                end_refiner=0
                start_refiner=len(unassigned_reactions.keys())
                refine_round=0
                while unassigned_reactions and start_refiner!=end_refiner:
                    start_refiner=len(unassigned_reactions.keys())
                    assigned_contigs = self.assign_contigs(unassigned_reactions,refine_round=refine_round)
                    self.add_to_bin_details(assigned_contigs,unassigned_reactions)
                    for contig_id in assigned_contigs:
                        unassigned_reactions.pop(contig_id)
                    end_refiner=len(unassigned_reactions.keys())
                    res[refine_round]=assigned_contigs
                    refine_round+=1
                    self.update_avg_vector_bin(assigned_contigs)
                #just removing refine_round key
                reformatted_res=self.reformat_refined_bins(res)
                self.export_refined_bins(reformatted_res)
                print('Binning Refinement took',time()-start_time,'seconds')
            else:
                print('No reactions to bin!')
            print('final RAM', (virtual_memory().total - virtual_memory().available) / (1024.0 ** 3))

        self.close_sql_connection()

        return res







if __name__ == '__main__':
    '''
    binny output for bins: cut -f1,2 H_S005_gsa_mapping_amber_hr_tax.tsv > bins.tsv

    '''

    depth_tsv = f'{f}contig_depth.txt'
    assembly_fasta = f"{f}assembly.fa"
    gff_annotation = f"{f}annotation_CDS_RNA_hmms_checkm.gff"
    bam_file = f"{f}sample.bam"
    #needs to be changed
    binned_contigs=f'{f}bins.tsv'
    annotations_tsv=f'{f}consensus_annotation.tsv'
    taxonomy_tsv=f'{f}taxonomy_lineage.tsv'

    binner=Binning_Refinement(assembly_fasta=assembly_fasta,gff_annotation=gff_annotation,bam_file=bam_file,depth_tsv=depth_tsv,binned_contigs=binned_contigs,annotations_tsv=annotations_tsv,taxonomy_tsv=taxonomy_tsv,output_folder=f)
    #l1,l2=(None, None, None, None, None, None, None), ('Bacteria', 'Bacteroidetes', None, None, None, None, None)
    #binner.calculate_taxa_distance(l1,l2)
    #print(1-distance.cosine([1, 0, 0], [1, 0, 0]))

    #binner.run_preprocessing()
    #binner.get_all_contig_keys()
    #binner.calculate_idf()
    binner.refine_bins()
    #benchmark_refiner(binned_contigs,f'{f}output.tsv',None)
