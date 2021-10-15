from scipy.stats import ttest_rel



folder='/small/'
refined_bins_tsv=f'{folder}output.tsv'
binner_bins_tsv=f'{folder}bins.tsv'
real_bins_tsv=f'{folder}real_bins.tsv'
test_scores_tsv=f'{folder}test_file_scores.txt'
assembly_path=f'{folder}assembly.fa'

def read_tsv(tsv_file,contig_col,bin_col):
    res={}
    with open(tsv_file) as file:
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n').split('\t')
            contig_id=line[contig_col]
            bin_id=line[bin_col]
            if bin_id not in res: res[bin_id]=set()
            res[bin_id].add(contig_id)
            line=file.readline()
    return res

def read_fasta_generator(fasta_path):
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

def load_assembly(assembly_path):
    res={}
    for query,seq in read_fasta_generator(assembly_path):
        res[query]=len(seq)
    return res

def get_contigs_length(contigs,assembly_info):
    res=0
    for contig in contigs:
        res+=assembly_info[contig]
    return res


def get_best_match(current_bin,real_bins,assembly_info):
    res={}
    for bin_id in real_bins:
        contigs_intersection=current_bin.intersection(real_bins[bin_id])
        contigs_intersection_len=get_contigs_length(contigs_intersection,assembly_info)
        real_bins_len=get_contigs_length(real_bins[bin_id],assembly_info)
        bin_score=contigs_intersection_len/real_bins_len
        res[bin_id]=bin_score
    best_bin=None
    best_bin_score=0
    for bin_id in res:
        if res[bin_id]>best_bin_score:
            best_bin=bin_id
            best_bin_score=res[bin_id]
    return best_bin

def get_stats(bins,real_bins,assembly_info,mapped_bins={}):
    all_assigned=[]
    if not mapped_bins:
        for bin_id in bins:
            bin_contigs=bins[bin_id]
            mapped_bin=get_best_match(bin_contigs,real_bins,assembly_info)
            all_assigned.append(mapped_bin)
            mapped_bins[bin_id]=mapped_bin
    if len(all_assigned)!=len(set(all_assigned)):
        print('multiple bins mapped to the same real bin')
    metrics_res={}
    for bin_id in mapped_bins:
        if bin_id in bins:
            real_bin_id=mapped_bins[bin_id]
            completeness,contamination,size_bin,size_real_bin=get_metrics_bin(bins[bin_id],real_bins[real_bin_id],assembly_info)
            metrics_res[bin_id]=completeness,contamination,size_bin,size_real_bin
    return metrics_res,mapped_bins


def get_metrics_bin(current_bin,real_bin,assembly_info):


    TP,TN,FP,FN=set(),set(),set(),set()
    for contig in current_bin:
        if contig in real_bin: TP.add(contig)
        elif contig not in real_bin: FP.add(contig)
    for contig in real_bin:
        if contig not in current_bin: FN.add(contig)
    TP=get_contigs_length(TP,assembly_info)
    FP=get_contigs_length(FP,assembly_info)
    FN=get_contigs_length(FN,assembly_info)
    completeness = TP/(TP+FN)
    purity=TP/(TP+FP)
    contamination= 1-purity
    completeness=round(completeness,4)
    contamination=round(contamination,4)
    return completeness,contamination,len(current_bin),len(real_bin)

def benchmark_refiner():
    assembly_info=load_assembly(assembly_path)
    refined_bins=read_tsv(refined_bins_tsv,contig_col=0,bin_col=1)
    binner_bins=read_tsv(binner_bins_tsv,contig_col=0,bin_col=1)
    for bin_id in refined_bins:
        refined_bins[bin_id].update(binner_bins[bin_id])
    real_bins=read_tsv(real_bins_tsv,contig_col=0,bin_col=1)
    metrics_binner,mapped_bins=get_stats(binner_bins,real_bins,assembly_info)
    metrics_refiner,_=get_stats(refined_bins,real_bins,assembly_info,mapped_bins=mapped_bins)

    average_contamination = [0, 0]
    all_contaminations=[[],[]]
    all_completeness=[[],[]]
    average_completeness = [0, 0]
    for i in metrics_refiner:
        completeness1,contamination1,size_bin1,size_real_bin1 = metrics_binner[i]
        completeness2,contamination2,size_bin2,size_real_bin2 = metrics_refiner[i]
        all_contaminations[0].append(contamination1)
        all_contaminations[1].append(contamination2)
        all_completeness[0].append(completeness1)
        all_completeness[1].append(completeness2)

        average_contamination[0] += 1
        average_contamination[1] += contamination2 - contamination1
        average_completeness[0] += 1
        average_completeness[1] += completeness2 - completeness1
        '''
        if completeness1<completeness2 and contamination1>contamination2:
            print('BETTER',completeness1,contamination1,completeness2,contamination2)
        else:
            if contamination1<contamination2:
                print('CONTAMINATION HIGHER',completeness1,contamination1,completeness2,contamination2)

            elif completeness1<completeness2:
                print('COMPLETENESS HIGHER',completeness1,contamination1,completeness2,contamination2)
        '''

        if completeness1!=completeness2 or contamination1!=contamination2:
            print('diff',completeness1, contamination1, completeness2, contamination2,i)
        else:
            print('same')
    if average_contamination[0]: print('AVERAGE CONTAMINATION',average_contamination[1]/average_contamination[0])
    if average_completeness[0]: print('AVERAGE COMPLETENESS',average_completeness[1]/average_completeness[0])
    print('#######')
    test_cont=ttest_rel(all_contaminations[0],all_contaminations[1])
    print('cont',test_cont)
    test_comp=ttest_rel(all_completeness[0],all_completeness[1])
    print('comp',test_comp)





def check_seed_results():
    assembly_info=load_assembly(assembly_path)
    binner_bins=read_tsv(binner_bins_tsv,contig_col=0,bin_col=1)
    real_bins=read_tsv(real_bins_tsv,contig_col=0,bin_col=1)
    metrics_binner,mapped_bins=get_stats(binner_bins,real_bins,assembly_info)
    for i in metrics_binner:
        completeness1,contamination1,size_bin1,size_real_bin1 = metrics_binner[i]
        print(i,mapped_bins[i],'\nbefore',completeness1,contamination1)
    '''
    test_scores=read_test_score(test_scores_tsv)
    for bin_id in metrics_binner:
        completeness1,contamination1,size_bin1,size_real_bin1 = metrics_binner[bin_id]
        if completeness1<0.7 and contamination1:
            print(bin_id, completeness1, contamination1, size_bin1, size_real_bin1)
            if bin_id in test_scores:
                for contig_id in test_scores[bin_id]:
                    candidate_score,bin_mcs_score,bin_dist_kmer_depth=test_scores[bin_id][contig_id]
                    print(contig_id,candidate_score,bin_mcs_score,bin_dist_kmer_depth)
            print('######')
    '''

def read_test_score(file_path):
    res={}
    with open(file_path) as file:
        for line in file:
            line=line.strip('\n').split()
            contig_id,bin_id,bin_mcs_score,bin_kegg_module_score,bin_dist_kmer_depth,bin_dist_taxa=line
            if bin_id not in res: res[bin_id]={}
            res[bin_id][contig_id]=float(bin_mcs_score),float(bin_kegg_module_score),float(bin_dist_kmer_depth),float(bin_dist_taxa)
    return res

def check_parameters():
    real_bins=read_tsv(real_bins_tsv,contig_col=0,bin_col=1)
    scores=read_test_score(test_scores_tsv)
    res={}
    for bin_id in real_bins:
        res[bin_id]={'mcs':[],'kegg':[],'kmer':[],'taxa':[]}
        for contig_id in real_bins[bin_id]:
            if contig_id in scores[bin_id]:
                bin_mcs_score, bin_kegg_module_score, bin_dist_kmer_depth, bin_dist_taxa = scores[bin_id][contig_id]
                if bin_mcs_score!=-1: res[bin_id]['mcs'].append(bin_mcs_score)
                if bin_kegg_module_score!=-1: res[bin_id]['kegg'].append(bin_kegg_module_score)
                if bin_dist_kmer_depth!=-1: res[bin_id]['kmer'].append(bin_dist_kmer_depth)
                if bin_dist_taxa!=-1: res[bin_id]['taxa'].append(bin_dist_taxa)
                print(bin_id,contig_id,bin_mcs_score)
        for k in ['mcs']:
        #for k in res[bin_id].keys():
            if len(res[bin_id][k]):
                avg_k=sum(res[bin_id][k])/len(res[bin_id][k])
                print(bin_id,k,avg_k)
    print('#######')
    all_res={'mcs':[],'kegg':[],'kmer':[],'taxa':[]}
    for bin_id in real_bins:
        for k in res[bin_id].keys():
            if len(res[bin_id][k]):
                avg_k=sum(res[bin_id][k])/len(res[bin_id][k])
                all_res[k].append(avg_k)
    for k in all_res:
        if len(all_res[k]):
            avg_k = sum(all_res[k]) / len(all_res[k])
            print(k,avg_k)




#check_seed_results()
#benchmark_refiner()
check_parameters()
