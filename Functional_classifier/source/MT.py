import statistics
import os
import pickle

mapping='D:/Data/Annotation_Output/mapping/'

def file_len(fname):
    try:
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    except: return 0


def save_metrics(pickle_path,quality_genomes):
    with open(pickle_path, 'wb') as handle:
        pickle.dump(quality_genomes, handle,protocol=4)

def load_metrics(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as handle:
            quality_genomes = pickle.load(handle)
            return quality_genomes


def read_summary(summary_file):
    res=[]
    c=0
    with open(mapping+summary_file) as file:
        line=file.readline()
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            date,mapped,unmapped,total_reads=line
            res.append([date,total_reads])
            line=file.readline()
    print(res)
    return res

def read_feature_counts(mapping_file):
    res=load_metrics(mapping+'mt_counts.pickle')
    if res: return res
    res={}
    with open(mapping+mapping_file) as file:
        line=file.readline()
        line=file.readline()
        feature_counts,sample,gene_id=None,None,None
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            gene_id,sample=line[0],line[1]
            if feature_counts:
                if sample not in res:
                    res[sample]={}
                res[sample][gene_id]=feature_counts
            feature_counts=line[4:]
            feature_counts=[int(i) for i in feature_counts]
            line=file.readline()
    save_metrics(mapping+'mt_counts.pickle',res)
    return res

#mapping_file='ALL.mt.annotation.featureCounts_final.txt'
#feature_counts=read_feature_counts(mapping_file)
#read_summary(summary_file='ALL.mt.annotation.featureCounts.txt.summary.txt')

def main():
    l=[]
    pickle_path = 'D:/Data/Annotation_Output/reference_proteomes/reference_all.pickle'
    intermediate_pickles='D:/Data/Annotation_Output/reference_proteomes/quality_proteomes/'
    res=load_metrics(pickle_path)
    for organism in res:
        N=sum(res[organism])
        precision=(res[organism][0]+res[organism][1])/N
        l.append(precision)
        print(organism,precision)
    print(statistics.mean(l))
    print(statistics.median(l))
    print(statistics.stdev(l))
    print(min(l))
    print(max(l))
main()