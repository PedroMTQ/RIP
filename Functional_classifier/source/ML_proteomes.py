import requests
from math import sqrt
import os
import re
#import plotly.graph_objects as go
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
#from kneed import KneeLocator
from sklearn.metrics import silhouette_score,davies_bouldin_score,classification_report
from sklearn import mixture as skmixture
from sklearn.preprocessing import StandardScaler
import statistics
import random
from sklearn import cluster as skcluster
from scipy.cluster import hierarchy
#import seaborn as sns; sns.set_theme()
from sklearn.decomposition import PCA
#import plotly.express as px
import sklearn.feature_selection as skfeature_select
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier,GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
#import plotly.figure_factory as ff
from scipy import stats

#https://www.uniprot.org/proteomes/?query=*&fil=reference%3Ayes+AND+taxonomy%3A%22Bacteria+%5B2%5D%22

linux_path = '/home/pedroq/Dropbox/'
windows_path = 'C:/Users/Pedro Queirós/Dropbox/'
if os.path.isdir(linux_path):
    path = linux_path
elif os.path.isdir(windows_path):
    path = windows_path
else:
    raise Exception

linux_path = '/home/pedroq/Python_projects/DRAX/'
#windows_path = 'C:/Users/Pedro Queirós/Documents/Python Projects/DRAX/'
windows_path = 'D:/'
if os.path.isdir(linux_path):
    drax_path = linux_path
elif os.path.isdir(windows_path):
    drax_path = windows_path
else:
    raise Exception

wanted_level=2

cluster_path='D:/Data/Annotation_Output/reference_proteomes/cluster.pickle'
ml_results='ml_results.pickle'
kegg_map='C:/Users/Pedro Queirós/Dropbox/PhD/Presentations/cet_2/kegg_map_fixed.txt'

def save_metrics(pickle_path,quality_genomes):
    with open(path+'/PhD/Presentations/cet_2/ML/'+pickle_path, 'wb') as handle:
        pickle.dump(quality_genomes, handle,protocol=4)

def load_metrics(pickle_path):
    if os.path.exists(path+'/PhD/Presentations/cet_2/ML/'+pickle_path):
        with open(path+'/PhD/Presentations/cet_2/ML/'+pickle_path, 'rb') as handle:
            quality_genomes = pickle.load(handle)
            return quality_genomes






def read_kegg_map():
    res={}
    level_1_pattern=re.compile('\d\.\s')
    level_2_pattern=re.compile('\d\.\d+\s')
    level_3_pattern=re.compile('\d{3,}:')
    with open(kegg_map) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            search_1=re.search(level_1_pattern,line)
            search_2=re.search(level_2_pattern,line)
            search_3=re.search(level_3_pattern,line)
            if search_1:
                line=line.split()
                level_1=line[0].strip('.')
                res[level_1]={'description':' '.join(line[1:])}
            elif search_2:
                line=line.split()
                level_2=line[0]
                res[level_1][level_2]={'description':' '.join(line[1:])}
            elif search_3:
                line=line.split(':')
                level_3=line[0]
                res[level_1][level_2][level_3]={'description':' '.join(line[1:])}
            line=file.readline()
    return res


def get_description(wanted_id,kegg_map):
    #need to remove Global and overview maps
    for level_1 in kegg_map:
        if level_1!='description':
            level_1_description=kegg_map[level_1]['description']
            if wanted_id==level_1:
                if wanted_level == 1:
                    return level_1_description
                elif not wanted_level:
                    return level_1_description
                else:
                    return
            for level_2 in kegg_map[level_1]:
                if level_2 != 'description':
                    level_2_description = kegg_map[level_1][level_2]['description']
                    if wanted_id == level_2:
                        if wanted_level == 1:
                            return level_1_description
                        elif wanted_level == 2:
                            return level_2_description
                        elif not wanted_level:
                            return level_2_description
                        else:
                            return None
                    for level_3 in kegg_map[level_1][level_2]:
                        if level_3 != 'description':
                            level_3_description = kegg_map[level_1][level_2][level_3]['description']
                            if wanted_id == level_3 and wanted_id in [
                                '01100',
                                '01110',
                                '01120',
                                '01200',
                                '01210',
                                '01212',
                                '01230',
                                '01220',
                            ]:
                                return level_3_description
                            #if wanted_id == level_3 and wanted_id in [
                            #    '00190',
                            #    '00195',
                            #    '00196',
                            #    '00710',
                            #    '00720',
                            #    '00680',
                            #    '00910',
                            #    '00920',
                            #]:
                            #    return level_3_description
                            elif wanted_id == level_3:
                                if wanted_level==1:
                                    return level_1_description
                                elif wanted_level==2:
                                    return level_2_description
                                elif wanted_level==3:
                                    return level_3_description
                                elif not wanted_level:
                                    return level_3_description
                                else:
                                    return


def get_frequency_kegg_map(sample_kegg_maps):
    freq={}
    for km in sample_kegg_maps:
        if km not in freq: freq[km]=0
        freq[km]+=1
    return freq


def get_seqs_count(target_sample):
    total_seqs=0
    with open(target_sample) as file:
        line=file.readline()
        while line:
            if line[0]=='>': total_seqs+=1
            line=file.readline()
    return total_seqs





def get_description_keys():


    kos= [
        'Carbohydrate metabolism',
        'Energy metabolism',
        'Lipid metabolism',
        'Nucleotide metabolism',
        'Amino acid metabolism',
        'Metabolism of cofactors and vitamins',
        #'Metabolism of terpenoids and polyketides',
        'Metabolism of secondary metabolites',
        'Xenobiotics biodegradation and metabolism',
        #'Chemical structure transformation maps',
            #'Transcription',
            #'Translation',
            #'Folding, sorting and degradation',
            #'Replication and repair',
            #'Membrane transport',
            #'Signal transduction',
            #'Signaling molecules and interaction',
            #'Transport and catabolism',
            #'Cell growth and death',
            #'Cellular community',
            #'Cell motility',
        #'Immune system',
        #'Endocrine system',
        #'Circulatory system',
        #'Digestive system',
        #'Excretory system',
        #'Nervous system',
        #'Sensory system',
        #'Development and regeneration',
        #'Aging',
        #'Environmental adaptation',

            #'Cancer',
            #'Immune disease',
            #'Neurodegenerative disease',
            #'Substance dependence',
            #'Cardiovascular disease',
            #'Endocrine and metabolic disease',
            #'Infectious disease',
            #'Drug resistance',

        #'Chronology',
        #'G protein-coupled receptors',
        #'Nuclear receptors',
        #'Ion channels',
        #'Transporters',
        #'Enzymes',
        #'Structure-based classification',
        #'Skeleton-based classification',
    ]
    kos=set(kos)
    c=0
    res={}
    for i in sorted(kos):
        res[i]=c
        c+=1
    return res

def get_description_keys_legacy():
    kos= [
        'Carbohydrate metabolism',
        'Energy metabolism',
        'Lipid metabolism',
        'Nucleotide metabolism',
        'Amino acid metabolism',
        'Metabolism of cofactors and vitamins',
        'Metabolism of terpenoids and polyketides',
        'Biosynthesis of secondary metabolites',
        'Xenobiotics biodegradation and metabolism',
        'Chemical structure transformation maps',
        'Transcription',
        'Translation',
        'Folding, sorting and degradation',
        'Replication and repair',
        'Membrane transport',
        'Signal transduction',
        'Signaling molecules and interaction',
        'Transport and catabolism',
        'Cell growth and death',
        'Cellular community',
        'Cell motility',
        'Immune system',
        'Endocrine system',
        'Circulatory system',
        'Digestive system',
        'Excretory system',
        'Nervous system',
        'Sensory system',
        'Development and regeneration',
        'Aging',
        'Environmental adaptation',
        'Cancer',
        'Immune disease',
        'Neurodegenerative disease',
        'Substance dependence',
        'Cardiovascular disease',
        'Endocrine and metabolic disease',
        'Infectious disease',
        'Drug resistance',
        'Chronology',
        'G protein-coupled receptors',
        'Nuclear receptors',
        'Ion channels',
        'Transporters',
        'Enzymes',
        'Structure-based classification',
        'Skeleton-based classification',
    ]
    kos=set(kos)
    c=0
    res={}
    for i in sorted(kos):
        res[i]=c
        c+=1
    return res



def vectorize_sample(sample,annotation):
    sample_size= get_seqs_count(sample)
    sample_ko_maps=[]
    with open(annotation) as sample_file:
        line=sample_file.readline()
        line=sample_file.readline()
        while line:
            line = line.strip('\n')
            line= line.split('\t')
            ko_maps=set([i.split(':map')[1] for i in line if 'kegg_map:' in i])
            sample_ko_maps.extend(ko_maps)
            line=sample_file.readline()
    freq = get_frequency_kegg_map(sample_ko_maps)
    freq_normalized={i:round(freq[i]/sample_size*100,3) for i in freq}
    #freq_normalized={i:round(freq_names[i]/sum(freq_names.values())*100,3) for i in freq_names}
    return freq_normalized

def convert_to_freq_names(freq):
    kegg_map = read_kegg_map()
    freq_names={}
    for km in freq:
        description=get_description(km,kegg_map)
        if description not in freq_names: freq_names[description]=0
        freq_names[description]+=freq[km]
    if None in freq_names:
        freq_names.pop(None)
    return freq_names

def get_map_names(vector):
    res={}
    description_keys = get_description_keys()
    description_keys = {description_keys[i]:i for i in description_keys}
    for dk in description_keys:
        if dk in vector:  res[description_keys[dk]]=vector[dk]
        else:  res[description_keys[dk]]=0
    return res

def get_map_ints(vector):
    res={}
    description_keys = get_description_keys()
    for dk in description_keys:
        if dk in vector:  res[description_keys[dk]]=vector[dk]
        else:  res[description_keys[dk]]=0
    return res



def get_vectors(samples):
    vectors_path='vectors_'+str(wanted_level)+'.pickle'
    res=load_metrics(vectors_path)
    if res:
        print('loaded vectors_path')
        return res
    print('reading vectors_path')
    vectors={}
    for s in samples:
        s_vector=vectorize_sample(s[0],s[1])
        converted_to_freq_names= convert_to_freq_names(s_vector)
        #here we add all the descriptions to the vector
        vectors[s[2]]=get_map_ints(converted_to_freq_names)
    save_metrics(vectors_path,vectors)
    return vectors

def generate_samples(proteomes_tsv):
    res=[]
    samples_path=''
    with open(proteomes_tsv) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            sample_name,ncbi_id=line[0],line[2]
            res.append([samples_path+sample_name+'/prodigal_proteins.faa',samples_path+sample_name+'/consensus_annotation.tsv',ncbi_id])
            line=file.readline()
    return res




def read_target_taxa_tsv():
    res={}
    with open(path+'PhD/Presentations/cet_2/ML/lao_taxa.tsv') as file:
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            sample,taxa_name,taxa_id=line
            res[sample]=[taxa_id,taxa_name]
            line = file.readline()
    return res

def read_proteomes_taxa_tsv():
    res={}
    with open(path+'PhD/Presentations/cet_2/ML/proteomes.tab') as file:
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            sample,taxa_name,taxa_id=line[0:3]
            taxa_name=taxa_name.split()[0:2]
            taxa_name=[i.strip('[]') for i in taxa_name if 'sp.' not in i and 'bacterium' not in i]
            taxa_name=' '.join(taxa_name)
            res[sample]=[taxa_id,taxa_name]
            line = file.readline()
    return res

def read_mantis_tsv(mantis_tsv_path):
    res={}
    with open(mantis_tsv_path) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            res[line[0]]=line[-1]
            line=file.readline()
    return res

def get_samples():
    samples_path='samples_'+str(wanted_level)+'.pickle'
    res=load_metrics(samples_path)
    if res:
        print('loaded samples')
        return res
    print('reading samples')
    res=[]
    for d in [
        drax_path+'Data/Annotation_Output/LAO_Isolates/',
        drax_path+'Data/Annotation_Output/Microthrix_parvicella/',
        drax_path+'Data/Annotation_Output/reference_proteomes/annotations/',
        drax_path+'Data/Annotation_Output/ben_must_2020-10-26T165726/',
              ]:
        list_dir = os.listdir(d)
        for inner_d in list_dir:
            if os.path.isdir(d+inner_d):
                if 'reference_proteomes' in d:
                    sample=d.replace('annotations','faa_folder')+inner_d+'.faa'
                    annot_path=d+inner_d+'/consensus_annotation.tsv'
                elif 'ben_must' in d:
                    sample=[d+inner_d+'/'+i for i in os.listdir(d+inner_d) if '.faa' in i][0]
                    annot_path=d+inner_d+'/consensus_annotation.tsv'
                else:
                    sample=d+inner_d+'/prodigal_proteins.faa'
                    annot_path=d+inner_d+'/output_mantis/consensus_annotation.tsv'
                    inner_d=inner_d.split('.')[0]
                res.append([sample,annot_path,inner_d])
    save_metrics(samples_path,res)
    return res





def get_organism_lineage(taxon_id, stdout_file=None):
    lineage_file_path = '/home/pedroq/Python_projects/DRAX/source/Pipelines/mantis/Resources/NCBI/' + 'taxidlineage.dmp'
    try:
        lineage_file = open(lineage_file_path, 'r')
    except:
        print_cyan(
            'Lineage dump is not present! If you\'d like to run taxonomic lineage annotation, please run < setup_databases >',
            flush=True, file=stdout_file)
        return []
    line = lineage_file.readline().strip('\n').replace('|', '')
    while line:
        line = line.split()
        if str(taxon_id) == str(line[0]):
            lineage_file.close()
            lineage=line[1:]
            lineage.append(taxon_id)
            return lineage
        line = lineage_file.readline().strip('\n').replace('|', '')
    lineage_file.close()
    return []

def get_taxa_ncbi(organism_name):
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + organism_name
    webpage = None
    c = 0
    while not webpage and c <= 10:
        req = requests.get(url)
        try:
            webpage = req.text
        except:
            c += 1
    taxa_id = re.search('<Id>\d+</Id>', webpage)
    if taxa_id: return re.search('\d+', taxa_id.group()).group()

def get_samples_taxa_ids():
    with open(path+'PhD/Presentations/cet_2/ML/taxa.tsv') as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            s,i=line[0],line[1]
            res=''
            if i:
                taxa=get_taxa_ncbi(i)
                if taxa:
                    res=s+'\t'+i+'\t'+taxa
                else:
                    res=s+'\t'+i+'\t'+''
            else:
                res = s + '\t' + '' + '\t' + ''
            print(res)
            line=file.readline()


def get_dataframe(samples):
    vectors = get_vectors(samples)
    description_keys=get_description_keys()
    description_keys={description_keys[i]:i for i in sorted(description_keys.keys())}
    data=[]
    keys=[]
    for v in vectors:
        temp=[vectors[v][i] for i in sorted(vectors[v].keys())]
        data.append(temp)
        keys.append(v)
    description_keys=get_description_keys()
    description_keys={description_keys[i]:i for i in description_keys}
    description_keys=[description_keys[i] for i in sorted(description_keys)]
    df=pd.DataFrame(np.array(data),index=keys,columns=description_keys)
    df=clean_dataframe(df)
    return df





def plot_distribution(df,ko_description,show_taxa=True):
    ref=[i for i in df.index if i in ref_taxa]
    ref_vals=df.loc[ref,ko_description]
    ko_description_list = df.loc[:,ko_description]
    target=[i for i in target_taxa]
    target_vals=df.loc[target,ko_description]
    mean=statistics.mean(ref_vals)
    mean_target=statistics.mean(target_vals)
    stdev=statistics.stdev(ref_vals)
    median=statistics.median(ref_vals)
    perc_90=np.percentile(ref_vals,90)
    perc_75=np.percentile(ref_vals,75)
    #print('mean',mean)
    #print('stdev',stdev)
    #print('median',median)
    #print('percentile 90',np.percentile(ko_description_list,90))

    #changing to taxa name
    if show_taxa:
        target=[target_taxa[i][1] for i in target_taxa]
        ref=[ref_taxa[i][1] for i in ref_taxa]

    above_mean=list(df[ko_description_list>=mean].loc[:,ko_description].index)
    above_median=list(df[ko_description_list>=median].loc[:,ko_description].index)
    below_mean=list(df[ko_description_list<mean].loc[:,ko_description].index)
    below_median=list(df[ko_description_list<median].loc[:,ko_description].index)
    #print('Samples above mean')
    #print(above_mean)
    #print([target_taxa[i][1] for i in above_mean if i in target_taxa])
    #print('Samples below mean')
    #print(below_mean)
    #print([target_taxa[i][1] for i in below_mean if i in target_taxa])
    # Group data together
    fig = go.Figure()
    fig.update_layout(title=ko_description)
    fig.add_trace(go.Scatter(
        x=target, y=target_vals,
        name='Target',
        mode='markers',
        marker_color='rgba(254, 180, 2, 1)',
    ))
    fig.add_trace(go.Scatter(x=ref, y=ref_vals,name='Reference',mode='markers',marker_color='rgba(0, 180, 0, 0.3)'))
    x_axis=['Neisseria','Alloscardovia macacae']
    fig.add_trace(go.Scatter(x=x_axis, y=[mean,mean],mode='lines',name='mean all'))
    fig.add_trace(go.Scatter(x=x_axis, y=[median,median],mode='lines',name='median all'))
    fig.add_trace(go.Scatter(x=x_axis,y=[perc_90,perc_90],mode='lines',name='perc_90'))
    fig.add_trace(go.Scatter(x=x_axis,y=[perc_75,perc_75],mode='lines',name='perc_75'))
    fig.add_trace(go.Scatter(x=x_axis,y=[mean_target,mean_target],mode='lines',name='mean_target'))
    fig.show()

def clean_dataframe(df):
    cleaned=[]
    #removing empty collumns
    c=0
    removed=[]
    var_threshold=0.1
    for ko in df:
        if not (df[ko] == 0).all():
            ko_mean=df[ko].mean()
            ko_stdev=df[ko].std()
            ko_var=df[ko].var()
            if ko_var<var_threshold:

                print('deleting 1',ko,ko_var)
                #del df[ko]
            else:
                standardized_ko=(df[ko]-ko_mean)/ko_stdev
                df[ko]=standardized_ko
        else:
            del df[ko]
    return df


def optimalK(data, clustering_function,nrefs=3, maxClusters=15):
    """
    Calculates KMeans optimal K using Gap Statistic from Tibshirani, Walther, Hastie
    Params:
        data: ndarry of shape (n_samples, n_features)
        nrefs: number of sample reference datasets to create
        maxClusters: Maximum number of clusters to test for
    Returns: (gaps, optimalK)
    """
    gaps = np.zeros((len(range(1, maxClusters)),))
    resultsdf = pd.DataFrame({'clusterCount': [], 'gap': []})
    for gap_index, k in enumerate(range(1, maxClusters)):
        refDisps = np.zeros(nrefs)
        # For n references, generate random sample and perform kmeans getting resulting dispersion of each loop
        for i in range(nrefs):
            # Create new random reference set
            randomReference = np.random.random_sample(size=data.shape)
            # Fit to it
            km = KMeans(k)
            km.fit(randomReference)
            refDisp = km.inertia_
            refDisps[i] = refDisp
        # Fit cluster to original data and create dispersion
        km = KMeans(k)
        km.fit(data)
        origDisp = km.inertia_
        # Calculate gap statistic
        gap = np.log(np.mean(refDisps)) - np.log(origDisp)
        # Assign this loop's gap statistic to gaps
        gaps[gap_index] = gap
        resultsdf = resultsdf.append({'clusterCount': k, 'gap': gap}, ignore_index=True)
    return  gaps.argmax() + 1

def get_optimal_cluster_number(df,clustering_function,random_state=None):
    silhouette_coefficients=[]
    n_kos=40
    algorithm=str(clustering_function).split('.')[-1].strip('\'>').lower()
    print('clustering method',algorithm)
    #print('Optimal cluster number gap static',optimalK(df,clustering_function,maxClusters=n_kos))
    if 'dbscan' in algorithm:
        clustering=clustering_function(eps=0.5, min_samples=5, metric='euclidean',
                                       metric_params=None, algorithm='auto', leaf_size=30, p=None,
               n_jobs=None).fit(df)
        labels=clustering.labels_
        if len(set(labels))>1:
            score = silhouette_score(df, labels)
            print("Silhouette Coefficient:", score)
            score = davies_bouldin_score(df, labels)
            print("davies bouldin score:", score)
        return
    elif \
        'affinity' in algorithm or \
        'meanshift' in algorithm:
        if random_state is not None:
            clustering=clustering_function(random_state=random_state).fit(df)
        else:
            clustering=clustering_function().fit(df)
        labels=clustering.labels_
        score = silhouette_score(df, labels)
        print("Silhouette Coefficient:", score)
        score = davies_bouldin_score(df, labels)
        print("davies bouldin score:", score)
        return
    sse = []
    silhouette_coefficients=[]
    davies_coefficients=[]

    if 'optics' in algorithm:
        for k in [0.5,1,2,3,4,5]:
            clustering = clustering_function(min_samples=k).fit(df)
            labels=clustering.labels_
            if len(set(labels)) > 1:
                score = silhouette_score(df, labels)
                silhouette_coefficients.append(score)
                score = davies_bouldin_score(df, labels)
                davies_coefficients.append(score)

    else:
        for k in range(2, n_kos):
            if 'gaussian' in algorithm:
                clustering=clustering_function(n_components=k).fit(df)
                labels=clustering.fit_predict(df)
            elif random_state is not None:
                clustering=clustering_function(n_clusters=k,random_state=random_state).fit(df)
                labels=clustering.labels_
            elif 'kmeans' in algorithm:
                clustering=clustering_function(n_clusters=k).fit(df)
                labels=clustering.labels_
                sse.append(clustering.inertia_)
            else:
                clustering=clustering_function(n_clusters=k).fit(df)
                labels=clustering.labels_
            score = silhouette_score(df, labels)
            print(k,score)
            silhouette_coefficients.append(score)
            score = davies_bouldin_score(df, labels)
            davies_coefficients.append(score)
        if 'kmeans' in algorithm and sse:
            kl = KneeLocator(range(2, n_kos), sse, curve = "convex", direction = "decreasing")
            print('optimal cluster number elbow',kl.elbow+2)
    print('optimal cluster number silhouette ',silhouette_coefficients.index(max(silhouette_coefficients))+2)
    print("Silhouette Coefficient:",max(silhouette_coefficients))
    print('optimal cluster number Davies Bouldin ',davies_coefficients.index(min(davies_coefficients))+2)
    print("Davies Bouldin Coefficient:",max(davies_coefficients))
    return None

def build_cluster(training_data,save_output=True,clustering_method='kmeans',test=False,n_cluster=None):
    '''
    silhouette and elbow gives the optimal number of clusters
    silhouette is -1 to 1, 1 is best
    empirical n of clusters should be sqrt(N/2) where N are the data points/samples
    davies_bouldin_score, 0 is best

    silhouette scores:
    agglomerative   0.166
    kmeans   0.193 , best 7 with silhouete, 15 with elbow
    affinity   0.125
    meanshit   0.323
    optics   -0.228
    gaussian mixture   0.176
    birch   0.174
    '''
    '''
        clustering method agglomerativeclustering
        optimal cluster number silhouette  2
        Silhouette Coefficient: 0.16569413785005205
        optimal cluster number Davies Bouldin  2
        Davies Bouldin Coefficient: 2.031553315832657
        clustering method kmeans
        optimal cluster number elbow 15
        optimal cluster number silhouette  7
        Silhouette Coefficient: 0.19300647036221974
        optimal cluster number Davies Bouldin  2
        Davies Bouldin Coefficient: 1.8547841833187495
        clustering method affinitypropagation
        Silhouette Coefficient: 0.12540753532128313
        davies bouldin score: 1.4186826638817278
        clustering method meanshift
        Silhouette Coefficient: 0.3227941869941836
        davies bouldin score: 0.8096767866668201
        clustering method optics
        optimal cluster number silhouette  2
        Silhouette Coefficient: -0.22790020766370345
        optimal cluster number Davies Bouldin  4
        Davies Bouldin Coefficient: 1.4867214015164418
        clustering method gaussianmixture
        optimal cluster number silhouette  2
        Silhouette Coefficient: 0.1761129009897275
        optimal cluster number Davies Bouldin  3
        Davies Bouldin Coefficient: 2.8107309946972006
        clustering method birch
        optimal cluster number silhouette  2
        Silhouette Coefficient: 0.17492752443593101
        optimal cluster number Davies Bouldin  2
        Davies Bouldin Coefficient: 2.0535503622172238
    
    '''

    cluster_path='cluster_'+str(wanted_level)+'.pickle'
    res=load_metrics(cluster_path)
    if res:
        print('loading cluster')
        return res
    else: print ('building cluster')

    clustering=None
    max_score=0
    if clustering_method=='agglomerative':
            clustering_function=skcluster.AgglomerativeClustering
            if not n_cluster:
                n_cluster=get_optimal_cluster_number(training_data,clustering_function=clustering_function)
            if not test:
                clustering = clustering_function(n_clusters=n_cluster).fit(training_data)
    elif clustering_method=='kmeans':
        for i in range(1000):
            clustering_function=skcluster.KMeans
            if not n_cluster:
                n_cluster=get_optimal_cluster_number(training_data,clustering_function=clustering_function)
            if not test:
                clustering=clustering_function(n_clusters=n_cluster,random_state=None).fit(training_data)
                cluster_score = silhouette_score(training_data, clustering.labels_)
                print(i,cluster_score)
                if cluster_score > max_score:
                    best_cluster = clustering
                    max_score=cluster_score
        if not test:
            clustering=best_cluster
    elif clustering_method=='affinity':
        clustering_function=skcluster.AffinityPropagation
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function,random_state=0)
        if not test:
            clustering = clustering_function(random_state=None).fit(training_data)
    elif clustering_method=='meanshift':
        clustering_function=skcluster.MeanShift
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function,)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method=='dbscan':
        clustering_function=skcluster.DBSCAN
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method=='optics':
        clustering_function=skcluster.OPTICS
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method=='gaussian':
        clustering_function=skmixture.GaussianMixture
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method=='bayesian_gaussian':
        clustering_function=skmixture.BayesianGaussianMixture
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method=='birch':
        clustering_function=skcluster.Birch
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    save_metrics(cluster_path,best_cluster)
    return clustering

def show_clusters(labels,df,show_taxa=False):
    if len(labels)!=df.shape[0]:
        print('different dimensions!')
        raise Exception
    fig_labels=[]
    cluster_by_taxa={}
    keys=df.index
    for l in range(len(labels)):
        if labels[l] not in cluster_by_taxa: cluster_by_taxa[labels[l]]=[]
        cluster_by_taxa[labels[l]].append(keys[l])

    for i in cluster_by_taxa:
        if show_taxa:
            print(i,[ref_taxa[j][1] if j in ref_taxa else target_taxa[j][1] for j in cluster_by_taxa[i]])
        else:
            print(i,[j for j in cluster_by_taxa[i]])

def get_cluster_coordinates(training_data,training_labels,test_labels,test_data,samples_to_plot=[],taxa_to_plot=[]):
    if len(training_labels)!=training_data.shape[0] or len(test_labels)!=test_data.shape[0]:
        print(len(training_labels),training_data.shape)
        print(len(test_labels),test_data.shape)
        print('different dimensions!')
        raise Exception
    res={}
    for ko_description in get_description_keys():
        try:
            ko_description_list_training = training_data.loc[:,ko_description]
            ko_description_list_test = test_data.loc[:,ko_description]
            cluster_values = {'target':{},'ref':{}}
            training_keys = training_data.index
            test_keys = test_data.index
            for l in range(len(test_labels)):
                if test_labels[l] not in cluster_values['target']: cluster_values['target'][test_labels[l]] = []
                cluster_values['target'][test_labels[l]].append(ko_description_list_test[l])
            for l in range(len(training_labels)):
                if training_labels[l] not in cluster_values['ref']: cluster_values['ref'][training_labels[l]] = []
                cluster_values['ref'][training_labels[l]].append(ko_description_list_training[l])
            cluster_median={'target':{},'ref':{}}
            for data_type in cluster_values:
                for i in cluster_values[data_type]:
                    if i not in cluster_median[data_type]: cluster_median[data_type][i]={'median':None,'stdev':None}
                    cluster_median[data_type][i]['median']=statistics.median(cluster_values[data_type][i])
                    try:
                        cluster_median[data_type][i]['stdev']=statistics.stdev(cluster_values[data_type][i])
                    except:
                        cluster_median[data_type][i]['stdev']=0
            res[ko_description]=cluster_median
        except: pass
    sample_coordinates=[]
    already_added=[]
    for sample in samples_to_plot:
        if sample not in already_added:
            if sample in training_data.index: sample_df= training_data.loc[sample,:]
            else: sample_df= test_data.loc[sample,:]
            sample_coordinates.append(sample_df)
            already_added.append(sample)
    taxa_to_sample_ref={}
    taxa_to_sample_target={}
    for s in ref_taxa:
        taxa = ref_taxa[s][1]
        if taxa:
            domain=taxa.split(' ')[0]
            if taxa not in taxa_to_sample_ref: taxa_to_sample_ref[taxa]=[]
            taxa_to_sample_ref[taxa].append(s)
            if domain not in taxa_to_sample_ref: taxa_to_sample_ref[domain]=[]
            taxa_to_sample_ref[domain].append(s)
    for s in target_taxa:
        taxa=target_taxa[s][1]
        if taxa:
            domain=taxa.split(' ')[0]
            if taxa not in taxa_to_sample_target: taxa_to_sample_target[taxa]=[]
            taxa_to_sample_target[taxa].append(s)
            if domain not in taxa_to_sample_target: taxa_to_sample_target[domain]=[]
            taxa_to_sample_target[domain].append(s)
    for taxa in taxa_to_plot:
        samples_taxa=taxa_to_sample_ref[taxa]
        samples_taxa.extend(taxa_to_sample_target[taxa])
        for sample in samples_taxa:
            if sample not in already_added:
                if sample in training_data.index: sample_df= training_data.loc[sample,:]
                elif sample in test_data.index: sample_df= test_data.loc[sample,:]
                else:
                    sample_df=None
                if sample_df is not None:
                    sample_coordinates.append(sample_df)
                    already_added.append(sample)
    return res,sample_coordinates

def plot_clusters(cluster_coordinates,show_error_bar=False,restrict_to_target=False,sample_coordinates=[]):
    fig = go.Figure()
    x=[]
    y={'target':{},'ref':{}}
    error_y={'target':{},'ref':{}}
    for ko_description in sorted(cluster_coordinates.keys()):
        x.append(ko_description)
        for data_type in cluster_coordinates[ko_description]:
            for cluster in cluster_coordinates[ko_description][data_type]:
                if cluster not in y[data_type]: y[data_type][cluster]=[]
                if cluster not in error_y[data_type]: error_y[data_type][cluster]=[]
                y[data_type][cluster].append(cluster_coordinates[ko_description][data_type][cluster]['median'])
                error_y[data_type][cluster].append(cluster_coordinates[ko_description][data_type][cluster]['stdev'])
    colors={0:'#000000',1:'#D20007',2:'#FFAA00',3:'#FFFB00',4:'#00FF1E',5:'#00FFFB',6:'#0800FF',7:'#D900FF'}
    colors=[colors[i] for i in colors]
    size={'ref':10,'target':15}
    for data_type in ['target','ref']:
        for cluster in y[data_type]:
            print(cluster,colors[cluster])
            if restrict_to_target:
                if cluster in y['target']:
                    if show_error_bar:
                        fig.add_trace(go.Scatter(x=x, y=y[data_type][cluster],
                                                 error_y=dict(type='data', array=error_y[data_type][cluster],
                                                              visible=True), mode='markers',
                                                 marker=dict(size=size[data_type], color=colors[cluster], ),
                                                 name=data_type + '_' + str(cluster)))
                    else:
                        fig.add_trace(go.Scatter(x=x, y=y[data_type][cluster], mode='markers',
                                                 marker=dict(size=size[data_type], color=colors[cluster], ),
                                                 name=data_type + '_' + str(cluster)))
            else:
                if show_error_bar:
                    fig.add_trace(go.Scatter(x=x, y=y[data_type][cluster],error_y=dict(type='data',array=error_y[data_type][cluster],visible=True),mode='markers',marker=dict(size=size[data_type],color=colors[cluster],),name=data_type+'_'+str(cluster)))
                else:
                    fig.add_trace(go.Scatter(x=x, y=y[data_type][cluster],mode='markers',marker=dict(size=size[data_type],color=colors[cluster],),name=data_type+'_'+str(cluster)))
    for sample in sample_coordinates:
        fig.add_trace(go.Scatter(x=sample.index,y=sample.values, mode='lines+markers',
                                 #marker=dict(size=size[data_type], colorscale=colors[data_type], ),
                                 name=sample.name))

    fig.show()


def plot_all_ko_distributions(training_data, training_labels, target_data, target_labels, cluster_to_ko,quantiles,add_ref=False,show_taxa=False,show_mix=False):
    for ko in get_description_keys():
        try:
            plot_distribution_clustered(training_data=training_data, training_labels=training_labels, target_data=target_data,
                                        target_labels=target_labels, ko_description=ko,cluster_to_ko=cluster_to_ko,quantiles=quantiles,add_ref=add_ref,show_taxa=show_taxa,show_mix=show_mix)
        except: pass



def check_clusters_quality(cluster,training_data):
    Y=cluster.labels_
    counter={}
    cluster_indexes={}
    for c in range(len(Y)):
        cluster_n=Y[c]
        if cluster_n not in cluster_indexes: cluster_indexes[cluster_n]=[]
        cluster_indexes[cluster_n].append(c)
    total_variance={}
    poor_clusters=set()
    for ko_description in get_description_keys():
        try:
            for c in cluster_indexes:
                df_indexes = [training_data.index[i] for i in cluster_indexes[c]]
                cluster_ko=training_data.loc[df_indexes,ko_description]
                cluster_ko_variance=cluster_ko.std()
                if cluster_ko_variance>1.5:
                    poor_clusters.add(c)
                    print(ko_description,c,cluster_ko_variance)
                if c not in total_variance: total_variance[c]=0
                total_variance[c]+=cluster_ko_variance
        except: pass
    for c in Y:
        if c not in counter: counter[c]=0
        counter[c]+=1
    print('cluster size',counter)
    print('poor clusters',poor_clusters)




def test_classifiers(training_labels,training_data):
    for i in range(5):
        X_train, X_test, y_train, y_test = train_test_split(training_data, training_labels, test_size = 0.3)
        n_features=training_data.shape[1]
        max_features=int(sqrt(n_features))+1
        clf = svm.SVC(C=0.5,random_state=0).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        clf = svm.SVC(C=0.75,random_state=0).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        clf = svm.SVC(random_state=0).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        #better for noisy data / noisy clusters
        clf = ExtraTreesClassifier(max_depth=None,min_samples_split = 2,max_features=max_features).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        clf = DecisionTreeClassifier(max_depth=None,min_samples_split = 2,max_features=max_features).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        clf = RandomForestClassifier(max_depth=None,min_samples_split = 2, random_state = 0,max_features=max_features).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        clf = GradientBoostingClassifier(learning_rate=1.0,max_depth = 1, random_state = 0).fit(X_train, y_train)
        score=clf.score(X_test, y_test)
        print(str(clf),score)
        #clf = MLPClassifier(solver='lbfgs', alpha=1e-5,hidden_layer_sizes = (5, 2), random_state = 0).fit(training_data, training_labels)
        #score=clf.score(X_test, y_test)
        #print(str(clf),score)


def build_classifier(training_labels,training_data):
    classifier_path='classifier_'+str(wanted_level)+'.pickle'
    res=load_metrics(classifier_path)
    if res:
        print('loading classifier')
        return res
    else: print ('building classifier')
    n_features = training_data.shape[1]
    max_features = int(sqrt(n_features)) + 1
    #clf = DecisionTreeClassifier(max_depth=None, min_samples_split=2).fit(training_data, training_labels)
    clf = ExtraTreesClassifier().fit(training_data, training_labels)
    #clf = ExtraTreesClassifier(max_depth=None, min_samples_split=2, random_state=0,max_features=max_features).fit(training_data, training_labels)
    #save_metrics(classifier_path,clf)
    return clf

def get_taxa_clusters(labels, df,threshold=1):
    if len(labels) != df.shape[0]:
        print('different dimensions!')
        raise Exception
    fig_labels = []
    taxa_by_cluster = {}
    keys = df.index
    for l in range(len(labels)):
        taxa=target_taxa[keys[l]][1] if keys[l] in target_taxa else ref_taxa[keys[l]][1]
        if taxa:
            if taxa not in taxa_by_cluster: taxa_by_cluster[taxa] = []
            taxa_by_cluster[taxa].append(labels[l])
    for taxa in taxa_by_cluster:
        taxa_by_cluster[taxa]=set(taxa_by_cluster[taxa])
    c=0
    for taxa in taxa_by_cluster:
        if len(taxa_by_cluster[taxa])>threshold:
            print(taxa, taxa_by_cluster[taxa])
            c+=1
    print('taxa incorrectly clustered',c)
    print('total taxa',len(taxa_by_cluster))
    print('taxa incorrectly clustered %',round(c/len(labels),5)*100)


def get_labels_df(df):
    target = [i for i in df.index if i in target_taxa]
    target_data = df.loc[target, :]
    training_data, training_labels, cluster_to_ko,ko_to_cluster,quantiles = get_training_labels(df)
    print(training_data.shape,len(training_labels))
    cluster_score = silhouette_score(training_data, training_labels)
    print('cluster score',cluster_score)
    fitted_classifier = build_classifier(training_labels, training_data)
    target_labels = fitted_classifier.predict(target_data)
    return training_data,training_labels,target_data,target_labels,cluster_to_ko,ko_to_cluster,quantiles

def plot_distribution_clustered(training_data,training_labels,target_data,target_labels,ko_description,cluster_to_ko,quantiles,add_ref=False,show_taxa=False,show_mix=False):
    colors={'Carbohydrate metabolism':'#000000',
            'Energy metabolism':'#D20007',
            'Lipid metabolism':'#FFAA00',
            'Nucleotide metabolism':'#793500',
            'Amino acid metabolism':'#008810',
            'Metabolism of cofactors and vitamins':'#7E7E7E',
            'Metabolism of secondary metabolites':'#0600BF',
            'Xenobiotics biodegradation and metabolism':'#85009C'}
    ko_to_cluster={i:cluster_to_ko[i] for i in cluster_to_ko}
    cluster_samples={}
    for i in range(len(target_labels)):

        if target_labels[i] not in cluster_samples: cluster_samples[target_labels[i]]=[]
        sample=target_data.index[i]
        cluster_samples[target_labels[i]].append(sample)
    ref_samples={}
    for i in range(len(training_labels)):
        if training_labels[i] not in ref_samples: ref_samples[training_labels[i]]=[]
        sample=training_data.index[i]
        ref_samples[training_labels[i]].append(sample)


    # Group data together
    fig = go.Figure()
    fig.update_layout(title=ko_description)#,paper_bgcolor='rgba(0,0,0,0)',plot_bgcolor='rgba(0,0,0,0)')
    x_axis=[]
    for cluster in cluster_samples:
        target_vals = list(target_data.loc[cluster_samples[cluster], ko_description].values)
        if show_taxa:    x_labels = [target_taxa[sample][1] for sample in cluster_samples[cluster]]
        else:            x_labels = cluster_samples[cluster]
        if ko_to_cluster[cluster] in colors:
            color=colors[ko_to_cluster[cluster]]
            x_axis.extend(x_labels)
            fig.add_trace(go.Scatter(x=x_labels, y=target_vals,name='Target enriched in '+cluster_to_ko[cluster],mode='markers',marker_color=color,marker={'size':10}))
        else:
            if show_mix:
                x_axis.extend(x_labels)
                fig.add_trace(go.Scatter(x=x_labels, y=target_vals,name='Target enriched in '+cluster_to_ko[cluster],mode='markers',marker={'size':7}))
    for cluster in ref_samples:
        target_vals = list(training_data.loc[ref_samples[cluster], ko_description].values)
        if show_taxa:    x_labels = [ref_taxa[sample][1] for sample in ref_samples[cluster]]
        else:            x_labels = ref_samples[cluster]
        if ko_to_cluster[cluster] in colors:
            print(cluster,len(x_labels))
            color=colors[ko_to_cluster[cluster]]
            if add_ref:
                x_axis.extend(x_labels)
                fig.add_trace(go.Scatter(x=x_labels, y=target_vals,name='Ref enriched in '+cluster_to_ko[cluster],mode='markers',marker_color=color,marker={'size':10}))
        else:
            if show_mix and add_ref:
                x_axis.extend(x_labels)
                fig.add_trace(go.Scatter(x=x_labels, y=target_vals,name='Ref enriched in '+cluster_to_ko[cluster],mode='markers',marker={'size':7}))
    x_axis=[x_axis[0],x_axis[-1]]
    fig.add_trace(go.Scatter(x=x_axis,y=[quantiles[ko_description],quantiles[ko_description]],mode='lines',line={'dash': 'dash', 'color': '#46B1C9'},name='Quantile 95'))
    fig.show()


def add_bens_taxa(target_taxa):
    bens_taxa = read_mantis_tsv(
        'C:/Users/Pedro Queirós/Documents/Python Projects/DRAX/Data/Annotation_Output/ben_must_2020-10-26T165726/ben_must.tsv')
    for bt in bens_taxa:
        if bens_taxa[bt] == '1304':
            target_taxa[bt] = [bens_taxa[bt], 'Streptococcus salivarius']
        elif bens_taxa[bt] == '1309':
            target_taxa[bt] = [bens_taxa[bt], 'Streptococcus mutans']
        elif bens_taxa[bt] == '482':
            target_taxa[bt] = [bens_taxa[bt], 'Neisseria']
    return target_taxa




def get_high_ko(df,ko_description):
    res=[]
    ref = [i for i in df.index if i in ref_taxa]
    ref_vals=df.loc[ref,ko_description]
    quantile=ref_vals.quantile(0.95)
    aboveq=df[df.loc[:,ko_description]>quantile].loc[:,ko_description]
    for s in aboveq.index:
        res.append(s)
    return res,quantile

def get_training_labels(df):
    res={}
    ref = [i for i in df.index if i in ref_taxa]
    training_df=df.loc[ref,:]
    quantiles={}
    for ko in get_description_keys():
        try:
            high_kos,quantile=get_high_ko(training_df,ko_description=ko)
            res[ko]=high_kos
            quantiles[ko]=quantile
        except: pass
    res_per_sample={}
    for ko in res:
        for s in res[ko]:
            if s not in res_per_sample: res_per_sample[s]=[]
            res_per_sample[s].append(ko)
    cluster_types={}
    merged_res_per_sample={}
    for s in res_per_sample:
        #if len(res_per_sample[s])==1:
        if True:
            cluster_type='_'.join(sorted(res_per_sample[s]))
            if cluster_type not in cluster_types: cluster_types[cluster_type]=[]
            cluster_types[cluster_type].append(s)
            merged_res_per_sample[s]=cluster_type
    all_samples=[]
    cluster_to_ko={}
    ko_to_cluster={}
    c=0
    for i in cluster_types:
        cluster_to_ko[c]=i
        ko_to_cluster[i]=c
        all_samples.extend(cluster_types[i])
        c+=1
    cluster_df=df.loc[all_samples,:]
    training_labels=[]
    for s in cluster_df.index:
        #training_labels.append(ko_to_cluster[res_per_sample[s][0]])
        training_labels.append(ko_to_cluster[merged_res_per_sample[s]])
    sample_int=int(list(cluster_df.index).index('UP000004699'))
    #removing outliers
    samples_to_remove=set()
    for cluster_type in cluster_to_ko:
        samples_to_remove.update(detect_outliers(cluster_df,training_labels,cluster_to_ko,cluster_type))
    indexes_to_remove=[]
    print('samples_to_remove',len(samples_to_remove),len(training_labels))
    for i in range(len(cluster_df.index)):
        if cluster_df.index[i] in samples_to_remove:
            indexes_to_remove.append(i)
    training_labels=np.delete(training_labels, indexes_to_remove)
    for sample in samples_to_remove:
        cluster_df=cluster_df.drop(sample)
    return cluster_df,training_labels,cluster_to_ko,ko_to_cluster,quantiles

def get_outliers_zscore(df):
    z_scores = stats.zscore(df)
    abs_z_scores = np.abs(z_scores)
    filtered_entries = (abs_z_scores < 3)
    samples_to_remove= [df.index[i] for i in range(len(filtered_entries)) if not filtered_entries[i]]
    return samples_to_remove

def get_outliers_modified_zscore(df):
    #placeholder
    med_col = df.median()
    med_abs_dev = (np.abs(df - med_col)).median()
    mod_z = 0.6745  * ((df - med_col) / med_abs_dev)
    filtered_entries = mod_z < 3.5
    samples_to_remove= [df.index[i] for i in range(len(filtered_entries)) if not filtered_entries[i]]
    return samples_to_remove


def get_outliers_iqr(df):
    Q1 = df.quantile(0.25)
    Q3 = df.quantile(0.75)
    IQR = Q3 - Q1
    filtered_entries = ~(df < (Q1 - 1.5 * IQR)) | (df > (Q3 + 1.5 * IQR))
    samples_to_remove= [df.index[i] for i in range(len(filtered_entries)) if not filtered_entries[i]]
    return samples_to_remove

def detect_outliers(training_data,training_labels,cluster_to_ko,cluster):
    ko_samples=[]
    for i in range(len(training_labels)):
        sample=training_data.index[i]
        if training_labels[i]==cluster:
            ko_samples.append(sample)
    ko_samples_df=training_data.loc[ko_samples,:]
    samples_to_remove=set()
    for col in ko_samples_df.columns:
        current_df=ko_samples_df.loc[:,col]
        #temp = get_outliers_zscore(current_df)
        #temp = get_outliers_modified_zscore(current_df)
        temp = get_outliers_iqr(current_df)
        samples_to_remove.update(temp)
    return samples_to_remove

def most_frequent(List):
    return max(set(List), key = List.count)

def get_labels_df_consensus(df):
    samples_labels={}
    for i in range(30):
        training_data, training_labels, target_data, target_labels, cluster_to_ko,ko_to_cluster, quantiles = get_labels_df(df)
        target_labels_ko=[cluster_to_ko[i] for i in target_labels]
        for i in range(len(target_labels)):
            sample=target_data.index[i]
            if sample not in samples_labels: samples_labels[sample]=[]
            samples_labels[sample].append(target_labels_ko[i])
    unstable_samples=[]
    for sample in samples_labels:
        if len(set(samples_labels[sample]))>1:
            unstable_samples.append(sample)
        samples_labels[sample]=most_frequent(samples_labels[sample])
    res=[]
    for i in range(len(target_labels)):
        sample = target_data.index[i]
        res.append(ko_to_cluster[samples_labels[sample]])
    return training_data, training_labels, target_data, res, cluster_to_ko, quantiles

def load_ml_results(df):
    res= load_metrics(ml_results)
    if res:
        return res
    training_data, training_labels, target_data, target_labels, cluster_to_ko, quantiles = get_labels_df_consensus(df)
    save_metrics(ml_results, (training_data, training_labels, target_data, target_labels, cluster_to_ko, quantiles))
    return training_data, training_labels, target_data, target_labels, cluster_to_ko, quantiles


target_taxa=read_target_taxa_tsv()
ref_taxa=read_proteomes_taxa_tsv()
#target_taxa={}
#target_taxa=add_bens_taxa(target_taxa)
proteomes_tsv='proteomes.tsv'
proteomes_path='proteomes.tab'

samples=get_samples()
df=get_dataframe(samples)

#df.T.to_csv(path+'/PhD/Presentations/cet_2/df.tsv', sep = '\t')

#plot_distribution(df,'Lipid metabolism',show_taxa=True)
#plot_all_ko_distributions(df)
#test_all_clustering_algorithms(df)
#test_cluster(df)
#main(df)
training_data, training_labels, target_data, target_labels,cluster_to_ko,quantiles=load_ml_results(df)
ko_description='Lipid metabolism'
plot_distribution_clustered(training_data, training_labels, target_data, target_labels, ko_description,cluster_to_ko,quantiles,show_taxa=False,add_ref=True,show_mix=False)
ko_description='Xenobiotics biodegradation and metabolism'
#plot_distribution_clustered(training_data, training_labels, target_data, target_labels, ko_description,cluster_to_ko,quantiles,show_taxa=False,add_ref=False,show_mix=False)

#get_training_labels(df)
#get_high_ko(df,ko_description)
#plot_distribution(df,ko_description='Lipid metabolism')
#plot_all_ko_distributions(training_data, training_labels, target_data, target_labels,cluster_to_ko,quantiles,add_ref=True,show_taxa=False,show_mix=False)
#test_classifiers(training_data=training_data,training_labels=training_labels)


#2 cluster score 0.09537799511608627
#3 cluster score 0.04835151124517616
#cluster score 0.03813917560058125
