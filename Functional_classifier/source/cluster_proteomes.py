from math import sqrt
import os
import re
import math
import pandas as pd
import numpy as np
import pickle
from sklearn.metrics import silhouette_score,davies_bouldin_score,classification_report
from sklearn import mixture as skmixture
from sklearn.preprocessing import StandardScaler
import statistics
import random
from sklearn import cluster as skcluster
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
from sklearn.decomposition import PCA
import sklearn.feature_selection as skfeature_select
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier,GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
from kneed import KneeLocator
import leidenalg
import igraph as ig

mash_genomes='D:/Data/Annotation_Output/reference_proteomes/mash_output_genomes/'
mash_proteomes='D:/Data/Annotation_Output/reference_proteomes/mash_proteomes'
pickle_path='D:/Data/Annotation_Output/reference_proteomes/distance_matrix.pickle'
cluster_path='D:/Data/Annotation_Output/reference_proteomes/cluster.pickle'
graph_path='D:/Data/Annotation_Output/reference_proteomes/distance_graph.gml'

def save_metrics(pickle_path,quality_genomes):
    with open(pickle_path, 'wb') as handle:
        pickle.dump(quality_genomes, handle,protocol=4)

def load_metrics(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as handle:
            quality_genomes = pickle.load(handle)
            return quality_genomes



def extract_mash(similarity_file,distance_matrix):
    with open(similarity_file) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            query=line[0].split('/')[-1].split('.')[0]
            reference=line[1].split('/')[-1].split('.')[0]
            similarity=round(float(line[2]),4)
            if query not in distance_matrix: distance_matrix[query]={}
            #if reference not in distance_matrix: distance_matrix[reference]={}
            #if query not in distance_matrix[reference]:
            distance_matrix[query][reference]=similarity
            line=file.readline()


def build_distance_matrix(mash_output_files):
    distance_matrix=load_metrics(pickle_path)
    if distance_matrix:
        return distance_matrix
    print('building distance matrix')
    distance_matrix={}
    for i in mash_output_files:
        extract_mash(i,distance_matrix)
    save_metrics(pickle_path,distance_matrix)
    return distance_matrix

def file_len(fname):
    try:
        with open(fname) as f:
            for i, l in enumerate(f):
                pass
        return i + 1
    except: return 0


def optimalK(data, clustering_function, nrefs=3, maxClusters=15):
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
    return gaps.argmax() + 1


def get_optimal_cluster_number(df, clustering_function, random_state=None):
    silhouette_coefficients = []
    n_kos = 20
    algorithm = str(clustering_function).split('.')[-1].strip('\'>').lower()
    print('clustering method', algorithm)
    # print('Optimal cluster number gap static',optimalK(df,clustering_function,maxClusters=n_kos))
    if 'dbscan' in algorithm:
        clustering = clustering_function(eps=0.5, min_samples=5, metric='euclidean',
                                         metric_params=None, algorithm='auto', leaf_size=30, p=None,
                                         n_jobs=None).fit(df)
        labels = clustering.labels_
        score = silhouette_score(df, labels)
        print("Silhouette Coefficient:", score)
        score = davies_bouldin_score(df, labels)
        print("davies bouldin score:", score)
        return
    elif \
            'affinity' in algorithm or \
                    'meanshift' in algorithm:
        if random_state is not None:
            clustering = clustering_function(random_state=random_state).fit(df)
        else:
            clustering = clustering_function().fit(df)
        labels = clustering.labels_
        score = silhouette_score(df, labels)
        print("Silhouette Coefficient:", score)
        score = davies_bouldin_score(df, labels)
        print("davies bouldin score:", score)
        return
    sse = []
    silhouette_coefficients = []
    davies_coefficients = []

    if 'optics' in algorithm:
        for k in [0.5, 1, 2, 3, 4, 5]:
            clustering = clustering_function(min_samples=k).fit(df)
            labels = clustering.labels_
            if len(set(labels)) > 1:
                score = silhouette_score(df, labels)
                silhouette_coefficients.append(score)
                score = davies_bouldin_score(df, labels)
                davies_coefficients.append(score)

    else:
        for k in range(2, n_kos):
            if 'gaussian' in algorithm:
                clustering = clustering_function(n_components=k).fit(df)
                labels = clustering.fit_predict(df)
            elif random_state is not None:
                clustering = clustering_function(n_clusters=k, random_state=random_state).fit(df)
                labels = clustering.labels_
            elif 'kmeans' in algorithm:
                clustering = clustering_function(n_clusters=k).fit(df)
                labels = clustering.labels_
                sse.append(clustering.inertia_)
            else:
                clustering = clustering_function(n_clusters=k).fit(df)
                labels = clustering.labels_
            score = silhouette_score(df, labels)
            print(k, score)
            silhouette_coefficients.append(score)
            score = davies_bouldin_score(df, labels)
            davies_coefficients.append(score)
        if 'kmeans' in algorithm and sse:
            kl = KneeLocator(range(2, n_kos), sse, curve="convex", direction="decreasing")
            print('optimal cluster number elbow', kl.elbow + 2)
    print('optimal cluster number silhouette ', silhouette_coefficients.index(max(silhouette_coefficients)) + 2)
    print("Silhouette Coefficient:", max(silhouette_coefficients))
    print('optimal cluster number Davies Bouldin ', davies_coefficients.index(min(davies_coefficients)) + 2)
    print("Davies Bouldin Coefficient:", min(davies_coefficients)+2)
    return None


def build_cluster(training_data, save_output=True, clustering_method='kmeans', test=False, n_cluster=None):
    res = load_metrics(cluster_path)
    if res:
        print('loading cluster')
        return res
    else:
        print('building cluster')

    clustering = None
    max_score = 0
    if clustering_method == 'agglomerative':
        clustering_function = skcluster.AgglomerativeClustering
        if not n_cluster:
            n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            #clustering = clustering_function(n_clusters=n_cluster).fit(training_data)
            clustering = clustering_function(n_clusters=n_cluster).fit(training_data)
            cluster_score = silhouette_score(training_data, clustering.labels_)
            print(cluster_score,len(clustering.labels_))


    elif clustering_method == 'kmeans':
        for i in range(1):
            clustering_function = skcluster.KMeans
            if not n_cluster:
                n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
            if not test:
                clustering = clustering_function(n_clusters=n_cluster, random_state=None).fit(training_data)
                cluster_score = silhouette_score(training_data, clustering.labels_)
                print(i, cluster_score)
                if cluster_score > max_score:
                    best_cluster = clustering
                    max_score = cluster_score
        if not test:
            clustering = best_cluster
    elif clustering_method == 'affinity':
        clustering_function = skcluster.AffinityPropagation
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function, random_state=0)
        if not test:
            clustering = clustering_function(random_state=None).fit(training_data)
    elif clustering_method == 'meanshift':
        clustering_function = skcluster.MeanShift
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function, )
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method == 'dbscan':
        clustering_function = skcluster.DBSCAN
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method == 'optics':
        clustering_function = skcluster.OPTICS
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method == 'gaussian':
        clustering_function = skmixture.GaussianMixture
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method == 'bayesian_gaussian':
        clustering_function = skmixture.BayesianGaussianMixture
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    elif clustering_method == 'birch':
        clustering_function = skcluster.Birch
        n_cluster = get_optimal_cluster_number(training_data, clustering_function=clustering_function)
        if not test:
            clustering = clustering_function().fit(training_data)
    save_metrics(cluster_path, (training_data.index,clustering))
    return training_data.index,clustering

def test_all_clustering_algorithms(df):
    for alg in [
        # kmeans (k=7 or 15) and meanshift are the best
         #'agglomerative',
         #'kmeans',
         #'gaussian',
         #'bayesian_gaussian',
         #'birch',
    ]:
        cluster= build_cluster(df, clustering_method=alg, test=True)






def generate_heatmap(distance_matrix):
    z=[]
    samples=list(distance_matrix.keys())
    for s in samples:
        temp=[distance_matrix[s][i] for i in samples]
        z.append(temp)
    fig = go.Figure(data=go.Heatmap(
        z=z,
        x=samples,
        y=samples,
        hoverongaps=False))
    fig.show()


def build_graph(distance_matrix):
    print('building graph')
    edges=[]
    edges_weights=[]
    samples=list(distance_matrix.keys())
    tuple_list=[]
    for s1 in samples:
        for s2 in samples:
            if s1!=s2:
                if distance_matrix[s1][s2]<0.05:
                    if (s1,s2) not in tuple_list and (s2,s1) not in tuple_list:
                        tuple_list.append((s1,s2))
    print('built tuple_list')
    g = ig.Graph.TupleList(tuple_list, weights=False)
    #ig.write(g,graph_path, format="gml")
    return g

def cluster_graph(distance_matrix,threshold=0.3):
    distance_matrix_graph = build_graph(distance_matrix)
    clusters=leidenalg.find_partition(distance_matrix_graph, leidenalg.ModularityVertexPartition)



#distance_matrix=build_distance_matrix([mash_genomes+i for i in os.listdir(mash_genomes)])
distance_matrix=build_distance_matrix([mash_proteomes])


#distance_matrix_graph=build_graph(distance_matrix)
#distance_matrix_graph=ig.load(graph_path)
#p=distance_matrix_graph.community()
#print('community_fastgreedy')
#print(p)

#part = leidenalg.find_partition(distance_matrix_graph, leidenalg.ModularityVertexPartition)
#print(part)
#for i in part:
#    print(i)
df=pd.DataFrame(distance_matrix)
print(df.shape)
#generate_heatmap(distance_matrix)
def count_labels(cluster):
    res={}
    for i in cluster.labels_:
        if i not in res: res[i]=0
        res[i]+=1
    print(res)

for i in range(20,50):
    print('building cluser',i)
    clustering_function = skcluster.AgglomerativeClustering
    clustering = clustering_function(n_clusters=i).fit(df)
    cluster_score = silhouette_score(df, clustering.labels_)
    print('cluster score',cluster_score)
    count_labels(clustering)
    print('######################################################')

