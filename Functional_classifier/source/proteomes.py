import requests
from math import sqrt
import os
import re
import plotly.graph_objects as go
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
from kneed import KneeLocator
from sklearn.metrics import silhouette_score,davies_bouldin_score,classification_report
from sklearn import mixture as skmixture
from sklearn.preprocessing import StandardScaler
import statistics
import random
from sklearn import cluster as skcluster
from scipy.cluster import hierarchy
import seaborn as sns; sns.set_theme()
from sklearn.decomposition import PCA
import plotly.express as px
import sklearn.feature_selection as skfeature_select
from sklearn.ensemble import ExtraTreesClassifier,RandomForestClassifier,GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn import svm
from sklearn.neural_network import MLPClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import train_test_split
import plotly.figure_factory as ff

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


def save_metrics(pickle_path,quality_genomes):
    with open(path+'/PhD/Presentations/cet_2/'+pickle_path, 'wb') as handle:
        pickle.dump(quality_genomes, handle,protocol=4)

def load_metrics(pickle_path):
    if os.path.exists(path+'/PhD/Presentations/cet_2/'+pickle_path):
        with open(path+'/PhD/Presentations/cet_2/'+pickle_path, 'rb') as handle:
            quality_genomes = pickle.load(handle)
            return quality_genomes

def get_proteome_ids(file_path):
    res={}
    with open(file_path) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            proteome_id,org_id,cpd=line[0],line[2],line[-1]
            if cpd=='Standard':
                res[proteome_id]=org_id
            line=file.readline()
    return res

def download_proteome(proteome_id,folder_path):
    list_dir= os.listdir(folder_path)
    if proteome_id+'.faa' not in list_dir:
        link='https://www.uniprot.org/uniprot/?query=proteome:'+proteome_id+'&format=fasta'
        r = requests.get(link)
        open(folder_path+proteome_id+'.faa', 'wb').write(r.content)

def download_all_proteomes(proteome_ids,proteomes_dump):
    c=1
    for p_id in proteome_ids:
        print('Downloading',p_id,str(c)+'/'+str(len(proteome_ids)))
        download_proteome(p_id,proteomes_dump)
        c+=1

def write_mantis_tsv(proteome_ids,proteomes_tsv):
    with open(proteomes_tsv,'w+') as file:
        for p_id in proteome_ids:
            file.write(p_id+'\t'+'path_here/'+p_id+'.faa\t'+proteome_ids[p_id]+'\n')



