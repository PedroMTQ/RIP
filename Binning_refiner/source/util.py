from sys import path as sys_path
if 'source' in sys_path[0]:
    from re import search as re_search
    p_to_add=sys_path[0][:re_search('Binning_refiner',sys_path[0]).span()[1]]
    if not sys_path[0].endswith('Binning_refiner') :  sys_path[0]=p_to_add

# External modules
import ast
import types
import difflib
import io
import os
import re
import sys
from html.parser import HTMLParser
from io import StringIO
from urllib.parse import quote_plus
from math import ceil
import urllib.request as request
import subprocess
import pickle

from contextlib import closing
import requests

from types import GeneratorType as generator
import csv
from sys import platform
from gzip import open as gzip_open
from pathlib import Path
from multiprocessing import cpu_count
from functools import wraps

if platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

current_path = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-1]
MAIN_FOLDER=[]
for i in current_path:
    if i=='Binning_refiner':
        MAIN_FOLDER.append(i)
        break
    MAIN_FOLDER.append(i)
MAIN_FOLDER = SPLITTER.join(MAIN_FOLDER) + SPLITTER
KEGG_MODULES=f'{MAIN_FOLDER}Data{SPLITTER}KEGG{SPLITTER}modules.pickle'


# to use pubchem API you can use requests.post
# inchi='InChI=1S/C15H12N2O2/c16-15(18)17-11-7-3-1-5-9(11)13-14(19-13)10-6-2-4-8-12(10)17/h1-8,13-14H,(H2,16,18)'
# r=requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON/',data={'inchi':inchi})
# print(r.text)

def save_metrics(pickle_path, to_pickle):
    with open(pickle_path, 'wb') as handle:
        pickle.dump(to_pickle, handle)


def load_metrics(pickle_path):
    if os.path.exists(pickle_path):
        with open(pickle_path, 'rb') as handle:
            pickled_results = pickle.load(handle)
            return pickled_results

def sql_cursor_restart(f):
    @wraps(f)
    def wrapper(self, *args, **kwargs):
        if hasattr(self,'cursor'):
            self.stop_sqlite_cursor()
            res = f(self, *args, **kwargs)
            self.restart_sqlite_cursor()
            return res
    return wrapper




def xstr(s):
    if not s:
        return ''
    if isinstance(s,generator): s=list(s)
    if isinstance(s,list):
        return ' , '.join(s)
    return str(s)

def is_ec(enz_id,required_level=3):
    if enz_id:
        ec_pattern = re.compile('^\d+\.\d+\.\d+(\.(-|\d+|([a-zA-Z]\d+)|[a-zA-Z]))?')
        if re.search(ec_pattern, enz_id):
            enz_id_copy=str(enz_id).replace('.-','')
            if len(enz_id_copy.split('.'))>=required_level:
                return True
    return False








def find_sign(string_to_search):
    #first search should be for one of these signs, which means directionality is set
    compiled_sign = re.compile('( = | <=> | <= | => | → | ← | ↔ | ⇄ | <-> | <--> | <- | -> | [=-]?&gt; | &lt;=&gt; | &lt;[=-]? )')
    sign = re.search(compiled_sign, string_to_search)
    if sign:
        return sign
    else:
        #the second search should be for unknown directionality
        #this way it avoids recognizing unknown compounds "?" as the direction sign of a reaciton
        interrogation_sign=re.compile('\?')
        sign = re.search(interrogation_sign, string_to_search)
        if sign: return sign
    return None

def uniform_sign(sign):
    both = re.compile('( = | <=> | ↔ | ⇄ | <-> | <--> | &lt;=&gt; )')
    right = re.compile('( => | → | -> | [=-]?&gt; )')
    left = re.compile('( <= | ← | <- | &lt;[=-]? )')
    if re.search(both,sign):  return ' <=> '
    if re.search(right,sign):  return ' => '
    if re.search(left,sign):  return ' <= '

def fix_html_sign(string_to_fix):
    res=str(string_to_fix)
    if re.search('[=-]?&gt;',res):
        res=res.replace(re.search('[=-]?&gt;',res).group(),'=>')
    if re.search('&lt;[=-]?',res):
        res=res.replace(re.search('&lt;[=-]?',res).group(),'=>')
    if re.search('&lt;[=-]?',res):
        res=string_to_fix.replace(re.search('&lt;=&gt;',res).group(),'<=>')
    if re.search('%2b',res):
        res=string_to_fix.replace(re.search('%2b',res).group(),'%2b')
    return res

def is_reversible(reaction_string):
    if '<=>' in reaction_string: return True
    else: return False

def get_id_or_str(rn_cpd):
    if not isinstance(rn_cpd,str): return rn_cpd.get_most_common_synonym()
    else: return rn_cpd




def test_match_possible_ids(possible_ids_1,possible_ids_2):
    if not possible_ids_1 or not possible_ids_2: return False
    if isinstance(possible_ids_1,str) or isinstance(possible_ids_1,int):
        temp_possible_ids_1= {possible_ids_1:1}
    elif isinstance(possible_ids_1,dict):
        temp_possible_ids_1=set(possible_ids_1.keys())
    else: temp_possible_ids_1=possible_ids_1


    if isinstance(possible_ids_2,str) or isinstance(possible_ids_2,int):
        temp_possible_ids_2= {possible_ids_2:1}
    elif isinstance(possible_ids_2,dict):
        temp_possible_ids_2=set(possible_ids_2.keys())
    else: temp_possible_ids_2=possible_ids_2

    for i in temp_possible_ids_1:
        for j in temp_possible_ids_2:
            if i and j:
                if i == j:
                    return True
    return False


def score_match_possible_ids(possible_ids_1,possible_ids_2,match_c=1,mismatch_c=-0.25):
    """
    just a scoring system. if one instance doesnt have an id we dont penalize
    if the ids match we benefit the match, if they dont we add a minor penalty.
    The penalty is lower since we assume databases already have some curation, thus the amount of errors shouldnt be too high
    """

    if not possible_ids_1 or not possible_ids_2: return 0
    if isinstance(possible_ids_1,set) or isinstance(possible_ids_1,list):
        temp_possible_ids_1= {i:1 for i in possible_ids_1}
    elif not isinstance(possible_ids_1,dict):
        temp_possible_ids_1= {possible_ids_1:1}
    else:
        temp_possible_ids_1=dict(possible_ids_1)


    if isinstance(possible_ids_2,set) or isinstance(possible_ids_2,list):
        temp_possible_ids_2= {i:1 for i in possible_ids_2}
    elif not isinstance(possible_ids_2, dict):
        temp_possible_ids_2 = {possible_ids_2: 1}
    else:
        temp_possible_ids_2=dict(possible_ids_2)
    for i in temp_possible_ids_1:
        for j in temp_possible_ids_2:
            if i and j:
                if i == j:
                    if not is_ec(i):
                        return match_c
    return mismatch_c

def unite_possible_ids(self_instance,instance_2,detail_type):
    if instance_2.get_detail(detail_type):
        possible_ids = instance_2.get_detail(detail_type,all_possible=True)
        if possible_ids:
            for possible_id in possible_ids:
                possible_id_count = possible_ids[possible_id]
                self_instance.set_detail(detail_type, {possible_id: possible_id_count})


#only used when we invoke a specific detail type
def list_has_common_items(iterable1,iterable2):
    if iterable1 and iterable2:
        for i in iterable1:
            for j in iterable2:
                if i.lower()==j.lower(): return True
    return False



def get_instance_type(instance_to_test):
    res=str(type(instance_to_test))
    res=res.split('.')[-1].replace('\'>','')
    return res


def regex_escape(regex_string):
    regex_characters=[
        "\\",
        ".",
        "+",
        "*",
        "?",
        "[",
        "^",
        "]",
        "$",
        "(",
        ")",
        "{",
        "}",
        "=",
        "!",
        "<",
        ">",
        "|",
        "\'",
        "\"",
        ":",
        "-"]
    l_regex_string=[i for i in regex_string]
    for i in range(len(l_regex_string)):
        if l_regex_string[i] in regex_characters:
            if l_regex_string[i]=='\\': replacer='\\\\'
            else: replacer='\\'
            l_regex_string[i]=replacer+l_regex_string[i]
    return ''.join(l_regex_string)


def uncompress_archive(source_filepath, extract_path=None, block_size=65536, remove_source=False, stdout_file=None):
    file_name = source_filepath.split(SPLITTER)[-1]
    dir_path = SPLITTER.join(source_filepath.split(SPLITTER)[0:-1])
    if not extract_path: extract_path = dir_path
    if '.tar' in file_name:
        unpack_archive(source_file=source_filepath, extract_dir=extract_path, remove_source=remove_source,
                       stdout_file=None)
    # only for files
    elif '.gz' in file_name:
        gunzip(source_filepath=source_filepath, dest_filepath=extract_path, block_size=block_size,
               remove_source=remove_source, stdout_file=stdout_file)
    elif '.zip' in file_name:
        unzip_archive(source_file=source_filepath, extract_dir=extract_path, remove_source=remove_source,
                      stdout_file=None)
    else:
        print('Incorrect format! ', source_filepath, flush=True, file=stdout_file)


# this unzips to the same directory!
def gunzip(source_filepath, dest_filepath=None, block_size=65536, remove_source=False, stdout_file=None):
    if not dest_filepath:
        dest_filepath = source_filepath.strip('.gz')
    if os.path.isdir(dest_filepath):
        file_name = source_filepath.split(SPLITTER)[-1].replace('.gz', '')
        dest_filepath = dest_filepath+SPLITTER + file_name
    print('Gunzipping ', source_filepath, 'to', dest_filepath, flush=True, file=stdout_file)
    with gzip_open(source_filepath, 'rb') as s_file, \
            open(dest_filepath, 'wb') as d_file:
        while True:
            block = s_file.read(block_size)
            if not block:
                break
            else:
                d_file.write(block)
        d_file.write(block)
    if remove_source: os.remove(source_filepath)


def unpack_archive(source_file, extract_dir, remove_source=False, stdout_file=None):
    print('Unpacking', source_file, 'to', extract_dir, flush=True, file=stdout_file)
    shutil.unpack_archive(source_file, extract_dir=extract_dir)
    if remove_source: os.remove(source_file)


def unzip_archive(source_file, extract_dir, remove_source=False, stdout_file=None):
    print('Unzipping', source_file, 'to', extract_dir, flush=True, file=stdout_file)
    with ZipFile(source_file, 'r') as zip_ref:
        zip_ref.extractall(extract_dir)
    if remove_source: os.remove(source_file)


def download_file_http(url, file_path, c):
    if c > 5:
        download_file_http_failsafe(url, file_path)
    else:
        with requests.get(url, stream=True) as r:
            with open(file_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)


# slower but safer
def download_file_http_failsafe(url, file_path):
    with requests.Session() as session:
        get = session.get(url, stream=True)
        if get.status_code == 200:
            with open(file_path, 'wb') as f:
                for chunk in get.iter_content(chunk_size=1024):
                    f.write(chunk)


def download_file_ftp(url, file_path):
    with closing(request.urlopen(url)) as r:
        with open(file_path, 'wb') as f:
            shutil.copyfileobj(r, f)


def download_file(url, output_folder='', stdout_file=None, retry_limit=10):
    file_path = output_folder + url.split('/')[-1]
    try:
        target_file = request.urlopen(url)
    except:
        print('Cannot download target url', url)
        return
    target_size = target_file.info()['Content-Length']
    transfer_encoding = target_file.info()['Transfer-Encoding']
    if target_size: target_size = int(target_size)
    if os.path.exists(file_path):
        if transfer_encoding == 'chunked':
            return
        elif os.stat(file_path).st_size == target_size:
            print('Not downloading from ' + url + ' since file was already found!', flush=True, file=stdout_file)
            return
        else:
            os.remove(file_path)
    print('Downloading from ' + url + '. The file will be kept in ' + output_folder, flush=True, file=stdout_file)
    c = 0
    while c <= retry_limit:
        if 'ftp' in url:
            try:
                download_file_ftp(url, file_path)
            except:
                try:
                    download_file_http(url, file_path, c)
                except: pass
        else:
            try:
                download_file_http(url, file_path, c)
            except:
                pass
        if transfer_encoding == 'chunked': return
        if os.stat(file_path).st_size == target_size: return
        c += 1
    print('Did not manage to download the following url correctly:\n' + url)
    raise Exception


def get_slurm_value(wanted_val, regex_pattern):
    res = None
    slurm_job_id = os.environ.get('SLURM_JOBID')
    if slurm_job_id:
        process = subprocess.run('sacct -j ' + str(slurm_job_id) + ' -o ' + wanted_val, shell=True,stdout=subprocess.PIPE)
        wanted = re.search(regex_pattern, str(process.stdout))
        if wanted: res = wanted.group()
    return res

def check_environment_cores():
    res = get_slurm_value('AllocCPUS', re.compile('\d+'))
    if res:
        if int(res):
            print('Cores allocated by slurm:', res)
            return int(res)
        else:
            res = cpu_count()
            print('Cores allocated:', res)
            return int(res)
    else:
        res = cpu_count()
        print('Cores allocated:', res)
        return int(res)

if __name__ == '__main__':
    print(MAIN_FOLDER)