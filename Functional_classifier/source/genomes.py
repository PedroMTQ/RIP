import os
import requests
import gzip
import shutil
import re
from difflib import SequenceMatcher

def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

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


def get_proteome_ids(file_path):
    res={}
    with open(file_path) as file:
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')

            proteome_id,organism,taxa_id,cpd=line[0],line[1],line[2],line[-1]
            strain_search=re.search('\(Strain: .*\)',organism)
            strain=None
            if strain_search: strain=strain_search.group().split('Strain:')[1].strip(')')
            if strain:
                if ' / ' in strain: strain=strain.split(' / ')
                else: strain=[strain]
            else: strain=[]
            strain=[i.strip() for i in strain]
            if cpd=='Standard':
                res[taxa_id]=[proteome_id,organism,strain]
            line=file.readline()
    return res

def download_and_unzip(link,folder_path,taxa_id):
    list_dir= os.listdir(folder_path)
    file_name = link.split('/')[-1]

    link=link.replace('ftp://','http://')
    if taxa_id+'.fasta' not in list_dir:
        r = requests.get(link)
        open(folder_path+file_name, 'wb').write(r.content)
        with gzip.open(folder_path+file_name, 'rb') as f_in:
            with open(folder_path+taxa_id+'.fasta', 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(folder_path+file_name)
    print(taxa_id)


def get_ftp_links(file_path,taxa_ids):
    res={}
    with open(file_path,encoding='utf8') as file:
        line=file.readline()
        line=file.readline()
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split('\t')
            taxa_id,genome_type,strain,ftp_link,assembly_level,genome_rep=line[5],line[4],line[8],line[-3],line[11],line[13]
            strain=strain.replace('strain=','').split(';')[0]
            if taxa_id in taxa_ids:
                assembly_id=ftp_link.split('/')[-1]
                if taxa_id not in res: res[taxa_id]={'strain':{}}
                link=ftp_link+'/'+assembly_id+'_genomic.fna.gz'
                res[taxa_id]['strain'][strain]={'genome_type':genome_type,'link':link,'assembly_level':assembly_level,'genome_rep':genome_rep}
                #downloading all genomes (there may be multiple genomes per taxon id- need to choose only representative genomes
            #else: print(assembly_level)
            line=file.readline()
    return res

def sort_strains(ftp_links,taxa_ids):
    res={}
    to_search=[]
    for i in taxa_ids:
        proteome_strains = set(taxa_ids[i][2])
        if i in ftp_links:
            genbank_strains = set()
            for ps in ftp_links[i]['strain']:
                if ftp_links[i]['strain'][ps]['assembly_level'] in ['Scaffold','Complete Genome']:
                    ps = ps.replace('type_strain:', '')
                    ps = ps.strip()
                    genbank_strains.add(ps)
            if genbank_strains:
                #when the proteome doesnt have a specific strain
                if not proteome_strains:
                    for ps in genbank_strains:
                        if ftp_links[i]['strain'][ps]['genome_type']=='representative genome':
                            res[i]=ftp_links[i]['strain'][ps]['link']
                            break
                    for ps in genbank_strains:
                        if ftp_links[i]['strain'][ps]['assembly_level']in ['Scaffold','Complete Genome']:
                            res[i]=ftp_links[i]['strain'][ps]['link']
                            break
                #when it does
                else:
                    #checking which strain intersect
                    strain_intersection=proteome_strains.intersection(genbank_strains)
                    if strain_intersection:
                        for ps in strain_intersection:
                            if ftp_links[i]['strain'][ps]['genome_type'] == 'representative genome':
                                res[i] = ftp_links[i]['strain'][ps]['link']
                                break
                        for ps in strain_intersection:
                            if ftp_links[i]['strain'][ps]['assembly_level'] in ['Scaffold','Complete Genome']:
                                res[i] = ftp_links[i]['strain'][ps]['link']
                                break
                    #if there is no intersection we apply a similarity between name
                    else:
                        strain_intersection=set()
                        for ps_proteome in proteome_strains:
                            for ps_genbank in genbank_strains:
                                similarity=similar(ps_proteome.lower(),ps_genbank.lower())
                                if similarity>0.5:
                                    strain_intersection.add(ps_genbank)
                        #we add strain names which are very similar
                        if strain_intersection:
                            for ps in strain_intersection:
                                if ftp_links[i]['strain'][ps]['genome_type'] == 'representative genome':
                                    res[i] = ftp_links[i]['strain'][ps]['link']
                                    break
                            for ps in strain_intersection:
                                if ftp_links[i]['strain'][ps]['assembly_level'] in ['Scaffold', 'Complete Genome']:
                                    res[i] = ftp_links[i]['strain'][ps]['link']
                                    break
                        #when there's no similar names in strains we just hope for the best
                        if not strain_intersection:
                            for ps in genbank_strains:
                                if ftp_links[i]['strain'][ps]['genome_type'] == 'representative genome':
                                    res[i] = ftp_links[i]['strain'][ps]['link']
                                    break
                            for ps in genbank_strains:
                                if ftp_links[i]['strain'][ps]['assembly_level'] in ['Scaffold', 'Complete Genome']:
                                    res[i] = ftp_links[i]['strain'][ps]['link']
                                    break

            else:
                genbank_strains = set()
                for ps in ftp_links[i]['strain']:
                    if ftp_links[i]['strain'][ps]['genome_type'] =='representative genome':
                        ps = ps.replace('type_strain:', '')
                        ps=ps.strip()
                        genbank_strains.add(ps)
                if genbank_strains:
                    for ps in genbank_strains:
                        res[i] = ftp_links[i]['strain'][ps]['link']
                        break
                else:
                    genbank_strains = set()
                    for ps in ftp_links[i]['strain']:
                        if ftp_links[i]['strain'][ps]['genome_rep'] == 'Full':
                            ps = ps.replace('type_strain:', '')
                            ps = ps.strip()
                            genbank_strains.add(ps)
                    if genbank_strains:
                        for ps in genbank_strains:
                            res[i] = ftp_links[i]['strain'][ps]['link']
                            break
                    else:
                        genbank_strains = set()
                        for ps in ftp_links[i]['strain']:
                            ps = ps.replace('type_strain:', '')
                            ps = ps.strip()
                            genbank_strains.add(ps)
                        if len(genbank_strains)==1:
                            ps=genbank_strains.pop()
                            res[i] = ftp_links[i]['strain'][ps]['link']
        else: to_search.append(i)
    #for taxa_id in res:
    #    download_and_unzip(res[taxa_id],'D:/Data/Annotation_Output/reference_proteomes/genomes/',taxa_id)

    for i in taxa_ids:
        print(itaxa_ids[i])

def rename_fastas(taxa_ids):
    for taxa_fasta in os.listdir(genomes_dir):
        taxa_id=taxa_fasta.split('.')[0]
        sample=taxa_ids[taxa_id][0]+'.fasta'
        os.rename(genomes_dir+taxa_fasta,genomes_dir+sample)

genomes_dir='D:/Data/Annotation_Output/reference_proteomes/genomes/'
proteomes_path='C:/Users/Pedro Queirós/Dropbox/PhD/Presentations/cet_2/ML/proteomes.tab'
taxa_ids=get_proteome_ids(proteomes_path)
rename_fastas(taxa_ids)
#assembly_path='D:/Data/Annotation_Output/reference_proteomes/assembly_summary_genbank.txt'
#ftp_links=get_ftp_links(assembly_path,taxa_ids)
#sort_strains(ftp_links,taxa_ids)



