import re

import plotly.graph_objects as go

def get_sample_size(sample_path):
    with open(sample_path) as file:
        line=file.readline()
        c=-1
        while line:
            line=file.readline()
            c+=1
    return c

def read_kegg_map():
    res={}
    level_1_pattern=re.compile('\d\.\s')
    level_2_pattern=re.compile('\d\.\d+\s')
    level_3_pattern=re.compile('\d{3,}:')
    with open('kegg_map.txt') as file:
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

def get_kegg_map_sample(sample):
    kegg_maps=[]
    with open(sample) as file:
        line=file.readline()
        while line:
            if 'kegg_map:' in line:
                line=line.strip('\n')
                line=line.split('\t')
                line_kegg_maps={i.split('map:')[1] for i in line if 'kegg_map:' in i}
                kegg_maps.extend(line_kegg_maps)
            line=file.readline()
    return kegg_maps

def get_description(wanted_id,kegg_map,wanted_level=2):
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
                            if wanted_id == level_3:
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


def get_frequency_kegg_map(sample,kegg_map,wanted_level=2):
    sample_kegg_maps = get_kegg_map_sample(sample)
    freq={}
    for km in sample_kegg_maps:
        if km not in freq: freq[km]=0
        freq[km]+=1
    freq_names={}
    for km in freq:
        description=get_description(km,kegg_map,wanted_level=wanted_level)
        if description not in freq_names: freq_names[description]=0
        freq_names[description]+=freq[km]
    #for km in freq_names: print(km,freq_names[km])
    freq_names.pop(None)
    return freq_names

def get_top_frequencies(frequency_table,top_n=10):
    sorted_frequency_table = sorted(frequency_table.items(), key=lambda item: item[1], reverse=True)
    total=sum(frequency_table.values())
    c=0
    res={}
    for k in sorted_frequency_table:
        description,frequency=k
        if c>top_n:
            break
        res[description]=round(100*frequency/total,2)
        c+=1
    res['others']=0
    for k in sorted_frequency_table:
        description,frequency=k
        if description not in res:
            res['others']+=round(100*frequency/total,2)
    return res

def get_wanted_frequencies(frequency_table,wanted_freqs,sample_total_kos):
    total=sum(frequency_table.values())
    total=sample_total_kos
    res={}
    res['others']=0
    for k in frequency_table:
        if k not in wanted_freqs:
            res['others']+=round(100*frequency_table[k]/total,2)
        else:
            if k not in res: res[k]=0
            res[k]+=round(100*frequency_table[k]/total,2)
    return res

def get_most_common_top(samples,top_n=10,wanted_level=2):
    kegg_map = read_kegg_map()
    sample_results=[]
    freq_top={}
    samples_total_ko=[]
    for s in samples:
        freq_names = get_frequency_kegg_map(s, kegg_map, wanted_level)
        samples_total_ko.append(sum(freq_names.values()))
        #print(s,freq_names)
        sample_results.append(freq_names)
        top_freqs = get_top_frequencies(freq_names,top_n=top_n)
        for tf in top_freqs:
            if tf not in ['others',]:#'Cellular community - prokaryotes']:
                if tf not in freq_top: freq_top[tf] = 0
                freq_top[tf]+=1
    sorted_frequency_table = sorted(freq_top.items(), key=lambda item: item[1], reverse=True)
    sorted_frequency_table=sorted_frequency_table[0:top_n+1]
    sorted_frequency_table=[i[0] for i in sorted_frequency_table]
    #print(sorted_frequency_table)
    res={}
    for s in range(len(sample_results)):
        top_wanted_freq=get_wanted_frequencies(sample_results[s],sorted_frequency_table,samples_total_ko[s])
        res[samples[s]]=top_wanted_freq
        print('sample: ',samples[s])
        #for sft in sorted_frequency_table:
        #    print(sft,top_wanted_freq[sft])
        print('others',top_wanted_freq['others'])
    return res



samples_oskar=[
    #'C:/Users/Pedro Queirós/Dropbox/PhD/Papers/must/Gut_MG_Kegg_sig_outkos_to_map.txt',
    'C:/Users/Pedro Queirós/Dropbox/PhD/Papers/must/Gut_MT_Kegg_sig_outkos_to_map.txt',
    #'C:/Users/Pedro Queirós/Dropbox/PhD/Papers/must/Oral_MG_Kegg_sig_outkos_to_map.txt',
    'C:/Users/Pedro Queirós/Dropbox/PhD/Papers/must/Oral_MT_Kegg_sig_outkos_to_map.txt',
]

samples_top_common=get_most_common_top(samples_oskar,wanted_level=2)
def plot_heatmap(samples_top_common):
    x=[]
    y=[]
    z=[]
    for sample in samples_top_common:
        sample_str=sample.split('/')[-1].split('.')[0]
        print(sample_str)
        y.append(sample_str)
        sorted_frequency_table=samples_top_common[sample]
        temp_z=[]
        for sft in sorted_frequency_table:
            if sft !='others':
                temp_z.append(sorted_frequency_table[sft])
                if sft not in x: x.append(sft)
        z.append(temp_z)
        #if 'others' not in x:x.append('others')


    fig = go.Figure(data=go.Heatmap(
                       z=z,
                       x=x,
                       y=y,
                       #hoverongaps = False,
                       colorscale='Cividis',
    ))
    fig.show()

plot_heatmap(samples_top_common)