


def read_kegg_kos(kegg_kos):
    #C  ID => kegg map
    #D  ID => kegg KO
    kos_in_map={}
    ko_to_map={}
    with open(kegg_kos) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')
            line=line.split()
            line_type=line[0]
            if len(line)>1:
                line_id=line[1]
                line_description=' '.join(line[2:])
                if line_type=='C':
                    kos_in_map[line_id]=set()
                    current_map=line_id
                elif line_type=='D':
                    kos_in_map[current_map].add(line_id)
                    ko_to_map[line_id]=current_map
            line=file.readline()
    return kos_in_map,ko_to_map




annot_path='/home/pedroq/Dropbox/PhD/Presentations/cet_2/microthrix_bio17.2020-10-05T204110/consensus_annotation.tsv'

def connect_sample_to_map(sample):
    kegg_kos = 'kegg_kos.txt'
    kos_in_map, ko_to_map = read_kegg_kos(kegg_kos)
    with open(sample) as sample_file:
        with open('/'.join(sample.split('/')[0:-1])+'/'+'kos_to_map.txt','w+') as output_file:
            line=sample_file.readline()
            line=sample_file.readline()
            while line:
                line = line.strip('\n')
                line= line.split('\t')
                sample_ko=None
                if sample_ko in ko_to_map:
                    output_file.write('kegg_map:'+str(ko_to_map[sample_ko])+'\n')
                else:
                    print('Unresolved KO:',sample_ko)
                line=sample_file.readline()



connect_sample_to_map(annot_path)