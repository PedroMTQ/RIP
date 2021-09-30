import os
import subprocess
from multiprocessing import Process, current_process,cpu_count, Manager,Pool
from time import sleep

def get_seqs(target_sample):
    res=[]
    with open(target_sample) as file:
        line=file.readline()
        while line:
            if line[0]=='>':
                line=line.strip(('\n'))
                line=line[1:]
                line=line.split('|')[1]
                res.append(line)
            line=file.readline()
    return res




def get_proteome_annotations():
    folder_proteomes='D:/Data/Annotation_Output/reference_proteomes/faa_folder/'
    #folder_ref_annotations='D:/Data/Annotation_Output/reference_proteomes/reference_annotations/'
    folder_ref_annotations='D:/Data/Annotation_Output/reference_proteomes/reference_annotations_all/'
    #all_proteins='D:/Data/Annotation_Output/reference_proteomes/reviewed_proteins.tsv'
    all_proteins='D:/Data/Annotation_Output/reference_proteomes/uniprot_references/all_proteins.tsv'
    res={}
    for proteome in os.listdir(folder_proteomes):
        seq_file=folder_proteomes+proteome
        out_file=folder_ref_annotations+proteome.replace('.faa','.tsv')
        seqs=get_seqs(seq_file)
        for seq in seqs:
            if seq not in res: res[seq]=[]
            res[seq].append(proteome)
    #print(res)
    print(len(res))
    #global_queue.append([all_proteins,seqs,out_file])
    #get_annotation_info_proteome(all_proteins,seqs,out_file)
    extract_reference(all_proteins,res)


def extract_reference(proteome_annot,seqs):
    folder_ref_annotations='D:/Data/Annotation_Output/reference_proteomes/reference_annotations_all/'
    with open(proteome_annot) as file:
        line=file.readline()
        while line and seqs:
            line=line.strip('\n')
            line=line.split('\t')
            entry=line[0]
            if entry in seqs:
                entry, protein_names, entry_name, ec_annotation, \
                go_ids, go_text, cazy_id, eggnog_id, kegg_ko, pfam_id, tigrfam_id, brenda_id, biocyc_id = line
                ec_annotation=set([i.strip() for i in ec_annotation.split(';') if i])
                go_ids=set([i.strip() for i in go_ids.split(';') if i])
                go_ids=[g.replace('GO:','').strip() for g in go_ids]
                protein_annotation=set([protein_names])
                #protein_annotation.add(entry_name)
                go_text=[i.strip() for i in go_text.split(';') if i]
                cazy_id=[i.strip() for i in cazy_id.split(';') if i]
                eggnog_id=[i.strip() for i in eggnog_id.split(';') if i]
                kegg_ko=[i.strip() for i in kegg_ko.split(';') if i]
                pfam_id=[i.strip() for i in pfam_id.split(';') if i]
                tigrfam_id=[i.strip() for i in tigrfam_id.split(';') if i]
                brenda_id=[i.strip() for i in brenda_id.split(';') if i]
                ec_annotation.update(brenda_id)
                protein_annotation.update(go_text)
                dict_annotation={
                                 'description':protein_annotation,
                                 'enzyme_ec':ec_annotation,
                                 'go':set(go_ids),
                                 'cazy':set(cazy_id),
                                 'eggnog':set(eggnog_id),
                                 'kegg_ko':set(kegg_ko),
                                 'pfam':set(pfam_id),
                                 'tigrfam':set(tigrfam_id),
                                 }
                line_list=[]
                for link_key in dict_annotation:
                    for link_str in dict_annotation[link_key]:
                        line_list.append(link_key + ':' + link_str)
                annotation_line='\t'.join(line_list)
                out_line=entry+'\t'+annotation_line+'\n'
                for proteome in seqs[entry]:
                    out_file_path = folder_ref_annotations + proteome.replace('.faa', '.tsv')
                    with open(out_file_path,'a+') as out_file:
                        out_file.write(out_line)
                seqs.pop(entry)

            line=file.readline()




if __name__ == '__main__':
    #generate_annotated_proteomes()
    get_proteome_annotations()
    #queue=generate_annotated_proteomes_queue()
    #processes_handler(queue,worker_get_annotation_info_proteome,cpu_count()-2)
