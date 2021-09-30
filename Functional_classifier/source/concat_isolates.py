import os

def rename_gene_call(file_path,sample):
    with open(file_path) as file:
        line=file.readline()
        while line:
            if '>' in line[0]:
                line=line.split(' ')[0]
                line=line.strip('\n')
                line+='\n'
                line=line[0]+sample+'_'+line[1:]
            yield line
            line=file.readline()

def read_samples(output_file,samples_fasta):
    current_dir=os.getcwd()
    print(current_dir)
    merged=open(current_dir+'/'+output_file,'w+')
    for i in os.listdir(current_dir):
        if 'Isolate_' in i:
            sample = i.split('/')[-1].split('.')[0]
            print(sample)
            if os.path.isdir(current_dir+'/'+i):
                isolate_dir=current_dir+'/'+i+'/'
                fasta_file=isolate_dir+samples_fasta
            else:
                fasta_file=current_dir+'/'+i
            print(fasta_file)
            lines=rename_gene_call(fasta_file,sample)
            for l in lines:
                merged.write(l)
    create_gff(merged)

def create_gff(merged_path):
    with open(merged_path) as file:
        line=file.readline()
        while line:
            line=line.strip('\n')

            line=file.readline()




read_samples('merged.faa','prodigal_proteins.faa')

