import os
import numpy as np
import sqlite3
from source.util import sql_cursor_restart
from time import sleep
from multiprocessing import Process


class SQLITE_Connector():
    def __init__(self):
        self.launch_sqlite_cursor()
        self.insert_step=50000
        self.fetch_step=500
        self.numpy_array_size=50000
        self.key_split='#_#'

    ###for multiprocessing
    @sql_cursor_restart
    def processes_handler(self, target_worker_function, add_sentinels=True):
        # os.getpid to add the master_pid
        if len(self.queue)<self.worker_count: worker_count=len(self.queue)
        else: worker_count=self.worker_count
        processes = [Process(target=target_worker_function, args=(self.queue, os.getpid(),)) for _ in range(worker_count)]
        # adding sentinel record since queue can be signaled as empty when its really not
        if add_sentinels:
            for _ in range(worker_count):   self.queue.append(None)
        for process in processes:
            process.start()
        # we could manage the processes memory here with a while cycle
        for process in processes:
            process.join()
            # exitcode 0 for sucessful exists
            if process.exitcode != 0:
                sleep(5)
                print('Ran into an issue, check the log for details. Exitting!')
                os._exit(1)

    def restart_sqlite_cursor(self):
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

    def stop_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def launch_sqlite_cursor(self):
        #this will probably need to be changed to an output_folder provided by the user
        self.db_file = f'{self.output_folder}refiner.db'
        print(self.db_file)
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()
        create_table_command = f'CREATE TABLE CONTIGS_VECTORS (' \
                            f'CONTIG TEXT UNIQUE,' \
                            f' VECTOR  TEXT )'
        self.cursor.execute(create_table_command)
        create_table_command = f'CREATE TABLE CONTIGS_DISTANCE (' \
                            f'CONTIG_KEY TEXT UNIQUE,' \
                            f'KMER_DEPTH_DISTANCE REAL,' \
                            f' TAXA_DISTANCE  REAL )'
        self.cursor.execute(create_table_command)

        create_table_command = f'CREATE TABLE REACTIONS (' \
                            f'REACTION_ID INTEGER,' \
                            f'REVERSIBILITY TEXT,' \
                            f'COMPOUNDS_LEFT TEXT,' \
                            f' COMPOUNDS_RIGHT  INTEGER )'
        self.cursor.execute(create_table_command)
        self.sqlite_connection.commit()


    def check_all_tables(self):
        self.cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
        all_tables = self.cursor.fetchall()
        print(all_tables)

    def close_sql_connection(self):
        if self.cursor:
            self.sqlite_connection.close()
            os.remove(self.db_file)

    def get_contig_key(self,contig1,contig2):
        contigs = sorted([contig1, contig2])
        contig_key = f'{contigs[0]}{self.key_split}{contigs[1]}'
        return contig_key

    def get_original_contigs(self,contig_key):
        contig1,contig2=contig_key.split(f'{self.key_split}')
        return contig1,contig2




    ###Saving and retrieving kmer_depth vectors

    def generate_inserts_vectors(self, contig_dict):
        list_keys = list(contig_dict.keys())
        list_values = list(contig_dict.values())
        step=self.insert_step
        for i in range(0, len(list_keys), step):
            yield zip(list_keys[i:i + step], list_values[i:i + step])


    def convert_vector_to_str(self,vector1):
        return '_'.join(str(i) for i in vector1)

    def convert_str_to_vector(self,str1):
        return np.fromstring(str1,sep='_')

    def store_kmer_depth_vectors(self,vector_dict):
        for contig_id in vector_dict:
            vector_dict[contig_id]=self.convert_vector_to_str(vector_dict[contig_id])
        generator_insert = self.generate_inserts_vectors(vector_dict)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO CONTIGS_VECTORS (CONTIG, VECTOR) values (?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def generate_fetch_keys(self, keys):
        keys=list(keys)
        step=self.fetch_step
        for i in range(0, len(keys), step):
            yield str(keys[i:i + step]).strip('[]')

    def fetch_kmer_depth_vectors(self,contigs):
        generator_fetch = self.generate_fetch_keys(contigs)
        res={}
        for chunk in generator_fetch:
            fetch_command = f"SELECT CONTIG,VECTOR FROM CONTIGS_VECTORS WHERE CONTIG in ({chunk})"
            for contig_id,vector in self.cursor.execute(fetch_command).fetchall():
                res[contig_id]=self.convert_str_to_vector(vector)
        return res


    def fetch_kmer_depth_vectors_np(self,contigs):
        generator_fetch = self.generate_fetch_keys(contigs)
        res=[]
        for chunk in generator_fetch:
            fetch_command = f"SELECT VECTOR FROM CONTIGS_VECTORS WHERE CONTIG in ({chunk})"
            for vector in self.cursor.execute(fetch_command).fetchall():
                res.append(self.convert_str_to_vector(vector[0]))
        return np.array(res)


    ###Saving and retrieving contig distance

    def get_all_contig_keys(self):
        res=set()
        fetch_command = f"SELECT CONTIG_KEY FROM CONTIGS_DISTANCE"
        for contig_key in self.cursor.execute(fetch_command).fetchall():
            contig_key=contig_key[0]
            res.add(contig_key)
        return res

    def generate_inserts_distance(self, contigs_keys, kmer_depth_distance,taxa_distance):
        step=self.insert_step
        for i in range(0, len(contigs_keys), step):
            yield zip(contigs_keys[i:i + step], kmer_depth_distance[i:i + step], taxa_distance[i:i + step])

    def store_distance(self,contigs_keys,kmer_depth_distance,taxa_distance):
        generator_insert = self.generate_inserts_distance(contigs_keys,kmer_depth_distance,taxa_distance)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO CONTIGS_DISTANCE (CONTIG_KEY, KMER_DEPTH_DISTANCE, TAXA_DISTANCE) values (?,?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()


    def fetch_distance(self,contigs_keys):
        res={}
        generator_fetch = self.generate_fetch_keys(contigs_keys)
        for chunk in generator_fetch:
            fetch_command = f"SELECT CONTIG_KEY,KMER_DEPTH_DISTANCE,TAXA_DISTANCE FROM CONTIGS_DISTANCE WHERE CONTIG_KEY in ({chunk})"
            for contig_key,kmer_depth_distance,taxa_distance in self.cursor.execute(fetch_command).fetchall():
                res[contig_key]=float(kmer_depth_distance),float(taxa_distance)
        return res


