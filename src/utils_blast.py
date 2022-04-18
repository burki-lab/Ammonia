# here we implement lines that help to document the fetching of 
# genomes as blast results.

from Bio import Entrez
from Bio import SeqIO
import csv
from datetime import datetime
import json
import os
import pandas as pd

from __init__ import dbex_utils

class BufferedGenomeLoader:
    '''Class that helps to smoothely download and track reference genomes from nuccore db.
    '''
    def __init__(self, csv_filepath, extend=False, force=False):
        '''Initialize a genome reference loader by a given path for the documenting csv file.
        '''
        self.FN = [
            "accession_id",
            "gb_taxid",
            "filepath",
            "date_fetched",
            "title",
            "genome_type",  # e.g. plastid
            "source_type",  # e.g. refseq
            "topology",  # e.g. circular,
            "organism",  # redundant but good for reader
            "length",  # in bp
            "last_update_on_db"  # could help with versions
        ]
        
        if os.path.exists(csv_filepath) and not extend and not force:
            raise RuntimeError(f"File {csv_filepath} already exists. Turn extend on.")
            
        self.filepath = csv_filepath
        
        if extend and force : raise RuntimeError("Conflict of 'force' (overwrite) and 'extend'.")
        
        if extend : self.flag = "a"
        else : self.flag = "w"
        
        self.accession_list = []
        return
    
    def build_writer(self, temp=False):
        '''Initialize the csv writer.
        '''
        if temp:
            self.tempfilepath = f"{self.filepath}.tmp"
            file = self.tempfile = open(self.tempfilepath, "w")
            write_header = True
        else:
            # we want a header for all new files.
            write_header = not os.path.exists(self.filepath) or self.flag=="w"
            file = self.file = open(self.filepath, self.flag)
            
        self.writer = csv.DictWriter(file, fieldnames=self.FN)
        if write_header : self.writer.writeheader()
        return
    
    def close(self, temp=False):
        '''Terminate writing process.
        '''
        if temp:
            self.tempfile.close()
        else:
            self.file.close()
            if hasattr(self, "tempfilepath") : os.remove(self.tempfilepath)
            # we change the flag to append the file if another round of genomes should be flushed one day...
            self.flag="a"
        del self.writer
        return
    
    def preload_genome(self, accession_id, filepath, genome_summary=None):
        '''Store info of a genome given its id and store metainfo to later fetch it.
        '''
        if not hasattr(self, "writer") : self.build_writer(temp=True)
        # we want to collectas much info as possible, so we need to catch some KeyErrors
        def try_key(dct, key):
            try : return dct[key]
            except KeyError : return None
        
        # collect the accession id and load summary if not given
        if accession_id in self.accession_list : return
        
        if genome_summary is None : genome_summary = get_genome_smry(accession_id)
        
        self.accession_list.append(accession_id)
        
        
        # assemble the info
        values = [
            accession_id,
            try_key(genome_summary, "taxid"),
            os.path.abspath(filepath),
            "PLACEHOLDER",
            try_key(genome_summary, "title"),
            try_key(genome_summary, "genome"),
            try_key(genome_summary, "sourcedb"),
            try_key(genome_summary, "topology"),
            try_key(genome_summary, "organism"),
            try_key(genome_summary, "slen"),
            try_key(genome_summary, "updatedate")
        ]
        
        csv_row = dict(zip(self.FN, values))
        self.writer.writerow(csv_row)
        return
    
    def flush(self, ret_type="fasta", db="nuccore"):
        '''Finally download all the genomes that were
        '''
        # terminate the genome collection
        self.close(temp=True)
        self.temp_df = pd.read_csv(self.tempfilepath, index_col=False)
        
        # at the same time we open the real csv writer
        self.build_writer()
        
        # download the accessions in batches of 1000
        while len(self.accession_list) > 0:
            batch = self.accession_list[:1000]
            self.accession_list = self.accession_list[1000:]
            query_list = ",".join(batch)
            kwargs = {"db": db,
                      "id": query_list,
                      "rettype": ret_type,
                      "retmode": "text"
                     }
            handle = dbex_utils.rescue_entrez(Entrez.efetch, **kwargs)
            fetching_time = str(datetime.today())
            # store the records
            self.save_and_log_records(handle, fetching_time, ret_type=ret_type)
        
        self.close()
        return
                
                
    def save_and_log_records(self, record_handle, fetching_time, ret_type="fasta"):
        '''Store the records given in a handle.
        '''
        from Bio import SeqIO
        for record in SeqIO.parse(record_handle, format=ret_type):
            rec_info = self.get_record_info(record)
            self.save_record(record, rec_info, ret_type=ret_type)
            self.log_record(record, rec_info, fetching_time)
        return
    
    
    def get_record_info(self, record):
        '''Obtain the information line of a given record.
        '''
        record_id = record.id.split(".")[0]
        rec_info_df = self.temp_df[self.temp_df["accession_id"]==record_id]
        return dict(rec_info_df.iloc[0, :])
    
    def save_record(self, record, rec_info_df, ret_type="fasta"):
        '''Store a record as fasta file.
        '''
        if ret_type=="fasta":
            with open(rec_info_df["filepath"], "w") as ffile:
                ffile.write(f">{rec_info_df['accession_id']}\n")
                ffile.write(f"{record.seq}\n")
        elif ret_type=="gb":
            SeqIO.write(record, rec_info_df["filepath"], "gb")
        return
    
    def log_record(self, record, rec_info_df, fetching_time):
        '''Document the record download.
        '''
        # fetching date should be updated
        csv_row = rec_info_df
        csv_row["date_fetched"] = fetching_time
        self.writer.writerow(csv_row)
        return
# end BufferedGenomeLoader


# helpers
def edo(query, method="esummary", db_from=None):
    '''Run esummary or elink on IDs.
    '''
    from __init__ import dbex_utils
    
    dbex_utils.authenticate()
    # we assemble a bunch of good kwargs
    kwargs = {"db": "nuccore",
          "id": query,
          "retmax": "10",
          "sort": "pdat",
          "rettype": "gb",
          "retmode": "json"
         }
    
    # 2 methods are enough for us
    if method == "esummary" : fct = Entrez.esummary
    elif method == "elink": 
        if db_from is None : kwargs["dbfrom"] = "protein"
        else : kwargs["dbfrom"] = db_from
        fct = Entrez.elink
    else : raise RuntimeError("Wrong method. Try elink or esummary.")
    
    handle = dbex_utils.rescue_entrez(fct, **kwargs)
    try : return json.load(handle)
    except json.JSONDecodeError: return None

def get_genome_smry(accession_id):
    '''Load the summary of a genome given its id.
    '''
    result = edo(accession_id, method="esummary")
    genome_summary = result["result"][result["result"]["uids"][0]]
    return genome_summary