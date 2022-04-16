# here we store functions that are used from multiple python
# kernelled ipynbs.

from Bio import Entrez
import os

# important functions for entrez downloading
def authenticate():
    """
    Log into NCBI account to fetch information faster.
    """
    Entrez.email = "FILL_IN"
    Entrez.api_key = "FILL_IN"
    return

def rescue_entrez(entrez_function, **kwargs):
    """
    Avoid HTTPErrors for Entrez requests.

    kwargs can be `db`, `id`, `dbfrom`, ... 
    """
    from http.client import IncompleteRead
    from urllib.error import HTTPError
    while True:
        try:
            handle = entrez_function(**kwargs)
        except (IncompleteRead, HTTPError) as err:
            print(f"We skipped: {err} with {entrez_function} and these arguments: {kwargs}.")
            continue
        break
    return handle

# nice functions for saving files #
def write_csv(filepath, entry_summary_dicts):
    """
    Store summary of all entries of downloading process in csv file.
    """
    import csv

    if check_file_existence(filepath): return False

    fieldnames = sorted(entry_summary_dicts[0].keys())
    with open(filepath, "w") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        [writer.writerow(entry) for entry in entry_summary_dicts]
    print(f"Metadata was successfully stored in {filepath}.")
    return True


def write_fasta(filepath, sequences_dict):
    """
    Store fasta-file of sequences with accession header in `filepath`.
    """
    if check_file_existence(filepath): return False
    with open(filepath, "w") as fastafile:
        [fastafile.write(f">{acc}\n"
                        f"{seq}\n") 
         for acc, seq in sequences_dict.items()]
    print(f"Sequence data was successfully stored in {filepath}.")
    return True


def write_log(filepath, taxon_list, dataentity_name,
              max_seq_per_taxon, seq_filepath, csv_filepath,
              database):
    """
    Store fasta-file of sequences with accession header in `filepath`.
    """
    from datetime import datetime
    NL_TAB = "\n\t"
    if check_file_existence(filepath): return False
    with open(filepath, "w") as logfile:
        logfile.write((
            f"Download of NCBI {database} data.\n"
            f"  Name of the {database}: {dataentity_name}\n"
            f"  Date: {datetime.now()}\n"
            f"  Given taxa as queries were:\n"
            f"\t{NL_TAB.join(taxon_list)}\n"
            f"  A maximum of {max_seq_per_taxon} RefSeq"
            " entries per taxon was downloaded.\n\n"
            f"Sequence data: {seq_filepath}\n"
            f"Metadata: {csv_filepath}\n"
        ))
    print(f"Download log was successfully stored in {filepath}.")
    return True


def check_file_existence(filepath):
    """
    Look up if `filepath` exists and return a bool value
    to prevent overwriting.
    """
    if os.path.exists(filepath):
        print(f"File {filepath} already exists. "
              "It will not be overwritten.")
        return True
    return False
