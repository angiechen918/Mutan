#!/usr/bin/env python3

import pandas as pd
import gzip
from Bio import SeqIO
from itertools import islice

def fastq_to_pandas(file_path:str):
    """
    Function to decompress a .fastq.gz file and form an initial pandas dataframe. 

    Input: str, file path of the .fastq.gz file to be decompressed
    Output: pandas dataframe with three columns. 
            'sequence_id' column, where each element is a str, 
            'sequence' column, where each element is a Bio.Seq.Seq
            'quality_scores' column, where each element is a list containing the per base numeric quality score

    """
    if file_path[-8:] != 'fastq.gz':
        raise TypeError('The input file should have a .fastq.gz format')

    with gzip.open(file_path, "rt") as handle:
        # Use SeqIO.parse() to read and parse the sequences
        record_list = []
        for record in SeqIO.parse(handle, "fastq"):
            # Access and process each sequence record as needed
            record_lib = {}
            record_lib['sequence_id'] = record.id
            record_lib['sequence'] = record.seq.upper()
            record_lib['quality_scores'] = record.letter_annotations["phred_quality"]
            record_list.append(record_lib)
    
    return pd.DataFrame(record_list)


def fasta_to_seq(file_path:str):
    """
    Function to read in a fasta file and return a sequence. 

    Input: str, file path of the .fasta file 
    Output: Bio.SeqRecord.SeqRecord object of the sequence
    """
    if file_path[-5:] != 'fasta':
        raise TypeError('The input file should have a .fasta format')
    
    return SeqIO.read(file_path,'fasta').upper()



def fastq_extract_seq(file_path:str, n=0):
    """
    Function to read in a fastq file and return the nth sequence in the file. 

    Input: str, file path of the .fasta file 
    Output: Bio.Seq.Seq object of the sequence
    """
    if file_path[-8:] != 'fastq.gz':
        raise TypeError('The input file should have a .fastq.gz format')
    
    with gzip.open(file_path, "rt") as handle:
        query_seqs = SeqIO.parse(handle, "fastq")
        query_seq = next(islice(query_seqs , n, None))
    
    return query_seq.upper()
