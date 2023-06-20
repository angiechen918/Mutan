#!/usr/bin/env python3

import pytest
import sys
sys.path.append('..')

import mutan.mutalign as mutalign
from Bio import SeqIO
from Bio import Seq
import pandas as pd

# Read the reference sequence
reference_file = "wtesta_p20m20.fasta"
reference_seq = SeqIO.read(reference_file, "fasta").seq

def test_align_DNA():
    
    # Create mismatches in the reference_seq
    read_seq = Seq.Seq(str(reference_seq)[:30]+'a'+str(reference_seq)[30:])
    read_seq = Seq.Seq(str(read_seq)[:54]+str(read_seq)[55:])
    lib = mutalign.align_DNA(read_seq, reference_seq)
    assert lib['deletions'] == [54] 
    assert lib['insertions'] == [29] 
    assert lib['substitutions'] == []
    
    
def test_align_DNA_write_to_csv():
    
    mutalign.align_DNA_write_to_csv('test.fastq.gz', 
                           'test.csv',
                           reference_seq,
                           set_num_iterations=False)
    
    test_df = pd.read_csv('test.csv')
    assert len(test_df) == 10
    
    mutalign.align_DNA_write_to_csv('test.fastq.gz', 
                           'test.csv',
                           reference_seq,
                           set_num_iterations=True,
                           num_iterations=5)
    
    test_df = pd.read_csv('test.csv')
    assert len(test_df) == 5
    