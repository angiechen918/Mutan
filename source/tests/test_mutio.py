#!/usr/bin/env python3

import pytest
import sys
sys.path.append('..')


import mutan.mutio as mutio

test_fastqgz_file = 'test.fastq.gz'
test_fasta_file = 'wtesta.fasta'

def test_fastq_to_pandas():
    """
    test function of fastq_to_pandas
    """
    
    df = mutio.fastq_to_pandas(test_fastqgz_file)
    row0 = df.iloc[0]
    
    assert df.columns[0] == 'sequence_id'
    assert df.columns[1] == 'sequence'
    assert df.columns[2] == 'quality_scores'
    assert type(row0['sequence_id']) == str
    assert type(row0['quality_scores']) == list
    assert type(row0['quality_scores'][0]) == int
    with pytest.raises(TypeError):
        mutio.fastq_to_pandas('falsetype.file')


def test_fasta_to_seq():
    """
    test function of fasta_to_seq
    """
    test_seq = mutio.fasta_to_seq(test_fasta_file).seq
    
    assert any([l in 'atgcksn'.upper() for l in test_seq ])
    with pytest.raises(TypeError):
        mutio.fasta_to_seq('falsetype.file')
        

def test_fastq_extract_seq():
    """
    test function of fastq_extract_seq
    """
    test_single = mutio.fastq_extract_seq(test_fastqgz_file,n=0)
    
    assert test_single.id == 'm54016_230501_202025/4194562/ccs'
    with pytest.raises(TypeError):
        mutio.fasta_to_seq('falsetype.file')


