#!/usr/bin/env python3

import pytest
import sys
sys.path.append('..')

import mutan.mutalign as mutalign
import mutan.mutio as mutio
import mutan.mutAA as mutAA
import mutan.mutlib as mutlib

target = mutio.fasta_to_seq('wtesta_p20m20.fasta').seq
file_path = 'test100.fastq.gz'
mutalign.align_DNA_write_to_csv(input_file=file_path,
                                output_file='test100.csv',
                                target=target,
                                set_num_iterations=True,
                                num_iterations=100)

df=mutAA.AA_analaysis_for_DNA_alignment(target=target,input_file='test100.csv', output_file='test_AA.csv')


correct_alignment_length=583
correct_AA_length=181

def test_library_composition():
    lib_comp = mutlib.library_composition(df, target)
    assert lib_comp == {'Good Mutants': 8,
                      'In-frame Stop Codon': 0,
                      'WT DNA': 0,
                      'Silence Mutants': 0,
                      'Rare Codons': 1}
