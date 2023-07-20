#!/usr/bin/env python3
import sys
sys.path.append('.')
import mutan.mutAA as mutAA

from Bio.Seq import Seq
import pandas as pd
import ast


def translate_target_aa(target):
    """
    Translates a DNA sequence to its AA sequence from the first start codon.
    
    Input:
    target -- str or Bio.Seq.Seq object, the WT DNA sequnece.
    
    Output:
    reference_AA -- Bio.Seq.Seq object, the WT AA sequnece.
    """
    reference_seq=Seq(target).upper()
    start_index = reference_seq.find('ATG')
    reference_AA = reference_seq[start_index:].translate(to_stop=True)
    
    return reference_AA

def translate_to_aa(row):
    """
    Alignment table processing function to calculate the AA translation of each aligned DNA read.
    
    Input:
    row -- an indexed row in a pandas dataframe that contains the 'sequence' column and the 'alignment_start' column, which is the string representation or Bio.Seq.Seq representation of an aligned read.
    
    Output:
    orf -- str, the AA translation of the aligned sequence.
    """
    sequence= row['sequence']
    seq = Seq(sequence).upper() 
    # start = row['alignment_start'] # To reduce the size of pandas dataframe, modified the alignment .csv file for the stored sequence to start only from the aligned region.
    start = 0
    start_codon = 'ATG'
    start_index = seq.find(start_codon, start = int(start))
    
    
    if start_index != -1:
        orf = seq[start_index:].translate(to_stop=True)
        return str(orf)
    else:
        return ''




def calculate_AA_substitutions(query_AA, reference_seq):
    """
    Alignment table processing function to calculate the AA substitutions of each aligned DNA read.
    
    Input:
    query_AA: str or Bio.Seq.Seq, single letter representation of the AA translation of an aligned read
    reference_seq: str or Bio.Seq.Seq, reference DNA sequence for each query_AA to compare to and calculate AA substitutions
    
    Output:
    str, the X000Y format of mutation representation separated by ', '. X is the amino acid at the 000 position of the WT amino acid, which is mutated to Y in the query AA sequence. 
    num_permutations_AA -- number of AA substitutions in a read
    """

    ## Calculate number of mutations:
    reference_AA = mutAA.translate_target_aa(reference_seq)

    #df_with_translation['AA_permutations'] = None
    permutations_AA=[]
    num_permutations_AA = []

    
    query_AA=str(query_AA)
    for i in range(min(len(reference_AA), len(query_AA))):
        if query_AA[i] != reference_AA[i]:
            permutations_AA.append(reference_AA[i]+str(i)+query_AA[i])
    num_permutations_AA=len(permutations_AA)
    
    return ', '.join(map(str, permutations_AA)), num_permutations_AA


# Document positions of AA substitutions
def mutation_position(permutations):
    """
    Alignment table processing function to calculate the positions of AA substitutions of each aligned DNA read.
    
    Input:
    permutations: str, the X000Y format of mutation representation separated by ', '. X is the amino acid at the 000 position of the WT amino acid, which is mutated to Y in the query AA sequence. 
    
    Output:
    str, the X000Y format of mutation representation separated by ', '. X is the amino acid at the 000 position of the WT amino acid, which is mutated to Y in the query AA sequence. 
    """
    
    if len(permutations) > 0:
        mut_list = permutations.split(',')
        pos_list = [int(mut.strip(' ')[1:-1]) for mut in mut_list]

        return ', '.join(map(str, pos_list))
    
    else: 
        return None
    
def count_rare_codons(row, rare_codons=['AGG', 'AGA', 'CTA', 'ATA', 'CGA', 'CCC']):
    """
    Alignment table processing function to calculate the number of rare codons in each aligned DNA read.
    
    Input:
    row -- an indexed row in a pandas dataframe that contains the 'sequence' column and the 'alignment_start' column, which is the string representation or Bio.Seq.Seq representation of an aligned read.
    rare_codons -- list of strs for codons to count.
    
    Output:
    int -- number of rare codons found in each sequence.
    
    """
    
    sequence= row['sequence']
    seq = Seq(sequence).upper()
    start = row['alignment_start']
    start_codon = 'ATG'
    start_index = seq.find(start_codon, start = int(start))
    
    if start_index != -1:
        found_codons = []
        for i in range(start_index, start_index + int(row['AA_length'])*3, 3):
            codon = str(sequence[i:i+3])
            if codon in rare_codons:
                found_codons.append(codon)
        return len(found_codons)
    else:
        return 0
    
def report_rare_codons(row, rare_codons=['AGG', 'AGA', 'CTA', 'ATA', 'CGA', 'CCC']):
    """
    Helper function to show the rare codons in an aligned DNA read.
    
    Input:
    row -- an indexed row in a pandas dataframe that contains the 'sequence' column and the 'alignment_start' column, which is the string representation or Bio.Seq.Seq representation of an aligned read.
    rare_codons -- list of strs for codons to report.
    
    Output:
    found_codons -- list of strs, str representations of rare codons found in the query sequence.
    locs -- list of ints, starting NT locations of rare codons found in the query sequence.
    
    """
    
    sequence= row['sequence']
    seq = Seq(sequence).upper()
    start = row['alignment_start']
    start_codon = 'ATG'
    start_index = seq.find(start_codon, start = int(start))
    
    if start_index != -1:
        found_codons = []
        locs=[]
        for i in range(start_index, start_index + int(row['AA_length'])*3, 3):
            codon = str(sequence[i:i+3])
            if codon in rare_codons:
                found_codons.append(codon)
                locs.append(i)
        return found_codons, locs
    else:
        return 0
    

def AA_analaysis_for_DNA_alignment(target, 
                                   df=None,
                                   output_file=None, 
                                   input_file=None, 
                                   rare_codons=['AGG', 'AGA', 'CTA', 'ATA', 'CGA', 'CCC']):
    """
    Function to post process the a pandas dataframe or .csv DNA alignment table and append columns of mutation information to the original table.
    
    Input: 
    target -- str or Bio.Seq.Seq object. WT DNA sequence for the library to align to.
    df -- pandas dataframe that contains the alignment information of reads in the library.
    output_file -- path and file name of the output csv to be saved.
    input_file -- path and file name of the input alignment csv to read from. If speficied, input_file overwrites the df input
    rare_codons -- list of strs, definition of rare codons in the current analysis.
    """
    if input_file is not None:
        df=pd.read_csv(input_file)
        
    # Apply the translation function to the 'Sequence' column and store the result in a new column 'AA'
    df.loc[:,'AA_translation'] = df.apply(translate_to_aa, axis=1)
    df.loc[:,'AA_length'] =  df['AA_translation'].str.len()
    df.loc[:,'num_rare_codons'] = df.apply(lambda x: count_rare_codons(x, rare_codons), axis=1)

    df.loc[:,'AA_substitutions'] = df['AA_translation'].apply(lambda x: calculate_AA_substitutions(x, target)[0])
    df.loc[:,'num_AA_substitutions'] =  df['AA_translation'].apply(lambda x: calculate_AA_substitutions(x, target)[1])


    df['substitutions'] = df['substitutions'].apply(lambda x: ', '.join(map(str, ast.literal_eval(x))))
    df['insertions'] = df['insertions'].apply(lambda x: ', '.join(map(str, ast.literal_eval(x))))
    df['deletions'] = df['deletions'].apply(lambda x: ', '.join(map(str, ast.literal_eval(x))))
    df['AA_substitution_positions'] = df['AA_substitutions'].apply(mutation_position)
    
    
    
    if output_file is not None: 
        df.to_csv(output_file)
    
    return df


