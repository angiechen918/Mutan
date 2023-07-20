#!/usr/bin/env python3
from Bio import Align
from Bio import SeqIO
import gzip
import pandas as pd
import csv
from Bio import pairwise2

def align_DNA(query, target,
              alignment_method="pairwise2",
              mode='local', 
              match_score=1, 
              mismatch_score=-2, 
              open_gap_score=-10, 
              extend_gap_score=-1) -> dict:
    """
    Perform local alignment of two DNA sequences.
    
    Input: query -- Bio.Seq.Seq object, the sequencing read to map to target sequence
           target -- Bio.Seq.Seq object, the reference sequence.
           alignment_method -- the method used to perform sequence alignment, "pairwise2" or "PairwiseAligner", "pairwise2" is the default
           mode -- parameter to pass to the biopython Bio.Align.PairwiseAligner object.
           match_score, mismatch_score, open_gap_score, extend_gap_score -- parameters to pass to the biopython Bio.Align.PairwiseAligner object.
           
    Output: alignment_lib -- A dictionary object containing the information to describe the DNA level mutations. 
              
    """
    
    # Standardize both sequences to upper cases.
    read_seq = query.upper()
    reference_seq = target.upper()
    
    if alignment_method=="PairwiseAligner":
    ### Alignment using PairwiseAligner
        # Create a PairwiseAligner object.
        aligner = Align.PairwiseAligner(mode='local', 
                                        match_score=match_score, 
                                        mismatch_score=mismatch_score, 
                                        open_gap_score=open_gap_score, 
                                        extend_gap_score=extend_gap_score)


        alignment = list(aligner.align(reference_seq, read_seq, strand='+'))[0]
        alignment_r = list(aligner.align(reference_seq, read_seq, strand='-'))[0]
        strand = 'forward'
        alignment_start=alignment.indices[1][alignment.indices[1] > 0].min()

        if alignment_r.score > alignment.score:
            strand = 'reverse'
            alignment = alignment_r
            read_seq = read_seq.reverse_complement()
            alignment_start = len(read_seq) - alignment.indices[1][alignment.indices[1] > 0].max() -1 ## the -1 is important. Otherwise the alignments on the reverse strands will be shifted by one digit to miss the one bp on the 5' end and gain one bp on the 3' end.

        reference = alignment[0] # str, with '-' for insertion
        query = alignment[1] # str, with '-' for deletion
        alignment_length = len(alignment[0]) # alignment[0] is the target sequence, with '-' for insertsion
    
    ### Alignment using PairwiseAligner done
    
    elif alignment_method=="pairwise2":
    ### Alignment using pairwise2 localms
    
        # perform alignment
        alignments = pairwise2.align.localms(read_seq, 
                                             reference_seq, 
                                             match_score, 
                                             mismatch_score, 
                                             open_gap_score, 
                                             extend_gap_score, 
                                             one_alignment_only=True)

        alignments_r = pairwise2.align.localms(read_seq.reverse_complement(),
                                               reference_seq, 
                                               match_score, 
                                               mismatch_score, 
                                               open_gap_score, 
                                               extend_gap_score, 
                                               one_alignment_only=True)
        strand = 'forward'

        if alignments_r[0][2] > alignments[0][2]:
            alignments = alignments_r
            read_seq = read_seq.reverse_complement()
            strand = 'reverse'

        # Get the alignment with the highest score (assuming unique alignment)
        best_alignment = alignments[0]

        # Get the starting position
        alignment_start = best_alignment[3]
        alignment_end = best_alignment[4]

        # Calculate alignment statistics
        reference = best_alignment[1].strip('-')
        # alignment_length = len(re.findall(r'ATG.*', reference)[0])
        alignment_length = alignment_end - alignment_start
        query = best_alignment[0][alignment_start: alignment_start+alignment_length]

    ### Alignment using pairwise2 localms done

    insertions = []
    deletions = []
    substitutions = []

    for j in range(len(query)):
        if query[j] == '-':
            deletions.append(j)
        elif reference[j] == '-':
            insertions.append(j)
        elif (query[j] != reference[j]):
            substitutions.append(reference[j]+str(j)+query[j])

    alignment_lib={"sequence": read_seq[alignment_start:alignment_start+alignment_length-len(deletions)],
                   "read_length": len(read_seq),
                   "alignment_start": alignment_start,
                   "alignment_length": alignment_length,
                   "insertions": insertions,
                   "num_insertions": len(insertions),
                   "deletions":deletions,
                   "num_deletions":len(deletions),
                   "substitutions": substitutions,
                   "num_substitutions":len(substitutions),
                   "strand":strand}
    
    return alignment_lib



def align_DNA_write_to_csv(input_file, 
                           output_file, 
                           target, 
                           alignment_method="pairwise2",
                           mode='local',
                           match_score=1, 
                           mismatch_score=-2, 
                           open_gap_score=-10, 
                           extend_gap_score=-1,
                           set_num_iterations=True,
                           num_iterations=500):
    """
    Function that reads in sequences from a .fastq.gz file, perform alignment and save each entry to .csv at a time. 
    
    input: input_file -- str, input file path
           output_file -- str, output file path
           target -- Bio.Seq.Seq object, the reference sequence, to be passed into mutalign.align_DNA()
           alignment_method -- the method used to perform sequence alignment, "pairwise2" or "PairwiseAligner", "pairwise2" is the default
           mode -- parameter to pass to the biopython Bio.Align.PairwiseAligner object, , to be passed into mutalign.align_DNA()
           match_score, mismatch_score, open_gap_score, extend_gap_score -- parameters to pass to the biopython Bio.Align.PairwiseAligner object, to be passed into mutalign.align_DNA()
           set_num_iterations -- bool, whether or not the user intends to set an upper limit on the number of entries to process, default is True. When set to false, the function processes the entire fastq.gz set.
           num_iterations -- int, to be set if set_num_iterations is True. Default is 500
    
    """
    with open(output_file, mode='w', newline='') as csv_file:
        
        # if user not intended to set the number of entries to align, the default is the iterate through the entire dataset
        if not set_num_iterations:
            with gzip.open(input_file, "rt") as handle:
                num_iterations = sum([1 for _ in SeqIO.parse(handle, 'fastq')])
                
        ct=0        
        with gzip.open(input_file, "rt") as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                alignment_result = align_DNA(record.seq, 
                                             target=target,
                                             alignment_method=alignment_method,
                                             mode="local", 
                                             match_score=match_score, 
                                             mismatch_score=mismatch_score, 
                                             open_gap_score=open_gap_score, 
                                             extend_gap_score=extend_gap_score)  # Assuming `target_sequence` is defined
                
                alignment_start=alignment_result["alignment_start"]
                alignment_length=alignment_result["alignment_length"]
                num_deletions=alignment_result["num_deletions"]
                # read_seq[alignment_start:alignment_start+alignment_length-num_deletions]
                
                alignment_result.update({'sequence_id':record.id})
                alignment_result.update({'quality':record.letter_annotations["phred_quality"][alignment_start:alignment_start+alignment_length-num_deletions]})
                
                if ct == 0: # first iteration writes the header
                    writer = csv.DictWriter(csv_file, fieldnames=alignment_result.keys())
                    writer.writeheader() # write the header for the output .csv file

                writer.writerow(alignment_result) # write one row per iteration
                
                ct+=1 # iterate until the set limit
                
                if ct%100 == 0:
                    print(f'{ct} iterations done')
                if ct >= num_iterations:
                    break



            
            

def align_all_DNA(file_path, target, exclude_short=False):
    """
    Extract every seq in a fastq.gz file, align them to a target seq and calculate 
    indels and substitutions. Save results to a pandas DataFrame
    
    input: filepath -- str, path of the fastq.gz file
           target -- Bio.Seq.Seq object to be used as the target for alignment
           exclude_short -- int if not False, the minimum number of bps in the read to allow the read to be included to the output df
    """
    if file_path[-8:] != 'fastq.gz':
        raise TypeError('The input file should have a .fastq.gz format')
        
    with gzip.open(file_path, "rt") as handle:
        # Use SeqIO.parse() to read and parse the sequences
        record_list = []
        for record in SeqIO.parse(handle, "fastq"):
            record_lib={}
            record_lib['sequence_id'] = record.id
            record_lib['sequence'] = record.seq
            # when exclude_short is specified, short reads are excluded from the df.
            if (exclude_short) and (len(record.seq)<exclude_short):
                continue
            # when exclude_short is not specified, or if the read is long, info included in the output df
            else:
                record_lib['quality_scores'] = record.letter_annotations["phred_quality"]
                record_lib.update(align_DNA(record.seq, target))
                record_list.append(record_lib)

    return pd.DataFrame(record_list)



