#!/usr/bin/env python3

# import sys
# sys.path.append('..')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
from Bio.Seq import Seq
# import mutan.mutAA as mutAA
import mutan.mutAA as mutAA



def effective_AAlib(df, target):
    """
    Filter out the rows in df that represent effective AA mutants satisfying the correct alignment length, the correct AA length and no NT insertions or deletions
    """
    correct_alignment_length=len(target)
    target_AA = mutAA.translate_target_aa(target)
    correct_AA_length=len(target_AA)
    
    df_effective_AAlib=df[(df['alignment_length'] == correct_alignment_length) & (df['AA_length'] == correct_AA_length) &(df['num_insertions']==0)&(df['num_deletions']==0)]

    return df_effective_AAlib


def effective_DNAlib(df, target):
    """
    Filter out the rows in df that represent effective AA mutants satisfying the correct alignment length, the correct AA length and no NT insertions or deletions
    """
    correct_alignment_length=len(target)
    df_effective_DNAlib=df[(df['alignment_length'] == correct_alignment_length) &(df['num_insertions']==0)&(df['num_deletions']==0)]

    return df_effective_DNAlib


def library_composition(df, target):
    """
    Function to calculate the library composition
    
    Input: df -- pandas Dataframe from the AA postprocessing of the alignment .csv
           target -- str or Bio.Seq.Seq, wildtype DNA sequence.
    Output: dict, a dictionary with library compositions, keys: Good Mutants, In-frame Stop Codon, WT DNA, Silence Mutants, Rare Codons
    """
    
    correct_alignment_length=len(target)
    target_AA = mutAA.translate_target_aa(target)
    correct_AA_length=len(target_AA)
    
    

    df_effective_DNAlib=effective_DNAlib(df, target)
    df_effective_AAlib=effective_AAlib(df, target)
    
    num_correct_alignment = len(df_effective_DNAlib)
    num_correct_AA_mutant = len(df_effective_AAlib)
    num_inframe_stop= num_correct_alignment - num_correct_AA_mutant
    num_wt_AA=len(df_effective_AAlib[df_effective_AAlib['num_AA_substitutions'] == 0])
    num_wt_DNA=len(df_effective_AAlib[df_effective_AAlib['num_substitutions'] == 0])
    num_silence_mut = num_wt_AA-num_wt_DNA

    num_rare_codons = (df_effective_AAlib['num_rare_codons']!=0).sum()

    num_good_mutants = num_correct_alignment -num_inframe_stop  - num_wt_AA - num_rare_codons

    library_composition={'Good Mutants':num_good_mutants,
                         'In-frame Stop Codon':num_inframe_stop,
                         'WT DNA': num_wt_DNA,
                         'Silence Mutants': num_silence_mut,
                         'Rare Codons': num_rare_codons}

    return library_composition



def plot_library_composition(lib_comp,  
                             figname=None,
                             colors=['C'+num for num in '30124'],
                             **kwargs):
    """
    Plot library composition into a pie chart
    
    input: lib_comp -- dict, a dictionary of composition key-count pairs
           figname -- str, if want to save figure, the path and name of the figure.
           colors -- list or other iterables that contains names of colors to pass to plt.pie.
    """
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    
    data = list(lib_comp.values())
    data_sum=sum(data)
   
    keys = [text+' '+str(round(val/data_sum*100,1))+'%' for text,val in zip(list(lib_comp.keys()),data)]
    
    

    colors = ['C'+num for num in '30124']

    wedges,texts=ax.pie(data, explode = 0.02*np.arange(len(keys)), colors=colors, **kwargs)

    ax.legend(wedges, keys,
              title="Library Composition",
              loc="center left",
              bbox_to_anchor=(1, 0, 0.5, 1))

    plt.tight_layout()
    
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
    plt.show()
    

    
def plot_mutations_per_read(df, 
                            interval=1, 
                            figname=None, 
                            marker = 'o',
                            markersize = 3,
                            **kwargs):
    """
    Plot function to generate the scattered plot of locations of mutations of each read in the library.
    
    Input:
    df -- Pandas DataFrame that contains the 'AA_substitution_positions' column
    interval -- One read is plotted every 'interval' reads in the dataframe. Default value is 1, where all reads in the dataframe are plotted.
    **kwargs -- Other key word arguments to pass to the plt.plot function.
    
    """

    # all_permutations = np.array([])
    fig, ax = plt.subplots(figsize=(15,10))
    # df = df.sort_values(by = 'Number of AA Substitutions', ascending = False)

    for i in range(0,len(df), interval):
        if df['AA_substitution_positions'].iloc[i] is not None:
            aa_permutation_positions = [int(pos.strip(' ')) for pos in df['AA_substitution_positions'].iloc[i].split(',')]
        
            ax.plot(i*np.ones(len(aa_permutation_positions)), aa_permutation_positions, ls='', marker=marker, markersize = markersize, **kwargs)
        
            ax.set_xlabel('Sequence Index', fontsize = 20)
            ax.set_ylabel('Positions of AA Substitutions', fontsize = 20)
            ax.tick_params(axis='x', labelsize=20)
            ax.tick_params(axis='y', labelsize=20)
            
    plt.tight_layout()
    
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        
    plt.show()

    
def slide_window(i, a, b, length):
    '''
    Helper function to calculate the indices of a [-a,+b] sliding window centered by indice i
    
    input: 
    i -- a digit index in a sequence
    a -- the number of neighbors before i to count in the window
    b -- the number of neighbors after i to count in the window
    length -- int length of the full sequence
    
    output:
    a list of indices, of [-a,+b] nearest neighbors
    '''
    if (i-a >= 0) and (i+b < length):
        return [n for n in range(i-a,i+b+1)]
    elif i-a<0:
        return [n for n in range(0,i+b+1)]
    elif i+b>= length:
        return [n for n in range(i-a,length)]
    
def nt_content(wtseq, ntbefore, ntafter):
    """
    Function to calculate the content of 4 types of NTs in the [-ntbefore,+ntafter] sliding window of a sequence
    
    input:
    wtseq -- str or Bio.Seq.Seq object. The DNA to analyze NT cotent for
    ntbefore -- the number of NTs to include in the sliding window before each NT posiiton
    ntafter -- the number of NTs to include in the sliding window after each NT posiiton
    
    """
    wtseq = str(wtseq).lower()
    seq_list = []
    for i in range(len(wtseq)):
        window_indices = slide_window(i, ntbefore, ntafter, len(wtseq))
        window_str = ''.join([wtseq[index] for index in window_indices])
        nt_window_dict = {}
        nt_window_dict['number of as'] = window_str.count('a')
        nt_window_dict['number of gs'] = window_str.count('g')
        nt_window_dict['number of cs'] = window_str.count('c')
        nt_window_dict['number of ts'] = window_str.count('t')
        seq_list.append(nt_window_dict)
    
    return pd.DataFrame(seq_list)
    

    
def plot_hist_NT_substitutions(df, 
                               figname=None,
                               kde=False, 
                               kde_kws={'gridsize': 100,'bw_method':0.25},
                               **kwargs
                               ):
    
    """
    Plot function to generate the histogram of NT substitutions per read in the library
    
    Input: 
    df -- Pandas DataFrame of alignment containing the 'substitutions' column
    figname -- str. The path and filename of the figure to be saved. Default is None, and no figure is saved
    kde -- bool, if a kernel density estimation is used. parameter to pass to the sns.histplot() function.
    kde_kws -- lib, parameter to pass to the sns.histplot() function to determine parameters for the kernel density estimation.
    **kwargs -- Other keyword arguments to pass to the sns.histplot() function.
    
    Output: 
    None
    """
    if df['num_substitutions'].max() == df['num_substitutions'].min():
        num_substitutions=df['num_substitutions'].max()
        raise ValueError(f'Every effective read has {num_substitutions} NT substitutions. No histogram can be plotted')
    
    fig,ax = plt.subplots()
    sns.histplot(df['num_substitutions'],  kde = kde, binwidth = 1, ax=ax, kde_kws=kde_kws, **kwargs)
    ax.set_xlabel('Number of Nucleotide Substitutions', fontsize = 15)
    ax.set_ylabel('Count', fontsize = 15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    
   
    
    if figname: 
        plt.savefig(figname, dpi=300, bbox_inches='tight')
    plt.show()
    
    return np.histogram(df['num_substitutions'], bins= int(df['num_substitutions'].max()-df['num_substitutions'].min()))

def plot_hist_AA_substitutions(df, 
                               figname=None,
                               kde=False, 
                               kde_kws={'gridsize': 100,'bw_method':0.25},
                               xlim=(0,30),
                               **kwargs):
    """
    Plot function to generate the histogram of NT substitutions per read in the library
    
    Input: 
    df -- Pandas DataFrame of alignment containing the 'substitutions' column
    figname -- str. The path and filename of the figure to be saved. Default is None, and no figure is saved
    kde -- bool, if a kernel density estimation is used. parameter to pass to the sns.histplot() function.
    kde_kws -- lib, parameter to pass to the sns.histplot() function to determine parameters for the kernel density estimation.
    **kwargs -- Other keyword arguments to pass to the sns.histplot() function.
    
    Output: 
    None
    """
    if df['num_AA_substitutions'].max() == df['num_AA_substitutions'].min():
        num_AA_substitutions=df['num_AA_substitutions'].max()
        raise ValueError(f'Every effective read has {num_AA_substitutions} AA substitutions. No histogram can be plotted')
        
    fig, ax = plt.subplots()
    sns.histplot(data = df['num_AA_substitutions'], 
                 binwidth = 1,
                 ax=ax, 
                 kde =kde, 
                 kde_kws=kde_kws,
                 **kwargs)
    
    ax.set_xlabel('Number of AA Substitutions', fontsize = 15)
    ax.set_ylabel('Count', fontsize = 15)
    ax.set_xlim(xlim)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
    return np.histogram(df['num_AA_substitutions'], bins= int(df['num_AA_substitutions'].max()-df['num_AA_substitutions'].min()))


def substitution_loc_bar(df):
    """
    Function to calculate the number of substitutions at each AA position. 
    
    Input:
    df -- pandas DataFrame that contains the 'AA_substitution_positions' column
    
    Output: 
    all_permutations -- np.array, the number of substitutions observed in the library.
    """
    
    all_permutations = np.array([])
    for i in range(len(df)):
        if df['AA_substitution_positions'].iloc[i] is not None:
            aa_permutation_positions = [int(pos.strip(' ')) for pos in df['AA_substitution_positions'].iloc[i].split(',')]
    
            all_permutations = np.hstack([all_permutations, aa_permutation_positions])
    return all_permutations
    

def plot_bar_substitutions_loc(df,
                               figname=None,
                               vline_positions=None,
                               GC_content=False,
                               wtseq=None,
                               ntbefore=None,
                               ntafter=None,
                               **kwargs):
    """
    Plot function to generate the bar plot to show the number of mutations at each AA position.
    
    Input:
    df -- pandas DataFrame that contains the 'AA_substitution_positions' column
    figname -- str, path and file name for the figure to be saved. Default value is none, where plotted figure is not saved.
    vline_positions -- bool, list or np.array, x axes where the user want vertical lines to be drawn. Default value is none, where no vertical lines are plotted
    GC_content -- bool, if true the GC content of the wtseq in a sliding window of [-ntbefore, ntafter] for each nucleotide is drawn at their corresponding AA location.
    wtseq -- str or Bio.Seq.Seq, the wildtype DNA sequnece of the target protein
    ntbefore -- number of nucleotides to count in a slideing window before the current NT index
    ntafter --number of nucleotides to count in a slideing window after the current NT index
    """
    
    
    all_permutations = substitution_loc_bar(df)

    fig, ax = plt.subplots(figsize=(14,3))
    sns.histplot(all_permutations, binwidth=1, ax = ax, **kwargs)
    ax.set_xlabel('Position of AA Substitution', fontsize = 15)
    ax.set_ylabel('Count', fontsize = 15)
    ax.tick_params(axis='x', labelsize=15)
    ax.tick_params(axis='y', labelsize=15)
    ax.set_xlim(0,181)
    
     #####
    if vline_positions is not None:
        ymin, ymax = ax.get_ylim()
        ax.vlines(x=vline_positions, ymin=0, ymax=1.1*ymax, colors='purple', ls='--', lw=2, label='10-AA windows')
        ax.legend()
    
    if GC_content is not False:
        ax1 = plt.twinx(ax)
        df_wtnt = nt_content(wtseq, ntbefore, ntafter)
        
        ax1.plot(np.arange(len(df_wtnt))/3, (df_wtnt['number of gs']+ df_wtnt['number of cs'])/(ntbefore+ntafter), label = f'GC content in {ntbefore+ ntafter}-NT sliding windows', color = 'orange')
        ax1.set_ylabel('GC content', fontsize = 15, rotation=270, labelpad=18)
        ax1.tick_params(axis='y', labelsize=15)
        ax1.legend()

    
    
    plt.legend()

    #####
    
    plt.tight_layout()
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches = 'tight')
    
    plt.show()
    
def descriptive_dict(df, column_name='num_AA_substitutions'):
    """
    Function to print the descriptive statistics of a library.
    
    Input:
    df -- the dataframe that contains post-alignment mutation information of reads in the library.
    column_name -- the name of the column in the dataframe to plot the descriptive statics of. Default value is 'num_AA_substitutions', which is the number of AA substitutions per read. 
    """
    
    return df[column_name].describe().to_dict()



    