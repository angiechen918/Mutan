#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from Bio.Seq import Seq
import mutan.mutAA as mutAA

def sort_AA_on_target_abundance(target):
    """
    helper function to sort AAs based on the abundance in target seq
    
    input: 
    target -- Bio.Seq.Seq object, target sequence to be sorted
    output:
    AAs_sorted -- str, 20 AAs sorted by abundance in target seq, in descending order
    reference_AA_ct -- 1d numpy array, the count of 20 AAs in the target sequence
    """
    
    reference_AA = mutAA.translate_target_aa(target)
    
    AAs='ARNDEQGHILKMFPSTWYVC'
    reference_AA_ct = np.zeros(len(AAs))
    for i,AA in enumerate(AAs):
        for j in range(len(reference_AA)):
            if str(reference_AA)[j] == AA:
                reference_AA_ct[i]+=1

    AAs_sort_index = np.argsort(reference_AA_ct)[::-1]
    reference_AA_ct = np.sort(reference_AA_ct)[::-1]
    AAs_sorted = ''.join([AAs[index] for index in AAs_sort_index])
    return AAs_sorted, reference_AA_ct

def matrix_mutto_loc(df, 
                     target, 
                     AA_sorted_on_wt_abundance=True, 
                     sorted_AAs='ARNDEQGHILKMFPSTWYVC'):
    """
    Function to calculate the mutate-to AA at each location of the peptide sequence.
    
    Input:
    df -- pandas DataFrame. AA analysis DataFrame containing the 'AA_substitutions' column
    target -- Bio.Seq.Seq object. WT target sequence 
    AA_sorted_on_wt_abundance -- Bool. Default value True. Whether choose to order the columns in the output DataFrame by the abundance of AAs in the WT sequence. If set to be False AAs will be sorted by the following order: 'ARNDEQGHILKMFPSTWYVC'. The AA_sorted_on_wt_abundance parameter overwrites the sorted_AAs parameter in this funciton.
    sorted_AAs -- If AA_sorted_on_wt_abundance is set to be False, provide a customized order of 20 AAs as columns in the output DataFrame. If AA_sorted_on_wt_abundance is set to be False and not specified with a customized value the default order is 'ARNDEQGHILKMFPSTWYVC'.
    
    Output: 
    df_matrix_mutto_loc -- pandas DataFrame with int dtype. Each column is one of 20 AAs. The columns are sorted by the abundance of AAs in the WT sequence or by an order specified by the user. Each column is a AA position in the WT sequence. The integer numbers in the dataframe specifies the number of mutations found in the library to the AA in the column index at the position of the row index.
    
    """
    

    AAs=sorted_AAs

    if AA_sorted_on_wt_abundance:
        AAs = sort_AA_on_target_abundance(target)[0]
    
    reference_AA = mutAA.translate_target_aa(target)
    
    # Initialize AA_list
    AA_list=[]
    for loc in range(len(reference_AA)):
        AA_dict = {}
        for AA in AAs:
            AA_dict[AA] = 0
        AA_list.append(AA_dict)
        


    for i in range(len(df)):
        if df['AA_substitutions'].iloc[i]!= '':
            mut_to_list = [(int(mut.strip(' ')[1:-1]), mut.strip(' ')[-1]) for mut in df['AA_substitutions'].iloc[i].split(',')]
            
            
            for _ in mut_to_list:
                AA_list[_[0]][_[1]]+=1
                
    df_matrix_mutto_loc = pd.DataFrame(AA_list)
                

    return df_matrix_mutto_loc


def plot_mut_loc_matrix(df_matrix_mutto_loc,
                        figname=None,
                        **kwargs):
    """
    Plot function to generate the mutation-location heatmap.
    
    Input: 
    df_matrix_mutto_loc -- pandas DataFrame calculated using the matrix_mutto function, where each column is one of the 20 AAs and each row is one of the AA positions in the WT seq.
    figname -- str. The path and filename to save the plotted figure. Default value is None, and the plotted function is not saved. 
    **kwargs -- keyword arguments to be passed to the ax.imshow function to specify plotting parameters.
    
    Output:
    None
    """

    fig,ax = plt.subplots(figsize = (20,5))
    
    ###########
    data = df_matrix_mutto_loc.to_numpy().transpose()

    value = 0

    masked_array = np.ma.masked_where(data == value, data)
    # Get the Spectral_r colormap
    spectral_r_cmap = mpl.cm.get_cmap("Spectral_r")

    # Create a new colormap based on the Spectral_r colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Spectral_r_copy",
                                                        spectral_r_cmap(np.linspace(0, 1, 256)))
    cmap.set_bad(color='black')
    
    gcf = ax.imshow(masked_array, cmap =cmap, **kwargs)
    cbar=plt.colorbar(gcf)
    cbar.ax.tick_params(labelsize=20)
    ###########

    #define y-unit to x-unit ratio
    ratio = 0.2

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

    #set tick labels
    ax.set_yticks(np.arange(0,20))
    ax.set_yticklabels([AA for AA in df_matrix_mutto_loc.columns.to_list()])
    
    
    #set plot labels
    ax.set_xlabel('AA position',fontsize=20)
    ax.set_ylabel('Mutate to',fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=13)
    ax.set_xlim(0,180)

    plt.tight_layout()
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        
        
def matrix_mutto_AA(df_matrix_mutto_loc, target):
    """
    Function to calculate the number of mutations from one of the 20 AAs to 19 other AAs in the library.
    
    Input:
    df_matrix_mutto_loc -- pandas DataFrame, calculated using the matrix_mutto function, where each column is one of the 20 AAs and each row is one of the AA positions in the WT seq
    target -- Bio.Seq.Seq 
    
    Output: df_matrix_mutto_AA -- Pandas DataFrame, each column is an AA sorted by the abundance of AAs in the WT seq or by a user-specified order. Each column represents one mutate-to AA. Rows are AAs in the same order as columns, but represent mutate-from AAs.  
    """
    
    reference_AA = mutAA.translate_target_aa(target)
    AA_to_AA = np.zeros([20,20])

    for i,AA in enumerate(df_matrix_mutto_loc.columns.to_list()):
        for j in range(len(reference_AA)):
            if reference_AA[j] == AA:
                AA_to_AA[i] += df_matrix_mutto_loc.to_numpy()[j]
                
    df_matrix_mutto_AA = pd.DataFrame(AA_to_AA,columns=df_matrix_mutto_loc.columns.to_list())
                
    return df_matrix_mutto_AA 



def plot_matrix_mutto_AA(df_matrix_mutto_AA, 
                         figname=None,
                         title=None, 
                         theoretical_mutto_AA=None,
                         target=None,
                         **kwargs):
    """
    Plot function to generate the AA-to-AA heatmap.
    
    Input:
    df_matrix_mutto_AA -- Pandas DataFrame calculated from the matrix_mutto_AA function, each column is an AA sorted by the abundance of AAs in the WT seq or by a user-specified order. Each column represents one mutate-to AA. Rows are AAs in the same order as columns, but represent mutate-from AAs. 
    title -- title of the figure to be printed.
    theoretical_mutto_AA -- pandas dataframe or numpy array to highlight theoretical_mutto_AA by single NT mutations using epPCR. Can also be specified as 'epPCR' for automatic calculation of highlights. Default value is None.
    target -- str or Bio.Seq.Seq. If theoretical_mutto_AA is specified as 'epPCR', the wt DNA seq has to be passed for automatic theoretical epPCR calculation.
    **kwargs -- keyword arguments to be passed to the ax.imshow function to specify plotting parameters.
    """
    
    fig,ax = plt.subplots(figsize = (10,8))
    
    ########### Set zero values to 0
    data = df_matrix_mutto_AA.to_numpy().transpose()

    value = 0

    masked_array = np.ma.masked_where(data == value, data)
    # Get the Spectral_r colormap
    spectral_r_cmap = mpl.cm.get_cmap("Spectral_r")

    # Create a new colormap based on the Spectral_r colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Spectral_r_copy",
                                                        spectral_r_cmap(np.linspace(0, 1, 256)))
    
    cmap.set_bad(color='black')
    #############


    gcf = ax.imshow(masked_array, cmap =cmap, **kwargs)

    cbar = plt.colorbar(gcf)
    cbar.set_label('Count', rotation=270,fontsize=15)
    cbar.ax.tick_params(labelsize=15)

    #define y-unit to x-unit ratio
    ratio = 1

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

    #set tick labels
    ax.set_yticks(np.arange(0,20))
    ax.set_yticklabels([AA for AA in df_matrix_mutto_AA.columns.to_list()],fontsize=15);

    ax.set_xticks(np.arange(0,20))
    ax.set_xticklabels([AA for AA in df_matrix_mutto_AA.columns.to_list()],fontsize=15);

    for i in range(20):
        for j in range(20):
            text = ax.text(i,j, round(df_matrix_mutto_AA.to_numpy()[i, j]),
                           ha="center", 
                           va="center", 
                           color="k",
                           fontsize=10)


    ax.set_ylabel('Mutate to',fontsize=20)
    ax.set_xlabel('Mutate from',fontsize=20)

    ax.set_title(title)
    
    if theoretical_mutto_AA is not None:
        if theoretical_mutto_AA == 'epPCR':
            if target is None:
                raise TypeError("Default value for target seq is None. For automatic calculation of epPCR mutations, a target seq in the str type or the Bio.Seq.Seq type has to be passed.")
            theoretical_mutto_AA = theoretical_epPCR_matrices_mutto(target)[1]
            
            
        if type(theoretical_mutto_AA) is pd.core.frame.DataFrame:
            theoretical_mutto_AA=theoretical_mutto_AA.to_numpy().astype(int)
            
        if type(theoretical_mutto_AA) is np.ndarray:
            highlight_theoretical_mutto_AA(theoretical_mutto_AA,ax=ax)
        
        else:
            raise TypeError("theoretical_mutto_AA should have a type of pd.DataFrame or np.ndarray, or should be specified as 'epPCR'")
            
    plt.tight_layout()
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        


def plot_matrix_mutto_AA_norm(df_matrix_mutto_AA, 
                         target, 
                         figname=None,
                         title=None, 
                         theoretical_mutto_AA=None,
                         **kwargs):
    """
    Function to plot the AA-to-AA mutation matrix, each location normalized by the count of AAs in the wt sequence.
    
    Input: 
    df_matrix_mutto_AA -- Pandas DataFrame calculated from the matrix_mutto_AA function, each column is an AA sorted by the abundance of AAs in the WT seq or by a user-specified order. Each column represents one mutate-to AA. Rows are AAs in the same order as columns, but represent mutate-from AAs. 
    target -- Target seq to normalize to. Note: this function can only be used correctly if the columns in the df_matrix_mutto_AA are sorted by abundance of AAs in the WT sequence.
    title -- title of the figure to be printed.
    theoretical_mutto_AA: pandas dataframe or numpy array to highlight theoretical_mutto_AA by single NT mutations using epPCR. Can also be specified as 'epPCR' for automatic calculation of highlights. Default value is None.
    **kwargs -- keyword arguments to be passed to the ax.imshow function to specify plotting parameters.
    """
    
    
    fig,ax = plt.subplots(figsize = (10,8))

    value = 0
    AA_to_AA = df_matrix_mutto_AA.to_numpy()
    reference_AA_ct = sort_AA_on_target_abundance(target)[1]
    data=(AA_to_AA[:-1,:]/(reference_AA_ct[:-1]).reshape(-1,1)).transpose()
    
    masked_array = np.ma.masked_where(data == value, data)
    # Get the Spectral_r colormap
    spectral_r_cmap = mpl.cm.get_cmap("Spectral_r")

    # Create a new colormap based on the Spectral_r colormap
    cmap = mpl.colors.LinearSegmentedColormap.from_list("Spectral_r_copy",
                                                        spectral_r_cmap(np.linspace(0, 1, 256)))
    
    cmap.set_bad(color='black')

    gcf = ax.imshow(masked_array, cmap =cmap, **kwargs) ## plot 'mutate from' on the x axis, 'mutate to' on the y axis
    cbar = plt.colorbar(gcf)
    cbar.set_label('Count per WT AA', rotation=270, fontsize=15, labelpad=18)
    cbar.ax.tick_params(labelsize=15)

    #define y-unit to x-unit ratio
    ratio = 1

    #get x and y limits
    x_left, x_right = ax.get_xlim()
    y_low, y_high = ax.get_ylim()

    #set aspect ratio
    ax.set_aspect(abs((x_right-x_left)/(y_low-y_high))*ratio)

    
    #set tick labels
    ax.set_yticks(np.arange(0,20))
    ax.set_yticklabels([AA for AA in df_matrix_mutto_AA.columns.to_list()],fontsize=15);

    ax.set_xticks(np.arange(0,19))
    ax.set_xticklabels([AA for AA in df_matrix_mutto_AA.columns.to_list()][:-1],fontsize=15);

    ax.set_ylabel('Mutate to',fontsize=20)
    ax.set_xlabel('Mutate from',fontsize=20)

    # Loop over data dimensions and create text annotations.
    for i in range(19):
        for j in range(20):
            text = ax.text(i,j, round(((AA_to_AA[:-1,:]/(reference_AA_ct[:-1].reshape(-1,1)))[i,j])),
                           ha="center", 
                           va="center", 
                           color="k",
                           fontsize=10)
    ax.set_title(title)
            
    if theoretical_mutto_AA is not None:
        if theoretical_mutto_AA == 'epPCR':
            theoretical_mutto_AA = theoretical_epPCR_matrices_mutto(target)[1]
            
            
        if type(theoretical_mutto_AA) is pd.core.frame.DataFrame:
            theoretical_mutto_AA=theoretical_mutto_AA.to_numpy().astype(int)
            
        if type(theoretical_mutto_AA) is np.ndarray:
            highlight_theoretical_mutto_AA(theoretical_mutto_AA,ax=ax)
        
        else:
            raise TypeError("theoretical_mutto_AA should have a type of pd.DataFrame or np.ndarray, or should be specified as 'epPCR'")

    plt.tight_layout()
    
    if figname:
        plt.savefig(figname, dpi=300, bbox_inches='tight')
        
        
def highlight_cell(x,y, ax=None, **kwargs):
    """
    Helper function to high light certain cells in a heatmap. This function is only called by highlight_theoretical_mutto_AA(theoretical_mutto_AA,ax)
    """
    
    rect = plt.Rectangle((x-.5, y-.5), 1,1, fill=False, **kwargs)
    ax = ax or plt.gca()
    ax.add_patch(rect)
    return rect

def highlight_theoretical_mutto_AA(theoretical_mutto_AA,ax,color='white',linewidth=3):
    """
    Helper function to highlight selected cells in heatmap with a given axis. 
    Input: 
    theoretical_mutto_AA -- numpy array or pandas dataframe. rows: Mutate from. Columns Mutate to.
    ax -- matplotlib axis object to add highlights to
    color -- color of the highlight, default is white
    linewidth -- linewidth of the highlight, default is 3.
    """
    
    num_mutate_from = np.shape(theoretical_mutto_AA)[0]
    num_mutate_to = np.shape(theoretical_mutto_AA)[1]


    for i in range(num_mutate_from):
        for j in range(num_mutate_to):
            if theoretical_mutto_AA[i,j]!=0:
                highlight_cell(i,j, ax=ax, color=color, linewidth=linewidth)
                
                
def possible_1mutations(original_codon):
    """
    Return the possible AA codon mutations when 1 nt mutation is allowed. Helper function to highlight theoretical epPCR cells
    
    input: str original 3 digit codon
    return: 1. list of strings, possible translations
            2. list of strings, possible mutated-to codons
    """
    original_codon = original_codon.upper()
    sing_muts = []
    sing_muts_AA = []
    for i,nt in enumerate(original_codon):
        for possible_nt in 'agct'.upper():
            sing_mut = original_codon[:i]+possible_nt+original_codon[i+1:]
            sing_muts.append(str(sing_mut))
            sing_muts_AA.append(str(Seq(sing_mut).translate()))
    
    sing_muts_AA=list(set(sing_muts_AA))
    sing_muts_AA.remove(Seq(original_codon).translate())
    sing_muts=list(set(sing_muts))
    sing_muts.remove(original_codon)
    
    return sing_muts_AA, sing_muts



def theoretical_epPCR_matrices_mutto(reference_seq, 
                                     AA_sorted_on_wt_abundance=True, 
                                     sorted_AAs='ARNDEQGHILKMFPSTWYVC'):
    """
    Helper function to calculate the theoretical matrix for AA mutations vs. location in epPCR.
    
    input: reference_seq -- str or Bio.Seq.Seq object
           AA_sorted_on_wt_abundance -- Bool. Default value True. Whether choose to order the columns in the output DataFrame by the abundance of AAs in the WT sequence. If set to be False AAs will be sorted by the following order: 'ARNDEQGHILKMFPSTWYVC'. The AA_sorted_on_wt_abundance parameter overwrites the sorted_AAs parameter in this funciton.
           sorted_AAs -- If AA_sorted_on_wt_abundance is set to be False, provide a customized order of 20 AAs as columns in the output DataFrame. If AA_sorted_on_wt_abundance is set to be False and not specified with a customized value the default order is 'ARNDEQGHILKMFPSTWYVC'.
    Output: 
    """
    
    AAs=sorted_AAs
    reference_seq = Seq(reference_seq)
    reference_AA = mutAA.translate_target_aa(reference_seq)
  
    if AA_sorted_on_wt_abundance:
        AAs=sort_AA_on_target_abundance(reference_seq)[0]
  
    
    # Initialize AA_list and then populate
    AA_list=[]
    for loc in range(len(reference_AA)):
        AA_dict = {}
        for AA in AAs:
            AA_dict[AA] = 0
        AA_list.append(AA_dict)
    
   
    # calculate possible mutations at each AA location
    start_codon = 'ATG'
    start_index = reference_seq.upper().find('ATG')
    for ntloc in range(start_index, start_index+3*len(reference_AA), 3):
        
        codon = reference_seq[ntloc:ntloc+3]
        possible_mutations=possible_1mutations(codon)[0]
        for AA in possible_mutations:
            
            if AA != '*': # The possible_1mutations(codon)[0] function returns stop codon info. 
                          # For plotting purpose this column is removed from AA_list before transformed into a df
                AA_list[(ntloc-start_index)//3][AA]=1
    
    df_theoretical_epPCR_matrix_mutto_loc = pd.DataFrame(AA_list)
    df_theoretical_epPCR_matrix_mutto_AA = matrix_mutto_AA(df_theoretical_epPCR_matrix_mutto_loc, reference_seq)


    return df_theoretical_epPCR_matrix_mutto_loc, df_theoretical_epPCR_matrix_mutto_AA


        
        
        



