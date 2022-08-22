#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 14:17:50 2022

@author: Jens Settelmeier
"""

import os
from Bio import SeqIO
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

def TDWrapper(fasta_file, output_name='out.txt', miss_cleavage=2, min_len=7, 
              max_len=30, write_with_t=False):
    """
    Parameters
    ----------
    fasta_file : fasta
        Fasta File with proteins.
    output_name : txt, optional
        file name from the ouput which contains the peptides. 
        The default is 'out.txt'.
    miss_cleavage : int, optional
        maximum considered missed cleavages. The default is 2.
    min_len : int, optional
        Minimum length of considered peptides. The default is 7.
    max_len : int, optional
        Maximum considered peptide length. The default is 40.

    Returns
    -------
    peptide_df : pandas frame
        returns da pandas frame with the peptides after digestion
        the proteins in the fasta file with trypsin.

    """
    if write_with_t is True:
        write_command = "t%s\n"
    else:
        write_command = "%s\n"
    # get the fasta file
    handle=SeqIO.parse(fasta_file,'fasta')
    peptide_list = list()
    # path things
    curr_path = os.getcwd()
    output_file = os.path.join(curr_path, output_name)
    with open(output_file,'w', encoding = 'utf-8') as f:
        for record in handle:
            proseq=str(record.seq)
            peptide_list=TrypsinDigestion(proseq,miss_cleavage)
            for peptide in peptide_list:
                f.write(write_command % (peptide))  
    # create the pandas frame with the unique peptides after digestion
    data=pd.read_csv(output_file, header = None, names = ['peptides'])
    data.drop_duplicates(inplace=True)
    tmp = data[data['peptides'].apply(len) >= min_len]
    peptide_df = tmp[tmp['peptides'].apply(len)<=max_len]
    return peptide_df
    

def TrypsinDigestion(proseq,miss_cleavage):
    '''
    original function from https://github.com/yafeng/trypsin
    modified by Jens Settelmeier
    
    Parameters
    ----------
    proseq : str
        protein sequence.
    miss_cleavage : int
        maximum number of missed cleavages.

    Returns
    -------
    peptides : list
        List of peptides after trypsin digestion of the protein from proseq.

    '''
    peptides=[]
    cut_sites=[0]
    for i in range(0,len(proseq)-1):
        if proseq[i]=='K' and proseq[i+1]!='P':
            cut_sites.append(i+1)
        elif proseq[i]=='R' and proseq[i+1]!='P':
            cut_sites.append(i+1)
    
    if cut_sites[-1]!=len(proseq):
        cut_sites.append(len(proseq))

    if len(cut_sites)>2:
        if  miss_cleavage==0:
            for j in range(0,len(cut_sites)-1):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])

        elif miss_cleavage==1:
            for j in range(0,len(cut_sites)-2):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
            
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])

        elif miss_cleavage==2:
            for j in range(0,len(cut_sites)-3):
                peptides.append(proseq[cut_sites[j]:cut_sites[j+1]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+2]])
                peptides.append(proseq[cut_sites[j]:cut_sites[j+3]])
            
            peptides.append(proseq[cut_sites[-3]:cut_sites[-2]])
            peptides.append(proseq[cut_sites[-3]:cut_sites[-1]])
            peptides.append(proseq[cut_sites[-2]:cut_sites[-1]])
    else: #there is no trypsin site in the protein sequence
        peptides.append(proseq)
    return peptides