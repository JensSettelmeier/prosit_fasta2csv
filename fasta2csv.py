#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 14:19:30 2022
This code converts a fasta or fas file to a csv file which is supported
by prosit on www.proteomicsdb.org/prosit 
It can be run in terminal, as script or imported for other applications
@author: Jens Settelmeier
jsettelmeier@ethz.ch
"""

import os
import argparse
import numpy as np
from insilicoDigestion import TDWrapper


def prosit_fasta2csv(fasta_file, collison_engery, charge_states):
    """

    Parameters
    ----------
    fasta_file : str
        "path/to/fasta".
    collison_engery : int
        The collision engery in keV as integer.
    charge_states : list of ints
        list of ints for charge states from which randomly is sampled.

    Returns
    -------
    peptide_df : pandas frame
        pandas data frame in the structue prosit expects a csv file.

    """
    peptide_df = TDWrapper(fasta_file)
    peptide_df.rename(columns={"peptides":"modified_sequence"},inplace=True)
    peptide_df['collision_energy'] = collison_engery
    random_charge_states = np.random.choice(charge_states,peptide_df.shape[0]) # note higher charge states are more unlikely - adjust prior probs. It can make sense to consider the same sequence with several charge states \ToDo
    peptide_df['precursor_charge'] = random_charge_states  
    peptide_df.set_index('modified_sequence', inplace=True)
    peptide_df.to_csv('peptides.csv')
    return peptide_df



def parse_args():
    """

    Returns
    -------
    args : TYPE
        Default argument parser for console execution.

    """
    parser = argparse.ArgumentParser(description='Convert a fasta file to a prosit compativle csv')
    parser.add_argument('--ff', '--fasta_file', type=str, default = os.path.join(os.getcwd(),'test_file.fasta'), help='path to the faster file: "path/to/file')
    parser.add_argument('--ce', '--collesion_energy', type=int, default = 28, help='integer between 10 to 50 identifying keV used for fragmentation')
    parser.add_argument('--cs', '--charge_states', type=list, default = [2], help='list of charge states from which randomly will be selected')
    args = parser.parse_args()
    return args

if __name__ == "__main__":
    args = parse_args()
    peptide_df = prosit_fasta2csv(args.ff, args.ce, args.cs) 
    print(peptide_df)