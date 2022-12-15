#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022

@author: sshanker
"""

import importlib
from colorama import Fore, Style

replace_msg = '0:do nothing; 1:replace; 2:delete; 3:delete if cannot replace; Default 0.'

def flag_validator(kw):
    procede_after_after_flag_check = 1 # dfault to run the code    
    # check openmm    
    if int(kw['minimize']) > 0:
        print (f'''{Fore.GREEN}
------------------------------------------------------------------
OpenMM minimization flag detected. This step takes more time than 
non-minimization calculations.'        

DECLARATION: 
a: Support for OpenMM Minimization is still under development. 
b: Currently, it supports docking with "-rmsd 0" flag. 
c: Non-standard amino acids can either be replaced by similar amino
   acids, or can be removed from the sequence.(From pdbfixer v1.7)
d: Non-standard amino-acids only in RECEPTOR is treated.  
e: No support for external parameter input for non-standard amino
   acids.{Style.RESET_ALL} 
   ''')        
        if importlib.find_loader('openmm') == None:  # checking presence of openmm
            print(f"{Fore.RED}OpenMM not found. Please install OpenMM or remove -nmin flag.{Style.RESET_ALL}")
            procede_after_after_flag_check = 0
    return procede_after_after_flag_check


def residue_support_validator(kw):    
    if not kw['minimize']:
        return 1
    
    from pdbfixer.pdbfixer import substitutions, proteinResidues
    from MolKit2 import Read
    import numpy as np
    
    substitutable_residues = list(substitutions.keys())
    # print(kw)
    rec = Read(kw['rec'])
    residues_in_pdb = rec._ag.select('name CA').getResnames()
    uniq_residues = np.unique(residues_in_pdb)
    check_non_standard_residues = [i for i in uniq_residues if proteinResidues.count(i)<1]
    
    treatment_options = [0] # default do nothing    
    if len(check_non_standard_residues) > 0:
        for ns_res in check_non_standard_residues:
            print(f"{Fore.GREEN}Non standard residue %s found" % ns_res)
            if substitutable_residues.count(ns_res):
                print (f'Residue %s can be substituted by %s.{Style.RESET_ALL}' % (ns_res, substitutions[ns_res]))
                treatment_options = [1,2,3] # replace; delete; delete if cannot replace                
            else:
                print(f'Residue %s CANNOT be substituted. Use remove option.{Style.RESET_ALL}' % ns_res)
                treatment_options = [2,3]
                
    if treatment_options.count(kw['omm_nst'])> 0:
        return 1
    else:
        print(f'{Fore.RED}Use -nst flag, available options:')
        print(replace_msg)
        print(f'Or, remove -nmin flag.{Style.RESET_ALL}')        
    return 0
        
    
def add_open_mm_flags(parser):
    parser.add_argument("-nmin", "--omm_nmin", type=int, default=0,
                       dest="minimize", help='Maximum number of poses for openmm minimization and energy calculations. Default is 0')    
    parser.add_argument("-ra", "--omm_rearrange", type=int, default=1,
                       dest="omm_rearrange", help='Option for rearranging rescored poses. Default is true.')    
    parser.add_argument("-omm_max_itr", "-omm_max_itr", type=int, default=5,
                       dest="omm_max_itr", help='Maximum steps for OpenMM minimization. Default is 5')    
    parser.add_argument("-omm_env", "-omm_environment", type=int, default=0,
                       dest="omm_environment", help='0 for vaccum; 1 for implicit. Default is 0')    
    parser.add_argument("-nst", "-omm_non_standard_res_treatment", type=int, default=0,
                       dest="omm_nst", help=replace_msg)    
    return parser