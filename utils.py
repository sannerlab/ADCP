#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022
@author: sshanker

This is a collection of functions to identify OpenMM support.
It is made as a separate file from openmmmethods.py to quick check requirements
without loading the main functions for openMM support.

* Checks availibility of openMM
* Identifies non-standand aa and suggests (finds) ways to treat them
* Parser flags for openMM parameters. 

"""

import importlib
from colorama import Fore, Style

replace_msg = ('this flag is used to specify the handling of non-standard \
amino acids when minimization is requested. When omitted, it defaults to \
"no-action" and the software will stop if non standard amino acids are \
found in the receptor and minimization is required. Alternative options \
could be: "replace" to replace non standard amino acid with similar standard AAs, \
"delete" to remove non standard amino acid if found, or "replace_delete" \
to replace if replaceable or delete')

def flag_validator(kw):
    procede_after_after_flag_check = 1 # dfault to run the code    
    # check openmm    
    if int(kw['minimize']) > 0:
        print (f'''{Fore.GREEN}
------------------------------------------------------------------
OpenMM minimization flag detected. This step takes more time than 
non-minimization calculations.'        

{Fore.BLUE}DECLARATION: 
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
    if not kw['minimize']:  # no need to validate residues if no NST
        return 1,''
    
    from pdbfixer.pdbfixer import substitutions, proteinResidues
    from MolKit2 import Read
    import numpy as np
    
    substitutable_residues = list(substitutions.keys())
    # print(kw)
    rec = Read(kw['rec'])
    residues_in_pdb = rec._ag.select('name CA').getResnames()
    uniq_residues = np.unique(residues_in_pdb) # for faster calculation
    list_of_identified_non_standard_residues = [i for i in uniq_residues if proteinResidues.count(i)<1]
    
    combined_treatment_options =['no-action']
    treatment_options_per_nsr = [['no-action']]*len(list_of_identified_non_standard_residues) # default do nothing   for all residues separately 
    if len(list_of_identified_non_standard_residues) > 0:
        for indx, ns_res in enumerate(list_of_identified_non_standard_residues):
            print(f"{Fore.GREEN}Non standard residue %s found" % ns_res)
            if substitutable_residues.count(ns_res):
                print (f'Residue %s can be substituted by %s.{Style.RESET_ALL}' % (ns_res, substitutions[ns_res]))
                treatment_options_per_nsr[indx] = ['replace','delete','replace_delete'] # replace; delete; delete if cannot replace                
            else:
                print(f'Residue %s CANNOT be substituted. Use remove option.{Style.RESET_ALL}' % ns_res)
                treatment_options_per_nsr[indx] = ['delete','replace_delete']
    else:
        return 1, rec # no need check further if no NST
     
    num_replaceable = len([0 for i in treatment_options_per_nsr if i.count('replace')>0])
    
    if  num_replaceable == len(treatment_options_per_nsr):  # when all residues are replacable
        print(f'{Fore.GREEN}All non-standard residues are replaceable. -nst "replace" (or "replace_delete") is suggested.{Style.RESET_ALL}')
        combined_treatment_options = ['replace','delete','replace_delete']
    elif num_replaceable == 0:                      # when none residues are replacable
        print(f'{Fore.GREEN}None of the non-standard residues are replaceable. -nst "delete" (or "replace_delete")  is suggested.{Style.RESET_ALL}')
        combined_treatment_options = ['delete','replace_delete']
    else:                                           # when some residues are replacable
        print(f'{Fore.GREEN}Some of the non-standard residues are replaceable. -nst "replace_delete" is suggested.{Style.RESET_ALL}')
        combined_treatment_options = ['delete','replace_delete']
    print("here:", kw['omm_nst'], combined_treatment_options, )
    if combined_treatment_options.count(kw['omm_nst'])> 0:
        return 1,'' #If possible -nst option is already provided
    else:
        print(f'{Fore.RED}Use -nst flag, available options:')
        print(replace_msg)
        print(f'Or, remove -nmin flag.{Style.RESET_ALL}')        
    return 0, rec #returning receptor to speed up calculation
        
    
def add_open_mm_flags(parser):
    parser.add_argument("-nmin", "--omm_nmin", type=int, default=0,
                       dest="minimize", help=( 'The -nmin option specifies the\
                      number of top ranking docking solutions to minimize with\
                      OpenMM after docking. If omitted this number defaults to\
                      0 meaning no solution is minimized. When solutions are\
                      minimized the reported solutions a ranked based on the\
                      minimized energy and by default ordered from best to worse\
                      based on the minimized complex energetic terms.'))    
    parser.add_argument("-dr", "--dockingRanking", action="store_false",
                       dest="dockingRanking", help=( 'When docking solutions are\
                       minimized, the -dr flag is used to prevent the re-ordering\
                       and report the minimized solutions based on docking score\
                       ordered from best to worse. Default is False.'))   
    parser.add_argument("-omm_max_itr", "-omm_max_itr", type=int, default=5,
                       dest="omm_max_itr", help='Maximum steps for OpenMM minimization. Default is 5')   
    parser.add_argument("-omm_env", "-omm_environment", type=str, default='vaccum',
                       dest="omm_environment", help='options: "vacuum" or "implicit". Default is "vacuum"')       
    parser.add_argument("-nst", "-omm_non_standard_res_treatment", type=str, default="no-action",
                       dest="omm_nst", help=replace_msg)  
    return parser