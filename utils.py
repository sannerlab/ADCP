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

replace_msg = ('''This flag is used to specify the handling of non-standard \
amino acids when minimization is requested. When omitted, the software will \
stop if non standard amino acids are found in the receptor and minimization \
is required. Alternatively if -fnst flag is given, non standard amino acids \
will be swapped with similar standard AAs using pdbfixer(v1.7), and mutate \
non-replaceables to "ALA".''')

def openmm_validator(kw):
    procede_after_flag_check = True # dfault to run the code    
    # check openmm    
    if int(kw['minimize']) > 0:
        print (f'''{Fore.GREEN}
------------------------------------------------------------------
OpenMM minimization flag detected. This step takes more time than 
non-minimization calculations.'        

{Fore.BLUE}DECLARATION: 
a: Support for OpenMM Minimization is still under development. 
b: Currently, it supports docking with "-rmsd 0" flag. 
c: Non-standard amino acids can either be replaced by similar amino \
acids using pdbfixer v1.7, or if pdbfixer does not identify the \
Non-standard amino acid, it can be replaced by ALA.
d: Non-standard amino-acids only in the RECEPTOR are treated.  
e: Currently no support for external parameter input for non-standard amino \
acids.
f: No cyclic peptides (flags "-cyc" or "-cys") are supported.{Style.RESET_ALL} 

   ''')        
        if importlib.find_loader('openmm') == None:  # checking presence of openmm
            print(f"{Fore.RED}OpenMM not found. Please install OpenMM or remove -nmin flag.{Style.RESET_ALL}")
            if procede_after_flag_check:
                procede_after_flag_check = False
            
    return procede_after_flag_check


def support_validator(kw):
    if not kw['minimize']:  # no need to validate residues if no NST
        return True,''
    
    all_flags_allowed = True
    detected_problems = []
    from pdbfixer.pdbfixer import (
        substitutions, proteinResidues,dnaResidues, rnaResidues)
    from MolKit2 import Read
    import numpy as np
    
    # check cyc flag 
    if kw['cyclic'] :
        detected_problems.append(f'{Fore.RED}"-cyc/--cyclic" flag detected. ' 
                                 f'OpenMM support for cyclic peptide is not ' 
                                 f'implemented yet.{Style.RESET_ALL}')
        if all_flags_allowed:
            all_flags_allowed = False
    
    # check cys flag
    if kw['cystein'] :
        detected_problems.append(f'{Fore.RED}"-cys/--cystein" flag detected. '
                                 f'OpenMM support for cyclic peptide through '
                                 f'disulfide bond is not implemented yet.{Style.RESET_ALL}')
        if all_flags_allowed:
            all_flags_allowed = False
    
    # print(kw)
    #check nst
    rec = Read(kw['rec'])
    residues_in_pdb = rec._ag.select('name CA').getResnames()
    uniq_residues = np.unique(residues_in_pdb) # for faster calculation
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues).union(['N','UNK','HOH'])
    list_of_identified_non_standard_residues = [i for i in uniq_residues if not i in keep]
    if len(list_of_identified_non_standard_residues) > 0:
        print(f"{Fore.GREEN}Non standard residue(s) detected.{Style.RESET_ALL}")
        replace_will_be_or_can_be = "will be"
        if not kw['omm_nst']: #If -fnst option is not provided
            replace_will_be_or_can_be = "can be"
            detected_problems.append(f'{Fore.RED}Use "-fnst" flag for non-standard amino acids:')
            detected_problems.append(replace_msg)
            detected_problems.append(f'Or, remove -nmin flag.{Style.RESET_ALL}') 
            if all_flags_allowed:
                all_flags_allowed = False  
         
        for indx, ns_res in enumerate(list_of_identified_non_standard_residues): 
            if ns_res in substitutions:
                print (f'{Fore.GREEN}Residue %s %s substituted with %s.{Style.RESET_ALL}' % (ns_res, replace_will_be_or_can_be,substitutions[ns_res]))
            else:
                print (f'{Fore.GREEN}Residue %s %s substituted with "ALA".{Style.RESET_ALL}' % (ns_res, replace_will_be_or_can_be))
    if not all_flags_allowed:
        print(f"{Fore.RED}Please resolve following issues to use openMM based ranking:{Style.RESET_ALL}")
        for msg in detected_problems:
            print(msg)
        #print(f"{Fore.RED}Exiting now{Style.RESET_ALL}")
       
    return all_flags_allowed, rec #returning receptor to speed up calculation
        

    
    
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
    parser.add_argument("-nitr", "--omm_max_itr", type=int, default=5,
                       dest="omm_max_itr", help='Maximum steps for OpenMM minimization. Default is 5')   
    parser.add_argument("-env", "--omm_environment", type=str, default='vaccum',
                       dest="omm_environment", help='options: "vacuum" or "implicit". Default is "vacuum"')       
    parser.add_argument("-fnst", "--fix_nst", action="store_true",
                       dest="omm_nst", help=replace_msg)  
    return parser