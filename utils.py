#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022
Last Update 5/23/23
@author: Sudhanshu Shanker

This is a collection of functions to identify OpenMM support.
It is made as a separate file from openmmmethods.py to quick check requirements
without loading the main functions for openMM support.

* Checks availability of openMM
* Identifies non-standand aa and suggests (finds) ways to treat them
* Parser flags for openMM parameters. 

"""

import importlib
#from colorama import Fore, Style

replace_msg = ("""This flag is used to specify the handling of non-standard \
amino acids when minimization is requested. When omitted, the software will \
stop if non standard amino acids are found in the receptor and minimization \
is required. Alternatively if -fnst flag is given, non standard amino acids \
will be swapped with similar standard AAs using pdbfixer(v1.7), and mutate \
non-replaceables to "ALA".""")

def openmm_validator(kw, myprint=print):
    # checks if openMM is installed or not
    procede_after_flag_check = True # default to run the code    
    # check if minimization is asked 
    if int(kw['minimize']) > 0:
        print ("""
------------------------------------------------------------------
OpenMM minimization flag detected. This step takes more time than 
non-minimization calculations.'        

DECLARATION (V1.1.0 build 3): 
a: Support for OpenMM Minimization is still under development. 
b: Currently, it supports docking with "-rmsd 0" flag. 
c: Non-standard amino acids can either be replaced by similar amino \
acids using pdbfixer v1.7, or if pdbfixer does not identify the \
Non-standard amino acid, it can be replaced by ALA.
d: Non-standard amino-acids only in the RECEPTOR are treated.  
e: Currently no support for external parameter input for non-standard amino \
acids. 
   """)        
        # check these packages:
        packages_to_check = ['openmm', 'parmed']
        
        for pkg in packages_to_check:             
            if importlib.find_loader(pkg) == None:  # checking presence of package
                myprint(("%s not found. Please install %s" +
                      " or remove -nmin flag.") % (pkg, pkg))
                if procede_after_flag_check:
                    procede_after_flag_check = False

            
    return procede_after_flag_check

    

def support_validator(kw,myprint=print):
    # this function validates the input options for openMM calculation
    if not kw['minimize']:  # no need to validate residues if no NST
        return True,''
    
    all_flags_allowed = True
    detected_problems = []
    from pdbfixer.pdbfixer import (
        substitutions, proteinResidues,dnaResidues, rnaResidues)
    from MolKit2 import Read
    import numpy as np
    
    # nmin
    if kw['minimize'] < 0:
        detected_problems.append('"-nmin" cannot be a negative integer.')
        if all_flags_allowed:
            all_flags_allowed = False
    #omm_max_itr        
    if kw['omm_max_itr'] < 0:
        detected_problems.append(f'"-nitr" cannot be a negative integer.')
        if all_flags_allowed:
            all_flags_allowed = False
    #ommenvironement      
    if not kw['omm_environment'] in ["vacuum" , "implicit"]:
        detected_problems.append(f'"-env" can be either "vacuum" or "implicit".')
        if all_flags_allowed:
            all_flags_allowed = False   
    # # check cyc flag 
    # if kw['cyclic'] :
    #     detected_problems.append(f'"-cyc/--cyclic" flag detected. ' 
    #                              f'OpenMM support for cyclic peptide is not ' 
    #                              f'implemented yet.')
    #     if all_flags_allowed:
    #         all_flags_allowed = False
    
    # check cys flag
    # if kw['cystein'] :
    #     detected_problems.append(f'"-cys/--cystein" flag detected. '
    #                               f'OpenMM support for cyclic peptide through '
    #                               f'disulfide bond is not implemented yet.')
    #     if all_flags_allowed:
    #         all_flags_allowed = False
    
    # check nst peptide
    if "<" in kw['sequence']:
        detected_problems.append(f'"Non-standard amino acid(s) detected in the PEPTIDE '
                                 f'sequence. OpenMM support for non-standard '
                                 f'amino acids in PEPTIDE is not implemented yet.')
        if all_flags_allowed:
            all_flags_allowed = False

    #check nst receptor
    rec = Read(kw['recpath'])
    residues_in_pdb = rec._ag.select('name CA').getResnames()
    uniq_residues = np.unique(residues_in_pdb) # for faster calculation
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues).union(['N','UNK','HOH'])
    list_of_identified_non_standard_residues = [i for i in uniq_residues if not i in keep]
    if len(list_of_identified_non_standard_residues) > 0:
        myprint("Non standard residue(s) detected.")
        replace_will_be_or_can_be = "will be"
        if not kw['omm_nst']: #If -fnst option is not provided
            replace_will_be_or_can_be = "can be"
            detected_problems.append(f'Use "-fnst" flag for non-standard amino acids:')
            detected_problems.append(replace_msg)
            detected_problems.append(f'Or, remove -nmin flag.') 
            if all_flags_allowed:
                all_flags_allowed = False  
         
        for indx, ns_res in enumerate(list_of_identified_non_standard_residues): 
            if ns_res in substitutions:
                myprint ('Residue %s %s substituted with %s.' % (ns_res, replace_will_be_or_can_be,substitutions[ns_res]))
            else:
                myprint ('Residue %s %s substituted with "ALA".' % (ns_res, replace_will_be_or_can_be))
    if not all_flags_allowed:
        myprint("Please resolve following issues to use openMM based ranking:")
        for msg in detected_problems:
            myprint(msg)
        #print("Exiting now")
       
    return all_flags_allowed, rec #returning receptor to speed up calculation
        

    
    
def add_open_mm_flags(parser):
    # Add flags for openMM based calculations
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
    parser.add_argument("-env", "--omm_environment", type=str, default='vacuum',
                       dest="omm_environment", help='options: "vacuum" or "implicit". Default is "vacuum"')       
    parser.add_argument("-fnst", "--fix_nst", action="store_true",
                       dest="omm_nst", help=replace_msg)  
    return parser
