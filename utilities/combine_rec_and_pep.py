#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 14:04:04 2023

@author: sshanker
"""
import os, sys
from MolKit2 import Read
from prody import writePDB
from ADCP.openMMmethods import (fix_my_pdb, makeChidNonPeptide, )
                                # make_disulfide_between_cys_in_last_chain, 
                                # cyclize_backbone_of_last_chain, 
                                # get_residue_templates_for_residues_with_same_graphs,
                                # AMBERFFXMLFILE,
                                # ffxml_path,
                                # GBN2IMPLICIT,DEFAULTSYSTEMFFXMLS)


def generate_complex (**kw):
    '''combines peptide from output_petide and receptor pdbqt file'''
    
    receptor_file = os.path.abspath(kw['recpdb'])
    peptide_file= os.path.abspath(kw['peppdb'] )
    out_dir= os.path.abspath(kw['outdir'])
    num_models=kw['peptidemodel']
    if num_models == 0:
        num_models = 1000000000000
        print("combining all peptide models")
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
        
    #receptor
    rec = Read(receptor_file)
    rec = rec._ag # receptor
    
    file_name_init = os.path.splitext(os.path.split(receptor_file)[1])[0]
    
    makeChidNonPeptide(rec) # removing numerical and 'Z' chains from receptor
    
    if not os.path.isfile(peptide_file):
        print("Clustered peptide file %s does not exist. OpenMM calculation terminated" % peptide_file)
        return
    
    pep_pdb = Read(peptide_file)
    num_mode_extract = min(pep_pdb._ag.numCoordsets(), num_models)   
    
    
    for f in range( num_mode_extract ):
        print("combining receptor and peptide model %5d" % f)
        pep_pdb._ag.setACSIndex(f)  
        pep_pdb._ag.setChids([x for x in 'Z'*pep_pdb._ag.numAtoms()]) # setting peptide chid to 'Z'
        
        #making complex for current peptide and receptor
        out_complex_name = os.path.join( os.path.join(out_dir , file_name_init+ "_recNpep_%d.pdb" % (f)) ) 
        combinedRecPep= rec + pep_pdb._ag
        writePDB(out_complex_name, combinedRecPep.select('not hydrogen'))   
        if not kw['nofix']:
            fix_my_pdb(out_complex_name, out_complex_name[:-4]+"_fixed.pdb", kw=kw)        
    
    
    
    
if __name__ == "__main__":
    print ("WARNING: THIS CODE IS UNDER DEVELOPMENT! IT MAY NOT WORK PROPERLY.")
    import argparse
    parser = argparse.ArgumentParser(description='Rec and peptide joiner', 
                                  usage="usage: python %(prog)s -r receptor -p peptide -o outdir")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="0.0.1" )
    parser.add_argument("-r", "--receptor",dest="recpdb",
                        help="receptor pdb or pdbqt file input")    
    parser.add_argument("-p", "--peptide",dest="peppdb",
                        help="peptide pdb or pdbqt file input")  
    

    parser.add_argument("-o", "--outdir",dest="outdir",
                        help="out directory for complexes")  
    

    parser.add_argument("-pmodel", "--peptidemodel", type=int, default=0,
                       dest="peptidemodel", help=( 'peptide model number(s) in a multi-model peptide pdb file to use in the complex. Default: all'))


    # parser.add_argument("-fnst", "--fix_nst", action="store_true",
    #                    dest="omm_nst", help=replace_msg)
    
    parser.add_argument("-nofix", "--DoNotFixModel",action="store_true",
                       dest="nofix", help="Do not fix model using openMM parameter sets. Default: False")
    
    
    parser.add_argument("-F", "--ffxmlset",dest="systemffxml", default= "sannerlabPLUSswiss",
                        help=("To specify the use of openmm parameter sets (swiss, sannerlab, or none) \
                        from 'ADCP/data/openMMff' directory. \
                        By default, 'sannerlab' is used as primary source of parameters. \
                        Parameter for only NSTs that are not identified by 'sannerlab' will be taken from \
                        'swiss'.'swiss' supports only non-terminal NST residues.\
                        option 'none' can be used to restrict use of any system default parameter sets. \
                        List of supported NSTs from both libraries is given in\
                        'AVAILABLE_PARAMETER.dat' file in 'ADCP/data/openMMff' \
                        directory."))
                        
    parser.add_argument("-f", "--userffxml",dest="userffxml",
                        help=("a ':' separated list of filenames(initials) to \
                        load openMM force-field data for NSTs. To support a new\
                        residue OpenMM requires three files  containing \
                        1) force-field data, 2) bond definitions, and 3) hydrogen\
                        definitions. User must provide all three files in the \
                        same directory with the same initial names, as for given input\
                        : '-f ./XXXX_ff.xml', the program will expect force-field file(./XXXX_ff.xml), \
                        bond definition file(./XXXX_residues.xml), and hydrogen definition file (./XXXX_hydrogen_def.xml) \
                        all three in the same directory."))
                        
                                         #OMM new line
    if len(sys.argv) == 1:
        parser.print_help()
        
    else:    
        # kw = vars(parser.parse_args())
        # runner = runADCP()        
        # runner(**kw)
        kwn, unknw = parser.parse_known_args()
        if len(unknw)> 0:
            print('Unknown options: ',end="")
            for u in unknw:
                if u.startswith("-"):
                    print(u, end=" ")
            print("")
        else:
            kw = vars(kwn)
            generate_complex(**kw)
    # generate_complex(receptor_file, peptide_file, out_dir)
