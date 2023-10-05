#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 14:30:50 2023

@author: sshanker
"""

# import os, sys
# parent_path = os.path.split(os.path.dirname(__file__))[0]
# sys.path.append(parent_path)
from ADCP.openMMmethods import (fix_my_pdb, makeChidNonPeptide, 
                                make_disulfide_between_cys_in_last_chain, 
                                cyclize_backbone_of_last_chain, 
                                get_residue_templates_for_residues_with_same_graphs,
                                AMBERFFXMLFILE,
                                ffxml_path,
                                GBN2IMPLICIT,DEFAULTSYSTEMFFXMLS)
from ADCP.utils import replace_msg

from prody import writePDB
from MolKit2 import Read
import tempfile

import os
from openmm.app import Modeller, PDBFile, ForceField, Simulation, HBonds
from openmm.unit import kilojoules_per_mole, nanometer
import parmed

SYSTEMTEMPDIR = tempfile.gettempdir()

def combine_rec_and_peptide_pdbs(receptorpdb=None, peptidepdb=None, peptide_model_number=0, out_complex_name=None):
    rec = Read(receptorpdb) 
    rec = rec._ag
    makeChidNonPeptide(rec)
    
    pep_pdb = Read(peptidepdb)
    pep_pdb._ag.setACSIndex(peptide_model_number) 
    pep_pdb._ag.setChids([x for x in 'Z'*pep_pdb._ag.numAtoms()]) # setting peptide chid to 'B'
    
    combinedRecPep= rec + pep_pdb._ag    
    writePDB(out_complex_name, combinedRecPep.select('not hydrogen')) 
    

def amber_inputs_from_pdb_files(**kw):#receptorpdb=None, peptidepdb=None, options =None, out_prefix='mol'):  
    receptorpdb = kw['recpdb']
    peptidepdb = kw['peppdb']
    out_prefix = kw['outprefix']
    env = kw['omm_environment']
    rec_plus_pep = kw['complex']
    
    if rec_plus_pep == None:    
        out_complex_pdb = tempfile.NamedTemporaryFile(suffix='.pdb')    
        out_complex_pdb_name = out_complex_pdb.name
        combine_rec_and_peptide_pdbs(receptorpdb, peptidepdb,out_complex_name=out_complex_pdb_name)    
        temp_file_fixed = out_complex_pdb_name[:-4]+"_fixed.pdb"
    else:
        out_complex_pdb_name = rec_plus_pep
        
    temp_file_fixed = out_complex_pdb_name[:-4]+"_fixed.pdb"    
    fixed_pdb = fix_my_pdb( out_complex_pdb_name,out = temp_file_fixed, NonstandardResidueTreatment=False,kw=kw)
    
    if 'out_complex_pdb' in locals():
        out_complex_pdb.close()
        
    pdb = PDBFile(fixed_pdb)
    topology = pdb.topology
    positions = pdb.positions
    residueTemplates = get_residue_templates_for_residues_with_same_graphs(topology)    
    
    ffxml_files = [AMBERFFXMLFILE]   
    
    if not 'none' in kw['systemffxml'].split(":"):            
        for ff in  kw['systemffxml'].split(":"):
            print('using', ff)
            ffxml_files.append(os.path.join(ffxml_path, ff+"_ff.xml"))
        

    if not kw['userffxml'] == None:
        user_ffxml_set = kw['userffxml'].split(":")
        for uffxml in user_ffxml_set:
            ffxml_files.append(uffxml)
            
    if env == 'implicit':
        ffxml_files.append(GBN2IMPLICIT)
        force_field = ForceField(*ffxml_files)      
        msg = "Using GBSA (gbn2) environment for energy minimization"
    else:
        force_field = ForceField(*ffxml_files)
        msg = "Using in-vacuo environment for energy minimization"

 
    system = force_field.createSystem(  topology, nonbondedCutoff=1*nanometer, 
                                      residueTemplates=residueTemplates)       
    structure = parmed.openmm.topsystem.load_topology( topology, system, positions)
    # print(pdb_str)
    structure.save(out_prefix+'.parm7', overwrite=True)
    structure.save(out_prefix+'.rst7' , format='rst7', overwrite=True)
    print("Amber input files are successfully created!")
    if os.path.isfile(temp_file_fixed):
        os.remove(temp_file_fixed)
        
        
def run_rec_and_pep_pdb_to_amber(**kw):
    amber_inputs_from_pdb_files(receptorpdb=kw['recpdb'], peptidepdb=kw['peppdb'],options=kw,out_prefix=kw['outprefix'])
   


import sys

if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser(description='ADCP file converter', 
                                  usage="usage: python %(prog)s -r receptor -p peptide -o parmOutPrefix")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="0.0.1" )
    parser.add_argument("-r", "--receptor",dest="recpdb",
                        help="receptor pdb or pdbqt file input")    
    parser.add_argument("-p", "--peptide",dest="peppdb",
                        help="peptide pdb or pdbqt file input")  
    
    parser.add_argument("-c", "--complex",dest="complex",
                        help="To specify a pdb file having both receptor and the peptide (complex) at the place of separate receptor and peptide pdbs.")  
    
    
    parser.add_argument("-o", "--outprefix",dest="outprefix",
                        help="outprefix for amber coorinate and toplogy files")  
    

    parser.add_argument("-pmodel", "--peptidemodel", type=int, default=0,
                       dest="peptidemodel", help=( 'peptide model number in a multi-model peptide pdb file to use in the complex'))

    parser.add_argument("-env", "--omm_environment", type=str, default='vacuum',
                       dest="omm_environment", help='options: "vacuum" or "implicit". Default is "vacuum"')
    parser.add_argument("-fnst", "--fix_nst", action="store_true",
                       dest="omm_nst", help=replace_msg)
    
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
            amber_inputs_from_pdb_files(**kw)