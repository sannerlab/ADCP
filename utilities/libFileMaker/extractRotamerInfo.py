#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 11:57:55 2023

@author: sshanker
"""


from Step1_libFile2Rotamer import run_extractRotamers

    
if __name__=='__main__':
    import sys, os
    import argparse
    parser = argparse.ArgumentParser(description='ADCPtools Rotamer data extractor from Dunbrack like lib files V0.1', 
                                  usage="usage: python %(prog)s -l libfiles -p pdbfiles -r rotamer_details.dat -i _random.lib -o outdir")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="0.1")
    parser.add_argument("-l", "--libFileDir",dest="libfileDir", required=True,
                        help="Path of directory containing Dunbrack like rotamer lib files for L variants of NSTs")   
    parser.add_argument("-i", "--libFileNameInit",dest="libfile_name_init", required=True,
                        help="Same pattern of lib file name for easy identification as for ASN_library.file, it should be _library.file") 
    parser.add_argument("-s", "--libFileSteriomer",dest="libfileStereomer", default='L',
                        help="L or D steriomer of lib file. Default is 'L'") 
    parser.add_argument("-p", "--pdbFileDir", dest="pdbfileDir", required=True,
                        help="Path of directory containing pdb files of NSTs organized in  D and L subdirectories")
    parser.add_argument("-r", "--rotamerDetailFile",dest="rotamerdetails", required=True,
                        help="Rotamer detail file with names for L and D format amino acids, base potential letters and, list of dihedrals in accordance with the respective rotamer lib file.")
    parser.add_argument("-o", "--outputfile",dest="tempOutputFileName",  required=True,
                        help="Output file path.")
  
    if len(sys.argv) == 1:
        parser.print_help()
        
    else:    
        kwn, unknw = parser.parse_known_args()
        if len(unknw)> 0:
            print('Unknown options: ',end="")
            for u in unknw:
                if u.startswith("-"):
                    print(u, end=" ")
            print("")
        else:
            kw = vars(kwn)
            kw["libfile_name_init"] = "%s%s" % ('%s',kw["libfile_name_init"])
            kw['nextfilename'] = "makeLibFile"
            run_extractRotamers(**kw)
