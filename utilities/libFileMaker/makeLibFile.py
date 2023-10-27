#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 13:34:16 2023

@author: sshanker
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 27 11:57:55 2023

@author: sshanker
"""


from Step2_ADCPlibGen import run_libGen

# def run_check(**kw):
#     run_allowed = True
#     for vals in ["tempOutputFileName", "pdbfileDir", "rotamerdetails", "outputLibFileName"]


    
if __name__=='__main__':
    import sys,os

    import argparse
    parser = argparse.ArgumentParser(description='ADCPtools Lib file maker from file with pre-extracted rotamer data V0.1', 
                                  usage="usage: python %(prog)s extractRotamerInfo -l libfiles -p pdbfiles -r rotamer_details.dat -i _random.lib -o outdir")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="0.1")
    parser.add_argument("-t", "--tempRotamerFile",dest="tempOutputFileName", required=True,
                        help="Temporary file with rotamer data file name extracted by extractRotamerInfo command") 
    parser.add_argument("-p", "--pdbFileDir", dest="pdbfileDir", required=True,
                        help="Path of directory containing pdb files of NSTs organized in  D and L subdirectories")
    parser.add_argument("-r", "--rotamerDetailFile",dest="rotamerdetails", required=True,
                        help="Rotamer detail file with names for L and D format amino acids, base potential letters and, list of dihedrals in accordance with the respective rotamer lib file.")
    parser.add_argument("-o", "--outputlibfile",dest="outputLibFileName", required=True,
                        help="Output library File name.")
    parser.add_argument("-s", "--libFileSteriomer",dest="libfileStereomer", default='L',
                        help="L or D steriomer of lib file. Default is 'L'") 
    parser.add_argument("-temp",dest="tempdir", default='tempdir', 
                        help="Directory for temp file generation. default is './tempdir'")
    
  
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
            kw['prevfilename'] = "extractRotamerInfo"
            run_libGen(**kw)
