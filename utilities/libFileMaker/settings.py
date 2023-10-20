#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:56:39 2023

@author: sudhanshu
"""

kw = {}
kw['data_dir'] = "./input_data/" 
kw['AAdihedralInfoFile'] = 'aa_dihedral_definition.dat'
kw['libfile_name_init'] = "%s_bbind_Gfeller.lib"  # so using only the name of AA lib file can be identified
kw['libfileStereomer'] = "L" # it could be D or L # tested for input of L type only
kw['outdir'] = "./output"
kw['outputFileName'] = "tempRotamer.py"
kw['tempdir'] = 'tempdir' # do not change
