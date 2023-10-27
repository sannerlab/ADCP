#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 16 16:56:39 2023

@author: sudhanshu
"""
import os


current_path = os.path.dirname(__file__)
kw = {}
kw['data_dir'] = os.path.join(current_path,"input_data") 
kw['rotamerdetails'] = os.path.join(kw['data_dir'],'sample_dihedral_definition.dat')
kw['libfile_name_init'] = "%s_bbind_Gfeller.lib"  # so using only the name of AA lib file can be identified
kw['libfileStereomer'] = "L" # it could be D or L # tested for input of L type only
kw['libfileDir'] = os.path.join(kw['data_dir'],"libfiles")
kw['pdbfileDir'] = os.path.join(kw['data_dir'],"pdbfiles")
kw['tempdir'] = os.path.join(current_path,"tempdir") 
kw['tempOutputFileName'] = os.path.join(kw['tempdir'],"tempRotamer.py")
kw['outdir'] = os.path.join(current_path,"output") 
kw['outputLibFileName'] = os.path.join(kw['outdir'] , "rotamer_swiss_nst.lib")