#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:38:09 2023

@author: sshanker
"""
from Step0_getDataFromSwisssidechain import download_and_extract_swisssidechain_data
from Step1_libFile2Rotamer import run_extractRotamers
from Step2_ADCPlibGen import run_libGen
from settings import kw
import os

def curate_swisssidechain_pdbfiles(dir_pdbfiles):
    '''Some pdb files have changed atom assignments compared to rotamer definition
    So we need to correct them manually'''
    
    #TFP4 and D4TM
    for fl in [['TFP4','L'],['D4TM','D']]:    
        tfp4_pdb = os.path.join(dir_pdbfiles,fl[1],'%s.pdb' % fl[0])
        print("Curating %s" % tfp4_pdb )
        fid = open(tfp4_pdb,'r')    
        data = fid.readlines()
        fid.close()
        new_data = []
        for line in data:
            if line.startswith('ATOM     12  CH '):
                line=line.replace('ATOM     12  CH ', 'ATOM     12  CH1')
            new_data.append(line)
            
        fid = open(os.path.join(dir_pdbfiles,fl[1],'%s.pdb' % fl[0]),'w+') 
        fid.writelines(new_data)
        fid.close()
            

download_success = download_and_extract_swisssidechain_data() # download required files
if download_success:  
    kw['AAdihedralInfoFile'] = os.path.join(kw['data_dir'],'swiss_dihedral_definition.dat')    
    kw['outdir'] =  os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','data', 'rotamers'))
    kw['pdbfileDir'] = os.path.join(kw['data_dir'],'pdbfiles')
    kw['libfileDir'] = os.path.join(kw['data_dir'],'libfiles')
    curate_swisssidechain_pdbfiles(kw['pdbfileDir'])    # clean some pdbfiles 
    kw['outputLibFileName'] = os.path.join(kw['outdir'] ,'swiss.lib')
    kw['tempOutputFileName'] = os.path.join(kw['tempdir'], "tempRotamer_swiss.py")
    
    run_extractRotamers(**kw)     # extract rotamers data from downloaded files
    run_libGen(**kw)              # generate ADCP acceptable coordinates and make final swiss.lib file

    
    
