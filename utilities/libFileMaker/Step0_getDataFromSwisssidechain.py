#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 11:30:52 2023

@author: sshanker
this program downloads required files from swisssidechain.ch to make swiss.lib 
rotamer file. swiss.lib is required for docking with non-standard amino-acids 
from swisssidchain.ch database. 
"""
import os, wget, zipfile, shutil

# we need L and D type PDBs and rotamer libraries for L type

def download_and_extract_swisssidechain_data():
    current_path = os.path.dirname(__file__)
    # make required temp and final output files
    urls = ['https://www.swisssidechain.ch/data/download/L_sidechain.zip', # for L PDB and rotamers
            'https://www.swisssidechain.ch/data/download/D_PDB.zip' # for D PDB
            ]
    
    
    
    temp_dir = os.path.join(current_path, 'tempdir')
    temp_extracted_files = os.path.join(temp_dir,'extracted')
    
    final_extracted_path = os.path.join(current_path,'input_data' )
    final_L_pdb_path = os.path.join(final_extracted_path,'pdbfiles','L')
    final_D_pdb_path = os.path.join(final_extracted_path,'pdbfiles','D')
    final_libfiles_path = os.path.join(final_extracted_path,'libfiles')
    
    for dir_ in [temp_dir, temp_extracted_files]:
        if not os.path.isdir(dir_):
            os.makedirs(dir_)
       
    download_success = True
    all_downloaded_files = []
    for url in urls:
        try:
            print("Downloading %s..." % url)
            filename = wget.download(url,out=temp_dir)   
            all_downloaded_files.append(filename)
            print(" Download successful!")
        except Exception as e:
            print(e)
            print("ERROR: %s could not be downloaded." % url)
            print(" Please check the download link and update it in the file '%s' if required."
                  % (os.path.dirname(__file__) ))
            download_success = False
        
        
    # if download successful make out_dirs     
    if download_success:    
        print("Making output directories...")
        for dir_ in [final_extracted_path, final_L_pdb_path, final_D_pdb_path,final_libfiles_path]:
            if not os.path.isdir(dir_):
                print("Making output directoy: %s" % dir_)
                os.makedirs(dir_)
                continue
            print("Output directoy: %s already exists!" % dir_)
       
        L_file_bundle_zip = all_downloaded_files[0]    
        print("Extracting required files from zip archive %s..." % L_file_bundle_zip)
        with zipfile.ZipFile(L_file_bundle_zip, 'r') as zip_ref:
                # if the trg file was renamed it might create a folder with a different name
                # so we first find out the name of the folder in which the maps are
                filenames = zip_ref.namelist()
                               
                for f in filenames:
                    if f.endswith('.pdb'):
                        zip_ref.extract(f,temp_extracted_files)
                        shutil.move(os.path.join(temp_extracted_files,f), os.path.join(final_L_pdb_path, os.path.split(f)[1]))
                        
                    elif f.endswith('bbind_Gfeller.lib.zip'):
                        zip_ref.extract(f,temp_extracted_files)
                        
                        with zipfile.ZipFile(os.path.join(temp_extracted_files, f),'r') as zip_ref2:
                            zip_ref2.extractall(final_libfiles_path)
                        
                        # os.remove(os.path.join(temp_extracted_files,f))
                    
        
        D_PDB_zip = all_downloaded_files[1] 
        print("Extracting required files from zip archive %s..." % D_PDB_zip)
        with zipfile.ZipFile(D_PDB_zip, 'r') as zip_ref:
            filenames = zip_ref.namelist()
            for f in filenames:
                if f.endswith('.pdb'):
                    zip_ref.extract(f,temp_extracted_files)
                    shutil.move(os.path.join(temp_extracted_files,f), os.path.join(final_D_pdb_path, os.path.split(f)[1]))
        print("All required files downloaded and extracted successfully!")
    #delete temp dirs  
    if not download_success:
        print("Downlodof all required data failed! Could not Proceed to next step!")
        
    print("Deleting temp files!")
    for f in all_downloaded_files:
        if os.path.isfile(f):
            os.remove(f)
    shutil.rmtree(temp_extracted_files)
    return download_success
    
            
if __name__ == "__main__":            
    download_and_extract_swisssidechain_data()