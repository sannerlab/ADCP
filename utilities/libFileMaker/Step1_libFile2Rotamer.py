#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  9 13:15:44 2023
@author: sshanker

(Step 1)
DESCRIPTION: This program reads rotamer lib file(s) (swiss) AA dihedral definition file "aa_and_dihedral_def.dat"
to generate a python file as "tempRotamers.py" that containes information of all rotamer dihedrals and possible dihedral combinations.


the "data_dir" directory should have sub-directoris "libfiles", "pdbfiles", 
   and the file "aa_dihedral_definition.dat"
   
format of files: 
   in libfiles dir:
       <AANAME>_bbind_Gfeller.lib # only for L type
       
   in pdbfile dir:
       L/
           <AANAME>.pdb
       D/ 
           <D-AANAME>.pdb
           
   data in sample_dihedral_definition.dat file:
    #AAName Base DTypeName DIHEDRALS
    004 F D004 N-CA-CB-CG1
    0A1 F D0A1 N-CA-CB-CG CA-CB-CG-CD1 CE1-CZ-OH-CH
    0AF W D0AF N-CA-CB-CG CA-CB-CG-CD1
    26P D D26P N-CA-CB-CG CA-CB-CG-CD CB-CG-CD-CE CG-CD-CE-CZ
    6CL D D6CL N-CA-CB-CG CA-CB-CG-CD CB-CG-CD-CE CG-CD-CE-CZ CD-CE-CZ-OH1
    HHK K DHHK N-CA-CB-CG CA-CB-CG-CD CB-CG-CD-CE CG-CD-CE-CZ CD-CE-CZ-CH CE-CZ-CH-NJ
          
"""
import numpy as np
import os


def read_independent_lib_file(libfile):    
    # Checked for bbind lib files obtained from www.swisssidechain.ch
    fid = open(libfile,'r')
    data = fid.readlines()
    fid.close()    
    
    l1_cheked = False # use line 1 to identify number of rotamers    
    dihedral_sets = []
    prob_set = []
    for i in data:
        if i.startswith("#"):  # comments
            continue
        if len(i.strip()) ==0: # empty lines
            continue
        ispl = i.split()
        
        if not l1_cheked:
            if len(ispl) < 20:
                count_dihedrals = ispl[1:5].count('1')
            elif len(ispl) > 20:                
                extra_idxes = (len(ispl) -19)//3  # size 19 is for 4 rotamer ; after that size increases by 3
                count_dihedrals = ispl[1:5+extra_idxes].count('1')
            l1_cheked = True
            
        chis = [float(v) for i,v in enumerate(ispl[-(count_dihedrals*2):]) if not i%2] ## reading each second value start from 0
        
        # np.array(ispl[-(count_dihedrals*2):][0::2]).astype(float)
        dihedral_sets.append(chis)  
        prob_v = 7
        if count_dihedrals> 4:
            prob_v=count_dihedrals+3        
        prob_set.append(float(ispl[prob_v]))    
        
    sort_arg = np.argsort(prob_set)[::-1]  # sorted by probability high to low
    prob_set = np.array(prob_set)[sort_arg]
    dihedral_sets = np.array(dihedral_sets)[sort_arg,:]        
    return dihedral_sets.tolist(), prob_set.tolist()


def getAAnameAndDihedralInfo(dihDefFile):    
    "Get information of dihedrals to rotate and AA name from a manually created file"
    fid = open(dihDefFile,'r')
    data = fid.readlines()
    fid.close()
    
    class my_dot_var:
        pass
    
    dihedral_info = {}    
    for line in data:
        if line.startswith("#"):
            continue  
        if len(line.strip()) < 1:
            continue
        line_before_comment = line.split("#")[0].split("//")[0]  # to igone comments in the same line
        line_spl = line_before_comment.split()        
        dihedral_list = [d.split("-") for d in line_spl[3:]]
        base_aa = line_spl[1]
        d_type_name = line_spl[2]
        
        dihedral_info[line_spl[0]] = my_dot_var()
        dihedral_info[line_spl[0]].dihedral_list = dihedral_list
        dihedral_info[line_spl[0]].d_type_name = d_type_name
        dihedral_info[line_spl[0]].base_aa = base_aa
        dihedral_info[line_spl[0]].l_type_name = line_spl[0]
             
    return dihedral_info

def run_extractRotamers(**kw):    
    AAdihedralInfoFile = kw['rotamerdetails']
    libfileDir = kw['libfileDir']
    pdbfileDir = kw['pdbfileDir']
    libfile_name_init = kw['libfile_name_init']
    input_lib_file_type = kw['libfileStereomer'] # We use L
    outfile_with_path = kw['tempOutputFileName']    
    
    if not outfile_with_path.endswith('.py'):
        outfile_with_path += ".py"
    
    temp_dir = os.path.split(outfile_with_path)[0]
    if 'nextfilename' in kw.keys():
        nextfilename = kw['nextfilename']
    else:
        nextfilename = 'Step2_ADCPlibGen.py'
    
    # outfile_with_path = os.path.join(temp_dir, outputfile)
    
    finalOutputDataCollector = getAAnameAndDihedralInfo( AAdihedralInfoFile)
    AAnames = list(finalOutputDataCollector.keys())
        
    if not (os.path.isdir(libfileDir) and 
            os.path.isdir(pdbfileDir)):
        print ('Required rotamer data not found! Please run "python Step0_getDataFromSwisssidechain.py" to \
               download it before running this step.')
        return
        
    
    
    for aa in AAnames:  ## it always should be L type name, D will be genearted by adding D infront    
        lib_file = os.path.join(libfileDir, libfile_name_init % aa)
        dihedral_angles, probs = read_independent_lib_file(lib_file)
        
        
        ## ALWAYS generate for L
        if input_lib_file_type == 'D':
            #sign should be changed
            dihedral_angles = (np.array(dihedral_angles)*-1).tolist()
            
        finalOutputDataCollector[aa].dihedral_angles = dihedral_angles
        finalOutputDataCollector[aa].probs = probs
    
    # # writing rotamer data
    if temp_dir:
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        else:
            if not os.path.isdir(temp_dir):
                print ("A file with the name %d exist. Please change the tempdir name.")
                return
    
    fout_id = open(outfile_with_path,"w+")   
    fout_id.write("#Data for L amino acids\n")    # it will be always for L type
    
    fout_id.write("chiAnglesDef = {\n")        
    for aa in AAnames:
        if len(finalOutputDataCollector[aa].dihedral_list)==0:
            continue    
        fout_id.write("  '%s':%s,\n" % (aa, str(finalOutputDataCollector[aa].dihedral_list)))
    fout_id.write("}\n")
    
    fout_id.write("chiValues = {\n")
    for aa in AAnames:
        if len(finalOutputDataCollector[aa].dihedral_angles)==0:
            continue 
        fout_id.write("  '%s':%s,\n" % (aa, str(finalOutputDataCollector[aa].dihedral_angles)))    
    fout_id.write("}\n")
    
    fout_id.write("chiValueProbs = {\n")
    for aa in AAnames:
        if len(finalOutputDataCollector[aa].probs)==0:
            continue 
        fout_id.write("  '%s':[" %aa)
        for pb in finalOutputDataCollector[aa].probs:
            fout_id.write("%8.6f, " % ( pb))
        fout_id.write("],\n")
    fout_id.write("}\n")
    
    fout_id.write("baseAA = {\n")
    for aa in AAnames:
        if len(finalOutputDataCollector[aa].base_aa)==0:
            continue 
        fout_id.write("  '%s':'%s',\n" % (aa, str(finalOutputDataCollector[aa].base_aa)))
    fout_id.write("}\n")
    
    fout_id.close()
    print("file '%s' is created successfully. Run '%s' to generate lib file." % (outfile_with_path, nextfilename))
    
if __name__ == "__main__":
    "Simple input type, modify it to use argument parser if needed"
    from settings import kw
    run_extractRotamers(**kw)
    
    

    
    
    
    
    
    
