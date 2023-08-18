#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022
@author: Sudhanshu Shanker

This file has collection of functions to evaluate input command and openMM 
support options given for adcp docking.

* Checks availability of openMM
* Identifies non-standand aa and suggests (finds) ways to treat them.
* Parser flags for openMM parameters.


Last Update 7/31/23 build 6
Build 6:
    Update 7/31/23 
    1: file objects to handle ffxml and rotamer library files.
    2: peptide sequence validations to avoid unknown errors by Crankite.
    3: Coarse potential validation to avoid runtime error for rotamers without default Coarse potentials.
    4: rotamer and ffxml(partially) lookup to find loaded rotamers/ffxml and suggest 
       use of specific rotamer/ffxml files.
    
    Update 8/8/23 
    1: Added openMM support for amino-acids with similar graphs e.g. (ALO & THR; ILE & IIL).
    2: Added support for terminal NSTs. using swiss_sannerlab.xml file.
    3: swiss_sannerlab.xml have AMBER compatible partial charges generated from local RED/pyRED server. 
    3: Added scheme to override default hydrogens.xml to avoid conflicts with some NST templates, like for ORN. 
    4: Solved atom-type errors for 2 letter element codes in rotamer files.
    
    Update 8/16/23
    1: Suggest correct NST name if case of NST wrong like for "af1" or 'AF1', it will suggest to use "Af1"
    2: Some other bug fixes.
    
    Update 8/16/23
    1: Support for post-docking openMm minimization added. ("-pdmin", "--postDockMinimize")
        ** required combinations 1) -pdmin with -nmin 2) if already minimized file in the directory -pdmin, -nmin and -O
    
    Remaining:
    1: ffxml parameters for all NSTs.
    2: clustering of similar atomtypes in ffxml file.
    
"""

import importlib
import os
import xml.etree.ElementTree as ET
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
OpenMM minimization flag detected. This step takes more time than non-minimization calculations.

DECLARATION (V1.1.0 build 6):
a: Support for OpenMM Minimization is still under development.
b: Currently, it supports docking with "-rmsd 0" flag.
c: Current version provides docking supports for peptides containing ~400 (L and D) NSTs \
and openMM support for 182 (182x2 for D and L) NSTs. 
d: Other unidentified non-standard amino acids can either be replaced by similar amino \
acids (pdbfixer v1.8), or if pdbfixer does not identify a non-standard amino \
acid, it can be replaced by ALA.
e: Currently no support for external parameter input for non-standard amino \
acids.
   """)
        # check availability of these packages:
        packages_to_check = ['openmm', 'parmed']

        for pkg in packages_to_check:
            if importlib.find_loader(pkg) == None:  # checking presence of package
                myprint(("%s not found. Please install %s" +
                      " or remove -nmin flag.") % (pkg, pkg))
                if procede_after_flag_check:
                    procede_after_flag_check = False

    return procede_after_flag_check

def loaded_ffxml_names(kw): # currently set to default # 
    """ 'swissaa' openmm parameters for non-terminal NSTs generated from 
    charges obtained from swissaa .rtp files. 
    
   'sannerlab' openmm parameters for terminal as well as intenal NSTs
    charges generated using local RED server with GAMESS. NME/ACE caps used
    for charge generation. N,H,C,O charges addjested as amber FFs as 
    N = -.4157; H = .2719; C = .5973 ; O = -.5679.    
    
    """
    #system_ffxml_sets = ['swiss'] #under
    system_ffxml_sets = ['sannerlab'] #under development 
    system_ffxml_sets = [kw['systemffxml']]    
    return system_ffxml_sets

class rotamerfile:
    ''' This class reads a rotamer lib file and provides various operations to 
    evaluate supports of the rotamer file
    '''
    def __init__(self, rotamerfile_in,myprint=print):
        '''Intiation'''
        self.rotamerfile_ = rotamerfile_in
        self.__aalist__ = 'ACDEFGHIKLMNPQRSTVWYO'
        self.residues = []        
        self.__read__()
        self.myprint = myprint
        self.get_all_residues()
        
    def __read__(self):
        '''Reads rotamer file and collects required details'''
        fid = open(self.rotamerfile_,'r')
        data = fid.readlines()
        fid.close()
        AA_dict ={}
        for line in data:
            if line.startswith('rotamer'):
                AA_dict[line.split()[1]] =  line.split()[2]
        
        self.AA_dict = AA_dict
        self.name = os.path.split(self.rotamerfile_)[1][:-4]
        
    def get_all_residues(self):
        '''returns names of all residues provided by the file'''
        if len(self.residues) > 0:
            return self.residues        
        self.residues = list(self.AA_dict.keys())
        self.residues.sort()
        return self.residues
    
    def get_case_insensitive_residues(self, in_res):
        '''returns residues names in capital letter and associated libraries'''
        upper_case_res = [r.upper() for r in self.residues]
        if not in_res.upper() in upper_case_res:
            return None        
        return self.residues[upper_case_res.index(in_res.upper())]
    
    def coarse_potential_for_residue(self,res_name):
        '''rerturns the coarse potential for specific residue'''
        return self.AA_dict[res_name]    
    
    def if_residue_exists(self,resname):
        '''Checks if residue exists in the rotamer file.'''
        if resname in self.residues:
            return True
        return False        
    
    def if_residue_and_course_combination_possible(self, res_name, coarse_p):
        '''Validates the combination of NST and coarse potential'''
        
        # checks if residue name valid
        if not self.if_residue_exists(res_name):
            self.myprint("Residue is not found in rotamer lib.")
            return False   
        
        # checks if corse potential is one letter long
        if not len(coarse_p)== 1:
            self.myprint("Coarse potential needs to be one-letter standard amino acid code!")
            return False
        
        #checks if coarse potential letter is in [amino acid letters + "O"]
        if not coarse_p.upper() in self.__aalist__:
            self.myprint("Coarse potential needs to be one-letter standard amino acid code!")
            return False
            
        # if it is not default as "O" any other amino acid coarse potential combination possible
        if not coarse_p.upper() == 'O':
            return True
        
        # if input is "O" (default) but rotamer file does not provide default
        # coarse potential ("O" at the place of any other AA letter), combination is not possible.
        if self.coarse_potential_for_residue(res_name).upper() == 'O':
            return False
        
        return True

        
class rotamerdata:
    ''' This class reads multiple rotamer files as rotamerfile object and provides
    required operations on multiple rotamer files'''

    def __init__(self,myprint=print):
        '''Initiation'''
        self.rotamerfiles = []
        self.rotamerobjects = []
        self.residues =[]
        self.myprint = myprint
        
    def readrotamerfile(self,rot_file):
        '''Reads rotamer file as rotamer file object'''
        rot_file_with_path = os.path.abspath(rot_file)
        
        self.rotamerfiles.append(rot_file_with_path)
        self.rotamerobjects.append(rotamerfile(rot_file_with_path,self.myprint))
        self.get_all_residues()
        
    def if_residue_exists(self,resname):
        '''checks if specific residue is present in loaded rotamer files'''
        for rot_file_object in self.rotamerobjects:
            if rot_file_object.if_residue_exists(resname):
                return True
        return False       
        
    def if_residue_and_course_combination_possible(self, res_name, coarse_p):
        '''evaluates if specific residue and coarse amino acid letter combination possible'''
        for rot_file_object in self.rotamerobjects:
            if rot_file_object.if_residue_and_course_combination_possible(res_name, coarse_p):
                return True        
        return False
    
    def get_all_residues(self):
        '''returns all available residues from all loaded rotamer files'''
        res =[]
        for rot_file_object in self.rotamerobjects:
            res = res + rot_file_object.residues            
        self.residues = list(set(res))
        self.residues.sort()
        
    def get_case_insensitive_residues(self,res_in):
        ''' returns all available residues from all loaded rotamer files in capital letter and associated rotamer libraries'''
        similar_res = []
        rotamer_file_names =[]
        for rot_file_object in self.rotamerobjects:
            similar_res_trial = rot_file_object.get_case_insensitive_residues(res_in)
            if not similar_res_trial == None:
                similar_res.append(similar_res_trial)   
                rotamer_file_names.append(rot_file_object.name)
        if len(similar_res) > 0:
            return similar_res, rotamer_file_names
        
        return None,None       
        
    def get_rotamerfile_name_for_residue(self, resname):   
        '''returns the names of rotamer files that provide the specific residue'''
        rotfile_names = []
        for rot_file_object in self.rotamerobjects:
            if rot_file_object.if_residue_exists(resname):
                rotfile_names.append(rot_file_object.name)        
        return rotfile_names
        
        
def currently_loaded_rotamer_data(kw,myprint=print):
    '''Provide rotamerdata object for currently loaded (provided by user)
    rotamer files'''  
    loaded_rotamers = rotamerdata(myprint)   
    
    # from system defaults
    if not kw['rotlibs']==None:
        sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data','rotamers')
        system_rot_libs = kw['rotlibs'].split(":")
        for filenm in system_rot_libs:
            rot_file = os.path.join(sys_rotamer_dir, filenm+".lib")
            loaded_rotamers.readrotamerfile(rot_file)
            
    # from user defined
    if not kw['userrotlibs']==None:
        sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data','rotamers')
        UserLibFiles = kw['userrotlibs'].split(":")
        for UlibFile in UserLibFiles:
            UlibPath = os.path.abspath(UlibFile)
            loaded_rotamers.readrotamerfile(UlibPath)
            
    return loaded_rotamers


def all_available_rotamers(kw,myprint=print):
    '''Provide rotamerdata object from currently loaded rotamer files and all other
    (not loaded but) system available rotamer files'''  
    loaded_rotamers = currently_loaded_rotamer_data(kw,myprint)    
    sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data','rotamers')    
    for libFile in os.listdir(sys_rotamer_dir):
        if not libFile.endswith(".lib"):
            continue
        if libFile in ['stdaa.lib']: # ignoring standard library
            continue    
        libfilepath = os.path.join(sys_rotamer_dir,libFile)
        if libfilepath in loaded_rotamers.rotamerfiles:
            continue        
        loaded_rotamers.readrotamerfile(libfilepath)        
    return loaded_rotamers
 

class ffxmlfile:
    '''Reads openMM ffxmlfile to evaluate availability of specific NST parameters'''
    def __init__(self,ffxmlFile, myprint = print):
        '''initiation'''
        self.ffxmlFile = ffxmlFile
        self.myprint = myprint
        self.residues =[]
        self.__read__()
        
    def __read__(self):
        '''Reading of file and required details'''
        tree = ET.parse(self.ffxmlFile)    
        root = tree.getroot()
        all_res = root.find('Residues').findall('Residue')
        for res in all_res:
            self.residues.append(res.get('name'))               
        self.residues.sort()
        self.name = os.path.split(self.ffxmlFile)[1][:-4]
        
    def get_all_residues(self):
        '''returns all resdiues available from the file'''
        if len(self.residues) > 0:
            return self.residues
        self._read__()
        return self.residues
        
    def if_residue_exists(self,resname):
        '''Checks if the residue is present in current ffxml'''
        if resname in self.residues:
            return True
        return False       
  

## FOR FFXML CURRENTLY ONLY BY DEFAULT SWISS_FF.XML IS LOADED
## USER CAN NOT PROVIDE A NEW FFXML FILE, BUT CAN USE -FNST FLAG TO
## RUN ADCP CALCULATIONS FOR FIXING UNKNOWN RSIDUES USING PDBFIXER

def currently_loaded_ffxml_data(kw,myprint=print):
    '''Provide rotamerdata object for currently loaded ffxml files'''  
    loaded_ffxmls = ffxmldata(myprint)   

    # currently implemented for default swiss only
    # from system defaults
    # if not kw['rotlibs']==None:
    if 1: # currently overriding for 'swiss' only
        sys_ffxml_dir = os.path.join( os.path.dirname(__file__), 'data/openMMff')
        # system_ffxml_libs = kw['rotlibs'].split(":")
        # currently overriding for 'swissaa' only
        system_ffxml_libs = loaded_ffxml_names(kw)
        
        for filenm in system_ffxml_libs:
            ffxml_file = os.path.join(sys_ffxml_dir, filenm+"_ff.xml")
            loaded_ffxmls.readffxmlfile(ffxml_file)
            
    # from user defined # will be implemented later
    # if not kw['userrotlibs']==None:
    #     sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data/rotamers')
    #     UserLibFiles = kw['userrotlibs'].split(":")
    #     for UlibFile in UserLibFiles:
    #         UlibPath = os.path.abspath(UlibFile)
    #         loaded_rotamers.readrotamerfile(UlibPath)
            
    return loaded_ffxmls

# currently not accepts openmm ff from users
def all_available_ffxml_data(kw,myprint=print):
    '''Provide rotamerdata object for currently loaded rotamer files and all other
    system available rotamer files'''  
        
    loaded_ffxmls = currently_loaded_ffxml_data(kw,myprint)      
    sys_ffxml_dir = os.path.join( os.path.dirname(__file__), 'data/openMMff')    
    for xmlFile in os.listdir(sys_ffxml_dir):
        if not xmlFile.endswith("_ff.xml"):
            continue        
        xmlFilePath = os.path.join(sys_ffxml_dir,xmlFile)
        if xmlFilePath in loaded_ffxmls.ffxmlfiles:
            continue        
        loaded_ffxmls.readffxmlfile(xmlFilePath)        
    return loaded_ffxmls
        
    
class ffxmldata:
    ''' This class reads multiple ffxml files as ffxmlfile object and provides
    required operations on multiple ffxml files'''

    def __init__(self,myprint=print):
        '''initiation'''
        self.ffxmlfiles = []
        self.ffxmlobjects = []
        self.residues =[]
        self.myprint = myprint
        
    def readffxmlfile(self,ffxml_file):
        '''Reading of the file and required details'''
        ffxml_file_with_path = os.path.abspath(ffxml_file)        
        self.ffxmlfiles.append(ffxml_file_with_path)
        self.ffxmlobjects.append(ffxmlfile(ffxml_file_with_path,self.myprint))
        self.get_all_residues()
        
    def if_residue_exists(self,resname):
        '''Evaulates if the residue is present in ffxml files'''
        for ffxml_file_object in self.ffxmlobjects:
            if ffxml_file_object.if_residue_exists(resname):
                return True
        return False       
   
    def get_all_residues(self):
        '''Retruns all residues from all read ffxml files'''
        res =[]
        for ffxml_file_object in self.ffxmlobjects:
            res = res + ffxml_file_object.residues            
        self.residues = list(set(res))
        self.residues.sort()
        
    def get_ffxmlfile_name_for_residue(self, resname):
        '''Returns names of files that provide specific residue'''
        ffxmlfile_names = []
        for ffxml_file_object in self.ffxmlobjects:
            if ffxml_file_object.if_residue_exists(resname):
                ffxmlfile_names.append(ffxml_file_object.name)        
        return ffxmlfile_names


def getNamesOfAAandNstsFromSequence(seq):
    '''Reads standard amino acids, NSTs, operations tokens to do something and
    and other unidentified characaters in the input sequence'''   
    NSTlist = []
    SAAlist = [] # standard amino acids
    aa=''
    d_tag = ''
    coarse_tag = ''
    readaa = False
    for pos, val in enumerate(seq):        
        if readaa:
            if val ==">":
                readaa = False
                NSTlist.append([aa+d_tag, coarse_tag])                
                continue
            if pos == len(seq)-1:
                SAAlist.append("->")
                continue                
            aa = aa + val
            continue        
        if val == "<":
            readaa=True
            aa= ''
            d_tag = ''
            coarse_tag = ''
            coarse_tag = seq[pos-1]
            if seq[pos-2] == '&':
                d_tag = "_D" 
            continue
        if val =='&':
            continue        
        if val.lower() == 'o':
            if seq[pos+1] == '<':
                SAAlist.append(val)
                continue
            else:
                SAAlist.append("-%s<" % val)
                # continue
        SAAlist.append(val)
    return SAAlist, NSTlist


def support_validator(kw,myprint=print):
    '''this function validates the input sequence and other options for ADCP run'''
    all_flags_allowed = True
    detected_problems = []
    
    ##>>>> STANDARD ADCP INPUT VALIDATOR START
    # check input peptide sequence 
    
    SAAseq, NSTlist = getNamesOfAAandNstsFromSequence(kw['sequence'])
    # print(SAAseq,NSTlist)
    bracket_error = False
    if '-<' in SAAseq:
        detected_problems.append('At least one "<" is missing in the peptide sequence.')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True
            
    if '-o<' in SAAseq:
        detected_problems.append('At least one "<" is missing after "o" in the peptide sequence.')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True
            
    if '-O<' in SAAseq:
        detected_problems.append('At least one "<" is missing after "O" in the peptide sequence.')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True
    
    if '->' in SAAseq:
        detected_problems.append('At least one ">" is missing in the peptide sequence.')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True
            
    if '&<' in kw['sequence']:    
        detected_problems.append('For D variants of NST use "&" followed by the coarse respresentation potential letter. e.g. "&o<OTYR>"')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True
    
    if '><' in kw['sequence'] or kw['sequence'].startswith('<'):    
        detected_problems.append('At least one coarse respresentation potential letter is missing in the peptide sequence.')                                 # 'NST example: "o<2AS>" or for D form "&o<2AS>"')
        if all_flags_allowed:
            all_flags_allowed = False
            bracket_error = True

    # check standard aa in sequence
    if bracket_error is False:
        saa = 'ACDEFGHIKLMNPQRSTVWY'
        unknownaa=[]
        for aa in SAAseq:
            if aa in ['o','O']:
                continue
            if not aa.upper() in saa:
                unknownaa.append(aa)     
        if len(unknownaa)>0:
            unknownaa = list(set(unknownaa))
            unknownaa.sort()
            detected_problems.append('Unknown amino acid letter(s) in input peptide sequence: "%s".' % '", "'.join(unknownaa))
            if all_flags_allowed:
                all_flags_allowed = False
                
        if len(SAAseq) < 5:  # as SAA contains the Coarse potential letter for NSTs, only len(SAAaeq) will be enough to know the sequence length
            if all_flags_allowed:
                detected_problems.append('Peptide sequence should be longer than 4 amino acids.')  
                all_flags_allowed = False
            

    # check errors and possible sources for NSTs
    if len(NSTlist)>0:  
        # print (NSTlist)     
        if not kw['postdockmin']:
            all_possible_rotamers = all_available_rotamers(kw,myprint)
            loaded_rotamers = currently_loaded_rotamer_data(kw,myprint)
            
            # checking NSTs
            for nst in NSTlist:      
               
                
                if not nst[0] in all_possible_rotamers.residues:
                    detected_problems.append('Undefined rotamers for NST: "%s" is not available from any default rotamer library files.' % nst[0])
                    detected_problems.append('To run ADCP for peptide with NST "%s", please provide an user-defined rotamer library using "-l" option.'
                                             % nst[0] )                
                    if all_flags_allowed:
                        all_flags_allowed = False
                        
                    if not nst[0] in all_possible_rotamers.residues:
                        possible_case_insentive_names,possible_lib_names = all_possible_rotamers.get_case_insensitive_residues(nst[0])
                        if not possible_case_insentive_names == None:
                            detected_problems.append("Possible residue(s) with similar name:")
                            for p_res,p_lib in zip(possible_case_insentive_names,possible_lib_names):
                                detected_problems.append(' "%s" from rotamer library "%s"'   %(p_res,p_lib))                    
                        
                else:
                    if not nst[0] in loaded_rotamers.residues:
                        name_of_required_rotamer_file = all_possible_rotamers.get_rotamerfile_name_for_residue(nst[0])                                        
                        detected_problems.append('Undefined rotamers for NST "%s". Rotamer library(ies) "%s" support(s) "%s" and '
                                                 % (nst[0], '", "'.join(name_of_required_rotamer_file), nst[0]) + 
                                                    'can be used with "-L" option.')                                             
                        if all_flags_allowed:
                            all_flags_allowed = False
                            
                    else:
                        if not loaded_rotamers.if_residue_and_course_combination_possible(nst[0], nst[1]):
                            if nst[1] in['o','O']:                            
                                rotamer_libs_providing_this_aa_params = loaded_rotamers.get_rotamerfile_name_for_residue(nst[0])
                                
                                detected_problems.append('NST: "%s" defined in library: "%s" does not provide a default coarse potential.' %(nst[0], ','.join(rotamer_libs_providing_this_aa_params)))
                                detected_problems.append('Do not use o<%s> instead replace "o" by a valid coarse potential.' %(nst[0]))
                                if all_flags_allowed:
                                    all_flags_allowed = False
                
    
        if not kw['minimize']:  # no need to check further if no minimization   
            if not all_flags_allowed: # print error messages and exit
                myprint("Please resolve following issues:")
                for msg in detected_problems:
                    myprint("* "+ msg)    
        
            return all_flags_allowed,''    
    
    ##<<<<< STANDARD ADCP INPUT VALIDATION ENDS
    
    
    # >>>>>> OPENMM RELATED ERROR CHECKS START    
    if 'swiss' in loaded_ffxml_names(kw):
        if kw['sequence'][0] =="<" or kw['sequence'][1] =="<" or kw['sequence'][-1] == ">" or (kw['sequence'][0] == "&" and kw['sequence'][2] =="<"):
            detected_problems.append('OpenMM calculation for terminal Non-standard amino acids are not supported with "swiss" ffxml set. Try "-F sannerlab"')
            #+
            #                        ' the ffxml set "swiss_sannerlab" can be used (under development) for terminal NSTs but currently'+
            #                        ' it does not support all NSTs listed on swisssidechain.ch. To use "swiss_sannerlab",')
            if all_flags_allowed:
                all_flags_allowed = False
                
    ## Terminal NST check<<<  
    if len(NSTlist)>0: 
        loaded_ffxmls = currently_loaded_ffxml_data(kw,myprint)
        # print(loaded_ffxmls.ffxmlfiles)
        all_possible_ffxmls = all_available_ffxml_data(kw,myprint) ## as currently user can not define so both are same
        
        # checking NSTs
        ## As for standard amino acids, amber uses same paramters for L and D type AAs. so nst[0].split("_") is used 
        ## to remove "_D" for D type res.
        for nst in NSTlist:
            if not nst[0].split("_")[0] in all_possible_ffxmls.residues: #nst[0].split("_") to remove "_D" for D type res
                if not kw['omm_nst']: #If -fnst option is not provided
                    detected_problems.append('Undefined openMM parameters for NST: "%s" is not available from any default ffxml files' % nst[0].split("_")[0] 
                                             + ' use -fnst option to mutate unknown NSTs to standard amino acids, or remove -nmin option'
                                             
                                             ) #+
                                             # '.To run ADCP for peptide with NST "%s", please provide an user-defined rotamer library using "-l" option'
                                             # % nst[0])                
                    if all_flags_allowed:
                        all_flags_allowed = False
                    
            else:
                if not nst[0].split("_")[0] in loaded_ffxmls.residues:
                    name_of_required_ffxmls_file = all_possible_ffxmls.get_ffxmlfile_name_for_residue(nst[0].split("_")[0])                                        
                    detected_problems.append('Undefined openMM parameters for NST "%s". ffxml library(ies) "%s" support(s) "%s" and '
                                             % (nst[0].split("_")[0] , nst[0].split("_")[0] ,'", "'.join(name_of_required_ffxmls_file))+ #+ 
                                                'can be used with "-F" option.')                                             
                    if all_flags_allowed:
                        all_flags_allowed = False
                        

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
        detected_problems.append('"-nitr" cannot be a negative integer.')
        if all_flags_allowed:
            all_flags_allowed = False
    #ommenvironement
    if not kw['omm_environment'] in ["vacuum" , "implicit"]:
        detected_problems.append('"-env" should be "vacuum" or "implicit".')
        if all_flags_allowed:
            all_flags_allowed = False
            

    #check nst receptor
    rec = Read(kw['recpath'])
    residues_in_pdb = rec._ag.select('name CA').getResnames()
    uniq_residues = np.unique(residues_in_pdb) # for faster calculation
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues).union(['N','UNK','HOH'])
    list_of_identified_non_standard_residues = [i for i in uniq_residues if not i in keep]
    if len(list_of_identified_non_standard_residues) > 0:
        myprint("Non standard residue(s) detected in receptor.")
           
        replace_will_be_or_can_be = "will be"
        if not kw['omm_nst']: #If -fnst option is not provided
            replace_will_be_or_can_be = "can be"
            detected_problems.append('Use "-fnst" option for non-standard amino acids:')
            detected_problems.append(replace_msg)
            detected_problems.append('Or, remove -nmin option.')
            if all_flags_allowed:
                all_flags_allowed = False

        for indx, ns_res in enumerate(list_of_identified_non_standard_residues):
            if ns_res in substitutions:
                myprint ('Receptor residue %s %s substituted with %s.' % (ns_res, replace_will_be_or_can_be,substitutions[ns_res]))
            else:
                myprint ('Receptor residue %s %s substituted with "ALA".' % (ns_res, replace_will_be_or_can_be))
    # all_flags_allowed = False  


    # check ffxml file
    system_ffxml_set = kw['systemffxml']     
    
    if not system_ffxml_set in ['swiss','sannerlab']:
        detected_problems.append('Unknown ffxmlset "%s". "-ffxmlset" should be "swiss" or "sannerlab".' % system_ffxml_set)
        if all_flags_allowed:
            all_flags_allowed = False
        
                
    if not all_flags_allowed:
        myprint("Please resolve following issues to use openMM based ranking:")
        for msg in detected_problems:
            myprint("* "+ msg) 
        #print("Exiting now")

    return all_flags_allowed, rec #returning receptor to speed up calculation



       
def evaluate_requirements_for_minimization(kw,myprint):
    '''To check all required options for post docking minimization'''
    return_val = True  
    
    ## if inputs are not supported, exit now
    if not support_validator(kw, myprint)[0]:     
        return  False  
    
    # as -nmin activates required minimization setting, it must be activated when postdoc minimization is asked
    if not kw['minimize']:                            
        print('Post docking minimization (-pdmin) should be used'+   #OMM new line
                     ' with minimize (-nmin) option') 
        return  False 

    jobName = kw['jobName']  
    omm_out_file = jobName + "_omm_rescored_out.pdb"
    
    # DO NOT WRITE OVER previously minimized file without -O flag, 
    if os.path.exists(omm_out_file) and not kw['overwriteFiles']:
        print("ERROR: openMM minimization output file %s exist! please use -O to force overwritting output files" % omm_out_file)
        return False                 
    
    # Look for required pdb and log file for output
    docked_file = '%s_out.pdb'% jobName
    if os.path.exists(docked_file):
        print('Docked pdb file "%s" for minimization found!' % docked_file)
    else:
        print('ERROR: Expected docked pdb file "%s" not found! Please check the availability of the file and if the sequence (-s), targe file (-T), and jobName (-o) inputs are same as used for docking.'  % docked_file)
        return_val = False
    
    summary_file = '%s_out.pdb'% jobName
    if os.path.exists('%s_summary.dlg'%jobName):
        print('Docking summary file "%s" found!' % summary_file)      
        
        ## check used commands are similar to current ones. rightnow I am only checking the input sequence (-s) and target file (-T)
        fid = open('%s_summary.dlg'%jobName,'r')
        data_lines = fid.readlines()
        fid.close()
        
        # previous command line options
        for dline in data_lines:
            if dline.startswith("MC search command"):
                prev_peptide = dline[dline.index("-t 2"): ].split()[2]  # peptide
                prev_trg_file = dline[dline.index("-T "): ].split()[1].split(os.sep)[-1] # target file name
                break
        
        # check peptide sequence
        if not prev_peptide.replace('"','') == kw['sequence']:
            print('ERROR: Peptide sequence used for docking: %s and currently provided peptide sequence: "%s" are different.' % (prev_peptide,kw['sequence'] ))
            return_val = False   
            
        # check target file
        currently_input_trg_file = kw['target'].split(os.sep)[-1].split(".")[0]
        if not prev_trg_file == currently_input_trg_file:
            print('ERROR: Target file used for docking: "%s" and currently provided targetfile: "%s" are different.' % (prev_trg_file, currently_input_trg_file))
            return_val = False   
    else:
        myprint('ERROR: Expected docking summary file "%s" not found! Please check the availability of the file and if the sequence (-s), targe file (-T), and jobName (-o) inputs are same as used for docking.' % summary_file)
        return_val = False
    return return_val


def extract_target_file(kw,myprint):    
    '''Extract rec from target file for openmm minimization'''
    if not os.path.exists(kw['target']):
        myprint("ERROR: no receptor files found, please provide a .trg file or path to inflated .trg")
        return False
    
    targetFile = kw['target']
    workingFolder = kw['jobName']
    
    if os.path.isdir(targetFile):
        target_folder = targetFile
        
    elif os.path.isfile(targetFile):
        # if target is zip file unzip and replace cmdline arguments
        import zipfile
        myprint( 'Inflating target file %s'%(targetFile))
        if not os.path.exists(workingFolder):
            os.mkdir(workingFolder)       
        
        with zipfile.ZipFile(targetFile, 'r') as zip_ref:
            # if the trg file was renamed it might create a folder with a different name
            # so we first find out the name of the folder in which the maps are
            filenames = zip_ref.namelist()
            folder = filenames[0].split(os.sep)[0]
            zip_ref.extract(os.path.join(folder, 'rigidReceptor.pdbqt' ), workingFolder) # we only need rigidReceptor for docking
        target_folder = os.path.join(workingFolder, folder)
            
    kw['recpath'] = os.path.join(target_folder, 'rigidReceptor.pdbqt')  
    return True
 

def add_open_mm_flags(parser):
    # Add flags for openMM based calculations
    parser.add_argument("-nmin", "--omm_nmin", type=int, default=0,
                       dest="minimize", help=( 'The -nmin option specifies the\
                      number of top ranking docking solutions to minimize with\
                      OpenMM after docking. If omitted this number defaults to\
                      0 meaning no solution is minimized. When solutions are\
                      minimized the reported solutions a ranked based on the\
                      minimized energy by default unless -dr is specified to \
                      retain the docking order before minimizationed. \
                      The solutions are ranked from best to worst according\
                      to the minimized complex energetic terms.'))
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
    
    parser.add_argument("-F", "--ffxmlset",dest="systemffxml", default= "swiss",
                        help=("To use openmm parameter (swiss or sannerlab) \
                        from 'ADCP/data/openMMff' directory (Default:'swiss'). \
                        parameter 'swiss' supports only non-terminal \
                        NST residues. For terminal NSTs use 'sannerlab' library\
                        List of supported NSTs for both libraries is given in\
                        'AVAILABLE_PARAMETER.dat' file in 'ADCP/data/openMMff' \
                        directory."))
                        
    parser.add_argument("-pdmin", "--postDockMinimize",dest="postdockmin",                         
                        action="store_true", default=False,
                        help=("To use openMM minimization on already docked\
                        poses. This option will stop ADCP to perform re-docking \
                        and activates openMM minimization of poses from \
                        already docked output file. If using this option, use target file (-T),\
                        jobName (-o), and sequence (-s) descriptions same as used for docking."))
                  
                      
    
    # THIS PART WILL BE ACTIVATED IN BUILD 7
    # parser.add_argument("-x", "--userffxml",dest="ffxmlibs",
    #                     help=("a ':' separated list of filenames(initials) to \
    #                     load openMM force-field data for NSTs. To support a new\
    #                     residue OpenMM requires three files  containing \
    #                     1) force-field data, 2) bond definitions, and 3) hydrogen\
    #                     definitions, user must provide all three files in the \
    #                     same directory with the same initial names  and the input\
    #                     can be given as: '-x ./XXXX' for ./XXXX_ff.xml, ./XXXX_bond_def.xml,\
    #                     ./XXXX_hydrogen_def.xml"))
                        
    return parser


# libNsts, ffxmlNsts = list_of_supported_NSTs()


