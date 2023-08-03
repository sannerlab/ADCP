#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 12 14:44:26 2022
@author: Sudhanshu Shanker

This is a collection of functions to identify OpenMM supports.
It is made as a separate file from openmmmethods.py to quick check requirements
without loading the main functions for openMM support.

* Checks availability of openMM
* Identifies non-standand aa and suggests (finds) ways to treat them.
* Parser flags for openMM parameters.


Last Update 7/31/23 build 6
Build 6:
    1: file objects to handle ffxml and rotamer library files.
    2: peptide sequence validations to avoid unknown errors by Crankite.
    3: Coarse potential validation to avoid runtime error for rotamers without default Coarse potentials.
    4: rotamer and ffxml(partially) lookup to find loaded rotamers/ffxml and suggest 
       use of specific rotamer/ffxml files.
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
OpenMM minimization flag detected. This step takes more time than
non-minimization calculations.'

DECLARATION (V1.1.0 build 6):
a: Support for OpenMM Minimization is still under development.
b: Currently, it supports docking with "-rmsd 0" flag.
c: Current version provides docking supports for peptides containing ~400 Non-standard amino acids. 
   other unknown non-standard amino acids can either be replaced by similar amino \
acids (pdbfixer v1.7), or if pdbfixer does not identify the non-standard amino \
acid, it can be replaced by ALA.
d: Treatment of non-standard amino-acids in PEPTIDE IS UNDER_DEVELOPMENT. (New from build 4)
e: Currently no support for external parameter input for non-standard amino \
acids.
   """)
        # check these packages:
        packages_to_check = ['openmm', 'parmed']

        for pkg in packages_to_check:
            if importlib.find_loader(pkg) == None:  # checking presence of package
                myprint(("%s not found. Please install %s" +
                      " or remove -nmin flag.") % (pkg, pkg))
                if procede_after_flag_check:
                    procede_after_flag_check = False

    return procede_after_flag_check

def loaded_ffxml_sets(): # current set to default # 
    """ 'swissaa' openmm parameters for non-terminal NSTs generated from 
    charges obtained from swissaa .rtp files. 
    
   'swiss_sannerlab' openmm parameters for terminal as well as intenal NSTs
    charges generated using local RED server with GAMESS. NME/ACE caps used
    for charge generation. N,H,C,O charges addjested as amber FFs as 
    N = -.4157; H = .2719; C = .5973 ; O = -.5679.    
    """

    system_ffxml_sets = ['swissaa']
    #system_ffxml_sets = ['swiss_sannerlab'] #under development 
    return system_ffxml_sets


class rotamerfile:
    ''' This class reads a rotamer lib file and provides various operations to 
    evaluate supports of the rotamer file
    '''
    def __init__(self, rotamerfile_in,myprint=print):
        self.rotamerfile_ = rotamerfile_in
        self.__aalist__ = 'ACDEFGHIKLMNPQRSTVWYO'
        self.residues = []        
        self.__read__()
        self.myprint = myprint
        self.get_all_residues()
        
    def __read__(self):
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
        if len(self.residues) > 0:
            return self.residues        
        self.residues = list(self.AA_dict.keys())
        self.residues.sort()
        return self.residues
    
    def coarse_potential_for_residue(self,res_name):
        return self.AA_dict[res_name]    
    
    def if_residue_exists(self,resname):
        if resname in self.residues:
            return True
        return False        
    
    def if_residue_and_course_combination_possible(self, res_name, coarse_p):
        if not self.if_residue_exists(res_name):
            self.myprint("Residue is not found in rotamer lib.")
            return False   
        
        if not len(coarse_p)== 1:
            self.myprint("Coarse potential needs to be one-letter standard amino acid code!")
            return False
        
        if not coarse_p.upper() in self.__aalist__:
            self.myprint("Coarse potential needs to be one-letter standard amino acid code!")
            return False
            
        if not coarse_p.upper() == 'O':
            return True
        
        if self.coarse_potential_for_residue(res_name).upper() == 'O':
            return False
        
        return True
        
        
class rotamerdata:
    ''' This class reads multiple rotamer files as rotamerfile object and provides
    required operations on multiple rotamer files'''

    def __init__(self,myprint=print):
        self.rotamerfiles = []
        self.rotamerobjects = []
        self.residues =[]
        self.myprint = myprint
        
    def readrotamerfile(self,rot_file):
        rot_file_with_path = os.path.abspath(rot_file)
        
        self.rotamerfiles.append(rot_file_with_path)
        self.rotamerobjects.append(rotamerfile(rot_file_with_path,self.myprint))
        self.get_all_residues()
        
    def if_residue_exists(self,resname):
        for rot_file_object in self.rotamerobjects:
            if rot_file_object.if_residue_exists(resname):
                return True
        return False       
        
    def if_residue_and_course_combination_possible(self, res_name, coarse_p):
        for rot_file_object in self.rotamerobjects:
            if rot_file_object.if_residue_and_course_combination_possible(res_name, coarse_p):
                return True        
        return False
    
    def get_all_residues(self):
        res =[]
        for rot_file_object in self.rotamerobjects:
            res = res + rot_file_object.residues            
        self.residues = list(set(res))
        self.residues.sort()
        
    def get_rotamerfile_name_for_residue(self, resname):
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
        sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data/rotamers')
        system_rot_libs = kw['rotlibs'].split(":")
        for filenm in system_rot_libs:
            rot_file = os.path.join(sys_rotamer_dir, filenm+".lib")
            loaded_rotamers.readrotamerfile(rot_file)
            
    # from user defined
    if not kw['userrotlibs']==None:
        sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data/rotamers')
        UserLibFiles = kw['userrotlibs'].split(":")
        for UlibFile in UserLibFiles:
            UlibPath = os.path.abspath(UlibFile)
            loaded_rotamers.readrotamerfile(UlibPath)
            
    return loaded_rotamers


def all_available_rotamers(kw,myprint=print):
    '''Provide rotamerdata object from currently loaded rotamer files and all other
    (not loaded but) system available rotamer files'''  
    loaded_rotamers = currently_loaded_rotamer_data(kw,myprint)    
    sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data/rotamers')    
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
        self.ffxmlFile = ffxmlFile
        self.myprint = myprint
        self.residues =[]
        self.__read__()
        
    def __read__(self):
        tree = ET.parse(self.ffxmlFile)    
        root = tree.getroot()
        all_res = root.find('Residues').findall('Residue')
        for res in all_res:
            self.residues.append(res.get('name'))               
        self.residues.sort()
        
    def get_all_residues(self):
        if len(self.residues) > 0:
            return self.residues
        self._read__()
        return self.residues
        
    def if_residue_exists(self,resname):
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
    if 1: # currently overriding for 'swissaa' only
        sys_ffxml_dir = os.path.join( os.path.dirname(__file__), 'data/openMMff')
        # system_ffxml_libs = kw['rotlibs'].split(":")
        # currently overriding for 'swissaa' only
        system_ffxml_libs = loaded_ffxml_sets()
        
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
# def all_available_ffxml_data(kw,myprint=print):
#     '''Provide rotamerdata object for currently loaded rotamer files and all other
#     system available rotamer files'''  
#     loaded_rotamers = currently_loaded_rotamer_data(kw,myprint)    
#     sys_rotamer_dir = os.path.join( os.path.dirname(__file__), 'data/rotamers')    
#     for libFile in os.listdir(sys_rotamer_dir):
#         if not libFile.endswith(".lib"):
#             continue
#         if libFile in ['stdaa.lib']: # ignoring standard library
#             continue    
#         libfilepath = os.path.join(sys_rotamer_dir,libFile)
#         if libfilepath in loaded_rotamers.rotamerfiles:
#             continue        
#         loaded_rotamers.readrotamerfile(libfilepath)        
#     return loaded_rotamers
        
    
class ffxmldata:
    ''' This class reads multiple ffxml files as ffxmlfile object and provides
    required operations on multiple ffxml files'''

    def __init__(self,myprint=print):
        self.ffxmlfiles = []
        self.ffxmlobjects = []
        self.residues =[]
        self.myprint = myprint
        
    def readffxmlfile(self,ffxml_file):
        ffxml_file_with_path = os.path.abspath(ffxml_file)
        
        self.ffxmlfiles.append(ffxml_file_with_path)
        self.ffxmlobjects.append(ffxmlfile(ffxml_file_with_path,self.myprint))
        self.get_all_residues()
        
    def if_residue_exists(self,resname):
        for ffxml_file_object in self.ffxmlobjects:
            if ffxml_file_object.if_residue_exists(resname):
                return True
        return False       
        
   
    def get_all_residues(self):
        res =[]
        for ffxml_file_object in self.ffxmlobjects:
            res = res + ffxml_file_object.residues            
        self.residues = list(set(res))
        self.residues.sort()
        
    def get_ffxmlfile_name_for_residue(self, resname):
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
        if val == 'o':
            if seq[pos+1] == '<':
                SAAlist.append(val)
                continue
            else:
                SAAlist.append("-o<")
                # continue
        SAAlist.append(val)
    return SAAlist, NSTlist


def support_validator(kw,myprint=print):
    # this function validates the input sequence and other options for ADCP run
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
            if aa in ['o']:
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
        all_possible_rotamers = all_available_rotamers(kw,myprint)
        loaded_rotamers = currently_loaded_rotamer_data(kw,myprint)
        
        # checking NSTs
        for nst in NSTlist:
            if not nst[0] in all_possible_rotamers.residues:
                detected_problems.append('Undefined rotamers for NST: "%s" is not available from any default rotamer library files' % nst[0] +
                                         '.To run ADCP for peptide with NST "%s", please provide an user-defined rotamer library using "-l" option'
                                         % nst[0] )                
                if all_flags_allowed:
                    all_flags_allowed = False
                    
            else:
                if not nst[0] in loaded_rotamers.residues:
                    name_of_required_rotamer_file = all_possible_rotamers.get_rotamerfile_name_for_residue(nst[0])                                        
                    detected_problems.append('Undefined rotamers for NST "%s". Rotamer library "%s" supports "%s" and '
                                             % (nst[0], '", "'.join(name_of_required_rotamer_file), nst[0]) + 
                                                'can be used with "-L" option.')                                             
                    if all_flags_allowed:
                        all_flags_allowed = False
                        
                else:
                    if not loaded_rotamers.if_residue_and_course_combination_possible(nst[0], nst[1]):
                        if nst[1] == 'o':
                            # FIXME: a better message would be "NST: %s defined in library: %s does not provide a default coarse potential. Do not use o<%s> instead replace o by a valid coarse potential" 
                            detected_problems.append('There is no default coarse potential defined for NST: "%s".' %(nst[0]))
                            detected_problems.append('Please provide a coarse potential amino-acid letter for calculation with NST: "%s".' %(nst[0]))
                            if all_flags_allowed:
                                all_flags_allowed = False
            

    if not kw['minimize']:  # no need to check further if no minimization   
        if not all_flags_allowed: # print error messages and exit
            myprint("Please resolve following issues:")
            for msg in detected_problems:
                myprint("* "+ msg)    
    
        return all_flags_allowed,''
    
    
    ##<<<<< STANDARD ADCP INPUT VALIDATOR ENDS
    
    
    # >>>>>> OPENMM RELATED ERROR CHECKS STARTS
        
    ## Terminal NSTs are not supported yet with OpenMM.Will be added in build 7>>>

    ## QUESTION: I thought you did the QM calculations to suport this case
    ##
    if loaded_ffxml_sets()[0] == 'swissaa':
        if kw['sequence'][0] =="<" or kw['sequence'][1] =="<" or kw['sequence'][-1] == ">":
            detected_problems.append('OpenMM calculation for terminal Non-standard amino acids are not supported with swiss ffxml.')
            #+
            #                        ' the ffxml set "swiss_sannerlab" can be used (under development) for terminal NSTs but currently'+
            #                        ' it does not support all NSTs listed on swisssidechain.ch. To use "swiss_sannerlab",')
            if all_flags_allowed:
                all_flags_allowed = False
                bracket_error = True
    ## Terminal NST check<<<
  
    if len(NSTlist)>0: 
        loaded_ffxmls = currently_loaded_ffxml_data(kw,myprint)
        # print(loaded_ffxmls.ffxmlfiles)
        all_possible_ffxmls = loaded_ffxmls ## as currently user can not define so both are same
        
        # checking NSTs
        for nst in NSTlist:
            if not nst[0] in all_possible_ffxmls.residues:
                if not kw['omm_nst']: #If -fnst option is not provided
                    ## QUESTION: can the user specify a user-defined  ffxml library file ?? This else branch seems to indicate this is possible
                    ##           if it is possible then this message should tell the use to specify his ffxml file and ideally point to some
                    ##           documentation explaining hot to create thus file
                    detected_problems.append('ffxml data for NST: "%s" is not available in default ffxml library files' % nst[0]
                                             + ' use -fnst option to mutate unknown NSTs to standard amino acids, or remove -nmin option'
                                             
                                             ) #+
                                             # '.To run ADCP for peptide with NST "%s", please provide an user-defined rotamer library using "-l" option'
                                             # % nst[0])                
                    if all_flags_allowed:
                        all_flags_allowed = False
                    
            else:
                if not nst[0] in loaded_ffxmls.residues:
                    name_of_required_ffxmls_file = all_possible_ffxmls.get_ffxmlfile_name_for_residue(nst[0])                                        
                    detected_problems.append('No ffxml file is provided for NST "%s". ffxml library that support "%s" is "%s" and '
                                             % (nst[0], nst[0],'", "'.join(name_of_required_ffxmls_file))) #+ 
                                                # 'can be used with "-L" option to perform ADCP docking.')                                             
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
                
    if not all_flags_allowed:
        myprint("Please resolve following issues to use openMM based ranking:")
        for msg in detected_problems:
            myprint("* "+ msg) 
        #print("Exiting now")

    return all_flags_allowed, rec #returning receptor to speed up calculation


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


