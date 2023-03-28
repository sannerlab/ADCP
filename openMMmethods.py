#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:29:31 2022
v2
@author: sshanker
"""
import os
import numpy as np
from colorama import Fore, Style

import prody
from prody.measure.contacts import findNeighbors
from prody import writePDB, writePDBStream
from MolKit2 import Read

from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer

'''
Steps:
  1: Identify receptor(s) and peptide chains.
  2: Identify interface residues and,
      Minimize peptide+ interacting receptor residues.
      Estimate E_complex.
  3: Split complex into receptor(s) and peptide and
      Estimate single point energy values for receptor(s) and petide, i.e.
      E_receptor and E_peptide
  4: Calculate different E metrics. and sort best to worst. Combine respective
      models in a single output file.
'''
      

def split_pdb_to_chain_A_and_Z(pdbfile):
    # It changes last chain (peptide) to "Z" 
    # After openMM minimization chain ids are changed back to consecutive letters.
    # like A,B,Z => A,B,C
    
    out_pdb_intial = pdbfile[:-4]    
    mol_object = Read(pdbfile)
    
    # handeling chain ids
    chids = mol_object._ag.getChids()    
    chains = list(set(chids))
    chains.sort()
    
    chids = ''.join(chids)
    chids = chids.replace(chains[-1],'Z')
    chids = [x for x in chids]
    mol_object._ag.setChids(chids)
    
    rec = mol_object._ag.select('not chid Z')
    pep = mol_object._ag.select('chid Z')  
    
    prody.writePDB(out_pdb_intial+"_A.pdb", rec)
    prody.writePDB(out_pdb_intial+"_Z.pdb", pep)
  
    
def identify_interface_residues(pdb_file, needs_b=0):
    mol_temp = Read(pdb_file)    
    rec = mol_temp._ag.select('not chid Z and not hydrogen')
    pep = mol_temp._ag.select('chid Z not hydrogen')    

    near_n = findNeighbors(rec, 5 , pep)
    interacting_chainA_res = []
    if needs_b ==1:
        interacting_chainZ_res = []
    for a1,a2,d in near_n:
        # import pdb ; pdb.set_trace()
        interacting_chainA_res.append(a1.getResindex())
        if needs_b ==1:
            interacting_chainZ_res.append(a2.getResindex())
    
    interacting_chainA_res = list(set(interacting_chainA_res))    
    interacting_chainA_res.sort()
    # import pdb; pdb.set_trace()
    res_list_prody=[]
    for i in interacting_chainA_res:
        res_list_prody.append(rec.select(' name CA and resindex %d' % i ).getResnames()[0])
    # print("prody ignore:",res_list_prody)
    
    if needs_b ==1:
        interacting_chainZ_res = list(set(interacting_chainZ_res))
        interacting_chainZ_res.sort()
        pep_list_prody=[]
        for i in interacting_chainZ_res:
            pep_list_prody.append(pep.select(' name CA and resindex %d' % i ).getResnames()[0])
        # print("prody ignore(B):",pep_list_prody)
        return interacting_chainA_res,interacting_chainZ_res
   
    # interacting_chainA_res=np.array(interacting_chainA_res)
    return interacting_chainA_res


def AddListForUnknownNSTsToALA(fixer):
    '''Do not run fixer.findNonstandardResidues() 
    after running this command. You can run 
    fixer.findNonstandardResidues() before this.
    
    Also, do not use fixer.removeHeterogens() before 
    this command otherwise, all NSTs will be deleted 
    from the pdbdata.
    '''
    
    from pdbfixer.pdbfixer import (
        substitutions, proteinResidues,dnaResidues, rnaResidues)
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues).union(['N','UNK','HOH'])
    for residue in fixer.topology.residues():
        if residue.name in keep:
            continue
        if not residue.name in substitutions:
            if not hasattr(fixer, 'nonstandardResidues'):
                fixer.nonstandardResidues =[]
            fixer.nonstandardResidues.append([residue, 'ALA'])


def fix_my_pdb(pdb_in,out=None, NonstandardResidueTreatment=False): 
    ''' Options for NonstandardResidueTreatment:
      false: do nothing; 
     -fnst: try to swap NSTs with similar amino acids, if cannot swap, mutate to ALA'''

      
    if out==None:
        pdb_out = "./fixed/" + pdb_in.split('/')[-1].split('.')[0]+'_fixed.pdb'
    else:
        pdb_out = out

    # Cleaning in Prody >>>        
    mol_tmp =Read(pdb_in)
    res_names = mol_tmp._ag.getResnames()
    atom_names = mol_tmp._ag.getNames()
    element_names = mol_tmp._ag.getElements()
    for indx, (i,j) in enumerate(zip(res_names, atom_names)):
        if i == 'LEU':
            if j == 'CD':
                atom_names[indx]='CD1'
                element_names[indx]='C'
            elif j == 'CG':
                element_names[indx]='C'                
        elif i == 'ASN':
            if j == '2HD':
                atom_names[indx]='2HD2'
                element_names[indx]='H'
        elif i == 'ARG':
            if j == '2HN1':
                atom_names[indx]='2HH1'
                element_names[indx]='H'
        elif i == 'VAL':
            if j == 'CD1':
                atom_names[indx]='CG1'
                element_names[indx]='C'
            elif j == 'CD2':
                atom_names[indx]='CG2'
                element_names[indx]='C'
       
    mol_tmp._ag.setNames(atom_names)
    mol_tmp._ag.setElements(element_names)
    prody.writePDB(pdb_out,mol_tmp._ag)
    # Cleaning in Prody done <<<<<
    
    # Cleaning in PDBFixer >>>>      
    fixer = PDBFixer(filename=pdb_out)
    
    
    if NonstandardResidueTreatment:
        fixer.findNonstandardResidues()
        AddListForUnknownNSTsToALA(fixer)
        fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()    
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.4)
    
    # import pdb; pdb.set_trace()
    with open( pdb_out , 'w') as outfile:
         PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile,
    keepIds=True)
    return pdb_out

  
def restrain_(system_, pdb, not_chain='Z', ignore_list=[]):
    # given chain will not be restrained;
    # ignore_list is for atoms not to restrain: e.g. interface residues
    # BUG resolved 3/24/2023
    # strong restrain
    restraint = CustomExternalForce("0.5 * k * ((x-x0)^2 + (y-y0)^2 + (z-z0)^2)")    
    restraint.addGlobalParameter('k', 20000.0*kilojoules_per_mole/nanometer*nanometer)    
    restraint.addPerParticleParameter('x0')
    restraint.addPerParticleParameter('y0')
    restraint.addPerParticleParameter('z0')
        
    # iteration over all atoms to select atoms for restrain
    for atom in pdb.topology.atoms():        
        # ignore atoms from ignore_list for restrain
        if len(ignore_list) > 0:
            if atom.residue.index in ignore_list:
                continue        
        # ignore hydrogen atoms
        if atom.element.name == 'hydrogen':
            continue

        # ignore atoms from not_chain chain.
        if atom.residue.chain.id == not_chain:
            continue
        
        # else apply restrain on other atoms
        restraint.addParticle(atom.index, pdb.positions[atom.index])
         
    system_.addForce(restraint)
    

def get_energy(pdbfile,mode='vacuum', verbose = 0):
    pdb_handle = PDBFile(pdbfile)
    
    if mode == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        msg="Using GBSA (gbn2) environment for energy calculation"    
    else:
        force_field = ForceField("amber99sb.xml")
        msg="Using In-Vacuo environment for energy calculation"
        
    if verbose == 1:
        print(msg)
    system = force_field.createSystem(pdb_handle.topology, nonbondedCutoff=1*nanometer, constraints=HBonds)
    restrain_(system, pdb_handle, 'All')
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    simulation = Simulation(pdb_handle.topology, system, integrator)
    simulation.context.setPositions(pdb_handle.positions)
    simulation.minimizeEnergy(maxIterations=1)
    state = simulation.context.getState(getEnergy=True)
    # import pdb; pdb.set_trace()        
    return state.getPotentialEnergy()._value


def openmm_minimize( pdb_str: str, env='implicit', verbose = 0, max_itr=5):
    """Minimize energy via openmm."""
    pdb = PDBFile(pdb_str)
    
    if env == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')        
        msg = "Using GBSA (gbn2) environment for energy minimization"
    else:
        force_field = ForceField("amber99sb.xml")
        msg = "Using in-vacuo environment for energy minimization"
    
    if verbose == 1:
        print('max itr is %d and env is %s' % (max_itr,env))
        print(msg)
    # integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    constraints = HBonds    
    system = force_field.createSystem(  pdb.topology, nonbondedCutoff=1*nanometer, 
                                      constraints=constraints)    
    ignore_list = identify_interface_residues(pdb_str)        
    restrain_(system, pdb, ignore_list=ignore_list)
    
    # platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
    simulation = Simulation( pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)      
    return_energy = {}
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    return_energy["einit"] = state.getPotentialEnergy()
    # ret["posinit"] = state.getPositions(asNumpy=True)
    simulation.minimizeEnergy(maxIterations=max_itr, tolerance=0.01)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    return_energy["efinal"] = state.getPotentialEnergy()
    # ret["pos"] = state.getPositions(asNumpy=True)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(simulation.topology, positions, open(pdb_str[:-4]+"_min.pdb", 'w'))
    # ret["min_pdb"] = _get_pdb_string(simulation.topology, state.getPositions())
    return return_energy,system


def makeChidNonPeptide(receptor):
    # Chid 'Z' is reserved for peptide. Will be mutated if present in receptor.
    use_letters_for_chain = 'ABCDEFGHIJKLMNOPQRSTUVWXY'  # Z is saved for peptide
    rec_chids_list = receptor.getChids()
    rec_chids_str = ''.join(rec_chids_list)
    
    rec_chains = list(set(rec_chids_list))       
    restore_chain_ids = [] # it can be used
    
    # removing letters with already present chains, 
    # so no other chain will be renamed to already present ones.
    for chain in rec_chains:
        use_letters_for_chain = use_letters_for_chain.replace(chain,'')  
        
    chain_letter_iterator = 0
    for ids in rec_chains:
        # chain letters should not be numeric
        if ids.isnumeric():
            curr_letter = use_letters_for_chain[chain_letter_iterator]
            chain_letter_iterator += 1
            rec_chids_str = rec_chids_str.replace(ids,curr_letter)
            restore_chain_ids.append([ curr_letter, ids ])
        
        # receptor chain letters should not be 'Z' as it is asigned for peptide
        elif ids == 'Z':
            curr_letter = use_letters_for_chain[chain_letter_iterator]
            chain_letter_iterator += 1
            rec_chids_str = rec_chids_str.replace(ids,curr_letter)
            restore_chain_ids.append([ curr_letter, ids ])
    rec_chids_list = [x for x in rec_chids_str]
    receptor.setChids(rec_chids_list)
    return restore_chain_ids
    

class ommTreatment:
    def __init__(self,name_proj, file_name_init, rec_data=None):
        print(f'{Fore.GREEN}Rescoring clustered poses using OpenMM ...{Style.RESET_ALL}')
        #omm defaults
        self.minimization_env = 'in-vacuo'
        self.minimization_steps = 100
        self.find_neighbor_cutoff = 5        
        #other settings
        self.omm_dir = "omm_dir"
        self.omm_proj_dir = self.omm_dir + "/" + name_proj.split(".")[0]
        self.file_name_init = file_name_init
        self.pdb_and_score=[]
        self.delete_at_end =[]
        self.rearranged_data_as_per_asked=[]
        self.output_file =''
        self.rearrangeposes = True
        self.combined_pep_file =''
        self.rec_data = rec_data
        self.CLEAN_AT_END = 1 # for debugging only
     
    def __call__(self, **kw):       
        # reading all flags
        self.read_settings_from_flags(kw)     
        
        # creating required directories        
        self.create_omm_dirs()    
        
        # reading receptor file
        if self.rec_data == None:
            if kw['rec']:
                rec = Read(kw['rec'])                         
        else:
                rec = self.rec_data[0]
        rec = rec._ag # receptor        
        makeChidNonPeptide(rec) # removing numerical and 'Z' chains from receptor
        
        # identifying peptide clustered file and initiating output file name        
        if not os.path.isfile(self.combined_pep_file):
            print(f"{Fore.RED}Clustered file %s does not exist. OpenMM calculation terminated.{Style.RESET_ALL}" % self.combined_pep_file)
            return
        
        # combining peptides and receptor; and scoring   
        print('2')
        pep_pdb = Read( self.combined_pep_file )
        num_mode_to_minimize = min(pep_pdb._ag.numCoordsets(), int(kw['minimize']))        
        print(f"{Fore.GREEN}From total %d models, minimizing top %d ...{Style.RESET_ALL}" %
              (pep_pdb._ag.numCoordsets(), num_mode_to_minimize))
        
        for f in range( num_mode_to_minimize ):
            print( f"{Fore.GREEN}\nWorking on #%d of %d models.{Style.RESET_ALL}" % (f+1, num_mode_to_minimize))
            
            #peptide data
            pep_pdb._ag.setACSIndex(f)  
            pep_pdb._ag.setChids([x for x in 'Z'*pep_pdb._ag.numAtoms()]) # setting peptide chid to 'B'
            
            #making complex for current peptide and receptor
            out_complex_name = self.omm_proj_dir + "/" + self.file_name_init + "_recNpep_%d.pdb" % (f)  
            combinedRecPep= rec + pep_pdb._ag
            writePDB(out_complex_name, combinedRecPep.select('not hydrogen'))
            
            # energy calculation
            enzs = self.estimate_energies_for_pdb( out_complex_name)
            print ( "E_Complex = %9.2f; E_Receptor = %9.2f; E_Peptide  = %9.2f" % (enzs[0],enzs[1],enzs[2]))
            print ("dE_Interaction = %9.2f; dE_Complex-Receptor = %9.2f" % (enzs[3], enzs[4]))
            
            # save required details
            self.make_post_calculation_file_lists(out_complex_name, enzs)        
        self.calculate_reranking_index()    
        #self.print_reranked_models()
        self.read_pdb_files_and_combine_using_given_index()
        if self.CLEAN_AT_END == 1:
            self.clean_temp_files()
        
    def make_post_calculation_file_lists(self, flnm, enzs):
        # for deleting at end
        self.delete_at_end.append(flnm)
        self.delete_at_end.append(flnm[:-4] + "_fixed.pdb")
        self.delete_at_end.append(flnm[:-4] + "_fixed_min_A.pdb")
        self.delete_at_end.append(flnm[:-4] + "_fixed_min_Z.pdb")
        self.delete_at_end.append(flnm[:-4] + "_fixed_min.pdb")
        # combining for final output
        if len(self.pdb_and_score)<1:
            self.pdb_and_score.append([0,flnm[:-4]+"_fixed_min.pdb", enzs])
        else:
            current_model = len(self.pdb_and_score)
            self.pdb_and_score.append([current_model,flnm[:-4]+"_fixed_min.pdb", enzs])
            
        # self.pdb_and_score.append([flnm[:-4]+"_fixed_min.pdb", enzs])
        
    def read_settings_from_flags(self,kw):
        self.rearrangeposes = bool(kw['dockingRanking'])
        self.combined_pep_file = self.file_name_init + "_" + kw['sequence'] + "_out.pdb"
        self.output_file = self.file_name_init + "_" + kw['sequence'] + "_omm_rescored_out.pdb"
        self.minimization_env = 'implicit'  if kw['omm_environment'] == 'implicit' else  'in-vacuo'
        self.minimization_steps = kw['omm_max_itr']
        self.NonstandardResidueTreatment = kw['omm_nst']
        print(f'{Fore.GREEN}OpenMM minimization settings: Environment="%s"; Max_itr=%d.{Style.RESET_ALL}'
              % (self.minimization_env, self.minimization_steps))
        
    def calculate_reranking_index(self):  
        '''Rearraging BUG removed'''
        metric_comp_minus_rec = []  
        ommrank_adcprank_data = []
        for i in self.pdb_and_score:
            metric_comp_minus_rec.append(i[2][-1])          
        reaarange_for_omm_ranking = np.argsort(metric_comp_minus_rec)
        
        for current_v, rearrage_idx in enumerate(reaarange_for_omm_ranking):
            ommrank_adcprank_data.append([current_v, self.pdb_and_score[rearrage_idx][0],
                                          self.pdb_and_score[rearrage_idx][1],
                                          self.pdb_and_score[rearrage_idx][2]])
        ommrank_adcprank_data = np.array(ommrank_adcprank_data)
        
            
        # print(self.omm_ranking,metric_comp_minus_rec)
        if self.rearrangeposes == True:            
            print (f"{Fore.GREEN}\nREARRANGING output poses using OpenMM energy{Style.RESET_ALL}")          
            self.rearranged_data_as_per_asked = ommrank_adcprank_data
            
        else:
            print (f"{Fore.GREEN}\nNOT REARRANGING output poses using OpenMM energy{Style.RESET_ALL}")
            self.rearranged_data_as_per_asked = ommrank_adcprank_data[ommrank_adcprank_data[:,1].argsort()]
        
         
    def create_omm_dirs(self):        
        if not os.path.exists(self.omm_dir):
            os.mkdir(self.omm_dir)            
        if not os.path.exists(self.omm_proj_dir):
            os.mkdir(self.omm_proj_dir)
            
    # def print_reranked_models(self):      
    #     print(f"{Fore.GREEN}Re-ranking by OpenMM{Style.RESET_ALL} (E_complex - E_receptor):")
    #     # print ("-------+---------+------------------+---------+")
    #     # print ("model |rank omm | rank by AD score |  E_Comp |")
    #     # print ("      |         |                  |  -E_Rec |")
    #     # print ("------+---------+------------------+---------+")  
        
        
    #     # for i, rank_v in enumerate(self.omm_ranking):
    #     #     print(" %8d %18d %8.1f" % (i+1, rank_v+1, self.pdb_and_score[rank_v][1][-1]))
            
    def clean_temp_files(self):
        for fl in self.delete_at_end:
            if os.path.isfile(fl):
                os.remove(fl)
        if os.path.isdir(self.omm_proj_dir):
            os.rmdir(self.omm_proj_dir)
        
    def read_pdb_files_and_combine_using_given_index(self):
        out_file = open(self.output_file,"w+")        

        print ("-------+------+------+------------+")
        print (" Model | Rank | Rank | E_Complex  |")
        print (" #     |OpenMM| ADCP |-E_Receptor |")
        # print ("       |      |      |  (kJ/mol)  |")
        print ("-------+------+------+------------+")  
        
            
        for model_num, arranged_line in enumerate(self.rearranged_data_as_per_asked):
            # file name and omm scores
            omm_rank, adcp_rank, flnm, scores = arranged_line
            
            # flnm = self.pdb_and_score[rank][0]
            # scores = self.pdb_and_score[rank][1]
            score_string_to_write1 = "E_Complex = %9.2f; E_Receptor = %9.2f; E_Peptide  = %9.2f\n" % (scores[0],scores[1],scores[2])
            score_string_to_write2 = "dE_Interaction = %9.2f; dE_Complex-Receptor = %9.2f\n" % (scores[3], scores[4])           
                        
            out_file.write("MODEL     %4d\n" % ( model_num+1 ))
            out_file.write("USER: OpenMM RESCORED SOLUTION %d\n" % (omm_rank+1))            
            out_file.write("USER: " + score_string_to_write1)
            out_file.write("USER: " + score_string_to_write2)
            flnm_data = Read(flnm)            
            writePDBStream(out_file, flnm_data._ag)                       
            out_file.write("ENDMDL\n")
            print(" %6d %6d %6d %10.1f" %(model_num+1 , omm_rank+1, adcp_rank+1, scores[4] ))
            
        out_file.close()
        print ("-------+------+------+------------+")
        
    def estimate_energies_for_pdb(self,pdb_fl, verbose=0):
        env=self.minimization_env
        if verbose == 1:
            print('Working on:',pdb_fl.split("/")[-1])
        fixed_pdb = fix_my_pdb(pdb_fl, pdb_fl[:-4] +"_fixed.pdb",NonstandardResidueTreatment=self.NonstandardResidueTreatment)
        if verbose == 1:
            print("Minimizing ...")
        _= openmm_minimize(fixed_pdb,env, max_itr=self.minimization_steps)
        minimized_pdb = fixed_pdb[:-4]+"_min.pdb"
        enzs=[]
        split_pdb_to_chain_A_and_Z(fixed_pdb[:-4]+"_min.pdb")
        for j in [minimized_pdb, minimized_pdb[:-4]+"_A.pdb",minimized_pdb[:-4]+"_Z.pdb" ]:
            enzs.append(get_energy(j,env))
        enzs.append(enzs[0] - enzs[1] -enzs[2])   # interaction energy
        enzs.append(enzs[0] - enzs[1])            # complex - protein energy
        e_str=''
        for i in enzs:
            e_str = e_str + "%10.2f" % i
        if verbose == 1:
            print(' E_Complex E_Receptr E_Peptide dE_Interc dE_ComRes')   
            print (e_str)
        return enzs
          
# In the next version, this class will be used, by either using tempfile for 
# temprory file writing or, streamIO for keeping temp data in memory.
        
# class buffer_file_methods:
#     def __init__(self):
#         self.name = 'buffer_file_methods'
        
#     def writePDB_prody_to_buffer( self, atoms, csets=None, autoext=True, **kwargs):
#         """Write *atoms* in PDB format to a file with name *filename* and return
#         *filename*.  If *filename* ends with :file:`.gz`, a compressed file will
#         be written."""
#         output_buffer = io.StringIO()
#         writePDBStream(output_buffer, atoms, csets, **kwargs)
#         return output_buffer
                
#     def Read_buffer(self, buffer_name):
     