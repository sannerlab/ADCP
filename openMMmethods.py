#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 13:29:31 2022
@author: Sudhanshu Shanker

Collection of methods to perform openMM based calculation of ADCP docked poses.

Last Update 7/26/23 build 6
Build 6:
    1: fix_my_pdb loads openMM parameters (bond definition and hydrogen definition)
    (currently only swiss) to support NSTs.
    2: fix_my_pdb returns topology and positions (also) rather than writing a pdb file for faster calculation.
    3: -fnst flag only mutates NSTs that is not present in loaded ffxml file (currently only swiss).
    4: 
"""
import os
import numpy as np
#from colorama import Fore, Style

import prody
from prody.measure.contacts import findNeighbors
from prody import writePDB, writePDBStream
from MolKit2 import Read
from utils import currently_loaded_ffxml_data, loaded_ffxml_sets

from openmm import *
from openmm.app import *
from openmm.unit import *
from pdbfixer import PDBFixer
import copy
import parmed


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

# Global Variables
ffxml_path = os.path.join(os.path.dirname(__file__), 'data/openMMff') 
current_ffxml_set = loaded_ffxml_sets()[0] # currently set to only one


class my_dot_variable:
    # a dot type variable for complex data handling
    pass


# PDB fixing and NST treatment 

def AddListForUnknownNSTsToALA(fixer,kw=None):
    '''Idenitifes unknow NSTs and makes list to replace them by Alanine
    
    Do not run fixer.findNonstandardResidues() 
    after running this command. You can run 
    fixer.findNonstandardResidues() before this.
    
    Also, do not use fixer.removeHeterogens() before 
    this command otherwise, all NSTs will be deleted 
    from the pdbdata.
    '''    
    if not kw == None:
        current_ffxmls = currently_loaded_ffxml_data(kw) # to ignore loaded NSTS to mutate
    
    from pdbfixer.pdbfixer import (
        substitutions, proteinResidues,dnaResidues, rnaResidues)
    keep = set(proteinResidues).union(dnaResidues).union(rnaResidues).union(['N','UNK','HOH'])
    for residue in fixer.topology.residues():
        if residue.name in keep:
            continue
        
        if not residue.name in substitutions:
            if not hasattr(fixer, 'nonstandardResidues'):
                fixer.nonstandardResidues =[]
                
            if not residue.name in current_ffxmls.residues:  # to ignore loaded NSTS to mutate
                fixer.nonstandardResidues.append([residue, 'ALA'])


def fix_my_pdb(pdb_in,out=None, NonstandardResidueTreatment=False, return_topology= False,kw=None): 
    '''Uses pdbfixer to solve problems and fix nsts; and prody to identify 
      required atom arrangements and fix pdb errors. 
      Options for NonstandardResidueTreatment:
      false: do nothing; 
      -fnst: try to swap NSTs with similar amino acids, if cannot swap, mutate to ALA'''
     

    if out==None:
        pdb_out = "./fixed/" + pdb_in.split('/')[-1].split('.')[0]+'_fixed.pdb'
    else:
        pdb_out = out
        
    if return_topology:
        pdb_out = './tmp.pdb'

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
    
    #    # adding required bond data for NSTs
    # fixer.topology.loadBondDefinitions(os.path.join(ffxml_path, 'swissaa_bond_def.xml')) ## default in B6
    fixer.topology.loadBondDefinitions(os.path.join(ffxml_path, '%s_bond_def.xml' % current_ffxml_set))

    # loading hydrogen info for NSTs
    modeller =Modeller(fixer.topology,fixer.positions)
    modeller.loadHydrogenDefinitions(os.path.join(ffxml_path, '%s_hydrogen_def.xml' % current_ffxml_set)) ## default in B6
    
    ## Replacement Treatment for unidentified NSTs
    if NonstandardResidueTreatment:
        fixer.findNonstandardResidues()
        AddListForUnknownNSTsToALA(fixer,kw)
        fixer.replaceNonstandardResidues()
    fixer.findMissingResidues()
    fixer.findMissingAtoms()    
    fixer.addMissingAtoms(seed=2112)
    fixer.addMissingHydrogens(7.4)
    
    
    ## AS pdbfixer uses its modeller topology to add hydrogens, it works well on N-terminal residues
    ## but it can not fix the C-terminal residues So we need to do it here for NSTs
    
    pep_chain = list(fixer.topology.chains())[-1] ## As last chain is peptide
    
    ## As we will face problem with last NST in peptide 
    # This problem can happen with the last NST in receptors too, but currently we are ignoring this
    last_res_in_pep = list(pep_chain.residues())[-1] 
    
    atomNames = list(atom.name for atom in last_res_in_pep.atoms())

    C_atom_idx = [atom.index for atom in last_res_in_pep.atoms() if atom.name == 'C'][0]
    O_atom_idx = [atom.index for atom in last_res_in_pep.atoms() if atom.name == 'O'][0]
    CA_atom_idx = [atom.index for atom in last_res_in_pep.atoms() if atom.name == 'CA'][0]      
    C_atom = list(last_res_in_pep.atoms())[atomNames.index('C')]

    if not 'OXT' in atomNames: # if OXT is not added
        from openmm.app.element import oxygen
        import openmm.unit as unit
    
        OXT_atom = fixer.topology.addAtom('OXT', oxygen, last_res_in_pep)
        fixer.topology.addBond(OXT_atom, C_atom)
        unit_v = fixer.positions.unit
        
        #       OXT
        #      /
        #CA---C
        #      \
        #       O
        # To get OXT position, we will use unit vectors C-->CA and C-->O, 
        # combining these vectors and inversing the direction will give unit vector for OXT position.
        
        d_c_ca = fixer.positions[CA_atom_idx].value_in_unit(unit_v) - fixer.positions[C_atom_idx].value_in_unit(unit_v)         
        v_c_ca = d_c_ca/unit.sqrt(unit.dot(d_c_ca, d_c_ca))        # unit vector C--> CA
                
        d_c_o = fixer.positions[O_atom_idx].value_in_unit(unit_v) - fixer.positions[C_atom_idx].value_in_unit(unit_v) 
        len_d_c_o = unit.sqrt(unit.dot(d_c_o, d_c_o))        
        v_c_o = d_c_o/ len_d_c_o # unit vector C--> O
        
        v = -(v_c_o + v_c_ca)/2   # vector in direction of OXT
        v =v/unit.sqrt(unit.dot(v, v)) # unit vector for OXT
        v *= len_d_c_o    # vector with C--O distance for OXT when C at origin
        
        c_val = fixer.positions[C_atom_idx].value_in_unit(unit_v) # Positon of 'C'
        v += c_val  # C dependent position of OXT        
        fixer.positions.append(v*unit_v) #As last residue and last atom, we do not need to check anything else
    
        
       
    if return_topology:
        return fixer.topology, fixer.positions
    
    # import pdb; pdb.set_trace()
    with open( pdb_out , 'w') as outfile:
         PDBFile.writeFile(fixer.topology, fixer.positions, file=outfile,
    keepIds=True)
    return pdb_out

  
## Interface identification and restrain     
def identify_interface_residues(pdb_file, needs_b=0):
    '''identifies interface residues between macromolecule(s) and peptide'''
    
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


def restrain_(system_, pdb, not_chain='Z', ignore_list=[]):
    '''To create a positional restrain for non-interface residue to avoid their
    movements during minimization'''
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

def cyclize_backbone_of_last_chain(topology, positions, myprint=print, display_info = False):
    '''To make bond between terminal residues of peptide if -cyc flag is given'''
    #cyclizes the backbone of last chain : peptide
    modeller = Modeller(topology, positions)
    nat = cat = None

    chains = list(topology.chains())
    # if not 'Z' in [i.id for i in chains]:
    #     print("NO peptide", [i.id for i in chains])
    #     return modeller
    
    residues = list(chains[-1].residues())
    # find N of first residue
    
    todelete = []
    for atom in residues[0].atoms():
        if atom.name =='N':
            nat = atom
        if atom.name =='H2':
            todelete.append(atom)
        if atom.name =='H3':
            todelete.append(atom)

    for atom in residues[-1].atoms():
        if atom.name=='C':
            cat = atom
        if atom.name =='OXT':
            todelete.append(atom)

    # make sure bond (nat, cat) does not exist
    # I beleive the bonds are built using templates, not based on distance
    # so this will always returns FALSE
     
    hasBond = False
    for b in residues[0].bonds():
        if (b.atom1==nat and b.atom2==cat) or (b.atom1==cat and
                                               b.atom2==nat):
            hasBond = True
            break

    # fixme if atoms are deleted before creating cyclic bond
    # it fails to build the system
    if not hasBond:
        if display_info:
            n_detail = nat.residue.name + nat.residue.id + "@" + nat.name #+"." +nat.residue.chain.id
            c_detail = cat.residue.name + cat.residue.id+ "@" + cat.name#+"." +cat.residue.chain.id
            myprint("Cyclizing peptide (last chain) by adding a bond between %s and %s."%(n_detail, c_detail))
        modeller.topology.addBond(nat, cat, 'single')
        # print ('deleting', todelete)
        modeller.delete(todelete)
    else:
        if display_info:
            myprint("cyclization bond was already present")
    
    return modeller

   
def make_disulfide_between_cys_in_last_chain(top, pos, myprint=print, debug = False):    
    '''Makes disulfide bridged for interacting cysteins.
    ;it is different than openmm default method as it useu 3.5A cutoff while
    openmm default is 3.0A'''
    # print("Inside check")
    class Cys_data:
        def __init__(self):
            self.HG = None
            self.SG = None
            self.SG_xyz = None
            self.index = None
        def __call__(self):
            return self.index, self.SG_xyz, self.SG, self.HG
            
    def isDisulfideBonded(atom, topology):
      for b in topology._bonds:
          if (atom in b and b[0].name == 'SG' and
              b[1].name == 'SG'):
              return True

      return False
    
    modeller = Modeller(top, pos)
    topology = modeller.topology
    positions = modeller.positions
    
    chains = list(topology.chains())
    # print("Inside cystein SSbond")
    
    # To make a disulfide link, openMM identifies CYS 
    # as CYX if HG atom is absent, therefore
    # delete HG to make CYS -> CYX    
        
    # Make a list of identified CYS in peptide
    # and all required data, if disulfide can be added    
    found_cys = [] #[itr_number, xyz of SG, and atom HG]
    SS_bond_cutoff = 3.5    # default openMM cutoff 3
    
    count_pos = 0
    for r in chains[-1].residues():
        if r.name == 'CYS':
            curr_cys_data = Cys_data()
            for a in r.atoms():
                # print(a.name)
                if a.name == 'SG':
                    p = positions[a.index]
                    xyz = np.array([p.x,p.y,p.z]) 
                    
                    curr_cys_data.index = count_pos
                    curr_cys_data.SG = a
                    curr_cys_data.SG_xyz = xyz

                    count_pos +=1
                elif a.name == 'HG':
                    curr_cys_data.HG = a 

            found_cys.append(curr_cys_data)
            # print(curr_cys_array)
    if positions.unit.get_symbol() == 'nm':
        nm_cutoff = SS_bond_cutoff/10
    else:
        nm_cutoff = SS_bond_cutoff
        
    
    # identify nearby CYS pairs which can make disulfide.    
    del_HG = [] # list of HG to delete to make CYS -> CYX (required for SSbond)
    make_bond_pair = []
    if len(found_cys) > 1:
        #print("more than 1 cysteins are present")    
        for cys1 in found_cys[:-1]:
            pos1, xyz1, SG1, HG1 = cys1()
            for cys2 in found_cys[pos1+1:]:
                pos2, xyz2, SG2, HG2 = cys2()
                dist_r1_r2 = sum((xyz1-xyz2)**2)**(1/2)  
                if dist_r1_r2 < nm_cutoff: # Check if CYS can make disulfide
                    hg1_detail = SG1.residue.name + SG1.residue.id + "@" + SG1.name #+"." +nat.residue.chain.id
                    hg2_detail = SG2.residue.name + SG2.residue.id+ "@" + SG2.name#+"." +cat.residue.chain.id
                    
                    if not ( isDisulfideBonded(SG1, topology) or  isDisulfideBonded(SG2, topology)):                    
                        myprint("Marking for disulfide link in peptide (last chain) between residue %s and residue %s." % 
                              (hg1_detail, hg2_detail))
    
                        make_bond_pair.append([SG1,SG2])
                    else: 
                        myprint("Disulfide present in peptide (last chain) between residue %s and residue %s." % 
                              (hg1_detail, hg2_detail))
                        
                    if not HG1 =='':
                        del_HG.append(HG1)
                    if not HG2 =='':
                        del_HG.append(HG2)
    
    if debug:
        #checking output
        PDBFile.writeFile(modeller.topology, modeller.positions, open('output1.pdb', 'w')) 
      
    for SGi,SGj in make_bond_pair:
        sg1_detail = SGi.residue.name + SGi.residue.id + "@" + SGi.name #+"." +nat.residue.chain.id
        sg2_detail = SGj.residue.name + SGj.residue.id+ "@" + SGj.name#+"." +cat.residue.chain.id
        myprint("Making disulfide link in peptide (last chain) between residue %s and residue %s." % 
              (sg1_detail, sg2_detail))
        modeller.topology.addBond(SGi, SGj, 'single')
    modeller.delete(del_HG) 
      
    if debug:
    #checking output
        PDBFile.writeFile(modeller.topology, modeller.positions, open('output2.pdb', 'w'))
    
    else:
        return modeller.topology, modeller.positions

 
def get_energy_on_modeller(modeller, mode='vacuum', verbose = 0):
    ''' Rather than using temp pdb files this method 
    takes modeller object as input to calculates potential energy.    
    '''
   
    if mode == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')
        msg="Using GBSA (gbn2) environment for energy calculation"    
    else:
        # force_field = ForceField("amber99sb.xml", os.path.join(ffxml_path, 'swissaa_ff.xml'))
        force_field = ForceField("amber99sb.xml", os.path.join(ffxml_path, '%s_ff.xml' % current_ffxml_set))
        msg="Using In-Vacuo environment for energy calculation"
    positions = modeller.positions
    topology = modeller.topology
    if verbose == 1:
        self.myprint(msg)
    system = force_field.createSystem( topology, nonbondedCutoff=1*nanometer, constraints=HBonds)
    # restrain_(system, pdb_handle, 'All')
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    simulation = Simulation( topology, system, integrator)
    simulation.context.setPositions(positions)
    simulation.minimizeEnergy(maxIterations=1)
    state = simulation.context.getState(getEnergy=True)
    # import pdb; pdb.set_trace()        
    return( state.getPotentialEnergy()._value)


def openmm_minimize( pdb_str: str, env='implicit', verbose = 0, max_itr=5, cyclize = False, SSbond=False, amberparminit=None, myprint=print):
    '''Minimize a protein-peptide complex using openMM'''
    #"""Minimize energy via openmm."""
    pdb = PDBFile(pdb_str)
    out_var = my_dot_variable()
    # for cyclic peptides
    if cyclize:
        modeller = cyclize_backbone_of_last_chain(pdb.topology, pdb.positions, myprint=myprint, display_info=True)
        topology = modeller.topology
        positions = modeller.positions
    else:
        topology = pdb.topology
        positions = pdb.positions
    
    # for disulfide bridges
    if SSbond:
        topology, positions = make_disulfide_between_cys_in_last_chain(topology, positions, myprint=myprint) 
    
    if env == 'implicit':
        force_field = ForceField("amber99sb.xml",'implicit/gbn2.xml')        
        msg = "Using GBSA (gbn2) environment for energy minimization"
    else:
        force_field = ForceField("amber99sb.xml", os.path.join(ffxml_path, '%s_ff.xml' % current_ffxml_set))
        msg = "Using in-vacuo environment for energy minimization"
    
    if verbose == 1:
        self.myprint('max itr is %d and env is %s' % (max_itr,env))
        self.myprint(msg)
    # integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    integrator = openmm.LangevinIntegrator(0, 0.01, 0.0)
    # constraints = HBonds    
    system = force_field.createSystem(  topology, nonbondedCutoff=1*nanometer)
                                      # constraints=constraints)    
    ignore_list = identify_interface_residues(pdb_str)        
    restrain_(system, pdb, ignore_list=ignore_list)
    
    
    # platform = openmm.Platform.getPlatformByName("CUDA" if use_gpu else "CPU")
    simulation = Simulation( topology, system, integrator)
    simulation.context.setPositions(positions)      
    return_energy = {}
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    return_energy["einit"] = state.getPotentialEnergy()
    # ret["posinit"] = state.getPositions(asNumpy=True)
    simulation.minimizeEnergy(maxIterations=max_itr, tolerance=0.01)
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    return_energy["efinal"] = state.getPotentialEnergy()
    # ret["pos"] = state.getPositions(asNumpy=True)
    positions = simulation.context.getState(getPositions=True).getPositions()
    PDBFile.writeFile(topology, positions, open(pdb_str[:-4]+"_min.pdb", 'w'))
    # ret["min_pdb"] = _get_pdb_string(simulation.topology, state.getPositions())
    
    out_var.minimized_energy = return_energy
    out_var.system = system
    out_var.topology = topology
    out_var.positions = positions
    out_var.env = env
    
    if not amberparminit == None:
        structure = parmed.openmm.topsystem.load_topology( topology, system, positions)
        # print(pdb_str)
        structure.save(amberparminit+'.parm7', overwrite=True)
        structure.save(amberparminit+'.rst7' , format='rst7', overwrite=True)

    # with open('system.xml', 'w') as output:
    #     output.write(XmlSerializer.serialize(system))

    return out_var


def return_comp_rec_pep_energies_from_omm_minimize_output(omm_min_in):
    '''for faster E_consensus and E_interface calulation it uses modeller object
    to calculate energy values'''
    
    topology = omm_min_in.topology
    positions = omm_min_in.positions
    env = omm_min_in.env
    
    chains = list(topology.chains())
    
    # complex single point energy
    modeller_comp = Modeller(topology, positions)
    comp_e = get_energy_on_modeller(modeller_comp, mode=env)        
    
    # receptor single point energy
    modeller_comp.delete([chains[-1]])
    rec_e = get_energy_on_modeller(modeller_comp, mode=env)
    
    # peptide single point energy
    modeller_pep = Modeller(topology, positions)
    modeller_pep.delete(chains[:-1])
    pep_e = get_energy_on_modeller(modeller_pep, mode=env)
    
    return [comp_e, rec_e, pep_e]


def makeChidNonPeptide(receptor):
    ''' to avoid chain-id related mis-calculations this code changes peptide chain
    to 'Z' (reserved). "Z" Will be mutated to another chain if present in receptor.'''
    
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
    

def read_and_write_pdb_data_to_fid(fid, pdbfile):
    '''as name says'''
    # for combining output pdbs 
    # prody gives error with modified residues like cys -> cyx
    fid2= open(pdbfile, 'r')
    data = fid2.readlines()
    data.remove('END\n')
    fid2.close()
    # print(data)
    fid.writelines(data)
    


class ommTreatment:
    #a class to perform all required operations for openmm calculations
    def __init__(self,target_file, jobName, myprint, rec_data=None):
        myprint('Rescoring clustered poses using OpenMM ..')
        #omm defaults
        self.minimization_env = 'in-vacuo'
        self.minimization_steps = 100
        self.find_neighbor_cutoff = 5        
        #other settings
        self.omm_temp_dir = "omm_temp_dir_" + jobName
        self.omm_proj_dir = self.omm_temp_dir + "/" + target_file.split(".")[0]
        self.amber_parm_out_dir = "%s_omm_amber_parm/" % jobName + target_file.split(".")[0]
        self.file_name_init = jobName
        self.pdb_and_score=[]
        self.delete_at_end =[]
        self.rearranged_data_as_per_asked=[]
        self.output_file =''
        self.rearrangeposes = True
        self.combined_pep_file =''
        self.rec_data = rec_data
        self.cyclize_by_backbone = False
        self.SSbond = False
        self.CLEAN_AT_END = True # for debugging only
        self.myprint = myprint
        
     
    def __call__(self, **kw):       
        # reading all flags
        self.kw = kw
        self.read_settings_from_flags(kw)     
        
        # creating required directories        
        self.create_omm_dirs()    
        
        # reading receptor file
        if self.rec_data == None:
            if kw['recpath']:
                rec = Read(kw['recpath'])                         
        else:
                rec = self.rec_data[0]
        rec = rec._ag # receptor
        makeChidNonPeptide(rec) # removing numerical and 'Z' chains from receptor
        
        # identifying peptide clustered file and initiating output file name        
        if not os.path.isfile(self.combined_pep_file):
            self.myprint("Clustered file %s does not exist. OpenMM calculation terminated" % self.combined_pep_file)
            return
        
        # combining peptides and receptor; and scoring   
        pep_pdb = Read( self.combined_pep_file )
        num_mode_to_minimize = min(pep_pdb._ag.numCoordsets(), int(kw['minimize']))        
        self.myprint("From total %d models, minimizing top %d ..." %
              (pep_pdb._ag.numCoordsets(), num_mode_to_minimize))
        
        for f in range( num_mode_to_minimize ):
            self.myprint( "\nWorking on #%d of %d models" % (f+1, num_mode_to_minimize))
            
            #peptide data
            pep_pdb._ag.setACSIndex(f)  
            pep_pdb._ag.setChids([x for x in 'Z'*pep_pdb._ag.numAtoms()]) # setting peptide chid to 'B'
            
            #making complex for current peptide and receptor
            out_complex_name = self.omm_proj_dir + "/" + self.file_name_init + "_recNpep_%d.pdb" % (f)  
            out_parm_dirname = self.amber_parm_out_dir + "/"+ self.file_name_init + "_%d" % f
            combinedRecPep= rec + pep_pdb._ag
            writePDB(out_complex_name, combinedRecPep.select('not hydrogen'))
            
            # energy calculation
            enzs = self.estimate_energies_for_pdb( out_complex_name, out_parm_dirname)
            self.myprint( "E_Complex = %9.2f; E_Receptor = %9.2f; E_Peptide  = %9.2f" % (enzs[0],enzs[1],enzs[2]))
            self.myprint("dE_Interaction = %9.2f; dE_Complex-Receptor = %9.2f" % (enzs[3], enzs[4]))
            
            # save required details
            self.make_post_calculation_file_lists(out_complex_name, enzs)        
        self.calculate_reranking_index()    
        #self.print_reranked_models()
        self.read_pdb_files_and_combine_using_given_index()
        if self.CLEAN_AT_END == True:
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
        # reading options
        self.rearrangeposes = bool(kw['dockingRanking'])
        # self.combined_pep_file = self.file_name_init + "_" + kw['sequence'] + "_out.pdb"
        self.combined_pep_file = self.file_name_init + "_out.pdb"
        # self.output_file = self.file_name_init + "_" + kw['sequence'] + "_omm_rescored_out.pdb"
        self.output_file = self.file_name_init + "_omm_rescored_out.pdb"
        self.minimization_env = 'implicit'  if kw['omm_environment'] == 'implicit' else  'in-vacuo'
        self.minimization_steps = kw['omm_max_itr']
        self.NonstandardResidueTreatment = kw['omm_nst']
        self.cyclize_by_backbone = kw['cyclic']
        self.SSbond = kw['cystein']
        self.myprint('OpenMM minimization settings: Environment="%s"; Max_itr=%d'
              % (self.minimization_env, self.minimization_steps))
        
        if self.cyclize_by_backbone:
            self.myprint('OpenMM minimization settings: Cyclizing by backbone')
        if self.SSbond:
            self.myprint('OpenMM minimization settings: Making Disulfide bonds if applicable')
            
        
    def calculate_reranking_index(self):  
        #for rearraging poses
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
        #ommrank_adcprank_data = np.array(ommrank_adcprank_data)
        adcp_ranks = np.array([ i[1] for i in ommrank_adcprank_data])
            
        # print(self.omm_ranking,metric_comp_minus_rec)
        if self.rearrangeposes == True:            
            self.myprint("\nREARRANGING output poses using OpenMM energy")          
            self.rearranged_data_as_per_asked = ommrank_adcprank_data
            
        else:
            self.myprint("\nNOT REARRANGING output poses using OpenMM energy")
            self.rearranged_data_as_per_asked = [ ommrank_adcprank_data[i] for i in adcp_ranks.argsort()]
        
         
    def create_omm_dirs(self):   
        # create required output directories
        if not os.path.exists(self.omm_temp_dir):
            os.mkdir(self.omm_temp_dir)            
        if not os.path.exists(self.omm_proj_dir):
            os.mkdir(self.omm_proj_dir)
        #for AMBER parameters
        if not os.path.exists(self.amber_parm_out_dir):
            os.makedirs(self.amber_parm_out_dir, exist_ok=True)
            
            
    # def print_reranked_models(self):      
    #     print(f"{Fore.GREEN}Re-ranking by OpenMM{Style.RESET_ALL} (E_complex - E_receptor):")
    #     # print ("-------+---------+------------------+---------+")
    #     # print ("model |rank omm | rank by AD score |  E_Comp |")
    #     # print ("      |         |                  |  -E_Rec |")
    #     # print ("------+---------+------------------+---------+")  
        
        
    #     # for i, rank_v in enumerate(self.omm_ranking):
    #     #     print(" %8d %18d %8.1f" % (i+1, rank_v+1, self.pdb_and_score[rank_v][1][-1]))
            
    def clean_temp_files(self):
        # works as function name says
        for fl in self.delete_at_end:
            if os.path.isfile(fl):
                os.remove(fl)
        if os.path.isdir(self.omm_proj_dir):
            os.rmdir(self.omm_proj_dir)
        if os.path.isdir(self.omm_temp_dir):
            os.rmdir(self.omm_temp_dir)
        
    def read_pdb_files_and_combine_using_given_index(self):
        # works as function name says
        out_file = open(self.output_file,"w+")        

        self.myprint("-------+------+------+------------+")
        self.myprint (" Model | Rank | Rank | E_Complex  |")
        self.myprint(" #     |OpenMM| ADCP |-E_Receptor |")
        # print ("       |      |      |  (kJ/mol)  |")
        self.myprint("-------+------+------+------------+")  
        
            
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
            # flnm_data = Read(flnm)            
            # writePDBStream(out_file, flnm_data._ag.select('not hydrogen')) 
            read_and_write_pdb_data_to_fid(out_file, flnm)
            out_file.write("ENDMDL\n")
            self.myprint(" %6d %6d %6d %10.1f" %(model_num+1 , omm_rank+1, adcp_rank+1, scores[4] ))
            
        out_file.close()
        self.myprint("-------+------+------+------------+")
        
    def estimate_energies_for_pdb(self,pdb_fl, amber_parm_init, verbose=0):
        # works as function name says
        env=self.minimization_env
        if verbose == 1:
            self.myprint('Working on:',pdb_fl.split("/")[-1])
        fixed_pdb = fix_my_pdb(pdb_fl, pdb_fl[:-4] +"_fixed.pdb",NonstandardResidueTreatment=self.NonstandardResidueTreatment,kw=self.kw)
        if verbose == 1:
            self.myprint("Minimizing ...")
        
        omm_min_list = openmm_minimize(fixed_pdb,env, max_itr=self.minimization_steps, 
                           cyclize=self.cyclize_by_backbone, SSbond=self.SSbond, amberparminit = amber_parm_init, myprint=self.myprint)
        
        enzs = return_comp_rec_pep_energies_from_omm_minimize_output(omm_min_list)
        
        enzs.append(enzs[0] - enzs[1] -enzs[2])   # interaction energy
        enzs.append(enzs[0] - enzs[1])            # complex - protein energy
        e_str=''
        for i in enzs:
            e_str = e_str + "%10.2f" % i
        if verbose == 1:
            self.myprint(' E_Complex E_Receptr E_Peptide dE_Interc dE_ComRes')   
            self.myprint (e_str)

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
     
