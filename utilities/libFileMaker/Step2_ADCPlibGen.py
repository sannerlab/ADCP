#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed May 10 12:08:31 2023
@author: sshanker


(Step 2)
DESCRIPTION: This program reads temp rotamer description file (tempdir/tempRotamer.py) 
created by Step1_libFile2Rotamer.py program
and creates output ADCP acceptabel lib file for L and D configurations of each AA        

"""

import prody
import numpy as np
from math import sqrt, pi
from MolKit2 import Read
from prody.measure import calcDihedral
import os,sys
from datetime import datetime
from mglutil.math.rotax import rotax
import openbabel
import importlib
degtorad = pi/180.


def normalizedVector(a,b):
    'To get normal of np vectors. (Author MS)'
    v = np.array(b)-np.array(a)
    norm = np.linalg.norm(v)
    if norm == 0: 
       return v
    return v / norm

def vectorProduct(a, b):
    'to get product of vectors. (Author MS)'
    return np.array([a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]])

def orient(mol):
    'to translate and rotate coordinates to origin of Molkit2 mol object.  (Author MS)'    
    Ns = mol._ag.select('name N')
    assert len(Ns)==1
    N = Ns.getCoords()[0]
    n = mol._ag[Ns.getIndices()[0]]
    # CA is C arbon bonded to N
    found = False
    for ca in n.iterBonded():
        if ca.getElement()=='C':
            found = True
            break
    if not found: raise
    CA = ca.getCoords()
    # CB is carbon bonded to CA and with O as a neighbor
    found = False
    for cb in ca.iterBonded():
        if cb.getElement()=='C':
            ok = True
            for n in cb.iterBonded():
                if n.getName()=='O':
                    ok = False
                    break
            if ok:
                found = True
                break
    if not found: raise
    CB = cb.getCoords()
    # translate to origin
    cc = mol._ag.getCoords() - CA
    
    # build rotation to align with axis N on the -x axis, and CB in the y=0 plane 
    v1 = normalizedVector(N, CA)
    v3 = normalizedVector(CA, CB);

    v2 = vectorProduct(v3, v1)
    no = 1. / sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
    v2 *= no
    
    v3 = vectorProduct(v1, v2)
    no = 1. / sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);
    v3 *= no

    # build transformation matrix
    mat = np.identity(3)
    mat[:3, 0] = v1
    mat[:3, 1] = v2
    mat[:3, 2] = v3

    # apply to current coordinates
    tc = np.dot(cc, mat)
    mol._ag.setCoords(tc)
    #prody.writePDB('test.pdb', mol._ag)


def getAtom(mol, name):
    'to get atom object from mol object'
    ats = mol.select('name %s'%name)
    assert len(ats)==1
    return mol._ag[ats.getIndices()[0]]


def getRotamerCoordsForAngles(mol, chiDef, chiAngles):
    'rotate coordinates for given chi definitions and angles.  (Author MS)'
    origCoords = mol._ag.getCoords().copy()
    newCoords = origCoords.copy()
    # create array of transformed coordinates (initially orig)
    hcoords = np.ones( (len(origCoords), 4), 'd')
    hcoords[:, :3] = origCoords
    cmat = None
    for angNum, angle in enumerate(chiAngles):
        at1 = getAtom(mol,chiDef[angNum][0])
        at2 = getAtom(mol,chiDef[angNum][1])
        at3 = getAtom(mol,chiDef[angNum][2])
        at4 = getAtom(mol,chiDef[angNum][3])
        currentAngle = calcDihedral(at1,at2,at3,at4)
        moved = mol.subTree(at2,at3) #moved.getIndices()[2:]

        angleDelta = angle - currentAngle
        if angle < 0:
            angle += 360.
            angleDelta += 360.
        mat = rotax(origCoords[at2.getIndex()][:3], origCoords[at3.getIndex()][:3],
                    angleDelta*degtorad, transpose=1)
        # add mat to cmat
        if cmat is None:
            cmat = mat
        else:
            cmat = np.dot(mat, cmat)

        # transform atoms effected by Chi x
        for j in moved.getIndices()[:]:  ## it shold be all not [2:]
            newCoords[j] = np.dot( [hcoords[j]], cmat)[:, :3]
    mol._ag.addCoordset(newCoords)
    return newCoords

def give_ADtype_idx(atm):
    'to get AutoDock AtomTypes'
    AD_atypes = [ "C", "N", "OA", "HD", "SA", "A", "NA", "H", "HS", "NS",
    		      "NX","OS","OX", "F",  "Mg", "MG", "P", "S", "SX", "Cl",
    		      "CL","Ca","Mn","MN","Fe","FE","Zn","ZN","Br","BR","I"]
    
    if type(atm) == str:
    
        if not atm in AD_atypes:
            print("not present in AD_atypes list:", atm)
            return -1
        else:
            return(AD_atypes.index(atm))
    elif type(atm) == int:
        if atm > len(AD_atypes):
            return 'UNK'
        else:
            return AD_atypes[atm]

def name_of_D_aa_for_given_L(datfile ): 
    'dictionary with L and D names'
    fid = open(datfile,'r')
    data = fid.readlines()
    fid.close()
    
    pdb_name_data = {}
    for i in data:
        if i.startswith("#"):
            continue   
        if len(i.strip()) < 1:
            continue
        i_spl = i.split()
        # dihedral_list = [d.split("-") for d in i_spl[3:]]
        # base_aa = i_spl[1]
        d_name = i_spl[2]
        pdb_name_data[i_spl[0]] = d_name       
    return pdb_name_data



def get_atom_types_for_pdb_file(molf, temp_dir):  
    'to get Atom types for writing lib file' 
    if temp_dir:
        if not os.path.isdir(temp_dir):
            os.mkdir(temp_dir)
    
    counter = 0
    while 1:
        counter += 1
        tmpfile = os.path.join(temp_dir, "tmp_%d.pdbqt" % counter)   
        if not os.path.exists(tmpfile):
            break
        
    fid = open(molf,'r')
    mol2data = fid.readlines()
    fid.close()

    mol_use = []    
    for i in mol2data:
        if i.startswith('ATOM') or i.startswith('HETATM'):
            # print(i)
            atp = i[12:16]
            if len(atp.strip()) ==4:
                if atp.startswith('H'):
                    atp = atp[-1]+atp[:-1]
            
            if i.strip()[-2:].upper().strip() in ['CL', 'BR', 'ZN', 'FE', 'MN', 'MG', 'CA']: ## for four letter NST, NSTname and element can join in single word
                i2=i[12:16]
                if i[12:16].startswith(' '):
                    i2 = i[12:16].strip().ljust(4,' ')                    
                mol_use.append([i[12:16].strip(),i[12:16], i2, np.array(i[30:54].split()).astype(float)])
                    
            else:            
                mol_use.append([i[12:16].strip(),i[12:16], atp, np.array(i[30:54].split()).astype(float)])
            
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("pdb", "pdbqt")
    
    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, molf)   # Open Babel will uncompress automatically
    obConversion.WriteFile(mol, tmpfile)
    
    fid = open(tmpfile,'r')
    data = fid.readlines()
    fid.close()
    os.remove(tmpfile)    
    
    pdbqt_use =[]
    for i in data:
        if i.startswith('ATOM'):
            # atom_data.append(i)
            pdbqt_use.append([i.split()[-1],np.array(i[31:54].split()).astype(float),i.split()[-2] ])   
    
    atom_and_type = {}    
    for anxyz in mol_use:
        an = anxyz[0]
        an_extra = anxyz[1:3]
        # print(an_extra)
        xyz1= np.array(anxyz[-1])
        
        for atxyz in pdbqt_use:
            at = atxyz[0]
            xyz2= np.array(atxyz[-2])
            # print(xyz1-xyz2)
            if np.sqrt(sum((xyz1-xyz2)**2)) < 0.1:
                atom_and_type[an]=[an_extra[0], an_extra[1],at, give_ADtype_idx(at), float(atxyz[-1]) ]
                break
                
    return atom_and_type     

def header_for_out_file(time_str, sterio_type):
    header ='''//
// %s
// Author Michel Sanner and Sudhanshu Shanker
//
// This file provides data defining full atoms amino acid sidechains in the
// following format:
// The first line defining rotamers for a a side chain named XXX starts with
// the keyword rotamer followed by the sidechain name, the name of the coarse
// potential for this side chain or "o" if undefined, the  number of rotamers
// and the number of side chain atoms. These values are space separated.
//
// example:
// rotamer ARG R 81 11
// 
// following the rotamer declaration line we expect to find a space separated
// list of rotamer probabilities which sum up to 1.0
//
// example for ARG: (0:C 1:N 3:H)
// 0.012345679 0.012345679 0.012345679 .... 0.012345679
//
// following the rotamer probabilitesline we expect to find a space separated
// list of atoms types
//
// example for ARG: (0:C 1:N 3:H)
// C  N  C  N  N  C  HD  HD  HD  HD  HD
//
// following the atom types line we expect to find a space separated
// list of atomic partial charges (i.e. floating point values)
//
// example for ARG: (0:C 1:N 3:H)
// 0.138 -0.227  0.023 -0.235 -0.235  0.665  0.177  0.174  0.174  0.174  0.174
// 
// following the atomic charges line we expect to find a list of strings
// providing atom names. NOTE: Each atom name has to be 5 characters matching
// the PDB atom name column i.e. " CA  " for Calpha and "CA   " for calcium
//
// example for ARG:
//  CD   NE   CG   NH1  NH2  CZ   HE  1HH1 2HH1 1HH2 2HH2 
// 
// following the atom names  line we expect multiple lines providing
// atomic coordinates for all rotameric conformations i.e.
// nbRot * nbAtoms * 3 floating point values that can be on an arbitrary
// number of lines. 
//
// empty lines and lines staring with "/" are ignored
//
// the rotamers are in this file are for the sidechains of non standard amino 
// acids listed on the swisssidechain website
// https://www.swisssidechain.ch/browse/family/table.php?family=all
// for %s non-standard amino acids
//
// the rotamers are generated by placing the CA is at origin,
// N on the -x axis, and CB in the y=0 using the script 
// canonicalFAAStructGen.py by Michel Sanner 2017
//

''' 
    return header % (time_str,sterio_type)

def createRotamersForPDBfile(pdb_file, l_name, temp_dir, rotamer_data, swap_sign=False, out_name =None ):
    rotamer_out_dir = os.path.join(temp_dir, 'rotamers')
    
    mol = Read(pdb_file)
    # import pdb; pdb.set_trace()
    elemDict = get_atom_types_for_pdb_file(pdb_file, temp_dir)
    orient(mol)
    atypes = []
    anames = []
    charges = []
    indices = []
    coords = []
    rotaDescr = {}  
    for a in mol._ag.select('not hydrogen'):
        if a.getName() in ['N', 'CA', 'C', 'O', 'O2', 'CB','OXT']: 
            continue
        
        if not a.getName() in elemDict.keys():
            print("check for name %s" % ( a.getName()))
            print(elemDict.keys())
            continue
        
        else:
            atypes.append(elemDict[a.getName()][2])
            # anames.append(a.getName())
            anames.append(elemDict[a.getName()][1])
        
        charges.append(elemDict[a.getName()][4])
        indices.append(a.getIndex())
        
    for chiList in rotamer_data.chiValues[l_name]:
        if swap_sign: ## D to L or L to D
            chiList = [-1*i for i in chiList]
            ## change resnames as ltype for proper mol writing
            resname_array = [l_name]*len(mol._ag)
            mol._ag.setResnames(resname_array)
            
        coords.append(getRotamerCoordsForAngles(mol, rotamer_data.chiAnglesDef[l_name], chiList)[indices])
    
    class my_dot_var:
        pass    
    rotaDescr[l_name] = my_dot_var()
    rotaDescr[l_name].num_rotamers = len(rotamer_data.chiValues[l_name])
    rotaDescr[l_name].num_atoms = len(atypes)
    rotaDescr[l_name].atypes = atypes
    rotaDescr[l_name].charges = charges
    rotaDescr[l_name].coords = coords
    rotaDescr[l_name].atmnames = anames
    rotaDescr[l_name].base_name = rotamer_data.baseAA[l_name]
    rotaDescr[l_name].probs = rotamer_data.chiValueProbs[l_name]
    
    # prody.writePDB(os.path.join(rotamer_out_dir ,'%s_rota.pdb'%out_name), mol._ag)  #Uncomment it for analysis of rotamers
    return rotaDescr
    

def write_rotamer_lib_file(out_lib_file_obj, rotaDescr,  Dname = None):
    'to write the final output lib file'
    D_write=""    
    out_lib_file_obj    
    for resname, v in rotaDescr.items():        
          # standard amino acid rotamer base
        if Dname:
            out_lib_file_obj.write('// D_AminoAcidName %s\n' % Dname)
            D_write = "_D"
        out_lib_file_obj.write("rotamer %s %s %d %d\n"%(resname+D_write, v.base_name, v.num_rotamers, v.num_atoms))
        
        for pb in v.probs:
            out_lib_file_obj.write("%8.6f "%(pb/100))
        out_lib_file_obj.write("\n")      
        
        for att in v.atypes:
            out_lib_file_obj.write("%s "%att)
        out_lib_file_obj.write("\n")

        for ch in v.charges:
            out_lib_file_obj.write("%5.3f "%ch)
        out_lib_file_obj.write("\n")
        
        for anm in v.atmnames:
            out_lib_file_obj.write("%s "%anm)
        out_lib_file_obj.write("\n")
        
        for i in range(v.num_rotamers):
            for j in range(v.num_atoms):
                out_lib_file_obj.write(" %7.3f  %7.3f  %7.3f\n"%tuple(v.coords[i][j]))
        out_lib_file_obj.write("\n")


def run_libGen(**kw):
    AAdihedralInfoFile = kw['rotamerdetails']
    rotamer_py_file_name = kw['tempOutputFileName']
    pdbfileDir = kw['pdbfileDir']
    temp_dir = kw["tempdir"]
    # rotamer_py_file_name = os.path.join(temp_dir, outputfile)
    out_lib_file_name = kw['outputLibFileName']
    aa_and_dihedral_def = name_of_D_aa_for_given_L(AAdihedralInfoFile)
    prev_file_name = "Step1_libFile2Rotamer.py"
    if 'prevfilename' in kw.keys():
        prev_file_name = kw['prevfilename' ]
            
    if not os.path.exists(rotamer_py_file_name):
        print("Extracted rotamer file %s not found" % rotamer_py_file_name)
        print("Run %s to generate required rotamer intermediate file" % prev_file_name)
        return
    
    out_dir= os.path.split(out_lib_file_name)[0]
    
    if temp_dir:
        if not os.path.exists(temp_dir):
            os.makedirs(temp_dir)
        else:
            if not os.path.isdir(temp_dir):
                print ("A file with the name %d exist. Please change the tempdir name.")
                return
    
    sys.path.append(os.path.split(rotamer_py_file_name)[0])
    rotamer_data = importlib.import_module(os.path.split(rotamer_py_file_name)[1][:-3])
    rotamer_out_dir = os.path.join(temp_dir, 'rotamers')
    
    if not os.path.exists(rotamer_out_dir): # for rotamers
        os.mkdir(rotamer_out_dir)
       
    if out_dir:
        if not os.path.exists(out_dir): # for libfile
            os.makedirs(out_dir)
    
    print('Using residue list from "%s" for lib file generation!' % AAdihedralInfoFile)
    ## First do it for L
    names = list(rotamer_data.chiAnglesDef.keys())#[:3] # you can change here to run for specific AAs. we took all   
    out_lib_file_obj = open(out_lib_file_name, 'w+')        
    time_now = datetime.now().strftime("%B %d %Y, %H:%M:%S" )
    out_lib_file_obj.writelines(header_for_out_file(time_now, "L"))
    out_lib_file_obj.close()
    for name in names:
        l_name = name
        print("Working on AA: %s"%( name))
        respective_pdb_file = os.path.join(pdbfileDir, 'L' ,'%s.pdb'% (name))
        rotaDescr = createRotamersForPDBfile(respective_pdb_file, l_name, temp_dir, rotamer_data,out_name=l_name)
        out_lib_file_obj = open(out_lib_file_name, 'a') 
        write_rotamer_lib_file(out_lib_file_obj, rotaDescr)
        out_lib_file_obj.close() ## Close to save changes
    
    ## Now do it for D
    time_now = datetime.now().strftime("%B %d %Y, %H:%M:%S" )
    out_lib_file_obj = open(out_lib_file_name, 'a') 
    out_lib_file_obj.writelines(header_for_out_file(time_now, "D"))
    out_lib_file_obj.close()
    for name in names:
        l_name = name
        d_name = aa_and_dihedral_def[l_name]
        print("Working on AA: %s"%( d_name))
        respective_pdb_file = os.path.join(pdbfileDir, 'D', '%s.pdb'% (d_name))
        rotaDescr = createRotamersForPDBfile(respective_pdb_file, l_name, temp_dir, rotamer_data, swap_sign=True,out_name=d_name)
        out_lib_file_obj = open(out_lib_file_name, 'a') 
        write_rotamer_lib_file(out_lib_file_obj, rotaDescr,Dname=d_name)
        out_lib_file_obj.close()
    print("ADCP acceptable lib file: '%s' generated successfully" % out_lib_file_name )
          
    
if __name__ == "__main__":
    "Simple input type, modify it to use argument parser if needed"
    from settings import kw
    run_libGen(**kw)
   
