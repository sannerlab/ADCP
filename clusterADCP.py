import os,time
import sys
import numpy
import argparse
import dockq
from MolKit2 import Read
from MolKit2.PDBresidueNames import AAnames
from mglutil.math.rmsd import RMSDCalculator
from ADFR.utils.cluster import clusterPoses
from prody import writePDB, writePDBStream, AtomMap
from MolKit2.AARotamer import AARotamer, CanonicalAARotamers, AARotamerMutator
from prody.measure.contacts import findNeighbors

## def isPalindrom(str):
##     reversed_string = str[::-1]
##     status=1
##     if(str!=reversed_string):
##         status=0
##     return status

class stringBuffer:
    def __init__(self):
        self.lines = []
    def write(self, line):
        if line.startswith('ATOM'):
            line = line[:67] + 'REC' + line[70:][:-1] # remove new line
        self.lines.append(line)

#Read the energy and Rotamer information for each structures in the pdb
def getEnergy(filename):
    f = open(filename)
    lines = f.readlines()
    f.close()
    if len(lines)==0:
        return []
    totalE=[]
    locE=[]
    extE=[]
    rotamers=[]
    for ln, l in enumerate(lines):
        if not l.startswith('Energy ='):
            continue
        enes=l.split()
        if enes[2] == 'totalE':
            if float(enes[3]) < -1000000 or float(enes[3]) > 1000000:
                return totalE,locE,extE,rotamers
            totalE.append(float(enes[3]))
            locE.append(float(enes[6]))
            #extE.append(float(enes[8]))

            extE.append(0.25*float(enes[3])+0.75*float(enes[8]))
        else:
            totalE.append(float(enes[2]))
            locE.append(float(enes[4]))
            extE.append(float(enes[5]))
        if enes[11] == 'Rotamers:':
            rotstring=''
            for rot in enes[12:]:
                rotstring += str(rot)+" "
            rotamers.append(rotstring)
    return totalE,locE,extE,rotamers

#build the side chain and write out the structure in pdb format.
def buildSC(clusterid,rotamers,receptor=None,outputname=None):
    prot = Read('tmp.pdb')
    mutator = AARotamerMutator()
    i = 0;
    rots = rotamers.split()
    for res in prot._ag.iterResidues():
        resname = res.getResname()
        if resname == 'ALA' or resname == 'GLY':
            i=i+1
            continue
        chid = res.getChid()
        resnum = res.getResnum()
        mutator.mutate(prot, chid, resnum, res.getResname())
        res = prot.select('chid %s and resnum %d'%(chid, resnum))
        rotamer = AARotamer(res, mutator.rotamer.angleDef,
                            mutator.rotamer.angleList, prot)
        res.setCoords(rotamer.getCoordsForRotamer(int(rots[i])))
        i=i+1
    #writePDB("%s_ranked_%i.pdb"%(outputname,clusterid+1),prot.select("not deleted and not hydrogen"))
    if receptor is not None:
        pairs = findResPairs(prot._ag.select('not hydrogen and not deleted'),receptor._ag.select("not hydrogen"))
        return prot, pairs, True
    else:
        return prot, True


def findResPairs(refAtoms,receptor):
    pairs = set()
    for pair in findNeighbors(refAtoms,5,receptor):
        pairs.add(str(pair[0].getResnum())+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
    return pairs

def clusterPosesInteraction(neighbors, order, cutOff=0.8):
    # cluster a set of solutions in the AutoDock fashion, i.e. use the best
    # solution as a seed and add all solution within cutOff RMSD to this cluster
    # then re-seed the algorithm with the best energy solution not yet clustered
    # until all solution given in "order" are clustered
    #
    # the coordinates are in coords (nsol, natoms, 3)
    # order id a list of indices into coords indicating the subset of
    # solutions to be clustered. These indices are expected to point to
    # solutions sorted by decreasing GA scores
    #
    remainder = order[:]
    clusters = []
    while len(remainder):
        seed = remainder[0]
        #print '%d left to cluster seed=%d'%(len(remainder), seed), remainder
        #import pdb;pdb.set_trace()
        seedPairs = neighbors[seed]
        cluster = [seed]
        notSelected = []
        for i in remainder:
            if i==seed:
                continue
            currPairs = neighbors[i]
            #print '   %d %f'%(i, rmsd)
            overlap = compareContact(seedPairs,currPairs,JC=True)
            if overlap>=cutOff:
                cluster.append(i)
            else:
                notSelected.append(i)
            #print overlap
        #print 'found cluster', cluster, seed
        clusters.append( cluster )
        remainder = notSelected
    return clusters

def compareContact(pairs1,pairs2,JC=False):
    # if JC is true, return the Jaccard coefficient, else, return fraction of contact
    # len1 is the reference
    len1 = len(pairs1)
    len2 = len(pairs2)
    overlap=len(frozenset(pairs1).intersection(pairs2))
    if len1 == 0:
        return 0
    if JC:
        return (overlap+0.0)/(len1+len2-overlap)
    else:
        return (overlap+0.0)/len1

def getPDBlines(atoms, cset):
    buf = stringBuffer()
    writePDBStream(buf, atoms, cset)
    return buf.lines
    # code below no longer needed as adcp C code outputs side chains atoms
    # order corretly within their residues
    ## resnums = list(set(atoms.getResnums()))
    ## atomsInNewOrder = []
    ## atind = {}
    ## for n, atom in enumerate(atoms):
    ##     atind[atom.getIndex()] = n
    ## for resnum in resnums:
    ##     if resnum < 0:
    ##         resatoms = atoms.select('resnum `%d`'%resnum)
    ##     else:
    ##         resatoms = atoms.select('resnum %d'%resnum)        
    ##     for atom in resatoms:
    ##         atomsInNewOrder.append(atom)

    ## lines = []
    ## for n, atom in enumerate(atomsInNewOrder):
    ##     if atom.hydrogen: continue
    ##     line = buf.lines[atind[atom.getIndex()]+1]
    ##     lines.append(line[:6] + '%5d'%(n+1) + line[11:]) # replace serial number
    ## return lines

class clusterADCP:

    def __call__(self, **kw):

        ## read the files containg all the docking poses to be clustered
        ## this file was created in runADCP.py by concatenating the trajectories
        ## of all MC runs
        syst = os.path.join(kw['workingFolder'],kw['input'],kw['input'])+".pdb"

        print = kw['print']
        
        models = Read(syst)
        his = models._ag.select('resname HIE HID')
        if his:
            his.setResnames(['HIS']*len(his))
            
        seq = ''
        for res in models._ag.iterResidues():
            seq += AAnames.get(res.getResname(), '?')

        ## palindrom =  isPalindrom(seq)

        if kw['rec']:
            rec = Read(os.path.join(kw['targetFolder'], kw['rec']))
            
        ## set cutoff used for clustering based on contacts (default)
        ## or rmsd is specified
        if kw['rmsd']:
            clusterCutoff = kw['rmsd']
            rmsd = True
        else:
            if kw['rec'] is None:
                raise RuntimeError('ERROR: contact-based clustering required "-rec receptor.pdbqt".\nAlternatively you can, use -rmsd cutoff to use rmsd to cluster poses in %s.'%syst)
            else:                
                natContact = True
                rmsd = False
                receptor = Read(os.path.join(kw['targetFolder'],kw['rec']))
                clusterCutoff = 0.8

        moldelAtoms = models._ag.select('name C N CA')
        modelAtomIndices = moldelAtoms.getIndices()

        hasRef = False
        if kw['ref']:
            hasRef = True
            ref = Read(kw['ref'])
            if ref._ag.numResidues() == models._ag.numResidues():
                ## if the pdbqt has a torsion tree the atoms and residues are not in the right order
                ## sort the refAtoms by increasing residue numbers. In addition resnums in the docked
                ## model are 1-n for an n-peptide, while the reference structure might have other resnums
                ## we need to build an atom map for ref to match the atoms in models
                refresnums = list(set(ref._ag.getResnums()))
                refresnums.sort()
                modelresnums = list(set(models._ag.getResnums()))

                # renumber resnums for ref molecule to match docked poses in models
                lookup = {}
                for rn,mn in zip(refresnums, modelresnums):
                    lookup[rn] = mn

                newResnums = [lookup[rn] for rn in ref._ag.getResnums()]
                ref._ag.setResnums(newResnums)

                # build mapping of model atoms to reference atoms for heavy atoms
                refAtomIndexForModelAtomIndex = {}
                print('Reference Atom matching model atom')
                ignore = []
                for atom in models._ag:
                    if atom.ishydrogen: continue
                    atRef = ref._ag.select('resnum %d name %s'%(atom.getResnum(), atom.getName()))
                    if atRef is None:
                        print(" %s %3d %4s has no match in docked model. This atom will be ignored"%(
                            atom.getResname(), atom.getResnum(), atom.getName()))
                        ignore.append(atom)
                    elif len(atRef)==1:
                        refAtomIndexForModelAtomIndex[atom.getIndex()] = atRef.getIndices()[0]
                        print(" %s %3d %4s -> %s %3d:%4s"%(atRef.getResnames()[0], atRef.getResnums()[0], atRef.getNames()[0],
                                                           atom.getResname(), atom.getResnum(), atom.getName()))
                    elif len(atRef)>1:
                        print(" %s %3d %4s has %d matches in reference molecule, reference will be ignored"%(
                            atom.getResname(), atom.getResnum(), atom.getName(),len(atRef)))
                        hasRef = False
                        break

                # now create tom mappings for computing RMSD and DockQ values
                ignoreStr = ''
                for atom in ignore:
                    ignoreStr += ' and not (resnum %s anme %s)'%(atom.getResnum(), atom.getResname())

                ## atoms lists for ADCP ref. RMSD calculations
                moldelAtoms = models._ag.select('name C N CA %s'%ignoreStr)
                modelAtomIndices = moldelAtoms.getIndices()
                refAtomsIndices = [refAtomIndexForModelAtomIndex[x] for x in modelAtomIndices]
                refAtomsMap = AtomMap(ref._ag, refAtomsIndices)

                ## atom lists for DockQ calculations
                moldelAtomsDQ = models._ag.select('name C N CA O %s'%ignoreStr)
                modelAtomIndicesDQ = moldelAtomsDQ.getIndices()
                refAtomsIndicesDQ = [refAtomIndexForModelAtomIndex[x] for x in modelAtomIndicesDQ]
                refAtomsMapDQ = AtomMap(ref._ag, refAtomsIndicesDQ)

                # instantiate the RMSD calculator for clustering
                rmsdCalc = RMSDCalculator()

            else:
                print("WARNING: reference molecule has %d residues while docked model has %d, reference will be ignored"%(
                    ref._ag.numResidues(), models._ag.numResidues()))
                hasRef = False

        # establish order of docked poses
        totE, locE, extE, allrotamers = getEnergy(syst)        

        order=[]
        scores=[]
        bestE = min(extE)
        for i,ext in enumerate(extE):
            if (ext-bestE)<20:
                order.append(i)
                scores.append(ext)

        oorder = numpy.argsort(scores)
        order = numpy.array(order)[oorder]

        # open file for saving docked poses that are reported
        solutionFile = open('%s_out.pdb'%(os.path.join(kw['workingFolder'], kw['input']),), 'w')
        
        # maximum number of reported docked poses 
        maxModes = int(kw['nmodes'])
        if rmsd: # clustering based on RMSD values between poses
            print("Clustering MC trajectories based in backbone N,CA,C atoms RMSD using cutoff: %f"%clusterCutoff)

            rmsdCalc = RMSDCalculator(models._ag._coords[0][modelAtomIndices])
            clusters = clusterPoses(models._ag._coords, order, rmsdCalc, clusterCutoff)

            print("mode |  affinity  | clust. | ref. |  5A cutoff |     10 A cutoff    |   CAPRI   |")
            print("     | (kcal/mol) | size   | rmsd | fnat !fnat |  iRMS  LRMS  DockQ |   class   |")
            print("-----+------------+--------+------+------------+--------------------+-----------+")

            for i, cl in enumerate(clusters):
                if i >= maxModes:
                    break;

                solutionFile.write("MODEL     %4d\n"%(i+1,))
                solutionFile.write("USER: ADCP SOLUTION %d\n"%(i+1,))
                solutionFile.write("USER: SCORE %7.3f\n"%( extE[cl[0]] * 0.59219))
                lines = getPDBlines(models.select("not deleted"), cl[0])
                [solutionFile.write(l+'\n') for l in lines]
                solutionFile.write("ENDMDL\n")
                    
                if hasRef:
                    rmsdCalc.setRefCoords(refAtomsMap.getCoords())
                    rmsdRef = rmsdCalc.computeRMSD(models._ag._coords[cl[0]][modelAtomIndices])

                    ## calculate DockQ
                    ## ****************************************************************
                    ## *                       DockQ                                  *
                    ## *   Scoring function for protein-protein docking models        *
                    ## *   Statistics on CAPRI data:                                  *
                    ## *    0.00 <= DockQ <  0.23 - Incorrect                         *
                    ## *    0.23 <= DockQ <  0.49 - Acceptable quality                *
                    ## *    0.49 <= DockQ <  0.80 - Medium quality                    *
                    ## *            DockQ >= 0.80 - High quality                      *
                    ## *   Reference: Sankar Basu and Bjorn Wallner, DockQ: A quality *
                    ## *   measure for protein-protein docking models, submitted      *
                    ## *                                                              *
                    ## *   For the record:                                            *
                    ## *   Definition of contact <5A (Fnat)                           *
                    ## *   Definition of interface <10A all heavy atoms (iRMS)        *
                    ## *   For comments, please email: bjorn.wallner@.liu.se          *
                    ## *                                                              *
                    ## ****************************************************************
                    moldelAtomsDQ.setACSIndex(cl[0])
                    result = dockq.Fnat(rec._ag, moldelAtomsDQ, rec._ag, refAtomsMapDQ, 5.0)
                    nbPoseResPairsNat5, nbPoseResPairsNonNat5, nbNatResPairs5, interface5 = result
                    if not nbNatResPairs5 == 0:
                        fnat = nbPoseResPairsNat5/float(nbNatResPairs5)
                    else:
                        fnat =0
                    
                    if not nbPoseResPairsNonNat5==0:                        
                        fnonnat = nbPoseResPairsNonNat5 / float(nbPoseResPairsNat5+nbPoseResPairsNonNat5)
                    else:
                        fnonnat = 0

                    result = dockq.Fnat(rec._ag, moldelAtomsDQ, rec._ag, refAtomsMapDQ, 10.0)
                    nbPoseResPairsNat10, nbPoseResPairsNonNat10, nbNatResPairs10, interface10 = result
                    iRMS, LRMS, DockQ = dockq.calcDockQ(moldelAtomsDQ, refAtomsMapDQ, fnat, interface10)

                    print("%4d  %11.1f   %6d  %5.1f   %4.2f  %4.2f  %6.3f %6.3f %6.3f  %10s"%(
                        i+1, extE[cl[0]] * 0.59219, len(cl), rmsdRef, fnat, fnonnat,
                        iRMS, LRMS, DockQ, dockq.capri_class_DockQ(DockQ)))
                else:
                    print("%4d  %11.1f %6d      NA     NA    NA      NA    NA     NA        NA"%(
                        i+1, extE[cl[0]] * 0.59219, len(cl)))

        elif natContact:
            # FIXME needs to be made to work again and add DockQ calculation (MS)
            print("Clustering MC trajectories based in contacts using cutoff: %f"%clusterCutoff)
            cset = models._ag._coords.shape[0]
            neighbors = []
            start_time = time.time()

            for i in range(cset):       
                pairs = set()
                if i not in order:
                    neighbors.append(pairs)
                    continue
                models._ag.setACSIndex(i)
                for pair in findNeighbors(models._ag.select("name CB or (name CA and resname GLY)"),8,receptor._ag.select("name CB or (name CA and resname GLY)")):
                    #import pdb;pdb.set_trace()
                    pairs.add(str(pair[0].getResnum())+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
                neighbors.append(pairs)
            print("finished calculating neighbors for %d poses with %3.1f seconds"%(cset,time.time()-start_time))
            clusters=clusterPosesInteraction(neighbors,order,clusterCutoff)

            if hasRef:
                refpairs = findResPairs(refAtomsMap,receptor._ag.select("not hydrogen"))
                #refpairs = set()
                #for pair in findNeighbors(refAtoms,5,receptor._ag.select("not hydrogen")):
                #    refpairs.add(str(pair[0].getResnum()-residDiff)+'_'+str(pair[1].getResnum())+'_'+str(pair[1].getChid()))
            print("mode |  affinity  | ref. | clust. | rmsd | energy | best |")
            print("     | (kcal/mol) | fnc  |  size  | stdv |  stdv  | run  |")
            print("-----+------------+------+--------+------+--------+------+")

            bestNCRef = 0.
            bestNCInd = -1.
            NCRef = 0.
            for i, cl in enumerate(clusters):
                if i >= maxModes:
                    break;

                solutionFile.write("MODEL %d\n"%(i+1,))
                solutionFile.write("USER: ADCP SOLUTION %d\n"%(i+1,))
                solutionFile.write("USER: SCORE %7.3f\n"%( extE[cl[0]] * 0.59219))
                solutionFile.write("USER: ROTAMERS %s\n"% allrotamers[cl[0]])
                lines = getPDBlines(models.select("not deleted"), cl[0])# and not hydrogen"))
                [solutionFile.write(l+'\n') for l in lines]
                solutionFile.write("ENDMDL\n")
                if hasRef:
                    currpairs = findResPairs(models._ag.select('not hydrogen and not deleted'),receptor._ag.select("not hydrogen"))
                    NCRef = compareContact(refpairs,currpairs)
                    if NCRef > bestNCRef:
                        bestNCInd = cl[0]
                        bestNCRef = NCRef
                ene0 = extE[cl[0]]
                ene = [ene0]        
                print("%4d  %11.1f  %7.1f  %6d      NA      NA    %03d "%(i+1, ene0 * 0.59219, NCRef, len(cl),cl[0]))
            if hasRef:
                writePDB("%s_best_%3.2f.pdb"%(kw['input'][:-4],bestNCRef),models._ag,bestNCInd)

        solutionFile.close()

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='cluster CrankPep output', usage="usage: %(prog)s -i out.pdb -r ref.pdb -c 2.5",version="%prog 0.1")
    parser.add_argument("-i", "--input",dest="input")
    parser.add_argument("-rec", "--rec",dest="rec")
    parser.add_argument("-nc", "--natContacts", type=float, default=0.8,
                       dest="nc", help='native contacts cutoff used in the clustering')
    parser.add_argument("-rmsd", "--rmsd", type=float,
                       dest="rmsd", help='backbone rmsd cutoff used in the clustering')
    parser.add_argument("-m", "--nmodes", type=int, default=100,
                       dest="nmodes", help='maximum number of reported docked poses')
    parser.add_argument("-ref", "--ref",dest="ref", help='reference peptide structure for calculating rmsd and fnc')
    parser.add_argument("-w", "--workingFolder", default=None,
                       dest="workingFolder", help='folder in which the target file is expanded and the MC runs happened. Will be deleted after the run finished unless -k is specified')
    kw = vars(parser.parse_args())
    runner = clusterADCP()
    runner(**kw)

