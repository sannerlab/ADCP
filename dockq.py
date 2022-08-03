###############################################################################
##
## This library is free software; you can redistribute it and/or
## modify it under the terms of the GNU Lesser General Public
## License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## 
## This library is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## Lesser General Public License for more details.
## 
## You should have received a copy of the GNU Lesser General Public
## License along with this library; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA
##
## (C) Copyrights Dr. Michel F. Sanner and TSRI 2019
##
################################################################################

#############################################################################
#
# Author: Michel F. SANNER
#
# Copyright: M. Sanner and TSRI 2022
#
#########################################################################
#
# relies on ProDy from MolKit2 to read and operate on molecules
# 
from MolKit2 import Read
from prody.measure.contacts import findNeighbors

## we assume atoms in pose and ref are ordered the same way

def Fnat(receptor, ligand, receptorRef, ligandRef, cutoff):
    """calculate the fraction of native and on native contact from the
    reference complex recapitulated by a given complex (receptor, ligand)
    for a given cutoff distance.

    return #natC, #nonNatC, #refC
    where:
         #natC is the number of contacting residue pairs correctly recapitulated 
         #nonNatC is the number of additional contacting residue pairs
         #refC the number of contacting pairs in the reference complex
    """
    intd = {}
    interface = []
    ## find pairs of residues with heavy atoms within cutoff Angrstroms
    ## in the reference complex
    #native = findNeighbors(receptorRef.select('not hydrogen'), cutoff, ligandRef.select('not hydrogen'))
    native = findNeighbors(receptorRef, cutoff, ligandRef)
    natResPairs = {}
    for a1,a2,d in native:
        key = '%d_%d'%(a1.getResnum(), a2.getResnum())
        # Resindex is sensitive to gaps and insertions
        #key = '%d_%d'%(a1.getResindex(), a2.getResindex()) 
        natResPairs[key] = True
        if intd.get(key, None) is None:
            interface.append((a1.getResnum(), a2.getResnum()))
            intd[key] = True
            
    ## find pairs of residues with heavy atoms within cutoff Angrstroms in the
    ## complex.
    #poseC = findNeighbors(receptor.select('not hydrogen'), cutoff, ligand.select('not hydrogen'))
    poseC = findNeighbors(receptor, cutoff, ligand)
    poseResPairsNat = {}
    poseResPairsNonNat = {}
    for a1,a2,d in poseC:
        key = '%d_%d'%(a1.getResnum(), a2.getResnum()) 
        # Resindex is sensitive to gaps and insertions
        #key = '%d_%d'%(a1.getResindex(), a2.getResindex()) 
        if natResPairs.get(key, None) is None:
            # this pair is NOT in the reference complexes so it is non native
            poseResPairsNonNat[key] = True
        else:
            # this pair is in the reference complexes so it is native
            poseResPairsNat[key] = True

    #print('Fnat %d %d %.5f'%(len(poseResPairsNat), len(natResPairs), len(poseResPairsNat)/len(natResPairs), ))
    #total = len(poseResPairsNonNat)+len(poseResPairsNat)
    #print('Fnonnat %d %d %.5f'%(len(poseResPairsNonNat), total, len(poseResPairsNonNat)/total))
    return len(poseResPairsNat), len(poseResPairsNonNat), len(natResPairs), interface

def capri_class_DockQ(DockQ):
    # function implemented by Bjorn Wallner see
    # https://github.com/bjornwallner/DockQ
    (c1,c2,c3)=(0.23,0.49,0.80)
    if(DockQ < c1):
        return 'Incorrect'
    elif(DockQ >= c1 and DockQ < c2):
        return 'Acceptable'
    elif(DockQ >= c2 and DockQ < c3):
        return 'Medium'
    elif(DockQ >= c3):
        return 'High'
    else:
        return 'Undef'

def calcDockQ(ligD, ligC, fnat, inter):
    # ligand interface backbone atoms in reference
    refIntBB = ligC.select('resnum %s and backbone'%(' '.join([str(x[1]) for x in inter])))

    # ligand interface backbone atoms in pose
    poseIntBB = ligD.select('resnum %s and backbone'%(' '.join([str(x[1]) for x in inter])))

    # compute iRMS
    from mglutil.math.rmsd import RMSDCalculator
    calc = RMSDCalculator(refIntBB.getCoords())
    iRMS = calc.computeRMSD(poseIntBB.getCoords())
    
    # compute LRMS
    calc.setRefCoords(ligC.select('bb').getCoords())
    LRMS = calc.computeRMSD(ligD.select('bb').getCoords())

    DockQ = (float(fnat) + 1/(1+(iRMS/1.5)*(iRMS/1.5)) + 1/(1+(LRMS/8.5)*(LRMS/8.5)))/3
    return iRMS, LRMS, DockQ
    
if __name__=='__main__':
    import sys
    recD = Read(sys.argv[1])
    ligD = Read(sys.argv[2])
    recC = Read(sys.argv[3])
    ligC = Read(sys.argv[4])
    cutoff = float(sys.argv[5])

    a, b, c, inter = Fnat(recD, ligD, recC, ligC, cutoff)
    fnat = a/c
    print('Fnat %d %d %.6f'%(a, c, fnat))
    print('Fnonnat %d %d %.6f'%(b, a+b, b/(a+b)))

    a, b, c, inter10 = Fnat(recD, ligD, recC, ligC, 2*cutoff)

    ## WARNING: for docking the receptor is fixed so we bypass aligning
    ## the receptor using the interface atoms
    iRMS, LRMS, DockQ = calcDockQ(ligD, ligC, fnat, inter10)
    print("iRMS %f"%iRMS)
    print("LRMS %f"%LRMS)
    print("DockQ %f %s"%(DockQ, capri_class_DockQ(DockQ)))
