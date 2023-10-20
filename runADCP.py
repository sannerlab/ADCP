################################################################################
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
# Copyright: M. Sanner and TSRI 2019
#
#########################################################################
#

import os, sys, numpy, shutil, random
from glob import glob
from time import time, sleep
from MolKit2 import Read
from utils import openmm_validator, add_open_mm_flags, support_validator,  BUILD_VERSION #OMM new line
#import os, sys, numpy, platform, datetime, tempfile, shutil, random, tarfile, pickle
# from glob import glob


def pepSeqStr(case, sequence):
    """this function will transform sequence into a peptide sequence for a
    helical initial confomation with upper case amino acid letters or extended
    with lower case letters but leaving non standard amino acids (<XXX> unchanged)"""
    newSeq = []
    inNST = False
    for char in sequence:
        if char is '<': inNST = True
        if inNST:
            newSeq.append(char)
        else:
            if case is 'upper':
                newSeq.append(char.upper())
            else:
                newSeq.append(char.lower())
        if char is '>': inNST = False
    return ''.join(newSeq)

class runADCP:

    def sanitize(self, name):
        name = name.replace('&', '\&')
        name = name.replace('<', '\<')
        name = name.replace('>', '\>')
        return name
    
    def myprint(self, txt, newline=True):
        if self.summaryFile:
            self.summaryFile.write(txt)
            if newline:
                self.summaryFile.write('\n')
            
        sys.stdout.write(txt)
        if newline:
            sys.stdout.write('\n')

    def myexit(self, error=True):
        if not self.keepWokingFolder and self.calcFolder:
            print("removing working folder %s"%self.calcFolder)
            try:
                shutil.rmtree(self.calcFolder)
            except OSError:
                pass        
        if self.summaryFile:
            self.summaryFile.close()
        if error:
            sys.exit()
            
    def __init__(self):

        import multiprocessing
        self.ncpu = multiprocessing.cpu_count()

        self.completedJobs = 0    
        self.numberOfJobs = 0
        self.outputBaseName = None # a folder with that name will be crated to store log files and ligands
                                   # a _summary.dlg file will be create with this name too
                                   # is specified using -o on the command line
        self.jobName = 'NoName'
        self.targetFile = None
        self.cleanup = True
        self.summaryFile = None
        self.workingFolder = None
        self.calcFolder = None

    def __call__(self, **kw):
        #
        # run ADFR GAs using the list of command line arguments from the sysargv list
        # 
        # to check different flag requirements (currently added openMM)

        import ADCP
        path = os.path.dirname(ADCP.__file__)

        import platform, subprocess
        system_info = platform.uname()
        _platform = system_info[0]

        self.keepWokingFolder = kw.pop('keepWorkingFolder')
        # find the binary
        if _platform == "Linux":
            binary = os.path.join(path, "CrankiteAD_Linux-x86_64_1.1")
        elif _platform == "Darwin":
            binary = os.path.join(path, "CrankiteAD_Darwin_1.1")
        else:
            binary = os.path.join(path, "CrankiteAD_Win-x86_64_1.1")

        assert binary is not None, "ERROR: binary for platform %s not found in path %s"%(_platform, path)
        assert os.path.exists(binary) is True, "ERROR: binary %s does not exist." % binary
        #print "ADCPpath", ADCPpath, "cwd:", os.getcwd(), "binary" , binary, "_platform", _platform
        if _platform == 'Windows':
            self.shell=False
        else:
            self.shell=True

        if openmm_validator(kw) == False:                                       #OMM new line
            return                                                              #OMM new line

        import subprocess, datetime
        dataDict = {}
        
        # make path to reference molecule absolute
        startFolder = os.path.abspath(os.getcwd())
        if kw['ref']:
            kw['ref'] = os.path.abspath(kw['ref'])
            if os.path.isfile(kw['ref']):
                try:
                    ref = Read(kw['ref'])
                except Exception as e:
                    print('ERROR: failed to load reference structure: %s, %s'%(kw['ref'], str(e)))

        seed = None
        rncpu= None
        nbRuns = 50
        numSteps = 2500000
        jobName = 'NoName'
        partition = 0
        skip = False

        self.targetFile = targetFile = kw.pop('target')
        if targetFile is None:
            print("ERROR: no receptor files found, please provide a .trg file or path to inflated .trg")
            self.myexit()

        # find receptor name
        if os.path.isdir(targetFile):
            l = targetFile.split(os.sep)
            # find lst non mepty string because '/a/b/c/'.split('/') gives ['a', 'b', 'c', '']
            i = len(l)-1
            while i>-1 and len(l[i])==0: i-=1
            targetName = l[i]
        else:
            targetName = os.path.splitext(os.path.basename(targetFile))[0]

        jobName = kw['jobName']
        if jobName is None:
            jobName = '%s_%s'%(targetName, self.sanitize(kw['sequence']))

        # create a working folder name after the jobname where we will perform calculations
        self.workingFolder = workingFolder = kw['workingFolder']

        # Post Docking Minimization step:                                       #OMM new line           
        if kw['postdockmin']:                                                   #OMM new line            
            from utils import ( evaluate_requirements_for_minimization,         #OMM new line
                               extract_target_file )                            #OMM new line
            kw['jobName'] = jobName                                             #OMM new line
            kw['target'] = targetFile                                           #OMM new line
            
            if not extract_target_file(kw, workingFolder, jobName):             #OMM new line
                self.myexit()                                                   #OMM new line
            
            # check everything needed is there for a safe & successful minimization
            if not evaluate_requirements_for_minimization(kw):                  #OMM new line              
                self.myexit()                                                   #OMM new line

            from openMMmethods import ommTreatment                              #OMM new line       
            runner_omm = ommTreatment(targetFile, kw['recpath'], workingFolder, jobName, PostDockMinCommands=sys.argv)          #OMM new line
            runner_omm(**kw)            
            self.myexit()                                                       #OMM new line
        
        if workingFolder == None:
                print('ERROR: Working folder is not defined. Working folder is required for temporary data generation. '+
                      'Use -w to define working folder.' )
                self.myexit()       
        
        if not os.path.exists(workingFolder):
            try:
                os.mkdir(workingFolder)
            except:
                raise RuntimeError("ERROR: uanble to create folder %s"%workingFolder)
                self.myexit()

        self.calcFolder = calcFolder = os.path.join(workingFolder, jobName)
        if not os.path.exists(calcFolder):
            os.mkdir(calcFolder)

        #check overwriting files          
        if os.path.exists('%s_out.pdb'%os.path.join(workingFolder, jobName)):
            if not kw['overwriteFiles']:
                print("ERROR: file %s_out.pdb exists! please use different jobName (-o) or use -O to force overwritting output files"%jobName)
                self.myexit()
        if os.path.exists('%s_summary.dlg'%os.path.join(workingFolder, jobName)):
            if not kw['overwriteFiles']:
                print("ERROR: file %s_summary.dlg exists! please use different jobName (-o) or use -O to force overwritting output files"%os.path.join(workingFolder, jobName))
                self.myexit()
        if os.path.exists(jobName):
            if not kw['overwriteFiles']:
                print("ERROR: file/folder %s exists! please use different jobName (-o) or use -O to force overwritting output files"%os.path.join(workingFolder, jobName))
                self.myexit()

        # setup targetFolder
        if not os.path.exists(targetFile):
            raise RuntimeError('ERROR: target file "%s" does not exist, fix -T option'%targetFile)

        if os.path.isdir(targetFile):
            target_folder = targetFile
            if not os.path.exists(os.path.join(targetFile, 'transpoints')):
                # create text version of translational points in target folder
                if os.path.exists(os.path.join(targetFile, 'translationPoints.npy')):
                    transPoints = numpy.load(os.path.join(targetFile, 'translationPoints.npy'))
                    TPfile = open(os.path.join(targetFile, 'transpoints'),'w')
                    TPfile.write('%s\n'%len(transPoints))
                    numpy.savetxt(TPfile,transPoints,fmt='%7.3f')
                    TPfile.close()

        elif os.path.isfile(targetFile):
            # if target is zip file unzip and replace cmdline arguments
            import zipfile
            self.myprint( 'Inflating target file %s'%(targetFile))
            with zipfile.ZipFile(targetFile, 'r') as zip_ref:
                # if the trg file was renamed it might create a folder with a different name
                # so we first find out the name of the folder in which the maps are
                filenames = zip_ref.namelist()
                folder = filenames[0].split(os.sep)[0]
                zip_ref.extractall(calcFolder)
            target_folder = os.path.join(calcFolder, folder)
            transPoints = numpy.load(os.path.join(target_folder, 'translationPoints.npy'))
            TPfile = open(os.path.join(target_folder, 'transpoints'),'w')
            TPfile.write('%s\n'%len(transPoints))
            numpy.savetxt(TPfile,transPoints,fmt='%7.3f')
            TPfile.close()

        # open summary file
        summaryFilename = '%s_summary.dlg'%os.path.join(self.workingFolder, jobName)
        if kw['reclusterOnly']:
            t0 = time()
            if not os.path.exists(summaryFilename):
                self.myprint('ERROR: summary file %s not found for reclustering'%summaryFilename)
                sys.exit(1)
            runFilenames = glob(os.path.join(calcFolder, 'run_*.out'))
            if len(runFilenames)==0:
                self.myprint('ERROR: run files run_%s.out files %s not found in folder %s for reclustering'%(summaryFilename, calcFolder))
                sys.exit(1)
            # read best target energies from run_%d.out files
            filenamesNums = [int(os.path.splitext(os.path.basename(fn))[0][4:]) for fn in runFilenames]
            order = numpy.argsort(filenamesNums)
            runEnergies = []
            for nn, i in enumerate(order):
                f = open(runFilenames[i])
                lines = f.readlines()
                f.close()
                for ln in lines:
                    if ln.startswith('best target energy'):
                        runEnergies.append(float(ln.rstrip().split()[3]))
                        #print(nn, runEnergies[nn], runFilenames[i])
                        break

            self.summaryFile = open('%s_summary.dlg'%os.path.join(self.workingFolder, jobName), 'a')

        else: # run docking
            self.summaryFile = open('%s_summary.dlg'%os.path.join(self.workingFolder, jobName), 'w')
            self.myprint("performing MC searches with: %s "%binary)
            if kw['ref']:
                self.myprint('reference srtucture: %s'%kw['ref'])
            self.myprint('target data from file: %s'%targetFile)
            self.myprint('job name: {0}, summary file {0}_summary.dlg, docked poses: {0}_out.pdb'.format(jobName))

            rncpu = kw.pop('maxCores')

            if rncpu is None:
                ncores = self.ncpu
                self.myprint( 'Detected %d cores, using %d cores.'%(self.ncpu, ncores))
            else:
                assert rncpu > 0, "ERROR: maxCores a positive number, got %d."%rncpu
                ncores = min(self.ncpu, rncpu)
                self.myprint( 'Detected %d cores, request %d cores, using %d cores.'%(self.ncpu, rncpu, ncores))

            if kw['nbRuns'] is not None:
                self.nbRuns = nbRuns = kw.pop('nbRuns')

            self.numberOfJobs = nbRuns
            self._jobStatus =  [None]*nbRuns        

            seed = kw.pop('seedValue')

            if seed is None:
                seed = str(random.randint(1,999999))

            # to check openmm supports for the residues in the receptor file
            kw['recpath'] = os.path.join(target_folder, 'rigidReceptor.pdbqt')      #OMM new line        
            if not support_validator(kw,self.myprint)[0]:                           #OMM new line
                self.myexit()                                                       #OMM new line
                return                                                              #OMM new line

            if not os.path.exists(os.path.join(target_folder, 'constrains')):
                fff = open(os.path.join(target_folder, 'constrains'),'w')
                fff.write('1\n')
                fff.close()

            self.dryRun = kw.pop('dryRun')

            seeds = [] # save seeds for runs to save in summary file
            # build cmdline args for adcp binary
            argv = ['cd %s; %s -t 2'%(calcFolder, binary)]

            if kw['sequence'] is None:
                if kw['input'] is None or kw['input'][-3:] != 'pdb':
                    self.myprint ("ERROR: no input for peptide found")
                    self.myexit()
                else:
                    argv.append('-f')
                    argv.append('%s'%kw['input'])
            else:
                if kw['partition'] is None:
                    argv.append('"%s"'%kw['sequence'])
                else:
                    partition = kw['partition']
                    argv.append('"%s"'%kw['sequence'])
                    if partition < 0:
                        partition = 0
                    elif partition > 100:
                        partition = 100

            if kw['rotlibs']:
                argv.append('-L %s'%kw['rotlibs'])
                self.myprint( 'using system rotamer libraries %s from %s'%(
                    kw['rotlibs'], os.path.join(os.path.dirname(binary), 'data', 'rotamers')))

            if kw['userrotlibs']:
                userrotlibswithpath = []
                for userrotlib in kw['userrotlibs'].split(':'):
                    userrotlibswithpath.append(os.path.abspath(userrotlib))
                    self.myprint( 'using user rotamer library %s'%userrotlib)
                print(":".join(userrotlibswithpath))
                argv.append('-l %s'% ":".join(userrotlibswithpath))      
                # argv.append('-l %s'%os.path.abspath(kw['userrotlibs']))
                # for word in kw['userrotlibs'].split(':'):
                #     self.myprint( 'using user rotamer library %s'%word)

            argv.append('-T %s'%os.path.abspath(target_folder))

            # set up the length for each run, 25M as default
            argv.append('-r')
            if kw['numSteps'] is not None:
                numSteps = kw['numSteps']        
            argv.append('1x%s'%numSteps)

            # set up other options for ADCP
            ADCPDefaultOptions = "-p Bias=NULL,external=5,constrains,1.0,1.0"
            if kw['cyclic']:
                ADCPDefaultOptions += ",external2=4,constrains,1.0,1.0"
            if kw['cystein']:
                ADCPDefaultOptions += ",SSbond=50,2.2,50,0.35"
            ADCPDefaultOptions += ",Opt=1,0.25,0.75,0.0"
            argv.append(ADCPDefaultOptions)

            # add arguments that will be set during the loop submitting jobs
            # for seed jubNum and outputName
            argv.extend(['-s', '-1', '-o', jobName,' '])
            jobs = {} # key will be process until process.poll() is not None (i.e. finished)

            t0 = time()
            runStatus = [None]*(nbRuns)

            runEnergies = [999.]*(nbRuns)
            # procToRun = {}
            # can't change dict during iteration so use 2 lists
            # outfiles = {}
            processes = [None]*(nbRuns)
            outfiles = [None]*(nbRuns)

            nbStart = 0 # number of started runs
            nbDone = 0 # number of completed runs
            # print(" ".join(argv))
            if seed == -1:
                self.myprint( 'Performing %d MC searches using %d evals each using a random seed.'%(self.nbRuns, numSteps))
            else:
                self.myprint( 'Performing %d MC searches using %d evals each using seed %g.'%(self.nbRuns, numSteps, seed))

            self.myprint( "Performing search (%d ADCP runs with %d steps each) ..."%
                          (nbRuns, numSteps))
            self.myprint ("0%   10   20   30   40   50   60   70   80   90   100%")
            self.myprint ("|----|----|----|----|----|----|----|----|----|----|")


            numHelix = nbRuns * partition / 100

            ## submit the first set of jobs
            for jobNum in range(1,min(nbRuns,ncores)+1):
                # overwrite seed
                if seed == -1:
                    argv[-4] = str(random.randint(1,999999))
                else:
                    argv[-4] = str(seed+jobNum-1)
                # overwrite jobNum
                argv[-2] = 'run_%d.pdb'%(jobNum)
                # since we set shell=False on Windows, we cannot redirect output to a file
                # like in the following commented out line. We will open a file "outfile"
                # and set stdout to the file handle.

                #argv[-1] = '> %s_%d.out 2>&1'%(jobName,jobNum)

                # overwrite the sequence if parition is found
                if partition > 0 and partition < 100 and kw['sequence'] is not None:
                    ## MS Aug 15 2023. Here we need to be careful because we can not
                    ## change the capitalization of NSTs
                    if jobNum <= numHelix:
                        #argv[1] = '"%s"'%kw['sequence'].upper()
                        argv[1] = pepSeqStr('upper', '"%s"'%kw['sequence'])
                    else:
                        #argv[1] = '"%s"'%kw['sequence'].lower()
                        argv[1] = pepSeqStr('lower', '"%s"'%kw['sequence'])
                if self.dryRun:
                    self.myprint ('/n*************** command ***************************\n')
                    self.myprint (' '.join(argv))
                    if kw['minimize']:                                                      #OMM new line
                        if openmm_validator(kw) == True:
                            from openMMmethods import openMMdryrunchecks
                            kw['target'] = targetFile   
                            openMMdryrunchecks(kw, self.myprint)
                    self.myexit()
                    sys.exit()
                elif jobNum==1:
                    command = ' '.join(argv)

                outfile =  open(os.path.join(calcFolder, 'run_%d.out'%(jobNum)), 'w')

                process = subprocess.Popen(' '.join(argv),
                                           stdout=outfile,#subprocess.PIPE , 
                                           stderr=outfile, #subprocess.PIPE, 
                                           bufsize = 1, shell=self.shell ,  cwd=os.getcwd())

                #procToRun[process] = jobNum-1
                #print('FUFU1', process, jobNum-1)
                #outfiles[jobNum] = outfile # process.stdout returns None. So save the file handle in the dictionary, so that we can close the file after the job finishes
                processes[jobNum-1]= process
                outfiles[jobNum-1] = outfile
                nbStart += 1
                #print("%3d %s %s"%(jobNum-1, process, outfile))

            # check for completion and start new runs until we are done
            while nbDone < nbRuns:
                # check running processes            
                #for proc, jnum in procToRun.items():
                for jnum, proc in enumerate(processes):
                    if proc is None: continue
                    #import pdb;pdb.set_trace()
                    if proc.poll() is not None: # process finished
                        if runStatus[jnum] is not None: continue
                        if proc.returncode != 0:
                            #print('FAFA1', proc, jnum, proc.returncode)
                            runStatus[jnum] = ('Error', '%s%04d'%(jobName, jnum+1))
                            error = '\n'.join(runStatus[jnum][1])
                            status = 'FAILED'
                            self.myprint( '%d ENDED WITH ERROR'%(jnum,))
                            print ('%d err'%jnum)
                        else:
                            #print('FAFA2', proc, jnum, proc.returncode)
                            status = 'OK'
                            error = ''
                            runStatus[jnum] = ('OKAY', '%s%04d'%(jobName, jnum+1))
                            #self.myprint( jnum, 'ENDED OK')
                            #print '%d ok'%jnum
                            #import pdb;pdb.set_trace()
                            #f = open('%s_%d.out'%(jobName,jnum+1))
                            f = outfiles[jnum] # the file should still be open
                            if not f.closed:
                                f.close()
                            f = open(os.path.join(calcFolder, 'run_%d.out'%(jnum+1)))
                            lines = f.readlines()
                            #print "Lines in %s_%d.out %d"%(jobName,jnum+1, len(lines))
                            f.close()
                            for ln in lines:
                                if ln.startswith('best target energy'):
                                    runEnergies[jnum] = float(ln.rstrip().split()[3])

                        nbDone += 1
                        self._jobStatus[jnum] = 2
                        self.completedJobs += 1
                        percent = float(self.completedJobs)/self.numberOfJobs
                        #print('FUGU percent %f %d %d %d %d'%(percent, self.completedJobs, self.numberOfJobs, int(50*percent), jnum))
                        sys.stdout.write('%s\r' % ('*'*int(50*percent)))
                        sys.stdout.flush()

                        if nbStart < nbRuns:
                            # start new one
                            jobNum += 1
                            if seed == -1:
                                argv[-4] = str(random.randint(1,999999))
                            else:
                                argv[-4] = str(seed+jobNum-1)
                            seeds.append(argv[-4])
                            # overwrite jobNum
                            argv[-2] = 'run_%d.pdb'%(jobNum)
                            #argv[-1] = '> %s_%d.out 2>&1'%(jobName,jobNum)
                            outfileName =  os.path.join(calcFolder, 'run_%d.out'%(jobNum))
                            # remove output file in case it exists
                            #try:
                            #    os.remove(argv[-1])
                            #except OSError:
                            #    pass
                            if os.path.exists(outfileName):
                                f = outfiles[jobNum-1]
                                if f and not f.closed:
                                    f.close()
                                os.remove(outfileName)
                            outfile = open(outfileName, "w")

                            # overwrite the sequence if parition is found
                            if partition > 0 and partition < 100 and kw['sequence'] is not None:
                                if jobNum <= numHelix:
                                    #argv[1] = '"%s"'%kw['sequence'].upper()
                                    argv[1] = pepSeqStr('upper', '"%s"'%kw['sequence'])
                                else:
                                    #argv[1] = '"%s"'%kw['sequence'].lower()
                                    argv[1] = pepSeqStr('lower', '"%s"'%kw['sequence'])
                            #process = subprocess.Popen(' '.join(argv),
                            #                           stdout=subprocess.PIPE , 
                            #                           stderr=subprocess.PIPE, 
                            #                           bufsize = 1, shell=self.shell, cwd=os.getcwd())
                            process = subprocess.Popen(' '.join(argv),
                                                       stdout=outfile, #subprocess.PIPE , 
                                                       stderr=outfile, #subprocess.PIPE, 
                                                       bufsize = 1, shell=self.shell, cwd=os.getcwd())
                            # print(process)
                            #procToRun[process] = jobNum-1
                            processes[jobNum-1]= process
                            outfiles[jobNum-1] = outfile
                            #print('FUFU2', process, jobNum-1)
                            #outfiles[jobNum] = outfile
                            nbStart += 1

                sleep(1)

            dt = time()-t0
            h,m,s = str(datetime.timedelta(seconds=dt)).split(':')
            self.myprint( 'Docking performed in %.2f seconds, i.e. %s hours %s minutes %s seconds '%(dt, h, m, s))

        sort_index = numpy.argsort(runEnergies)

        # discard crazy low energies that are sometimes generated by
        # CrankiteAD, not sure why
        Elow = -100000
        Ehigh = -Elow
        validLowestE = -Elow
        for nn,i in enumerate(sort_index):
            if runEnergies[i] > Elow:
                validLowestE = runEnergies[i]
                break
    
        if validLowestE==-Elow:
            self.myprint("no run obtained a reasonable low energy solution.\n Stopping calculation without removing calcualtion folder %d"%calcFolder)
            sys.exit(1)

        self.myprint('bestEnergies %s '%(str(runEnergies)))
        self.myprint('bestEnergy in run %d %f (%d)'%(sort_index[nn]+1, validLowestE, nn))

        #try:
        #    from MolKit2.AARotamer import AARotamer, CanonicalAARotamers, AARotamerMutator
        #except ImportError:
        #    #write out energy for top 5 solutions and exit  # MS why do that if import failed ??? makes no sense          
        #    for i in range(min(5,nbRuns)):
        #        self.myprint('No. %d energy found is %3.1f kcal/mol at %s_%d.pdb '%(i+1, runEnergies[sort_index[i]]*0.59219, jobName, sort_index[i]+1))
        #    self.myexit()

        self.myprint('Analyzing results ....')
        kw['print'] = self.myprint
        kw['workingFolder'] = workingFolder
        kw['input'] = jobName
        kw['rec'] = 'rigidReceptor.pdbqt'
        kw['targetFolder'] = os.path.abspath(target_folder)

        ## concatenate all runs into a single pdb file
        ntrajsaved = 0
        self.myprint('concatenating trajectories for run with best energy < %f'%(validLowestE + 20))
        with open(os.path.join(calcFolder, "%s.pdb"%jobName), 'wb') as outFile:
            for i in range(nbRuns):
                if runEnergies[i] >= validLowestE and runEnergies[sort_index[0]] + 20:
                    ntrajsaved += 1
                    with open(os.path.join(calcFolder, "run_%d.pdb"%(i+1)), 'rb') as com:
                        shutil.copyfileobj(com, outFile)
        self.myprint('concatenated %d trajectories'%ntrajsaved)

        from clusterADCP import clusterADCP
        runner = clusterADCP()
        runner(**kw)

        dt = time()-t0
        if kw['reclusterOnly']:
            h,m,s = str(datetime.timedelta(seconds=dt)).split(':')
            self.myprint( 'Reclustering completed %.2f seconds, i.e. %s hours %s minutes %s seconds '%(dt, h, m, s))
            self.keepWokingFolder = True
        else:
            self.myprint( 'Calculations completed %.2f seconds, i.e. %s hours %s minutes %s seconds '%(dt, h, m, s))
            self.myprint('MC search command: %s'%command)
            self.myprint('seed: %s'%str(seeds))
            
        if kw['minimize']:                                                      #OMM new line
            self.myprint('Minimizing docked poses ....')                        #OMM new line
            from openMMmethods import ommTreatment                              #OMM new line
            if self.summaryFile:
                self.summaryFile.close()                                        #OMM new line
            runner_omm = ommTreatment(targetFile, kw['recpath'], workingFolder, jobName)          #OMM new line
            runner_omm(**kw)                                                    #OMM new line

        self.myexit(error = False)




if __name__=='__main__':

    #from ADFR.utils.runADFR import runADFR
    #from ADFR.utils.optParser import ArgParser
    import argparse
    parser = argparse.ArgumentParser(description='AutoDock CrankPep V1.1.%d' % BUILD_VERSION, 
                                  usage="usage: python %(prog)s -s GaRyMiChEL -T rec.trg -w WorkFolder -o output")
                                  # version="%prog 0.1")
    parser.add_argument('--version', action='version', version="1.1.%d" % BUILD_VERSION )
    parser.add_argument("-s", "--sequence",dest="sequence",
                        help="initialize peptide from sequence, lower case for coil and UPPER case for helix")
    parser.add_argument("-p", "--partition",dest="partition",type=int,
                        help="partition for starting from a mixture of helix/coil conformation, percentage(helix)=partition/100 note this option will overwrite the CaSe in sequence")
    parser.add_argument("-i", "--input",dest="input",
                        help="use conformation from pdb file as input")
    parser.add_argument("-T", "--target",dest="target",
                        help="a .trg file created by AGFR or the path to the inflated file")
    parser.add_argument("-L", "--rotlibs",dest="rotlibs",
                        help="a ':' separated list of filenames to load rotamers from ADCP/data/rotamers eg")
    parser.add_argument("-l", "--userRotlibs",dest="userrotlibs",
                        help="a ':' separated list of filenames with path and extension to load user-specified rotamer, eg. -L ./myrots.lib")
    parser.add_argument("-n", "--numSteps", type=int,
                             default=2500000, dest="numSteps", help='max step for one replica')
    parser.add_argument("-N", "--nbRuns", type=int,
                             default=50, dest="nbRuns", help='number of replicas')
    parser.add_argument("-c", "--maxCores", type=int, dest="maxCores")
    parser.add_argument("-k", "--keepWorkingFolder", action="store_true", default=False,
                        help="prevents the deletion of the working Folder")
    parser.add_argument("-o", "--jobName",dest="jobName",
                        help="root name for the _out.pdb docked poses and _summary.dlg log file.\
                        when ommited it is build as [basename(targetFile)]_[peptideSequence]")
    parser.add_argument(
            '-y', "--dryRun", dest="dryRun", action="store_true",
            default=False,
            help="print the first adcp command line and exit")
    parser.add_argument(
            '-cyc', "--cyclic", dest="cyclic", action="store_true",
            default=False,
            help="option for cyclic peptide through backbone")
    parser.add_argument(
            '-cys', "--cystein", dest="cystein", action="store_true",
            default=False,
            help="option for cyclic peptide through CYS-S-S-CYS")
    parser.add_argument(
            '-O', "--overwriteFiles", dest="overwriteFiles",
            action="store_true", default=False,
            help="overwrite existing output files silently")
    parser.add_argument(
            '-S', "--seed", dest="seedValue", type=int, default=-1,
            help="seed for random number generator")
    parser.add_argument("-nc", "--natContacts", type=float,
            dest="nc", help='native contacts cutoff used in the clustering')
    parser.add_argument("-rmsd", "--rmsd", type=float,
            dest="rmsd", help='backbone rmsd cutoff used in the clustering')
    parser.add_argument("-ref", "--ref",dest="ref", help='reference peptide structure for calculating rmsd and fnc')
    parser.add_argument("-m", "--nmodes", type=int, default=100,
                       dest="nmodes", help='maximum number of reported docked poses')
    parser.add_argument("-w", "--workingFolder", default=None,
                       dest="workingFolder", help='folder in which the target file is expanded and the MC runs happen. Will be deleted after the run finished unless -k is specified')
    parser.add_argument("-r", "--reclusterOnly", default=False,  dest="reclusterOnly",
                        action="store_true", help="only perform reclustering, assuming the calcualtion folder still exists")

    # openMM flag support
    parser = add_open_mm_flags(parser)                                          #OMM new line
    if len(sys.argv) == 1:
        parser.print_help()
        
    else:    
        # kw = vars(parser.parse_args())
        # runner = runADCP()        
        # runner(**kw)
        kwn, unknw = parser.parse_known_args()
        if len(unknw)> 0:
            print('Unknown options: ',end="")
            for u in unknw:
                if u.startswith("-"):
                    print(u, end=" ")
            print("")
            if "-t" in unknw:
                print('Option "-t" is deprecated and replaced by "-T".')
        else:
            kw = vars(kwn)
            runner = runADCP()
            runner(**kw)
