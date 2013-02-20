'''
    NRGsuite: PyMOL molecular tools interface
    Copyright (C) 2011 Gaudreault, F., Morin, E. & Najmanovich, R.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''


'''
@title: FlexAID - Simulation.py

@summary: Class that handle the flexAID simulation.

@contain: dictAdjAtom, dictDisAngDih, CreateTempPDB, progressBarHandler
          getVarAtoms, buildcc

@organization: Najmanovich Research Group
@creation date:  Sept. 24, 2010
'''

from collections import defaultdict
from pymol import cmd
from subprocess import Popen, PIPE, STDOUT

import math, os, time, re
import sys
import shutil
import threading
import Color
import Geometry
import UpdateScreen


# Start the simulation with FlexAID
class Start(threading.Thread):

    def __init__(self, top, commandline):

        #print "New instance of Start Class"

        threading.Thread.__init__(self)

        self.commandline = commandline
        print(self.commandline)
        
        self.top = top
        self.FlexAID = self.top.top

        self.FlexAID.ProcessError = False
        self.FlexAID.ProcessRunning = True

        self.start()

    # Start FlexAID on a side thread
    def run(self):        
        
        print("FlexAID starting thread has begun.")
 
        try:
            logfile = open(self.top.Manage.LOGFILE, "w")
            
            if self.FlexAID.OSid == 'WIN':
                self.FlexAID.Run = Popen(self.commandline, shell=False, stdout=logfile, stderr=STDOUT)
            else:
                self.FlexAID.Run = Popen(self.commandline, shell=True, stdout=logfile, stderr=STDOUT)
            
            self.FlexAID.Run.wait()
            
            if self.FlexAID.Run.returncode != 0:
                self.FlexAID.ProcessError = True
                
        except IOError:
            print('  FATAL ERROR: Could not open logfile for FlexAID.')
            self.FlexAID.ProcessError = True
            
        except:
            print('  FATAL ERROR: Could not run the executable FlexAID.')
            print('  Make sure you downloaded NRGsuite for the right platform.')
            self.FlexAID.ProcessError = True
            
        else:
            logfile.close()
        
        self.FlexAID.Run = None
        self.FlexAID.ProcessRunning = False

        print("FlexAID starting thread has ended.")
        

class Parse(threading.Thread):
#class Parse:
    
    # 100 msec
    INTERVAL = 0.10
    
    # 1 minute timeout
    TIMEOUT = INTERVAL * 600

    def __init__(self, top, queue):
        
        threading.Thread.__init__(self)
        
        #print "New instance of Parse Class"

        self.top = top
        self.FlexAID = self.top.top
        self.queue = queue

        self.READ = self.top.Manage.READ
        self.UPDATE = self.top.Manage.UPDATE
        self.LOGFILE = self.top.Manage.LOGFILE
        self.LOGFILETMP = self.top.Manage.LOGFILETMP
        self.ParseFile = self.LOGFILE

        self.ReferencePath = self.FlexAID.IOFile.ReferencePath.get()
        self.listSideChain = self.FlexAID.Config1.Vars.TargetFlex.listSideChain
        self.BindingSiteDisplay = self.FlexAID.Config1.BindingSiteDisplay
                
        self.FlexStatus = self.FlexAID.Config2.FlexStatus.get()
        self.dictFlexBonds = self.FlexAID.Config2.Vars.dictFlexBonds
        self.RngOpt = self.FlexAID.Config1.RngOpt.get()
        self.NbTotalGen = int(self.FlexAID.GAParam.NbGen.get())
        self.DefaultDisplay = self.top.SimDefDisplay.get()

        self.Translation = self.FlexAID.Config2.IntTranslation.get()
        self.Rotation = self.FlexAID.Config2.IntRotation.get()

        self.NbGen = ' / ' + str(self.FlexAID.GAParam.NbGen.get())
        self.FloatNbGen = float(self.FlexAID.GAParam.NbGen.get())        
        self.NbGenFreq = int(self.FlexAID.GAParam.NbGenFreq.get())      # Draw every XX generation
        self.NbTopChrom = int(self.FlexAID.GAParam.NbTopChrom.get())    # Number of Lines READ per Generation

        self.Error = False
        self.ErrorMsg = '*FlexAID ERROR: An unexpected error occured.'
            
        self.Generation = -1
        self.Best = ''
        self.TOP = -1
        
        self.State = 1
        self.CurrentState = 0

        self.Ori =  [ 0.0, 0.0, 0.0 ]    # Origin coordinate
        self.OriX = [ 0.0, 0.0, 0.0 ]   # Origin coordinate with X+1
        self.OriY = [ 0.0, 0.0, 0.0 ]   # Origin coordinate with Y+1

        self.ParseGA = False

        self.nRead = dict()
        self.ListAtom = list()
        self.dictSideChainNRot = dict()
        self.dictSideChainRotamers = dict()
        self.GridVertex = dict()
        self.FixedAngle = dict()
        
        # References
        self.dictSimData = self.top.dictSimData

        self.ReferenceLines = self.top.Manage.ReferenceLines
        self.VarAtoms = self.top.Manage.VarAtoms
        self.RecAtom = self.top.Manage.RecAtom
        self.DisAngDih = self.top.Manage.DisAngDih
        self.dictCoordRef = self.top.Manage.dictCoordRef
        self.listTmpPDB = self.top.Manage.listTmpPDB

        self.nbAtoms = len(self.DisAngDih)
                
        self.top.ProcessParsing = True

        self.start()
        
    '''
    @summary: SUBROUTINE run: Start the simulation
    '''    
    def run(self):
        
        print("FlexAID parsing thread has begun.")

        # Set the auto_zoom to off
        cmd.set("auto_zoom", 0)
        cmd.delete("TOP_*__")
        cmd.delete("RESULT_*")
        cmd.refresh()
        cmd.frame(1)
        
        self.queue.put(lambda: self.top.InitStatus())
        self.queue.put(lambda: self.top.progressBarHandler(0, self.NbTotalGen))
        
        while self.FlexAID.Run is not None and self.FlexAID.Run.poll() is None:
            
            time.sleep(self.INTERVAL)
            
            # Once copied cannot go change file (safe-protection)
            ParseFile = self.ParseFile 
            
            if self.top.SimStatus.get() == 'Paused.':
                continue
                
            elif self.CopyRead(ParseFile):
                self.Error = True
                self.ErrorMsg = '*NRGsuite ERROR: Could not successfully copy/read temporary files.'
                break
                
            if self.nRead.get(ParseFile):
                # Resume a file that was not read completely
                self.Lines = self.Lines[self.nRead[ParseFile]:]
            else:
                self.nRead[ParseFile] = 0
            
            
            for Line in self.Lines:
                
                # check line completion
                m = re.search('\n', Line)
                if m:
                    self.nRead[ParseFile] = self.nRead[ParseFile] + 1
                else:
                    # EOF signal - will resume from here next time
                    break
                
                m = re.match("Grid\[(\d+)\]=", Line)
                if m:
                    index = int(m.group(1))
                
                    strcoor = Line[(Line.find('=')+1):]
                    self.GridVertex[index] = [  float(strcoor[0:8].strip()), 
                                                float(strcoor[8:16].strip()),
                                                float(strcoor[16:24].strip())]

                    continue


                m = re.match("(\s*(\d+) \()", Line)
                if m:
                    #print Line
                    
                    '''
                    Generation:   0
                    best by energy
                        0 (   16.819   128.976    15.591   154.488   171.496  -137.480 )  value= -406.907 fitnes=  100.000
                        1 (   17.286   -58.110   -43.937  -148.819   -26.929   -83.622 )  value= -494.421 fitnes=   99.000
                        2 (   19.386    -1.417    63.780   120.472   -15.591    29.764 )  value= -498.632 fitnes=   98.000
                        3 (   22.186   -38.268    72.283  -106.299   157.323    52.441 )  value= -515.372 fitnes=   97.000
                        4 (   22.653    -1.417    38.268  -109.134    89.291  -126.142 )  value= -533.335 fitnes=   96.000
                    best by fitnes
                        0 (   16.819   128.976    15.591   154.488   171.496  -137.480 )  value= -406.907 fitnes=  100.000
                        1 (   17.286   -58.110   -43.937  -148.819   -26.929   -83.622 )  value= -494.421 fitnes=   99.000
                        2 (   19.386    -1.417    63.780   120.472   -15.591    29.764 )  value= -498.632 fitnes=   98.000
                        3 (   22.186   -38.268    72.283  -106.299   157.323    52.441 )  value= -515.372 fitnes=   97.000
                        4 (   22.653    -1.417    38.268  -109.134    89.291  -126.142 )  value= -533.335 fitnes=   96.000
                    '''

                    self.TOP = int(m.group(2))

                    # Find starting index where to parse column values
                    colNo = len(m.group(1))

                    if self.Best == 'energy' and self.Generation != -1 and self.TOP != -1:

                        # Reading the values calculated for the generation
                        if (self.Generation % self.NbGenFreq) == 0 or self.Generation == self.NbTotalGen:

                            ID = str(self.Generation) + '.' + str(self.TOP)
                            #print("updating " + ID)

                            #print Line
                            Update = UpdateScreen.UpdateScreen( self, ID, colNo, Line, self.CurrentState, self.TOP, 
                                                                self.Translation, self.Rotation )
                            
                            if (self.TOP+1) == self.NbTopChrom:
                                self.State = self.CurrentState

                                # Update energy/fitness table
                                self.queue.put(lambda: self.top.update_DataList())

                                self.top.Refresh_LigDisplay()
                                self.top.Refresh_CartoonDisplay()
                        
                        # Ready to read another file
                        if (self.TOP+1) == self.NbTopChrom:
                        
                            self.nRead[ParseFile] = 0
                                
                            if self.Generation == self.NbTotalGen:
                                self.ParseFile = self.LOGFILE
                                
                            else:
                                if self.Remove_UPDATE():
                                    self.Error = True
                                    self.ErrorMsg = '*NRGsuite ERROR: Could not successfully copy/read temp files.'
                                    break
                                                                    
                    continue

                m = re.match("Generation:\s*(\d+)\s+", Line)
                if m:
                    #print Line
                    
                    self.Generation = int(m.group(1))
                    #print("Generation " + str(self.Generation))
                    self.CurrentState = self.State + 1

                    #print "will update progressbar"
                    self.queue.put(lambda: self.top.progressBarHandler(self.Generation, self.NbTotalGen))

                    continue

                m = re.match("best by (\w+)\s+", Line)
                if m:
                    #print Line

                    self.Best = m.group(1)
                    #print "Best by " + self.Best
                    continue
            
                m = re.match("clustering all individuals", Line)
                if m:
                    #print Line
                    self.queue.put(lambda: self.top.ClusterStatus())

                    continue

                m = re.match("Rotamer for", Line)
                if m:                        
                    self.AddRotamerFromLine(Line)
                
                    continue

                m = re.match("Done.", Line)
                if m:                    
                    continue

                m = re.match("shiftval=", Line)
                if m and self.FlexStatus != '':

                    # Shift values for fixed dihedrals between pair-triplets of atoms
                    fields = Line.split()
                    MergeAtomsAB = fields[1] + fields[2]
                    self.FixedAngle[MergeAtomsAB] = fields[5]                   
                                            
                    continue

                m = re.match("lout\[\d+\]=\s*(\d+)\s+", Line)
                if m:                        
                    self.ListAtom.append(int(m.group(1)))
                    if len(self.ListAtom) == self.nbAtoms:
                        # Order the FLEDIH based on the atoms occurences
                        self.OrderFledih()
                
                    continue

                m = re.match("the protein center of coordinates is:\s+(\S+)\s+(\S+)\s+(\S+)\s+", Line)
                if m:
                    self.Ori[0] = float(m.group(1))
                    self.Ori[1] = float(m.group(2))
                    self.Ori[2] = float(m.group(3))

                    self.OriX[0] = self.Ori[0] + 1.0  # X
                    self.OriX[1] = self.Ori[1]        # Y
                    self.OriX[2] = self.Ori[2]        # Z
                
                    self.OriY[0] = self.Ori[0]        # X
                    self.OriY[1] = self.Ori[1] + 1.0  # Y
                    self.OriY[2] = self.Ori[2]        # Z
                
                    continue            
                
                m = re.match("SIGMA_SHARE", Line)
                if m:                            
                    # send ready to simulate signal  
                    print('  Signal sent to start simulation')
                    self.queue.put(lambda: self.top.RunStatus())

                    #self.ParseGA = True
                    self.ParseFile = self.UPDATE

                    continue
            
                # track errors
                m = re.match("ERROR", Line)
                if m:
                    Line = Line.rstrip('\n')
                    self.ErrorMsg = '*FlexAID ' + Line
                    
                    try:
                        self.FlexAID.Run.terminate()
                    except:
                        pass
                    
                    break

            if self.Error:
                break
        
        if self.FlexAID.ProcessError or self.Error:
            self.queue.put(lambda: self.top.ErrorStatus(self.ErrorMsg))
        else:
            self.queue.put(lambda: self.top.SuccessStatus())
        
        # Put back the auto_zoom to on
        cmd.set("auto_zoom", -1)
        cmd.disable("TOP_*__")
        cmd.refresh()
        cmd.enable("RESULT_*")
        cmd.refresh()
        cmd.frame(1)

        self.top.ProcessParsing = False

        print("FlexAID parsing thread has ended.")

    '''
    @summary: SUBROUTINE OrderFledih: Order the FLEDIH atoms number based on
                                      the ligand construction (lout) if REQUIRED!                
    '''
    def OrderFledih(self):

        # Flexible bond(s) selected?
        # LOOK for ALL the EXTRA column (1 par Flexbond SELECTED)
        if self.FlexStatus != '':
            
            # Find occurrence order
            tot = len(self.ListAtom)
            
            # Get all the keys in the Dictionary (ALL the Atoms)
            order = self.dictFlexBonds.keys()
            order.sort()
        
            for k in order:
        
                if (self.dictFlexBonds[k][2] > 1):                                                        
                                                                
                    # Get the Atoms LIST
                    AtList = list()
                    AtList = self.dictFlexBonds[k][3:]
                    
                    AtListSorted = list()
                    
                    counter = 0
                    totalAtom = self.dictFlexBonds[k][2]                    
                            
                    for an in range(0, tot):
                        NoAtom = str(self.ListAtom[an])
                        
                        for noAt in range(0, totalAtom):                            
                        
                            if NoAtom == AtList[noAt]:
                                AtListSorted.append(NoAtom)
                                AtList[noAt] = 0
                                counter += 1
                                break
                            
                        if counter == totalAtom:
                            break

                    elemTot = len(self.dictFlexBonds[k])

                    for elem in range(3, elemTot):
                        self.dictFlexBonds[k][elem] = AtListSorted[elem-3]            
        

    '''
    @summary: SUBROUTINE AddRotamerFromLine: Adds a rotamer from the output of FlexAID
    '''
    def AddRotamerFromLine(self, Line):
        
        words = Line.split()
        
        #  0      1    2     3      4           5        6        7
        #Rotamer for GLU81- with dihedrals  -87.647  -80.859   27.236            

        if self.dictSideChainNRot.get(words[2],''):
            for i in range(5,len(words)):
                self.dictSideChainRotamers[words[2]].append(float(words[i].strip()))
        else:
            self.dictSideChainNRot[words[2]] = 0
            self.dictSideChainRotamers[words[2]] = [float(words[5].strip())]
            for i in range(6,len(words)):
                self.dictSideChainRotamers[words[2]].append(float(words[i].strip()))
        
        self.dictSideChainNRot[words[2]] += 1
        
        #print "%s now has %d rotamer(s)" % (words[2], self.dictSideChainNRot[words[2]])
    
    '''
    @summary: SUBROUTINE: Remove_UPDATE: Tries to remove the .update to be able to generate another one
    '''   
    def Remove_UPDATE(self):
    
        TIME = 0
        while TIME < self.TIMEOUT:
            try:
                os.remove(self.UPDATE)
                break
                
            except OSError:
                pass
            
            time.sleep(self.INTERVAL)

            TIME = TIME + self.INTERVAL
            
        if TIME >= self.TIMEOUT:
            return 1
            
        return 0
    
    '''
    @summary: SUBROUTINE: CopyRead_UPDATE: Tries to copy then read the .read (log.txt OR .update) file
    '''   
    def CopyRead(self, ParseFile):
    
        TIME = 0
        while TIME < self.TIMEOUT:
            
            try:
                shutil.copy(ParseFile, self.READ)

                readhandle = open(self.READ, 'r')
                self.Lines = readhandle.readlines()
                readhandle.close()                
                break
                
            except OSError:
                pass
            
            except IOError:
                pass
            
            time.sleep(self.INTERVAL)

            TIME = TIME + self.INTERVAL
            
        if TIME >= self.TIMEOUT:
            return 1
            
        return 0

