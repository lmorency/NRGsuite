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
import threading
import Color
import Geometry
import UpdateScreen

# Start the simulation with FlexAID      
class Start(threading.Thread):

    def __init__(self, top,commandline):

        #print "New instance of Start Class"

        threading.Thread.__init__(self)

        self.commandline = commandline
        print self.commandline

        self.top = top
        self.FlexAID = self.top.top

        self.top.ProcessError = False
        self.start()

    # Start FlexAID on a side thread
    def run(self):        
        
        #print "FlexAID starting thread has begun."
        
        self.FlexAID.ProcessRunning = True

        try:
            if self.FlexAID.OSid == 'WIN':
                self.FlexAID.Run = Popen(self.commandline, shell=False, bufsize=1, stdout=PIPE, stderr=STDOUT)
            else:
                self.FlexAID.Run = Popen(self.commandline, shell=True,  bufsize=1, stdout=PIPE, stderr=STDOUT)
                
            self.FlexAID.Run.wait()
        except:

            self.FlexAID.DisplayMessage('   Fatal error: Could not run the executable FlexAID', 1)
            self.FlexAID.DisplayMessage('   Make sure you downloaded NRGsuite for the right platform', 1)
            self.top.ProcessError = True

        #print "FlexAID starting thread has ended."
        
        self.FlexAID.Run = None
        self.FlexAID.ProcessRunning = False
        


class Parse(threading.Thread):
#class Parse():

    def __init__(self, top):
        threading.Thread.__init__(self)
        
        #print "New instance of Parse Class"

        self.top = top
        self.FlexAID = self.top.top

        self.FlexStatus = self.FlexAID.Config2.FlexStatus.get()
        self.OSid = self.FlexAID.OSid
        self.BindingSiteDisplay = self.FlexAID.Config1.BindingSiteDisplay

        self.LOGFILE = self.top.Manage.LOGFILE
                
        self.NbTopChrom = int(self.FlexAID.GAParam.NbTopChrom.get())
        
        self.dictFlexBonds = self.FlexAID.Config2.dictFlexBonds

        self.RngOpt = self.FlexAID.Config1.RngOpt.get()  # Possibility: GLOBAL, LOCCEN, LOCCLF
        
        self.NbTotalGen = int(self.FlexAID.GAParam.NbGen.get())
        self.DefaultDisplay = self.top.SimDefDisplay.get()
        
        self.LOGFILE = self.top.Manage.LOGFILE     
        self.LOGFILETMP = self.top.Manage.LOGFILE + '.tmp'

        self.Generation = -1
        self.Best = ''
        self.TOP = -1
        self.GridVertex = {}

        self.Translation = self.FlexAID.Config2.IntTranslation.get()
        self.Rotation = self.FlexAID.Config2.IntRotation.get()
        
        #self.Updating = 0

        self.State = 0
        self.CurrentState = 0

        self.PrevTextValue = 0
        self.StrCount = '   0'
        self.NbGen = ' / ' + str(self.FlexAID.GAParam.NbGen.get())
        self.FloatNbGen = float(self.FlexAID.GAParam.NbGen.get())
        self.ParseGA = False
        self.FixedAngle = {}

        self.ReferencePath = self.FlexAID.IOFile.ReferencePath

        self.listSideChain = self.FlexAID.Config1.TargetFlex.listSideChain
        self.dictSideChainNRot = {}
        self.dictSideChainRotamers = {}

        self.dictSimData = self.top.dictSimData

        self.ModuloGEN = int(self.FlexAID.GAParam.NbGenFreq.get())     # Draw every XX generation        
        self.NBLineGEN = int(self.FlexAID.GAParam.NbTopChrom.get())    # Number of Lines READ per Generation        

        self.listTmpPDB = list()
        self.dictCoordRef = dict()

        self.CreateTempPDB(self.FlexAID.FlexAIDSimulationProject_Dir)

        # Set the Colors
        self.NBCOLOR = Color.NBCOLOR
        self.PymolColorList = Color.GetHeatColorList(self.NBLineGEN, False)

        #Creation of the dictionaries
        self.RecAtom = self.dictRecAtom(self.top.Manage.INPFlexAIDRunSimulationProject_Dir)        # 3 Neighbors that lead to the middle atoms
        self.DisAngDih = self.dictDisAngDih(self.top.Manage.ICFlexAIDRunSimulationProject_Dir)     # Distance, Angle, Dihedral angle

        self.Ori = [0.0, 0.0, 0.0]    # Origin coordinate
        self.OriX = [0.0, 0.0, 0.0]   # Origin coordinate with X+1
        self.OriY = [0.0, 0.0, 0.0]   # Origin coordinate with Y+1
        
        self.ListAtom = []           # List of the atoms of the ligand

        # In order, atoms that need their values to be modified
        self.VarAtoms = self.getVarAtoms(self.top.Manage.INPFlexAIDRunSimulationProject_Dir)    # 1st = 3 positions, 2nd = 2 positions, 3rd = 1 position
         
        self.start()
        
    '''
    @summary: SUBROUTINE run: Start the simulation
    '''    
    def run(self):
    #def start(self):
        
        #print "FlexAID parsing thread has begun."

        # 20 msec
        INTERVAL = 0.05
        TIMEOUT = INTERVAL * 10

        # Number of ligand atoms
        nbAtoms = len(self.DisAngDih)

        # Always print docking results in frame range (1...Xn)
        self.State = 1

        try:
            file = open(self.ReferencePath,'r')
            self.ReferenceFile = file.readlines()
            file.close()

            if self.FlexAID.Config2.UseReference.get():
                self.dictCoordRef = self.Get_CoordRef()

        except:
            print "Could not read ligand PDB File"
            return

    
        print "Waiting for the simulation to start..."
        while not self.FlexAID.ProcessRunning:
            time.sleep(INTERVAL)
        
        # Wait for the simulation to start...
        while TIMEOUT:
            if not self.FlexAID.Run is None:
                #print "FlexAID is running..."
                break
            elif self.top.ProcessError:
                print "An error occured while trying to run FlexAID"
                return
                        
            time.sleep(INTERVAL)
            TIMEOUT -= INTERVAL

        
        self.top.progressBarHandler(0,self.NbTotalGen)
        self.top.Init_Table()
        
        # Set the auto_zoom to off
        cmd.set("auto_zoom", 0)
        cmd.delete("TOP_*__")
        cmd.frame(1)

        self.top.InitStatus()
        
        while self.FlexAID.Run is not None: # and self.FlexAID.Run.poll() is None:

            while (1):
                # Parsing output
                try:
                    Line = self.FlexAID.Run.stdout.readline()
                    #print Line
                except:
                    break

                # stop reading from stdout buffer
                if Line == '':
                    break
                
                if self.ParseGA:
                             
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

                    m = re.match("(\s*(\d+) \()", Line)
                    if m:
                        self.TOP = int(m.group(2))

                        # Find starting index where to parse column values
                        colNo = len(m.group(1))

                        if self.Best == 'energy' and self.Generation != -1 and self.TOP != -1:

                            # Reading the values calculated for the generation
                            if (self.Generation % self.ModuloGEN) == 0 or self.Generation == self.NbTotalGen:

                                #self.Updating += 1

                                ID = str(self.Generation) + '.' + str(self.TOP)
                                #print("updating " + ID)

                                #print Line
                                Update = UpdateScreen.UpdateScreen( self, ID, colNo, Line, self.CurrentState, self.TOP, 
                                                                    self.Translation, self.Rotation )

                                    
                                #Update.join()

                                #while(self.Updating):
                                    #print "is already updating..."
                                #    time.sleep(INTERVAL)
                                
                                    
                                #    while(self.Updating):
                                #        time.sleep(INTERVAL)

                                if (self.TOP+1) == self.NBLineGEN:
                                    self.State = self.CurrentState

                                    try:
                                        # append range object to next state
                                        cmd.create(self.BindingSiteDisplay, self.BindingSiteDisplay, 1, self.CurrentState)
                                    except:
                                        self.FlexAID.DisplayMessage("ERROR: Could not display binding-site. Object no longer exists", 1)

                                    # Update energy/fitness table
                                    self.top.update_DataList()

                        continue


                    m = re.match("best by (\w+)\s+", Line)
                    if m:
                        self.Best = m.group(1)
                        #print "Best by " + self.Best
                        continue

                    m = re.match("Generation:\s*(\d+)\s+", Line)
                    if m:

                        self.Generation = int(m.group(1))
                        #print("Generation " + str(self.Generation))
                        self.CurrentState = self.State + 1

                        #print "will update progressbar"
                        # ProgressionBar Handler
                        try:
                            self.top.progressBarHandler(self.Generation, self.NbTotalGen)
                        except:
                            pass

                        continue
                    
                    m = re.match("clustering all individuals", Line)
                    if m:
                        self.top.ClusterStatus()

                    m = re.match("Done", Line)
                    if m:
                        self.top.SuccessStatus()
                        
                else:
    
                    # track errors
                    if Line.startswith('ERROR'):
                        self.FlexAID.DisplayMessage(str("A critical error occured\n" + Line), 1)
                        self.top.ErrorStatus()
                        try:
                            self.FlexAID.Run.terminate()
                        except:
                            pass                            
                        break
                    

                    m = re.match("Rotamer for", Line)
                    if m:                        
                        self.AddRotamerFromLine(Line)
                        

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

                        if len(self.ListAtom) == nbAtoms:
                        
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

                    m = re.match("Grid\[(\d+)\]=", Line)
                    if m:
                        index = int(m.group(1))
                        
                        strcoor = Line[(Line.find('=')+1):]
                        self.GridVertex[index] = [float(strcoor[0:8]), float(strcoor[8:16]), float(strcoor[16:24])]
                        continue
                        
                    m = re.match("SIGMA_SHARE", Line)
                    if m:
                        self.ParseGA = True
                        self.top.RunStatus()
                        continue
                                
            time.sleep(INTERVAL)


        #print "FlexAID parsing thread has ended."
        self.FlexAID.ProcessRunning = False

        # Empty rotamers data
        self.dictSideChainNRot.clear()
        self.dictSideChainRotamers.clear()

        # Re-enable buttons
        self.top.Btn_Start.config(state='normal') 
        self.top.Btn_PauseResume.config(state='disabled')
        self.top.Btn_Stop.config(state='disabled')
        self.top.Btn_Abort.config(state='disabled')

        # Put back the auto_zoom to on
        cmd.set("auto_zoom", -1)


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
    @summary: SUBROUTINE getVarAtoms: Get the atoms that need their values to be modified                  
    @return: listAtVar - list
    '''
    def getVarAtoms(self, inpPath):
        
        NoAtom3 = 0
        NoAtom2 = 0
        NoAtom1 = 0

        file = open(inpPath)
        inpFile = file.readlines()
        file.close()
        
        # Creation of a dictionary containing the 3 neighbors of each atoms
        for line in inpFile:
            if line.startswith('HETTYP'):
                noAtom = int(line[7:11])
                if int(line[32:36]) == 0:
                    if int(line[27:31]) == 0:
                        if int(line[22:26]) == 0:
                            NoAtom3 = noAtom
                        else:
                            NoAtom2 = noAtom
                    else:
                        NoAtom1 = noAtom
                        
                
        return [NoAtom3, NoAtom2, NoAtom1]
    

    '''
    @summary: SUBROUTINE Get_CoordRef: Get the list of initial coordinates of the reference
    '''
    def Get_CoordRef(self):

        dictCoordRef = {}

        for Line in self.ReferenceFile:
            
            if Line.startswith('HETATM'):
                index = int(Line[6:11].strip())
                CoordX = float(Line[30:38].strip())
                CoordY = float(Line[38:46].strip())
                CoordZ = float(Line[46:54].strip())

                dictCoordRef[index] = [ CoordX, CoordY, CoordZ ]

        return dictCoordRef
    
    '''
    @summary: SUBROUTINE dictRecAtom: Create a dictionary containing the neighbors of
              each atoms for the reconstruction.                  
    @param inpPath: Path to the inp file
    @return: RecAtom - dictionary
    '''
    def dictRecAtom(self, inpPath):
        RecAtom = {}

        file = open(inpPath)
        inpFile = file.readlines()
        file.close()
        
        # Creation of a dictionary containing the 3 neighbors of each atoms
        for line in inpFile:
            if line.startswith('HETTYP'):
                noLine = int(line[7:11])
                RecAtom[noLine] = [int(line[22:26]), int(line[27:31]), int(line[32:36])]
                
        return RecAtom


    def printDictDisAngDih(self):
        
        print("******************************************")
        for k, v in self.DisAngDih.items():
            print("ATOM: " + str(k) + "  " + str(self.DisAngDih[k][0]) + " " + str(self.DisAngDih[k][1]) + " " + str(self.DisAngDih[k][2]))

    '''
    @summary: SUBROUTINE dictRecAtom: Create a dictionary containing the neighbors of
              each atoms.                  
    @param icPath: Path to the ic file
    @return: DisAngDih - dictionary
    '''
    def dictDisAngDih(self, icPath):

        DisAngDih = {}        
                
        file = open(icPath)
        icFile = file.readlines()
        file.close()
        
        # Creation of a dictionary containing the 3 neighbors of each atoms
        for line in icFile:
            if line[0:6] != 'REFPCG':
                noLine = int(line[1:5])
                DisAngDih[noLine] = [float(line[7:15]), float(line[16:24]), float(line[25:33])]
                    
        return DisAngDih
    
    '''
    @summary: SUBROUTINE CreateTempPDB: Create X copy of the ligand PDB file                  
    '''   
    def CreateTempPDB(self, Path):
        
        # Create the custom path for the temporary pdb files
        for i in range (self.NBLineGEN + 1):
            self.listTmpPDB.append(os.path.join(Path,'LIG' + str(i) + '.pdb'))
            #print Path + "/LIG" + str(i) + '.pdb'
