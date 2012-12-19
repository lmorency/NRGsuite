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
from pymol import cmd

#import threading
import Geometry
import Constants


#class UpdateScreen(threading.Thread):
class UpdateScreen():

    def __init__(self, top, ID, colNo, Line, State, TOP):
        #threading.Thread.__init__(self)

        self.top = top

        # unique ID of complex
        self.ID = ID
        #print self.ID + " initialized!"

        # TOPn to edit
        self.TOP = TOP

        # Starting index of line
    	self.colNo = colNo

        # input line to parse
        self.Line = Line

        # State on which updating is happening
        self.State = State

        self.dictSimData = self.top.dictSimData
        self.dictFlexBonds = self.top.dictFlexBonds
        self.dictCoordRef = self.top.dictCoordRef

        self.LigandName = self.top.FlexAID.IOFile.LigandName.get()
        self.ProtName = self.top.FlexAID.IOFile.ProtName.get()

        self.LigandObj = self.LigandName + '_' + str(self.TOP+1)
        self.ProteinObj = self.ProtName + '_' + str(self.TOP+1)
        
        # Selections of molecules (ligand/side-chain/protein)
        self.selProtein = '(' + self.ProteinObj + ' & present)'
        self.selLigand = '(' + self.LigandObj + ' & present)'
        self.selSideChains = ''

        # start thread
        self.start()


    #def run(self):
    def start(self):
    
        self.Update()
        #self.top.Updating -= 1


    # Updates the PyMOL interface
    def Update(self):
        
        try:
            # Copy the initial protein (Frame 1) into the working state
            cmd.create(self.ProteinObj, self.ProtName, 1, self.State)

            # Display the last frame
            cmd.frame(self.State)
            #print "Switched to frame " + str(self.State)

        except:
            self.CriticalError("Object " + str(self.ProtName) + " no longer exists")
        

        if self.UpdateLigandAnchorPoint():
            self.selSideChains = self.UpdateSideChainConformations()

        if self.UpdateLigandFlexibility() and self.WriteOutLigand() and self.EditView() and self.UpdateDataList():
            try:
                # delete temporary protein PDB file
                cmd.delete(self.ProteinObj)
            except:
                self.CriticalError("Object " + str(self.ProteinObj) + " no longer exists")


    '''=========================================================================
       UpdateDataList: Updates the table containing energy/fitness values
    ========================================================================='''
    def UpdateDataList(self):

        try: 

            self.dictSimData[self.TOP+1] = []

            #Get index position of energy column
            self.colNo = self.Line.rfind("value=")
            
            if self.colNo != -1:
                #Update the energy of the dictionary
                self.dictSimData[self.TOP+1].append(self.Line[self.colNo+6:self.colNo+15].strip())     # Energy
            else:
                self.dictSimData[self.TOP+1].append('N/A')


            #Get index position of fitnes column
            self.colNo = self.Line.rfind("fitnes=")
            if self.colNo != -1:
                #Update the fitnes of the dictionary
                self.dictSimData[self.TOP+1].append(self.Line[self.colNo+7:self.colNo+16].strip())    # Fitness
            else:
                self.dictSimData[self.TOP+1].append('N/A')

            # RMSD of ligand
            RMSD = Geometry.rmsd(self.dictCoord, self.dictCoordRef)

            if RMSD != 'N/A':
                RMSD = '%.3f' % RMSD

            self.dictSimData[self.TOP+1].append(RMSD)


        except:
            self.CriticalError("Error while updating data list")
            return 0

        return 1

    '''=========================================================================
       EditView: Edit the visual aspects of the PyMOL interface
    ========================================================================='''
    def EditView(self):
        
        #Display the ligand with the new coordinates
        #print "will load " + self.LigandObj + " in state " + str(self.State) 
        try:
            cmd.load(self.top.listTmpPDB[self.TOP+1], self.LigandObj, self.State)

            # No flexible side-chain(s)
            if self.selSideChains == '':
                selString = self.selLigand + " or " + self.selProtein
            else:
                selString = self.selSideChains + " or " + self.selLigand + " or " + self.selProtein                         
            #print selString


            # Create solution object
            # Object contains the whole protein-ligand complex
            SolutionObj = "TOP_" + str(self.TOP+1) + "__"
            cmd.create(SolutionObj, selString, self.State, self.State)

            # Color ligand of solution TOP
            cmd.color(self.top.PymolColorList[self.TOP], "(resn LIG & " + SolutionObj + " & present)")
           
            # Color side-chains of solution TOP
            if self.selSideChains != '':
                cmd.color(self.top.PymolColorList[self.TOP], self.selSideChains.replace(self.ProteinObj,SolutionObj))

            #cmd.show(self.DefaultDisplay, "(resn LIG & sol_*__ & present)")
            cmd.show(self.top.DefaultDisplay, "(resn LIG & " + SolutionObj + " & present)")

            if self.selSideChains != '':
                #cmd.show("sticks", self.selSideChains.replace(self.protName__,"sol_*__"))
                cmd.show("sticks", self.selSideChains.replace(self.ProteinObj,SolutionObj))

            cmd.delete(self.LigandObj)                        

        except:
            self.CriticalError("Error while editing view")
            return 0

        return 1

    '''=========================================================================
       WriteOutLigand: Ouputs PDB ligand file with new coordinates
    ========================================================================='''
    def WriteOutLigand(self):

        try:
        
            #Get the new coordinates of the ligand
            self.dictCoord = {}    
            self.dictCoord = Geometry.buildcc(self.top.ListAtom,self.top.RecAtom,self.top.DisAngDih,self.top.Ori)
        
            #Replace the coordinate in pdb file with the new one
            #print "writing to " + self.top.listTmpPDB[self.TOP+1]
            text_file = open(self.top.listTmpPDB[self.TOP+1], 'w')
            
            for pdbLine in self.top.ReferenceFile:
                
                type = pdbLine[0:6].strip()
                if (type == 'HETATM') or (type == 'ATOM'):

                    NoAtom = int(pdbLine[7:11])
                    #print NoAtom
                    
                    atomX = float(self.dictCoord[NoAtom][0])
                    atomY = float(self.dictCoord[NoAtom][1])
                    atomZ = float(self.dictCoord[NoAtom][2])
                    
                    tmpLine = pdbLine[0:30]
                    tmpLine += '%8.3f' % atomX     # The atom X coordinate
                    tmpLine += '%8.3f' % atomY     # The atom Y coordinate
                    tmpLine += '%8.3f' % atomZ     # The atom Z coordinate
                    tmpLine += pdbLine[54:]
                    
                    text_file.write(tmpLine)

            else:                                        
                text_file.write(pdbLine)                                    

            text_file.close()                               


        except IOError:
            self.CriticalError("Error while writing PDB ligand file.")
            return 0

        return 1

    '''=========================================================================
       UpdateLigandFlexibility: Updates the dihedral angles of the ligand
    ========================================================================='''
    def UpdateLigandFlexibility(self):

        # Flexible bond(s) selected?
        # LOOK for ALL the EXTRA column (1 par Flexbond SELECTED)
        try:
            if self.top.FlexStatus != '':

                for k in sorted(iter(self.dictFlexBonds)):

                    # k: Represent the FLEDIH number 
                    # Which is the key of the dictionary

                    # ==================================================================

                    # If the dihangle between atoms NEED to be modified
                    if self.dictFlexBonds[k][0] == 1: 

                        # Deplacement in the LogFile Column
                        ColValue = float(self.Line[self.colNo:self.colNo+9])
                        self.colNo = self.colNo + 10

                        # Is there ONLY 1 atom that define the flexible bond
                        if self.dictFlexBonds[k][2] == 1:

                            self.top.DisAngDih[int(self.dictFlexBonds[k][3])][2] = ColValue

                        # Is there MULTIPLE atoms that define the flexible bond
                        elif self.dictFlexBonds[k][2] > 1:                                                        

                            # Multiple possibilities, need to mix up the atoms number
                            # Example: [1 ,2 ,3] will give [1], [1, 2], [2, 3]

                            # SET the 1st ATOM Dihedral Angle...
                            self.top.DisAngDih[int(self.dictFlexBonds[k][3])][2] = ColValue                                               

                            for flexA in range(1, self.dictFlexBonds[k][2]):                 

                                ATflex_A = self.dictFlexBonds[k][flexA + 2]
                                ATflex_B = self.dictFlexBonds[k][flexA + 3]                                                                 

                                #print("*****************************************")
                                ATmerge = ATflex_A + ATflex_B

                                # SEARCH for the 2 atoms in FLEDIH that will be present
                                # in SHIFTVAL

                                # Be sure the key exist before calling the value
                                if self.top.FixedAngle.has_key(ATmerge):

                                    #print("1st Atom: " + str(ATflex_A))  
                                    #print("2nd Atom: " + str(ATflex_B))
                                    #print("ShiftVal: " + str(ATmerge))                                                        

                                    # Get the constant angle IN SHIFTVAL
                                    ConstAngle = float(self.top.FixedAngle[ATmerge])                                                                     

                                    # ADD the constant Angle WITH the Column value in the LOGFILE
                                    ColValue = ColValue + ConstAngle

                                    # SET the 2nd ATOM Dihedral Angle...
                                    self.top.DisAngDih[int(ATflex_B)][2] = ColValue                                                                   
        except:
            self.CriticalError("Error while updating ligand flexibility")
            return 0

        return 1

    '''=========================================================================
       UpdateLigandAnchorPoint: Updates the position of the ligand on the screen
    ========================================================================='''
    def UpdateLigandAnchorPoint(self):

        try: 
            # Alter atom X,Y,Z of ligand according to GPA positionning columns
            if self.top.RngTYPE == 'LOCCLF':

                index = int(float(self.Line[self.colNo:self.colNo+9]))
                self.colNo += 10

                coordX = float(self.top.GridLines[index][30:38].strip())     # The atom X coordinate
                coordY = float(self.top.GridLines[index][39:46].strip())     # The atom Y coordinate
                coordZ = float(self.top.GridLines[index][47:54].strip())     # The atom Z coordinate

                pointA = [coordX, coordY, coordZ]
                pointB = [self.top.OriX[0], self.top.OriX[1], self.top.OriX[2]]
                pointC = [self.top.Ori[0], self.top.Ori[1], self.top.Ori[2]]
                pointD = [self.top.OriY[0], self.top.OriY[1], self.top.OriY[2]]

                self.top.DisAngDih[self.top.VarAtoms[0]][0] = Geometry.distance(pointA, pointB)                             
                self.top.DisAngDih[self.top.VarAtoms[0]][1] = Geometry.angle(pointA, pointB, pointC)
                self.top.DisAngDih[self.top.VarAtoms[0]][2] = Geometry.dihedralAngle(pointA, pointB, pointC, pointD)

                self.top.DisAngDih[self.top.VarAtoms[1]][1] = float(self.Line[self.colNo:self.colNo+9])
                self.top.DisAngDih[self.top.VarAtoms[1]][2] = float(self.Line[self.colNo+10:self.colNo+19])
                self.top.DisAngDih[self.top.VarAtoms[2]][2] = float(self.Line[self.colNo+20:self.colNo+29])

                self.colNo += 30

            else:

                self.top.DisAngDih[self.top.VarAtoms[0]][0] = float(self.Line[self.colNo:self.colNo+9])
                self.top.DisAngDih[self.top.VarAtoms[0]][1] = float(self.Line[self.colNo+10:self.colNo+19])
                self.top.DisAngDih[self.top.VarAtoms[0]][2] = float(self.Line[self.colNo+20:self.colNo+29])
                self.top.DisAngDih[self.top.VarAtoms[1]][1] = float(self.Line[self.colNo+30:self.colNo+39])
                self.top.DisAngDih[self.top.VarAtoms[1]][2] = float(self.Line[self.colNo+40:self.colNo+49])
                self.top.DisAngDih[self.top.VarAtoms[2]][2] = float(self.Line[self.colNo+50:self.colNo+59])

                self.colNo += 60

                #print("*******************************************")
                #print("Line Log: " + str(self.Line))
                #print("ATOM: " + str(self.top.VarAtoms[0]) + "  " + str(self.top.DisAngDih[self.top.VarAtoms[0]][0]) + " " +
                #      str(self.top.DisAngDih[self.top.VarAtoms[0]][1]) + " " + str(self.top.DisAngDih[self.top.VarAtoms[0]][2]))
                #print("ATOM: " + str(self.top.VarAtoms[1]) + "           " +
                #      str(self.top.DisAngDih[self.top.VarAtoms[1]][1]) + " " + str(self.top.DisAngDih[self.top.VarAtoms[1]][2]))
                #print("ATOM: " + str(self.top.VarAtoms[2]) + "                       " + str(self.top.DisAngDih[self.top.VarAtoms[2]][2]))


                # In the dictionary: key = ConnNo, value = 0/1 NOT Selected / Selected (to be flexible),
                #                                  value = 0/1 NOT Forced / Forced, 
                #                                  value = Number of HETTYP (atom no)
                #                                  value = List of ATOMS no (from FLEDIH)

        except:
            self.CriticalError("Error while updating ligand anchor point")
            return 0

        return 1

    '''=========================================================================
      .UpdateSideChainConformations: Update side-chain dihedral angles using rotamer library
    ========================================================================='''
    def UpdateSideChainConformations(self):
        
        try:
            # temporary sel. var
            strSelectSC  = ''

            # Loop through Flexible side-chains
            for j in range(0,len(self.top.listFlexSideChain)):

                #print "Setting dihedrals for " + self.top.listFlexSideChain[j]

                # Were any rotamers accepted for this side-chain
                if self.top.listFlexSideChainNRot[j] > 0:

                    #print "List of possible rotamers:" + str(self.top.listFlexSideChainNRot[j])

                    # Residu Name
                    Res = self.top.listFlexSideChain[j][0:3]
                    Num = self.top.listFlexSideChain[j][3:len(self.top.listFlexSideChain[j])-1]
                    Chn = self.top.listFlexSideChain[j][len(self.top.listFlexSideChain[j])-1:len(self.top.listFlexSideChain[j])]

                    strSelectSC += "(resn " + Res + " & resi " + Num
                    if Chn != '-':
                        strSelectSC += " & chain " + Chn

                    strSelectSC += " & ! name C+O+N " + " & " + self.ProteinObj + " & present) or "

                    # Get Integer value from GA.
                    IntVal = int(float(str(self.Line[self.colNo:(self.colNo+9)]).strip()) + 0.5)

                    if IntVal >= 1: # 0 is the default PDB side-chain conf.

                        # Starting index in list of rotamers for Residu Name
                        StartIndex = int(self.top.dictRotamersIndex[self.top.listFlexSideChain[j]][IntVal-1]*Constants.nFlexBonds[Res])

                        # Get List of Dihedrals to rebuild
                        for k in range(0,Constants.nFlexBonds[Res]):

                            # Set dihedrals for side-chain
                            cmd.set_dihedral(self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+0]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+1]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+2]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+3]),
                                             self.top.dictRotamers[Res][StartIndex+k], self.State)

                            
                    # Get next column starting index
                    self.colNo = self.colNo + 10


            # Side-chain selection string - remove last 4 chars
            if strSelectSC != '':
                strSelectSC = strSelectSC[:len(strSelectSC)-4]

            return strSelectSC

        except:

            self.CriticalError("Error while updating side-chain conformations")
    
    '''=========================================================================
      Get_AtomString: Retrives the PyMOL atom selection string
    ========================================================================='''
    def Get_AtomString(self, R, N, C, Atom):

        AtomString  = "resn " + R + " & resi " + N
        if C != '-':
            AtomString += " & chain " + C
            
        AtomString += " & name " + Atom
        AtomString += " & " + self.ProteinObj

        return AtomString

    '''=========================================================================
      CriticalError: Abort simulation if there was a fatal error
    ========================================================================='''
    def CriticalError(self, text):

        self.top.FlexAID.DisplayMessage("CRITICAL ERROR: " + text, 1)
    
        #Create the .abort file
        abort_file = open(self.top.top.Manage.ABORT, 'w')
        abort_file.close()

        