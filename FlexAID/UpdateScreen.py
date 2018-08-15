'''
    NRGsuite: PyMOL molecular tools interface
    Copyright (C) 2011 Gaudreault, F., Morency, LP. & Najmanovich, R.

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

from __future__ import print_function

from pymol import cmd
from pymol import util

#import threading
import Geometry
import Constants


class UpdateScreen(object):

    def __init__(self, top, ID, colNo, Line, State, TOP, Translation, Rotation):
        #threading.Thread.__init__(self)

        self.top = top
        
        self.Translation = Translation
        self.Rotation = Rotation
        
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
        # self.State = State
        self.State = 1
        
        self.dictFlexBonds = self.top.dictFlexBonds
        
        self.LigandName = self.top.LigandName
        self.TargetName = self.top.TargetName
        
        self.LigandObj = self.LigandName + '_' + str(self.TOP+1)
        self.TargetObj = self.TargetName + '_' + str(self.TOP+1)
        
        # Selections of molecules (ligand/side-chain/protein)
        self.selTarget = '(' + self.TargetObj + ' & present)'
        self.selLigand = '(' + self.LigandObj + ' & present)'
        self.selSideChains = ''

        # start thread
        self.start()
    
    def start(self):
    
        self.Update()
    
    # Updates the PyMOL interface
    def Update(self):
        
        try:
            # Copy the initial protein (Frame 1) into the working state
            cmd.create(self.TargetObj, self.TargetName, 1, self.State)
            cmd.refresh()

            # Display the last frame
            cmd.frame(self.State)
            #print "Switched to frame " + str(self.State)

        except:
            self.CriticalError("Object " + str(self.TargetName) + " no longer exists")
        
        
        if not self.UpdateLigandAnchorPoint() and not self.UpdateLigandFlexibility():
        
            self.selSideChains = self.UpdateSideChainConformations()
            
            if self.WriteOutLigand() or self.EditView() or \
               self.top.UpdateDataList(self.Line, self.TOP, self.top.Reference, self.dictCoord):
                self.Delete_Object()

        self.Delete_Object()
        
        return
        
    '''=========================================================================
       Delete_Object: deleted the updating object if an error occurs
    ========================================================================='''
    def Delete_Object(self):

        try:
            # delete temporary protein PDB file
            cmd.delete(self.TargetObj)
            cmd.refresh()

            cmd.delete(self.LigandObj)
            cmd.refresh()

        except:
            self.CriticalError("Object " + str(self.TargetObj) + " no longer exists")

    '''=========================================================================
       EditView: Edit the visual aspects of the PyMOL interface
    ========================================================================='''
    def EditView(self):
        
        #Display the ligand with the new coordinates
        #print "will load " + self.LigandObj + " in state " + str(self.State) 
        try:
            cmd.load(self.top.listTmpPDB[self.TOP+1], self.LigandObj, self.State)
            cmd.refresh()

            # No flexible side-chain(s)
            if self.selSideChains == '':
                selString = self.selLigand + " or " + self.selTarget
            else:
                selString = self.selSideChains + " or " + self.selLigand + " or " + self.selTarget                         
            #print selString

            # Create solution object
            # Object contains the whole protein-ligand complex
            SolutionObj = "TOP_" + str(self.TOP+1) + "__"
            cmd.create(SolutionObj, selString, self.State, self.State)
            cmd.refresh()

            # Color ligand of solution TOP
            Selection = "(resn LIG & " + SolutionObj + " & present)"
            cmd.color(self.top.top.PymolColorList[self.TOP], Selection)
            util.cnc(Selection)
            cmd.refresh()
           
            # Color side-chains of solution TOP
            if self.selSideChains != '':
                selSC = self.selSideChains.replace(self.TargetObj,SolutionObj)
                cmd.color(self.top.top.PymolColorList[self.TOP], selSC)
                util.cnc(selSC)
                cmd.refresh()

            #cmd.show(self.DefaultDisplay, "(resn LIG & sol_*__ & present)")
            cmd.show(self.top.DefaultDisplay, "(resn LIG & " + SolutionObj + " & present)")
            cmd.refresh()

            if self.selSideChains != '':
                #cmd.show("sticks", self.selSideChains.replace(self.TargetName__,"sol_*__"))
                cmd.show("sticks", self.selSideChains.replace(self.TargetObj,SolutionObj))
                cmd.refresh()

        except:
            self.CriticalError("Could not refresh the visual display")
            return 1

        return 0

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
            
            for pdbLine in self.top.ReferenceLines:
                
                type = pdbLine[0:6].strip()
                if type == 'HETATM' or type == 'ATOM':

                    NoAtom = int(pdbLine[6:11])
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
            self.CriticalError("Could not write PDB ligand file.")
            return 1

        return 0

    '''=========================================================================
       UpdateLigandFlexibility: Updates the dihedral angles of the ligand
    ========================================================================='''
    def UpdateLigandFlexibility(self):
        
        # Flexible bond(s) selected?
        try:
            if self.top.FlexStatus != '':

                for k in sorted(self.dictFlexBonds.keys()):

                    # k: Represent the FLEDIH number 
                    # Which is the key of the dictionary

                    # ==================================================================

                    # If the dihangle between atoms NEED to be modified
                    if self.dictFlexBonds[k][0] == 1:

                        # Read and move to next column
                        ColValue = float(self.Line[self.colNo:self.colNo+10])
                        self.colNo = self.colNo + 11

                        # Is there ONLY 1 atom that define the flexible bond
                        if self.dictFlexBonds[k][2] == 1:

                            self.top.DisAngDih[int(self.dictFlexBonds[k][3])][2] = ColValue

                        # Is there MULTIPLE atoms that define the flexible bond
                        elif self.dictFlexBonds[k][2] > 1:

                            #print "How many defines this flexible bonds", self.dictFlexBonds[k][2]
                            # Multiple possibilities, need to mix up the atoms number
                            # Example: [1 ,2 ,3] will give [1], [1, 2], [2, 3]

                            # SET the 1st ATOM Dihedral Angle...
                            self.top.DisAngDih[int(self.dictFlexBonds[k][3])][2] = ColValue

                            for flexA in range(1, self.dictFlexBonds[k][2]):
                                
                                ATflex_A = self.dictFlexBonds[k][flexA + 2]
                                ATflex_B = self.dictFlexBonds[k][flexA + 3]

                                ATmerge = ''
                                ATmergeAB = ATflex_A + ATflex_B
                                ATmergeBA = ATflex_B + ATflex_A
                                
                                # Be sure the key exist before calling the value
                                if ATmergeAB in self.top.FixedAngle:
                                    ATmerge = ATmergeAB
                                    Factor = 1
                                elif ATmergeBA in self.top.FixedAngle:
                                    ATmerge = ATmergeBA
                                    Factor = -1
                                
                                if ATmerge:
                                    # Get the constant angle IN SHIFTVAL
                                    ConstAngle = float(self.top.FixedAngle[ATmerge])
                                    #print "ConstAngle", ConstAngle
                                    
                                    # ADD the constant Angle WITH the Column value in the LOGFILE
                                    ColValue = ColValue + float(Factor)*ConstAngle
                                    #print "ColValue", ColValue

                                    # SET the 2nd ATOM Dihedral Angle...
                                    self.top.DisAngDih[int(ATflex_B)][2] = ColValue
        except:
            self.CriticalError("Could not update ligand flexibility")
            return 1

        return 0

    '''=========================================================================
       UpdateLigandAnchorPoint: Updates the position of the ligand on the screen
    ========================================================================='''
    def UpdateLigandAnchorPoint(self):

        try: 
                                    
            if self.Translation:
                index = int(float(self.Line[self.colNo:self.colNo+10]))
                
                coordX = self.top.GridVertex[index][0]     # The atom X coordinate
                coordY = self.top.GridVertex[index][1]     # The atom Y coordinate
                coordZ = self.top.GridVertex[index][2]     # The atom Z coordinate
                
                pointA = [coordX, coordY, coordZ]
                pointB = [self.top.OriX[0], self.top.OriX[1], self.top.OriX[2]]
                pointC = [self.top.Ori[0], self.top.Ori[1], self.top.Ori[2]]
                pointD = [self.top.OriY[0], self.top.OriY[1], self.top.OriY[2]]

                self.top.DisAngDih[self.top.VarAtoms[0]][0] = Geometry.distance(pointA, pointB)                             
                self.top.DisAngDih[self.top.VarAtoms[0]][1] = Geometry.angle(pointA, pointB, pointC)
                self.top.DisAngDih[self.top.VarAtoms[0]][2] = Geometry.dihedralAngle(pointA, pointB, pointC, pointD)

                self.colNo += 11
            
            if self.Rotation:
                self.top.DisAngDih[self.top.VarAtoms[1]][1] = float(self.Line[self.colNo:self.colNo+10])
                self.top.DisAngDih[self.top.VarAtoms[1]][2] = float(self.Line[self.colNo+11:self.colNo+21])
                self.top.DisAngDih[self.top.VarAtoms[2]][2] = float(self.Line[self.colNo+22:self.colNo+32])
                self.colNo += 33

        except:
            self.CriticalError(" Could not update ligand anchor point")
            return 1

        return 0

    '''=========================================================================
      .UpdateSideChainConformations: Update side-chain dihedral angles using rotamer library
    ========================================================================='''
    def UpdateSideChainConformations(self):
        
        try:
            # temporary sel. var
            strSelectSC  = ''

            # Loop through Flexible side-chains
            for residue in self.top.listSideChain:

                #print "Setting dihedrals for " + residue
                
                # Were any rotamers accepted for this side-chain
                if self.top.dictSideChainNRot.get(residue,''):

                    #print "List of possible rotamers:" + str(self.top.dictSideChainNRot[residue])

                    # Residu Name
                    Res = residue[0:3]
                    Num = residue[3:len(residue)-1]
                    Chn = residue[len(residue)-1:len(residue)]

                    strSelectSC += "(resn " + Res + " & resi " + Num
                    if Chn != '-':
                        strSelectSC += " & chain " + Chn
                    else:
                        strSelectSC += " & chain ''"                    
                    
                    strSelectSC += " & ! name C+O+N " + " & " + self.TargetObj + " & present) or "

                    # Get Integer value from GA.
                    IntVal = int(float(self.Line[self.colNo:(self.colNo+10)].strip()) + 0.5)
                    nFlex = Constants.nFlexBonds[Res]
                    
                    #print("IntVal", str(IntVal))
                    #print("nFlex", str(nFlex))
                    #print strSelectSC
                    
                    if IntVal > 0: # 0 is the default PDB side-chain conf.
        
                        # Get List of Dihedrals to rebuild
                        for k in range(0,nFlex):
                        
                            '''
                            print "for k=" + str(k) + " for residue=" + str(residue)
                            print self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+0])
                            print self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+1])
                            print self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+2])
                            print self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+3])
                            print self.top.dictSideChainRotamers[residue]
                            print self.top.dictSideChainRotamers[residue][(IntVal-1)*nFlex+k]
                            '''
                            
                            # Set dihedrals for side-chain
                            cmd.set_dihedral(self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+0]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+1]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+2]),
                                             self.Get_AtomString(Res,Num,Chn,Constants.setDihedrals[Res][4*k+3]),
                                             self.top.dictSideChainRotamers[residue][(IntVal-1)*nFlex+k], self.State)
                            cmd.refresh()

                            
                    # Get next column starting index
                    self.colNo = self.colNo + 11


            # Side-chain selection string - remove last 4 chars
            if strSelectSC != '':
                strSelectSC = strSelectSC[:len(strSelectSC)-4]

            return strSelectSC

        except:

            self.CriticalError("Could not update side-chain conformations")
    
    '''=========================================================================
      Get_AtomString: Retrives the PyMOL atom selection string
    ========================================================================='''
    def Get_AtomString(self, R, N, C, Atom):

        AtomString  = "resn " + R + " & resi " + N
        if C != '-':
            AtomString += " & chain " + C
        else:
            AtomString += " & chain ''"
        
        AtomString += " & name " + Atom
        AtomString += " & " + self.TargetObj

        return AtomString

    '''=========================================================================
      CriticalError: Abort simulation if there was a fatal error
    ========================================================================='''
    def CriticalError(self, text):

        print("  CRITICAL ERROR: " + text)
    
        #Create the .abort file
        abort_file = open(self.top.top.Manage.ABORT, 'w')
        abort_file.close()

        
