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
@title: FlexAID - FlexBonds.py

@summary: Permit to select the flexible bonds desired of the ligand.

@organization: Najmanovich Research Group
@creation date:  Sept. 16, 2010
'''

from pymol.wizard import Wizard
from pymol import cmd
from pymol.cgo import *
from pymol import util

import pymol
import General_cmd
import Geometry

class flexbond(Wizard):

    LigDisplay = 'FLEXIBLE_LIGAND__'
    PossFlexDisplay = 'POSS_FLEX_BONDS__'
    SelFlexDisplay = 'SELECTED_BONDS__'
    AtomDisplay = 'HIGHLIGHT_ATOM__'
    
    Translation = [1000,1000,1000]
    
    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top):
        
        #print "New instance of flexbond Class.\n"

        Wizard.__init__(self)

        self.top = top
        self.FlexAID = self.top.top
        self.FlexAID.WizardError = False
        
        self.RefLigand = self.FlexAID.IOFile.ReferencePath.get()

        self.View = cmd.get_view()
        self.State = cmd.get_state()
        
        self.pick_count = 0
        
        self.panelForceBond = 'Adv.: Force bond OFF'
        self.Force = False

        self.point1 = list()
        self.point2 = list()

        self.atom1 = list()
        self.atom2 = list()

    #=======================================================================
    ''' Executes the first steps of the Wizard'''
    #=======================================================================    
    def Start(self):

        self.ErrorCode = 1        

        try:
            #self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
            self.selection_mode = cmd.get("mouse_selection_mode")
            cmd.set("mouse_selection_mode", 0) # set selection mode to atomic

            # Mask objects
            self.exc = [self.LigDisplay, self.PossFlexDisplay, self.SelFlexDisplay]
            General_cmd.mask_Objects(self.exc)

            self.ErrorCode = 0

        except:
            self.top.DisplayMessage("  ERROR: Could not start the Flexible Bonds wizard", 1)
            self.top.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        # Display the ligand from the PDB file
        if self.DisplayLigand():
            self.top.DisplayMessage("  ERROR: Could not display the ligand", 1)
            self.top.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        # Display all Possible Flexible Bonds
        if self.show_FlexibleBonds():
            self.top.DisplayMessage("  ERROR: Could not display the flexible bonds", 1)
            self.top.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        # Display all Selected Flexible Bonds
        if self.show_SelectedBonds():
            self.top.DisplayMessage("  ERROR: Could not display the selected flexible bonds", 1)
            self.top.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return
        
        # remove any possible selection before selecting atoms
        cmd.deselect()
              

    #=======================================================================
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self):
        
        try:
            General_cmd.unmask_Objects(self.exc)
            cmd.set('mouse_selection_mode', self.selection_mode)

            cmd.delete(self.LigDisplay)
            cmd.delete(self.SelFlexDisplay)
            cmd.delete(self.PossFlexDisplay)         
            cmd.delete(self.AtomDisplay)         

            cmd.deselect()
        except:
            pass
      
        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True
        
        self.top.FlexBondsRunning(False)
        self.FlexAID.ActiveWizard = None

        cmd.set_wizard()
        cmd.set_view(self.View)

    #=======================================================================
    ''' Displays the ligand to be modified '''
    #=======================================================================    
    def DisplayLigand(self):
        
        try:
            cmd.load(self.RefLigand, self.LigDisplay, state=self.State)
            #print self.RefLigand
            
            # Display the atoms spheres
            cmd.show('spheres', self.LigDisplay)
            cmd.alter(self.LigDisplay,'vdw=0.25')
            cmd.rebuild()
        
            util.cbag(self.LigDisplay)
            cmd.translate(self.Translation,self.LigDisplay)
            cmd.zoom(self.LigDisplay)

        except:
            self.ErrorCode = 1

        return self.ErrorCode

    #=======================================================================
    ''' Reset all the variables '''
    #=======================================================================       
    def reset(self):
                
        Del_Down2 = 0

        #Unselect ALL the flexible bonds
        for index in iter(self.top.dictFlexBonds):
            self.top.dictFlexBonds[index][0] = 0 
            
            if self.top.dictFlexBonds[index][1]:
                if Del_Down2 == 0 or index < Del_Down2:
                    Del_Down2 = index
            
        # Remove forced bond(s) from dict
        if Del_Down2:
            for i in range(Del_Down2, len(self.top.dictFlexBonds)+1):
                del self.top.dictFlexBonds[i]
    
        self.show_SelectedBonds()

    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):
                
        Status = ''

        # Have to store the User Flexible bond selection...
        nBonds = 0
        nPoss = 0
        nBondsForced = 0
        for index in iter(self.top.dictFlexBonds):
            if self.top.dictFlexBonds[index][0]:
                #if self.top.dictFlexBonds[index][1]:
                #    nBondsForced += 1
                nBonds += 1
            if not self.top.dictFlexBonds[index][1]:
                nPoss += 1

        self.FlexAID.WizardResult = nBonds
        
        self.Quit_Wizard()
  
    #=======================================================================
    ''' Button Force Clicked: Force the selection of a NONE recognized bond '''
    #=======================================================================    
    def btn_Force(self):
        
        if self.Force:
            self.Force = False
            self.panelForceBond = 'Adv.: Force Bond OFF'
        else:
            #self.top.DisplayMessage('  FORCE the bond selection...',0)
            self.Force = True
            self.panelForceBond = 'Adv.: Force bond ON!'
        
        cmd.refresh_wizard()
        
          
     #=======================================================================
    ''' Show the possible flexible bonds '''
    #=======================================================================        
    def show_FlexibleBonds(self):

        point1 = list()
        point2 = list()

        try:
            View = cmd.get_view()

            PossFlexBonds = []

            for index in iter(self.top.dictFlexBonds):

                point1 = []
                point2 = []

                #print self.top.dictFlexBonds[index]
                #print self.top.dictNeighbours[self.top.dictFlexBonds[index][3]]

                # if bond is not Forced
                if not self.top.dictFlexBonds[index][1]:
                    # Get coordinates of 1st and 2nd neighbours
                    if int(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][0]) != 0 and \
                       int(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][1]) != 0 and \
                       int(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][2]) != 0:

                        if self.get_Coords(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][0], point1) or \
                           self.get_Coords(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][1], point2):
                            return 1

                    PossFlexBonds.extend(self.highlight_Possible(point1, point2))
                    #print PossFlexBonds
                    
            cmd.load_cgo(PossFlexBonds, self.PossFlexDisplay, state=self.State)   
            cmd.set_view(View)

        except:
            self.ErrorCode = 1

        return self.ErrorCode

     #=======================================================================
    ''' Show the selected flexible bonds '''
    #=======================================================================        
    def show_SelectedBonds(self):

        point1 = list()
        point2 = list()

        try:
            View = cmd.get_view()
            cmd.delete(self.SelFlexDisplay)            

            SelFlexBonds = []

            print self.top.dictFlexBonds
            for index in iter(self.top.dictFlexBonds):

                point1 = []
                point2 = []

                print self.top.dictFlexBonds[index]
                
                # if bond is flexible
                if self.top.dictFlexBonds[index][0]:
                    # Get coordinates of 1st and 2nd neighbours
                    if self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][0] != 0 and \
                       self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][1] != 0 and \
                       self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][2] != 0:

                        if self.get_Coords(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][0], point1) or \
                           self.get_Coords(self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][1], point2):
                            return 1

                    print "point1", point1
                    print "point2", point2
                    
                    SelFlexBonds.extend(self.highlight_Selected(point1, point2))
                    
            cmd.load_cgo(SelFlexBonds, self.SelFlexDisplay, state=self.State)   
            cmd.set_view(View)

        except:
            self.ErrorCode = 1

        return self.ErrorCode

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def btn_SelectAllBonds(self):
        
        for index in iter(self.top.dictFlexBonds):
            if not self.top.dictFlexBonds[index][1]:
                # Get coordinates of 1st and 2nd neighbours
                if self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][0] != 0 and \
                   self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][1] != 0 and \
                   self.top.dictNeighbours[self.top.dictFlexBonds[index][3]][2] != 0:

                    self.top.dictFlexBonds[index][0] = 1

        self.show_SelectedBonds()
        
    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
         [ 1, '* Flexible Bond Options *',''],
         [ 2, 'Select all flexible bonds','cmd.get_wizard().btn_SelectAllBonds()'],
         [ 2, 'Clear flexible bonds','cmd.get_wizard().reset()'],         
         [ 2, 'Show atom names','cmd.get_wizard().show_AtomsName()'],
         [ 2, 'Show atom IDs','cmd.get_wizard().show_AtomsNumber()'],
         [ 2, 'Hide labels','cmd.get_wizard().hide_Labels()'],
         #[ 2, self.panelForceBond,'cmd.get_wizard().btn_Force()'],        
         [ 2, 'Done','cmd.get_wizard().btn_Done()'],         
         ]
    

    #=======================================================================
    ''' Display the index of each atom of the ligand on the Pymol interface '''
    #=======================================================================
    def show_AtomsNumber(self):
        
        cmd.label(self.LigDisplay, "\"%d\" % ID")        
     
    #=======================================================================   
    ''' Display the name of each atom of the ligand on the Pymol interface '''
    #=======================================================================
    def show_AtomsName(self):
              
        cmd.label(self.LigDisplay, "\"%s\" % name")
    
    #=======================================================================   
    ''' Hide labels on the ligand '''
    #=======================================================================
    def hide_Labels(self):
              
        cmd.hide('labels', self.LigDisplay)

    #=======================================================================   
    ''' Highlight a selected flex bond '''
    #=======================================================================    
    def highlight_Selected(self, p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2

        # Draw white cylinder
        return [ CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2),
                 0.08, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000 ]

    #=======================================================================   
    ''' Highlight a possible flex bond '''
    #=======================================================================    
    def highlight_Possible(self, p1, p2):
        x1, y1, z1 = p1
        x2, y2, z2 = p2

        # Draw orange cylinder
        return [ CYLINDER, float(x1), float(y1), float(z1), float(x2), float(y2), float(z2), 
                 0.06, 0.800, 0.300, 0.000, 0.800, 0.300, 0.000 ]

    #=======================================================================   
    ''' Highlight atom upon clicking '''
    #=======================================================================    
    def highlight_Atom(self, atom):

        try:
            cmd.pseudoatom(self.AtomDisplay, pos=atom[5:], vdw=0.30, color='white')
            cmd.hide('nonbonded', self.AtomDisplay)
            cmd.show('spheres', self.AtomDisplay)
            cmd.mask(self.AtomDisplay)
        except:
            self.top.DisplayMessage("Failed to highlight atom upon selecting atom", 1)
            return
            
    #=======================================================================   
    ''' get x,y,z coordinates of an atom '''
    #=======================================================================    
    def get_Coords(self, atom_number, point):

        found = False

        try:
            atom_number = int(atom_number)

            atoms = cmd.get_model(self.LigDisplay, state=cmd.get_state())
            for at in atoms.atom:
                
                if General_cmd.get_ID(at.index, self.LigDisplay) == atom_number:
                    point.extend([at.coord[0],at.coord[1],at.coord[2]])
                    found = True
                    break

            if not found:
                return 1

        except:
            return 1

        return 0

    #=======================================================================   
    ''' gets atom information (coordinates and index)'''
    #=======================================================================    
    def get_Atom(self, name):

        info = []

        try: 
            list = cmd.index(name + " & " + self.LigDisplay)

            if len(list) > 0:
                atoms = cmd.get_model(name, state=cmd.get_state())
                for at in atoms.atom:
                    info.extend([ at.index, at.resn, at.resi, at.chain, at.name,
                                  at.coord[0], at.coord[1], at.coord[2] ])                  
            else:
                self.top.DisplayMessage("You must click in the object " + self.LigDisplay, 1)
                return info
                
            cmd.deselect()

        except:
            self.top.DisplayMessage("Error while retrieving atom info", 1)
            self.Quit_Wizard()
        
        return info
 
    #=======================================================================   
    ''' update the prompt '''
    #=======================================================================    
    def pickNextAtom(self):

        self.pick_count += 1

        # necessary to force update of the prompt
        cmd.refresh_wizard()

    #=======================================================================   
    ''' Checks if a force bond was already added '''
    #=======================================================================    
    def is_Already_Forced(self, atom_id):
        
        for index in iter(self.top.dictFlexBonds):
            if self.top.dictFlexBonds[index][1]:
                if atom_id in self.top.dictFlexBonds[index][3:]:
                    return index

        return 0

    #=======================================================================   
    ''' Forces a new bond in the dictionary'''
    #=======================================================================    
    def Force_NewBond(self, atom_id1, atom_id2):

        list = []
        list.extend([1, 1, 0])

        for atom in iter(self.top.dictNeighbours):
            
            if (atom_id1 == int(self.top.dictNeighbours[atom][0]) and \
                atom_id2 == int(self.top.dictNeighbours[atom][1])) or \
               (atom_id2 == int(self.top.dictNeighbours[atom][0]) and \
                atom_id1 == int(self.top.dictNeighbours[atom][1])):

                list[2] += 1
                list.append(atom)

        if not self.is_Already_Forced(list[3]):      
            Index = len(self.top.dictFlexBonds) + 1
            self.top.dictFlexBonds[Index] = list
        else:
            self.top.dictFlexBonds[self.is_Already_Forced(list[3])][0] = 0
            
    #=======================================================================   
    ''' Check if bond can be flexible using the neighbours'''
    #=======================================================================    
    def is_Definable(self, atom_id1, atom_id2):

        for atom in iter(self.top.dictNeighbours):
            
            if (atom_id1 == int(self.top.dictNeighbours[atom][0]) and \
                atom_id2 == int(self.top.dictNeighbours[atom][1])) or \
               (atom_id2 == int(self.top.dictNeighbours[atom][0]) and \
                atom_id1 == int(self.top.dictNeighbours[atom][1])):

                return 1

        return 0

    #=======================================================================   
    ''' Check bond can only be selected if Forced '''
    #=======================================================================    
    def is_Flexible(self, atom_id1, atom_id2):
        
        if self.is_Definable(atom_id1, atom_id2):
            
            for index in iter(self.top.dictFlexBonds):
                
                for i in range(3,3 + int(self.top.dictFlexBonds[index][2])):
                    if (atom_id1 == int(self.top.dictNeighbours[self.top.dictFlexBonds[index][i]][0]) and \
                        atom_id2 == int(self.top.dictNeighbours[self.top.dictFlexBonds[index][i]][1])) or \
                       (atom_id2 == int(self.top.dictNeighbours[self.top.dictFlexBonds[index][i]][0]) and \
                        atom_id1 == int(self.top.dictNeighbours[self.top.dictFlexBonds[index][i]][1])):

                        return index

        else:
            return -1

        return 0
            
    #=======================================================================   
    ''' Check if the bond is Valid '''
    #=======================================================================    
    def is_Valid(self, atom1, atom2):

        # 2 hydrogens
        if atom1[4].strip()[0] == 'H' and atom2[4].strip()[0] == 'H':
            return 0
        # Distance is too far
        elif Geometry.distance(atom1[5:],atom2[5:]) > 2.0:
            return 0

        self.FlexIndex = self.is_Flexible(General_cmd.get_ID(atom1[0], self.LigDisplay), 
                                          General_cmd.get_ID(atom2[0], self.LigDisplay))

        # The bond cannot be defined
        if self.FlexIndex == -1:
            return 0
        # The bond can only be defined using force
        elif self.FlexIndex == 0 and not self.Force:
            return 0

        return 1

    #=======================================================================
    ''' Pick an atom, then display the information related to it'''
    #=======================================================================
    def do_select(self, name):

        if self.pick_count < 1:
            self.atom1 = self.get_Atom(name)
            
            if len(self.atom1) > 0:
                if int(self.atom1[2]) == self.FlexAID.IOFile.ResSeq.get() and \
                        self.atom1[1] == 'LIG':
                    self.pickNextAtom()
                    self.highlight_Atom(self.atom1)
                else:
                    self.top.DisplayMessage("You can only select atoms in the object " + self.LigDisplay, 2)
            else:
                self.top.DisplayMessage("No atom could be selected in " + self.LigDisplay, 2)

        else:
            self.atom2 = self.get_Atom(name)

            if len(self.atom2) > 0:
                if int(self.atom2[2]) == self.FlexAID.IOFile.ResSeq.get() and \
                        self.atom2[1] == 'LIG':
                    self.pickNextAtom()
                    self.highlight_Atom(self.atom2)

                    if not self.is_Valid(self.atom1, self.atom2):
                        self.top.DisplayMessage("The selected bond cannot be defined as flexible.", 2)
                    else:
                        # Already existing Possible Bond
                        if self.FlexIndex != 0:
                            if self.top.dictFlexBonds[self.FlexIndex][0]:
                                self.top.dictFlexBonds[self.FlexIndex][0] = 0
                            else:
                                self.top.dictFlexBonds[self.FlexIndex][0] = 1

                        # Force new bond
                        else:
                            self.Force_NewBond(General_cmd.get_ID(self.atom1[0], self.LigDisplay), 
                                               General_cmd.get_ID(self.atom2[0], self.LigDisplay))

                    self.show_SelectedBonds()
                    self.pick_count = 0

                    cmd.refresh_wizard()
                    cmd.delete(self.AtomDisplay)                    

                else:
                    self.top.DisplayMessage("No atom could be selected in " + self.LigDisplay, 2)
            
            else:
                self.top.DisplayMessage("You can only select atoms in the object " + self.LigDisplay, 2)
                
    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):

        if self.pick_count == 0:
            return ['Please click on the FIRST atom that defines the flexible bond']
        elif self.pick_count == 1:
            return ['Please click on the SECOND atom that defines the flexible bond']
        elif self.pick_count == 2:
            return ['Analyzing bond...']
