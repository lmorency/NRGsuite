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
@title: FlexAID - Constraint.py

@summary: Permits the addition of covalent/interaction constraints between atoms

@organization: Najmanovich Research Group
@creation date:  Sept. 16, 2010
'''

from Tkinter import *
from pymol.wizard import Wizard
from pymol import cmd
from pymol.cgo import *
from pymol import util

import pymol
import General
import General_cmd
import Geometry

class constraint(Wizard):

    CYLINDER_WIDTH = 0.06
    
    DefaultConsDist = '2.5'
    LigDisplay = 'CONSTRAINT_LIGAND__'
    MiddleDisplay = 'DISTANCE_OBJECT__'
    MiddleRadius = 0.05
    
    CONSTRAINT = 'CONS_'
    ACTIVE = 'ACTIVE_CONS__'

    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top):
              
        Wizard.__init__(self)

        self.top = top
        self.FlexAID = self.top.top

        self.RefLigand = self.FlexAID.IOFile.ReferencePath.get()

        self.ActiveCons = self.top.ActiveCons
        self.dictConstraints = self.top.Vars.dictConstraints
        
        self.pick_count = 0
        self.ErrorStatus = [ "The active constraint is shown as a solid white line.",
                             "Use the scaler in the interface to edit the interaction distance." ]

        self.View = cmd.get_view()
        self.State = cmd.get_state()
        self.auto_zoom = cmd.get("auto_zoom")

    #=======================================================================
    ''' Executes the first steps of the Wizard'''
    #=======================================================================    
    def Start(self):
    
        self.ErrorCode = 1

        try:
            self.selection_mode = cmd.get("mouse_selection_mode")
            cmd.set("mouse_selection_mode", 0) # set selection mode to atomic

            # Mask objects
            self.exc = [ self.FlexAID.IOFile.ProtName.get() ]
            General_cmd.mask_Objects(self.exc)

            self.ErrorCode = 0

        except:
            self.FlexAID.DisplayMessage("  ERROR: Could not start the constraints wizard", 1)
            self.FlexAID.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        if self.Display_Ligand():
            self.FlexAID.DisplayMessage("  ERROR: Could not display the ligand", 1)
            self.FlexAID.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        if self.refresh_display():
            self.FlexAID.DisplayMessage("  ERROR: Could not display the constraints", 1)
            self.FlexAID.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return
        
        # remove any possible selection before selecting atoms
        cmd.deselect()
              
    #=======================================================================
    ''' Displays the ligand to be modified '''
    #=======================================================================    
    def Display_Ligand(self):
        
        try:
            cmd.set("auto_zoom", 0)
            
            cmd.load(self.RefLigand, self.LigDisplay, state=self.State)
            cmd.refresh()

            # Display the atoms spheres
            cmd.show('spheres', self.LigDisplay)
            cmd.refresh()

            cmd.alter(self.LigDisplay,'vdw=0.25')
            cmd.rebuild(self.LigDisplay)
            
            util.cbag(self.LigDisplay)
            cmd.refresh()

            cmd.zoom(self.LigDisplay)            
            cmd.refresh()

        except:
            self.ErrorCode = 1

        cmd.set("auto_zoom", self.auto_zoom)
        
        return self.ErrorCode
        
    #=======================================================================
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self):
        
        try:

            General_cmd.unmask_Objects(self.exc)
            cmd.set("mouse_selection_mode", self.selection_mode)
            
            cmd.delete(self.ACTIVE)
            cmd.refresh()

            cmd.delete(self.MiddleDisplay)
            cmd.refresh()

            cmd.delete(self.CONSTRAINT + '*')
            cmd.refresh()

            cmd.delete(self.LigDisplay)
            cmd.refresh()

            cmd.deselect()
            cmd.unpick()
            
        except:
            pass

        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True

        self.top.ConsRunning(False)
        self.FlexAID.ActiveWizard = None
        
        cmd.set_view(self.View)
        cmd.set_wizard()

    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):
    
        self.FlexAID.WizardResult = len(self.dictConstraints)

        self.Quit_Wizard()

    #=======================================================================
    ''' Button Clear Clicked: Delete all constraints '''
    #=======================================================================       
    def btn_Clear(self):

        self.dictConstraints.clear()
        self.ActiveCons.set('')
        
    #=======================================================================
    ''' Resets to pick the first atom '''
    #=======================================================================
    def btn_Reset(self):

        self.pick_count = 0

        cmd.refresh_wizard()

    #=======================================================================   
    ''' gets atom information (coordinates and index)'''
    #=======================================================================    
    def get_Atom(self, name):

        info = []

        try: 
            list = cmd.index(name)
            listlen = len(list)

            if listlen == 1:
                atoms = cmd.get_model(name, state=self.State)
                for at in atoms.atom:
                    info.extend([ at.index, at.resn, at.resi, at.chain, at.name,
                                  at.coord[0], at.coord[1], at.coord[2] ])
            elif listlen > 1:
                self.FlexAID.DisplayMessage("  ERROR: Multiple atoms were selected",1)
                return info

            elif listlen == 0:
                self.FlexAID.DisplayMessage("  ERROR: No atoms could be selected on the current state",1)
                return info
                
            cmd.deselect()
            
        except:
            self.FlexAID.DisplayMessage("  ERROR: Could not retrieve atom info", 1)
            self.Quit_Wizard()
        
        return info

    #=======================================================================
    ''' Pick next atom in PyMOL Viewer '''
    #=======================================================================
    def pickNextAtom(self):

        self.pick_count += 1

        cmd.refresh_wizard()

    #=======================================================================
    ''' Executed upon selection in PyMOL '''
    #=======================================================================
    def do_select(self,name):

        if self.pick_count == 0:
            self.atom1 = self.get_Atom(name)
            
            if len(self.atom1) > 0:
                self.pickNextAtom()
            
        else:
            self.atom2 = self.get_Atom(name)
            
            if len(self.atom2) > 0:
                if self.AnalyzeConstraint(self.atom1,self.atom2):
                    self.ErrorStatus = [ "An unexpected error occured. Try again." ]

            self.pick_count = 0
            cmd.refresh_wizard()

    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):
        
        if self.pick_count == 0:
            return self.ErrorStatus + \
                    [ "Please click on the FIRST atom that defines the constraint." ]
        elif self.pick_count == 1:
            return [ "Please click on the SECOND atom that defines the constraint." ]
        
    #=======================================================================
    ''' Analyze if constraint can be added '''
    #=======================================================================
    def AnalyzeConstraint(self, atom1, atom2):
        
        try:
            #        0     1      2     3     4
            #atom1 [907, 'BTN', '300', 'A', 'O3', 14.692000389099121, -0.44600000977516174, -8.2049999237060547] + name
            #atom2 [234, 'SER', '45', 'A', 'OG', 14.543999671936035, 0.81300002336502075, -12.081999778747559]
            
            if atom1[3] == '': atom1[3] = '-'
            if atom2[3] == '': atom2[3] = '-'

            # example : '#907(O3) BTN300A'
            # example : '#234(OG) SER45A'

            leftsel = 'resn ' + atom1[1] + ' & resi ' + atom1[2]
            if atom1[3] != '-':
                leftsel += ' & chain ' + atom1[3]
                
            obj = self.Get_Prot_or_Lig(leftsel)
            if not obj:
                self.ErrorStatus = [ "You can only select atoms from the protein or the constraint ligand objects. Try again." ]
                
            else:
                leftsel += obj
                
                rightsel = 'resn ' + atom2[1] + ' & resi ' + atom2[2]
                if atom2[3] != '-':
                    rightsel += ' & chain ' + atom2[3]
                    
                obj = self.Get_Prot_or_Lig(rightsel)
                if not obj:
                    self.ErrorStatus = [ "You can only select atoms from the protein or the constraint ligand objects. Try again." ]
                    
                else:
                    rightsel += obj


                    leftkey  = '#' + str(General_cmd.get_ID(atom1[0], leftsel))
                    while len(atom1[1]) != 3:
                        atom1[1] = atom1[1] + '-'
                    leftkey += ' ' + atom1[1] + atom1[2] + atom1[3]

                    rightkey  = '#' + str(General_cmd.get_ID(atom2[0], rightsel))
                    while len(atom2[1]) != 3:
                        atom2[1] = atom2[1] + '-'
                    rightkey += ' ' + atom2[1] + atom2[2] + atom2[3]
                    
                    key = leftkey + ' :: ' + rightkey
                    altkey = rightkey + ' :: ' + leftkey
                    
                    if leftkey == rightkey:
                        self.ErrorStatus = [ "Atoms selected are the same. Try again." ]

                    elif len( [ cons for cons in self.dictConstraints if self.dictConstraints[cons][2] == key ] ) or \
                         len( [ cons for cons in self.dictConstraints if self.dictConstraints[cons][2] == altkey ] ):
                         
                        self.ErrorStatus = [ "The selected constraint already exists. Try again." ]

                    else:
                    
                        pointA = [ atom1[5], atom1[6], atom1[7] ]
                        pointB = [ atom2[5], atom2[6], atom2[7] ]
                        
                        dist = Geometry.distance(pointA, pointB)
                        if dist > self.top.MAX_DIST_CONSTRAINT:
                            dist = self.DefaultConsDist
                        
                        middle = Geometry.middle(pointA, pointB)
                        
                        # Create new entry in dictionary
                        distobj = self.get_ConstraintID()
                        
                        self.dictConstraints[distobj] = [   leftkey, rightkey,
                                                            key,                                 # constraint name
                                                            pointA,                              # coordinates of atom 1
                                                            pointB,                              # coordinates of atom 2
                                                            dist,                                # value of slider
                                                            middle                               # pseudoatom location
                                                        ]
                                                        
                        if not self.ActiveCons.get():
                            self.ActiveCons.set(distobj)
                        else:
                            self.refresh_display()

        except:            
            return 1
        
        return 0

    #=======================================================================
    ''' Gets the right object '''
    #=======================================================================
    def Get_Prot_or_Lig(self, sele):

        if cmd.count_atoms(sele + ' & ' + self.FlexAID.IOFile.ProtName.get()):
            return ' & ' + self.FlexAID.IOFile.ProtName.get()
        elif cmd.count_atoms(sele + ' & ' + self.LigDisplay):
            return ' & ' + self.LigDisplay
        else:
            return ''
    
    #=======================================================================
    ''' Highlights the active constraint '''
    #=======================================================================
    def highlight_Active(self):
        
        Error = 0

        
        try:
            # Cleaning
            cmd.delete(self.ACTIVE)
            cmd.refresh()
        except:
            pass
        
        try:
            cmd.set("auto_zoom", 0)
            
            for key in self.dictConstraints.keys():
                if key == self.ActiveCons.get():
                    
                    cons = self.dictConstraints[key]
                    
                    Highlight = []
                    Highlight.extend([ CYLINDER,
                                       float(cons[3][0]), float(cons[3][1]), float(cons[3][2]), 
                                       float(cons[4][0]), float(cons[4][1]), float(cons[4][2]),
                                       self.CYLINDER_WIDTH + 0.02,
                                       1.000, 1.000, 1.000,
                                       1.000, 1.000, 1.000 ])

                    cmd.load_cgo(Highlight, self.ACTIVE, state=self.State)
                    cmd.refresh()
                    
                    self.refresh_distance()
                    
                    break            
    
        except:
            Error = 1

        cmd.set("auto_zoom", self.auto_zoom)

        return Error

    #=======================================================================
    ''' Gets a unique constraint identifier '''
    #=======================================================================
    def get_ConstraintID(self):
        
        ID = 1
        CONSTRAINT_ID = self.CONSTRAINT + str(ID) + '__'
        
        while CONSTRAINT_ID in self.dictConstraints:
            ID = ID + 1
            CONSTRAINT_ID = self.CONSTRAINT + str(ID) + '__'
        
        return CONSTRAINT_ID
        
    #=======================================================================
    ''' Move up to the next active constraint '''
    #=======================================================================
    def Next_Active(self):
            
        Active = ''
        First = ''
        Next = False
        
        for key in sorted(self.dictConstraints.keys(), key=str.lower):
            
            if Next:
                Active = key
                break
            
            if not First:
                First = key
            
            if key == self.ActiveCons.get():
                Next = True
            
        if not Active:
            Active = First

        self.ActiveCons.set(Active)
    
    #=======================================================================
    ''' Move up to the previous active constraint '''
    #=======================================================================
    def Previous_Active(self):
            
        Active = ''
        First = ''
        Next = False
        
        for key in reversed(sorted(self.dictConstraints.keys(), key=str.lower)):
            
            if Next:
                Active = key
                break
            
            if not First:
                First = key
            
            if key == self.ActiveCons.get():
                Next = True
            
        if not Active:
            Active = First
            
        self.ActiveCons.set(Active)
    

    #=======================================================================
        ''' deletes the active constraint '''
    #=======================================================================
    def delete(self):

        if self.ActiveCons.get():
            del self.dictConstraints[self.ActiveCons.get()]
            self.Next_Active()

    #=======================================================================
        ''' refreshes the display of the constraints '''
    #=======================================================================
    def refresh_display(self):
        
        try:
            cmd.delete(self.CONSTRAINT + '*')
            cmd.refresh()
        except:
            pass
            
        try:
            # Create if not exists constraint
            for key in self.dictConstraints.keys():
                distobj = self.dictConstraints[key][2]
                if not General_cmd.object_Exists(distobj) and self.create_cons(key):
                    return 1
            
            self.highlight_Active()
            
        except:
            self.ErrorCode = 1
        
        return self.ErrorCode
        
    #=======================================================================
        ''' refreshes the distance label of the distance object '''
    #=======================================================================
    def refresh_distance(self):

        try:
            cmd.delete(self.MiddleDisplay)
            cmd.refresh()
        except:
            pass
        
        try:
            cmd.set("auto_zoom", 0)
            
            Active = self.ActiveCons.get()
            
            cmd.pseudoatom(self.MiddleDisplay,
                           pos=self.dictConstraints[Active][6],
                           vdw=self.MiddleRadius,
                           state=self.State)
            cmd.refresh()
            
            cmd.hide('nonbonded', self.MiddleDisplay)
            cmd.refresh()
            
            dist = '%.2f' % self.dictConstraints[Active][5]
            cmd.label(self.MiddleDisplay, str(dist))
            cmd.refresh()
            
        except:
            pass
        
        cmd.set("auto_zoom", self.auto_zoom)
        
    #=======================================================================
    ''' Shows the constraint using the distance object '''
    #=======================================================================
    def create_cons(self, key):

        Error = 0
        
        try:
        
            cmd.set("auto_zoom", 0)
            
            cons = self.dictConstraints[key]

            Highlight = []
            Highlight.extend([ CYLINDER,
                               float(cons[3][0]), float(cons[3][1]), float(cons[3][2]),
                               float(cons[4][0]), float(cons[4][1]), float(cons[4][2]),
                               self.CYLINDER_WIDTH,
                               0.800, 0.300, 0.000,
                               0.800, 0.300, 0.000 ])

            cmd.load_cgo(Highlight, key)
            cmd.refresh()
            
        except:
            Error = 1

        cmd.set("auto_zoom", self.auto_zoom)

        return Error

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
         [ 1, '* Constraint Options *',''],
         [ 2, 'Go to next active','cmd.get_wizard().Next_Active()'],
         [ 2, 'Go to previous active','cmd.get_wizard().Previous_Active()'],
         [ 2, 'Delete active','cmd.get_wizard().delete()'],
         [ 2, 'Clear constraints', 'cmd.get_wizard().btn_Clear()'],
         [ 2, 'Reset', 'cmd.get_wizard().btn_Reset()'],
         [ 2, 'Done','cmd.get_wizard().btn_Done()'],         
         ]
