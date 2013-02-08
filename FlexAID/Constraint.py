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

import pymol
import General
import General_cmd

class constraint(Wizard):

    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top):
              
        Wizard.__init__(self)

        self.top = top
        self.FlexAID = self.top.top

        self.DisplayVar = BooleanVar()

        self.Active = self.top.ActiveCons
        self.dictConstraints = self.top.Vars.dictConstraints
        self.Scale = self.top.sclConsDist
        self.ScaleVar = self.top.ConsDist
        
        self.DefaultConsDist = '1.5'
        self.CONSTRAINT = 'CONSTRAINT_'
        self.ACTIVE = 'ACTIVE_CONS__'

        self.pick_count = 0
        self.ErrorStatus = [ "The active constraint is shown as a solid white line.",
                             "Use the scaler in the interface to edit the interaction distance." ]

        self.View = cmd.get_view()
        self.State = cmd.get_state()
        
        # Trigger controls state
        self.Active.set(self.Active.get())
        
        cmd.rebuild()

    #=======================================================================
    ''' Executes the first steps of the Wizard'''
    #=======================================================================    
    def Start(self):
    
        self.ErrorCode = 1

        try:
            self.selection_mode = cmd.get("mouse_selection_mode")
            cmd.set("mouse_selection_mode", 0) # set selection mode to atomic

            # Mask objects
            self.exc = [ self.FlexAID.IOFile.ProtName.get(), self.FlexAID.IOFile.LigandName.get() ]
            General_cmd.mask_Objects(self.exc)

            self.ErrorCode = 0

        except:
            self.FlexAID.DisplayMessage("  ERROR: Could not start the constraints wizard", 1)
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
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self):
        
        try:

            General_cmd.unmask_Objects(self.exc)
            cmd.set("mouse_selection_mode", self.selection_mode)
            
            cmd.delete(self.ACTIVE)
            cmd.delete(self.CONSTRAINT + '*')

            cmd.deselect()
        except:
            pass

        cmd.refresh()

        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True

        self.top.ConsRunning(False)
        self.FlexAID.ActiveWizard = None
        
        cmd.set_wizard()
        cmd.set_view(self.View)

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
        self.Active.set('')
        
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
            #atom1 [907, 'BTN', '300', 'A', 'O3', 14.692000389099121, -0.44600000977516174, -8.2049999237060547]
            #atom2 [234, 'SER', '45', 'A', 'OG', 14.543999671936035, 0.81300002336502075, -12.081999778747559]
            
            if atom1[3] == '': atom1[3] = '-'
            if atom2[3] == '': atom2[3] = '-'

            # example : '#907(O3) BTN300A'
            # example : '#234(OG) SER45A'

            leftsel = 'resn ' + atom1[1] + ' & resi ' + atom1[2]
            if atom1[3] != '-':
                leftsel += ' & chain ' + atom1[3]
            
            rightsel = 'resn ' + atom2[1] + ' & resi ' + atom2[2]
            if atom2[3] != '-':
                rightsel += ' & chain ' + atom2[3]

            leftkey = '#' + str(General_cmd.get_ID(atom1[0], leftsel)) + ' ' + atom1[1] + atom1[2] + atom1[3]
            rightkey = '#' + str(General_cmd.get_ID(atom2[0], rightsel)) + ' ' + atom2[1] + atom2[2] + atom2[3]
            key = leftkey + ' :: ' + rightkey
            altkey = rightkey + ' :: ' + leftkey
            
            if leftkey == rightkey:
                self.ErrorStatus = [ "Atoms selected are the same. Try again." ]

            elif key in self.dictConstraints or altkey in self.dictConstraints:
                self.ErrorStatus = [ "The selected constraint already exists. Try again." ]

            else:

                # Create new entry in dictionary
                distobj = self.CONSTRAINT + str(len(self.dictConstraints) + 1) + '__'
                self.dictConstraints[key] = [ leftkey, rightkey, 
                                              distobj,                             # constraint name
                                              [ atom1[5], atom1[6], atom1[7] ],    # coordinates of atom 1
                                              [ atom2[5], atom2[6], atom2[7] ],    # coordinates of atom 2
                                              self.DefaultConsDist ]               # value of slider
                
                if not self.Active.get():
                    self.Active.set(key)
                else:
                    self.refresh_display()

        except:            
            return 1

        return 0


    #=======================================================================
    ''' Highlights the active constraint '''
    #=======================================================================
    def highlight_Active(self):
        
        try:
            View = cmd.get_view()
            
            for key in self.dictConstraints:
                if key == self.Active.get():
                    
                    cons = self.dictConstraints[key]
                    
                    Highlight = []
                    Highlight.extend([ CYLINDER, 
                                       float(cons[3][0]), float(cons[3][1]), float(cons[3][2]), 
                                       float(cons[4][0]), float(cons[4][1]), float(cons[4][2]),
                                       0.05,
                                       1.000, 1.000, 1.000,
                                       1.000, 1.000, 1.000 ])

                    cmd.load_cgo(Highlight, self.ACTIVE)
                    
                    break

            cmd.set_view(View)
                
        except:
            return 1

        return 0

    #=======================================================================
    ''' Move up to the next active constraint '''
    #=======================================================================
    def Next_Active(self):
            
        Active = ''
        First = ''
        Next = False
        
        for key in sorted(self.dictConstraints, key=str.lower):
            
            if Next:
                Active = key
                break
            
            if not First:
                First = key
            
            if key == self.Active:
                Next = True
            
        if not Active:
            Active = First

        self.Active.set(Active)

        return 0
    
    #=======================================================================
    ''' Move up to the previous active constraint '''
    #=======================================================================
    def Previous_Active(self):
            
        Active = ''
        First = ''
        Next = False
        
        for key in reversed(sorted(self.dictConstraints, key=str.lower)):
            
            if Next:
                Active = key
                break
            
            if not First:
                First = key
            
            if key == self.Active:
                Next = True
            
        if not Active:
            Active = First

        self.Active.set(Active)

        return 0
    

    #=======================================================================
        ''' deletes the active constraint '''
    #=======================================================================
    def delete(self):

        if self.Active.get():
            del self.dictConstraints[self.Active.get()]
            self.Next_Active()

    #=======================================================================
        ''' refreshes the display of the constraints '''
    #=======================================================================
    def refresh_display(self):
        
        try:
            # Cleaning
            cmd.delete(self.ACTIVE)
            cmd.delete(self.CONSTRAINT + '*')
        except:
            pass
            
        try:
            # Create if not exists constraint
            for key in self.dictConstraints:
                distobj = self.dictConstraints[key][2]
                if not General_cmd.object_Exists(distobj) and self.create_cons(key):
                    return 1
            
            self.highlight_Active()
            cmd.refresh()
            
        except:
            self.ErrorCode = 1
        
        return self.ErrorCode
        
    #=======================================================================
    ''' Shows the constraint using the distance object '''
    #=======================================================================
    def create_cons(self, key):

        try:
            atom1 = self.top.parse_cons(self.dictConstraints[key][0])
            atom2 = self.top.parse_cons(self.dictConstraints[key][1])

            sel1 = "id " + str(atom1[0]) + " & resn " + atom1[1] + " & resi " + atom1[2]
            if atom1[3] != '-':
                sel1 += " & chain " + atom1[3]

            sel2 = "id " + str(atom2[0]) + " & resn " + atom2[1] + " & resi " + atom2[2]
            if atom2[3] != '-':
                sel2 += " & chain " + atom2[3]

            cmd.distance(self.dictConstraints[key][2], sel1, sel2)
            cmd.hide('labels', self.dictConstraints[key][2])
            cmd.color('white', self.dictConstraints[key][2])
            
        except:
            return 1

        return 0

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
