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

class constraint(Wizard):

    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top, type):
              
        Wizard.__init__(self)

        self.top = top
        self.FlexAID = self.top.top

        self.DisplayVar = BooleanVar()

        if type == 'Covalent':
            self.dictConstraints = self.top.dictCovConstraints
            self.optConstraints = self.top.optCovCons
            self.ConsSelection = self.top.CovConsSelection
            self.Scale = self.top.sclCovDist
            self.ScaleVar = self.top.CovDist
            self.Check = None
            self.CheckVar = IntVar()   # dummy Var
            self.DefaultValue = self.top.DefaultCovDist
            self.CONSTRAINT = 'COVALENT_CONS__'
        else:
            self.dictConstraints = self.top.dictIntConstraints
            self.optConstraints = self.top.optIntCons
            self.ConsSelection = self.top.IntConsSelection
            self.Scale = self.top.sclIntFactor
            self.ScaleVar = self.top.IntFactor
            self.Check = self.top.chkIntForced
            self.CheckVar = self.top.IntForced
            self.DefaultValue = self.top.DefaultIntFactor
            self.CONSTRAINT = 'INTERACT_CONS__'

        self.DisplayVarTrace = self.DisplayVar.trace('w', self.DisplayVar_Toggle)
        self.ScaleVarTrace = self.ScaleVar.trace('w', self.ScaleVar_Toggle)
        self.CheckVarTrace = self.CheckVar.trace('w', self.CheckVar_Toggle)
        self.ConsSelectionTrace = self.ConsSelection.trace('w', self.ConsSelection_Toggle)

        self.HIGHLIGHT = self.top.HIGHLIGHT

        self.DisplayVar.set(True)
        self.ErrorStatus = ''

        self.pick_count = 0

        self.View = cmd.get_view()
        cmd.rebuild()        


    #=======================================================================
    ''' Toggle between displaying or not constraints '''
    #=======================================================================    
    def DisplayVar_Toggle(self, *args):

        if self.DisplayVar.get():
            self.DisplayText = 'Hide all'
        else:
            self.DisplayText = 'Show all'

        self.refresh_display()

    #=======================================================================
    ''' Highlight the right selection '''
    #=======================================================================    
    def ConsSelection_Toggle(self, *args):

        self.refresh_display()

    #=======================================================================
    ''' Sets the value of slider of selected key '''
    #=======================================================================    
    def ScaleVar_Toggle(self, *args):

        self.dictConstraints[self.ConsSelection.get()][5] = self.ScaleVar.get()

    #=======================================================================
    ''' Sets the value of check of selected key '''
    #=======================================================================    
    def CheckVar_Toggle(self, *args):

        self.dictConstraints[self.ConsSelection.get()][6] = self.CheckVar.get()
        
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
            self.FlexAID.DisplayMessage("  ERROR: Could not start the Constraint wizard", 1)
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
                
            cmd.delete(self.CONSTRAINT + '*')
            cmd.delete(self.HIGHLIGHT)

            cmd.deselect()
        except:
            pass

        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True

        self.top.ConsRunning(self.dictConstraints, self.optConstraints, self.ConsSelection, False)

        self.ScaleVar.trace_vdelete('w',self.ScaleVarTrace)
        self.CheckVar.trace_vdelete('w',self.CheckVarTrace)
        self.ConsSelection.trace_vdelete('w',self.ConsSelectionTrace)

        self.FlexAID.ActiveWizard = None
        
        cmd.set_wizard()
        cmd.set_view(self.View)

    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):

        self.Quit_Wizard()
        
    #=======================================================================
    ''' Button Display: Toggles between on/off displaying constraints '''
    #=======================================================================       
    def btn_display(self):

        if self.DisplayVar.get():
            self.DisplayVar.set(False)
        else:
            self.DisplayVar.set(True)

        cmd.refresh_wizard()

    #=======================================================================
    ''' Button Reset Clicked: Reset atom selection '''
    #=======================================================================       
    def btn_reset(self):

        self.ErrorStatus = ''
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
                atoms = cmd.get_model(name, state=cmd.get_state())
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
                    self.ErrorStatus = "An unexpected error occured.\n"
            
                self.btn_reset()
     
    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):

        if self.pick_count == 0:
            return [self.ErrorStatus + 
                    'Please click on the FIRST atom that defines the constraint.']
        elif self.pick_count == 1:
            return ['Please click on the SECOND atom that defines the constraint.']
        

    #=======================================================================
    ''' Counts the total number of constraints '''
    #=======================================================================
    def count_cons(self):
        
        tot = 0

        for key in self.dictConstraints:
            for item in self.dictConstraints[key]:
                tot = tot + 1
        
        return tot

    #=======================================================================
    ''' Analyze if constraint can be added '''
    #=======================================================================
    def AnalyzeConstraint(self, atom1, atom2):
        
        try:
            #        0     1      2     3     4
            #atom1 [907, 'BTN', '300', 'A', 'O3', 14.692000389099121, -0.44600000977516174, -8.2049999237060547]
            #atom2 [234, 'SER', '45', 'A', 'OG', 14.543999671936035, 0.81300002336502075, -12.081999778747559]

            #print atom1
            #print atom2

            if atom1[3] == '': atom1[3] = '-'
            if atom2[3] == '': atom2[3] = '-'

            # example : '#907(O3) BTN300A'
            # example : '#234(OG) SER45A'
            #leftkey = '#' + str(atom1[0]) + '(' + atom1[4] + ') ' + atom1[1] + atom1[2] + atom1[3]
            #rightkey = '#' + str(atom2[0]) + '(' + atom2[4] + ') ' + atom2[1] + atom2[2] + atom2[3]        

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

            #print leftsel
            #print rightsel
            #print leftkey
            #print rightkey

            if leftkey == rightkey:
                self.ErrorStatus = "Atoms selected are the same. Try again.\n"

            elif key in self.dictConstraints or altkey in self.dictConstraints:
                self.ErrorStatus = "Constraint already exists. Try again.\n"

            else:

                # Add item in drop-down-list
                self.optConstraints["menu"].add_command(label=key, command=lambda temp = key: self.optConstraints.setvar(self.optConstraints.cget("textvariable"), value = temp))
                self.ConsSelection.set(key)

                # Create new entry in dictionary
                self.dictConstraints[key] = [ leftkey, rightkey, 
                                              self.CONSTRAINT + str(len(self.dictConstraints) + 1),  # constraint name
                                              [ atom1[5], atom1[6], atom1[7] ],    # coordinates of atom 1
                                              [ atom2[5], atom2[6], atom2[7] ],    # coordinates of atom 2
                                              0.25,                                 # value of slider
                                              0 ]                                  # checkbox value
                
                # Set default value of slider
                self.ScaleVar.set(self.DefaultValue)

                # Set default value of checkbox
                if self.Check != None:
                    self.CheckVar.set(0)

                # Create object
                if self.create_cons(key):
                    return 1
                
                # Highlight newly created object
                self.highlight_Active(self.ConsSelection)

        except:
            return 1

        return 0


    #=======================================================================
    ''' Disable all constraints of type i '''
    #=======================================================================
    def highlight_Active(self, selection):
        
        try:

            View = cmd.get_view()

            cmd.delete(self.HIGHLIGHT)
            #cmd.color('orange', self.CONSTRAINT + '*')
            
            self.Scale.config(state='disabled')

            if self.Check != None:
                self.Check.config(state='disabled')
            
            for key in self.dictConstraints:

                if key == selection.get():
                    ref = self.dictConstraints[key]
                    # Draw white cylinder

                    Highlight = []
                    Highlight.extend([ CYLINDER, 
                                       float(ref[3][0]), float(ref[3][1]), float(ref[3][2]), 
                                       float(ref[4][0]), float(ref[4][1]), float(ref[4][2]),
                                       0.05, 
                                       1.000, 1.000, 1.000, 
                                       1.000, 1.000, 1.000 ])

                    cmd.load_cgo(Highlight, self.HIGHLIGHT)
                    cmd.refresh()

                    self.Scale.config(state='normal')
                    self.ScaleVar.set(ref[5])

                    if self.Check != None:
                        self.Check.config(state='normal')
                        self.CheckVar.set(ref[6])
            
                    cmd.set_view(View)
                    return 0
                
            # Default values of blank item ''
            self.ScaleVar.set(0.25)
            if self.Check != None:
                self.CheckVar.set(0)

        except:
            return 1

        return 0

    #=======================================================================
    ''' Deletes the active constraint '''
    #=======================================================================
    def delete(self):
        
        self.top.Del_Constraint(self.dictConstraints,self.optConstraints,self.ConsSelection)

    #=======================================================================
        ''' refreshes the display of the constraints '''
    #=======================================================================
    def refresh_display(self):
        
        try:
            
            # Cleaning
            cmd.delete(self.HIGHLIGHT)
            cmd.delete(self.CONSTRAINT + '*')

            # Create if not exists constraint
            for key in self.dictConstraints:
                obj = self.dictConstraints[key][2]
                if not General_cmd.object_Exists(obj) and self.create_cons(key):
                    return 1

                if not self.DisplayVar.get():
                    cmd.disable(obj)
                else:
                    cmd.enable(obj)

            self.highlight_Active(self.ConsSelection)

        except:
            self.ErrorCode = 1
        
        return self.ErrorCode
        
    #=======================================================================
    ''' Disable all constraints of type i '''
    #=======================================================================
    def hide_all(self):
        
        try:
            #cmd.disable(self.HIGHLIGHT)
            cmd.disable(self.CONSTRAINT + '*')
        except:
            return 1
        
        return 0

    #=======================================================================
    ''' Loop through contraints of type i and show '''
    ''' if not exists, create new object           '''
    #=======================================================================
    def show_all(self):
        
        try:
            
            cmd.delete(self.CONSTRAINT + '*')

            # Create if not exists
            for key in self.dictConstraints:
                if not General_cmd.object_Exists(self.dictConstraints[key][2]) and self.create_cons(key):
                    return 1

            # Enable all constraints
            cmd.enable(self.CONSTRAINT + '*')

            self.highlight_Active(self.ConsSelection)

        except:
            return 1
        
        return 0
                    
    #=======================================================================
    ''' Shows the constraint using the distance object '''
    #=======================================================================
    def create_cons(self, key):

        try:
            atom1 = parse_cons(self.dictConstraints[key][0])
            atom2 = parse_cons(self.dictConstraints[key][1])

            sel1 = "id " + str(atom1[0]) + " & resn " + atom1[1] + " & resi " + atom1[2]
            if atom1[3] != '-':
                sel1 += " & chain " + atom1[3]

            sel2 = "id " + str(atom2[0]) + " & resn " + atom2[1] + " & resi " + atom2[2]
            if atom2[3] != '-':
                sel2 += " & chain " + atom2[3]

            #print sel1
            #print sel2

            cmd.distance(self.dictConstraints[key][2], sel1, sel2)
            cmd.hide('labels', self.dictConstraints[key][2])
            #cmd.color('orange', self.dictConstraints[key][2])
            
        except:
            return 1

        return 0

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
         [ 1, '* Constraint Options *',''],
         #[ 2, 'Show all','cmd.get_wizard().show_all()'],
         #[ 2, 'Hide all','cmd.get_wizard().hide_all()'],
         [ 2, self.DisplayText,'cmd.get_wizard().btn_display()'],
         [ 2, 'Delete active','cmd.get_wizard().delete()'],
         [ 2, 'Reset selection','cmd.get_wizard().btn_reset()'],         
         [ 2, 'Done','cmd.get_wizard().btn_Done()'],         
         ]

#=======================================================================
''' Parse constraints and returns an array with each element '''
#=======================================================================
def parse_cons(key):

    l = list()

    #l.append(key[1:key.find('(')])
    l.append(key[1:key.find(' ')])

    #st = key.find(')')+2
    st = key.find(' ') + 1
    rnc = key[st:]

    l.append(rnc[0:3])
    l.append(rnc[3:len(rnc)-1])
    l.append(rnc[len(rnc)-1:len(rnc)])

    return l
