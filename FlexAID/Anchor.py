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

from Tkinter import *
from pymol.wizard import Wizard
from pymol import cmd
from pymol.cgo import *
from pymol import util

import pymol
import General_cmd

class anchor(Wizard):
    
    LigDisplay = 'ANCHOR_LIGAND__'
    AtomDisplay = 'ANCHOR_ATOM__'

    Translation = [1000,1000,1000]
    
    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top, LigandPath, AnchorAtom):
        
        #print("New instance of anchor Class.\n")

        Wizard.__init__(self)

        self.top = top
        self.FlexAID = self.top.top
        self.FlexAID.WizardError = False
        
        self.LigandPath = LigandPath
        
        self.PrevAnchorAtom = AnchorAtom
        self.AnchorAtom = self.PrevAnchorAtom
        
        self.View = cmd.get_view()
        self.State = cmd.get_state()
        self.auto_zoom = cmd.get("auto_zoom")

        self.ErrorStatus =  []
        
    #=======================================================================
    ''' Executes the first steps of the Wizard'''
    #=======================================================================    
    def Start(self):

        self.ErrorCode = 1
        
        try:
            self.selection_mode = cmd.get("mouse_selection_mode")
            cmd.set("mouse_selection_mode", 0) # set selection mode to atomic

            # Mask objects
            self.exc = [self.LigDisplay]
            General_cmd.mask_Objects(self.exc)

            self.ErrorCode = 0
            
        except:
            self.FlexAID.DisplayMessage("  ERROR: Could not start the Anchor Atom wizard", 1)
            self.FlexAID.DisplayMessage("         The wizard will abort prematurely", 1)
            self.Quit_Wizard()
            return

        # Display the ligand from the PDB file
        if self.DisplayLigand() or self.RefreshDisplay():
            self.FlexAID.DisplayMessage("  ERROR: Could not display the ligand", 1)
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
            cmd.set('mouse_selection_mode', self.selection_mode)

            cmd.delete(self.LigDisplay)
            cmd.refresh()

            cmd.delete(self.AtomDisplay)
            cmd.refresh()
            
            cmd.deselect()
            
        except:
            pass
      
        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True
        
        self.FlexAID.WizardResult = self.AnchorAtom
        self.top.AnchorRunning(False)
        self.FlexAID.ActiveWizard = None

        cmd.set_wizard()
        cmd.set_view(self.View)
        
    #=======================================================================
    ''' Displays the ligand to change anchor '''
    #=======================================================================    
    def DisplayLigand(self):
        
        try:
            cmd.set("auto_zoom", 0)
            cmd.load(self.LigandPath, self.LigDisplay, state=self.State)
            cmd.refresh()

            cmd.translate(self.Translation,self.LigDisplay)
            cmd.refresh()

            cmd.zoom(self.LigDisplay)
            cmd.refresh()
            
        except:
            self.ErrorCode = 1

        cmd.set("auto_zoom", self.auto_zoom)
        
        return self.ErrorCode

    #=======================================================================
    ''' Refreshes the ligand to update anchor atom '''
    #=======================================================================    
    def RefreshDisplay(self):
        
        try:
            # Display the atoms spheres
            cmd.show('spheres', self.LigDisplay)
            cmd.refresh()

            cmd.alter(self.LigDisplay,'vdw=0.25')
            cmd.rebuild(self.LigDisplay)

            util.cbag(self.LigDisplay)
            cmd.refresh()
            
            if self.AnchorAtom != -1:
                AtomSel = self.LigDisplay + ' & id ' + str(self.AnchorAtom)
                cmd.color('white',AtomSel)
                cmd.refresh()

                cmd.alter(AtomSel ,'vdw=0.30')
                cmd.rebuild(AtomSel)
                        
        except:
            self.ErrorCode = 1

        return self.ErrorCode

    #=======================================================================
    ''' Clears the atom that will serve as anchor '''
    #=======================================================================    
    def clear(self):
        
        self.AnchorAtom = -1
        self.RefreshDisplay()
        
    #=======================================================================
    ''' Resets the atom to atom set before running wizard '''
    #=======================================================================    
    def reset(self):

        self.AnchorAtom = self.PrevAnchorAtom
        self.RefreshDisplay()

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
            [ 1, '* Anchor atom options *',''],
            [ 2, 'Clear anchor atom','cmd.get_wizard().clear()'],
            [ 2, 'Reset anchor atom','cmd.get_wizard().reset()'],                  
            [ 2, 'Show atom names','cmd.get_wizard().show_AtomsName()'],
            [ 2, 'Show atom IDs','cmd.get_wizard().show_AtomsNumber()'],
            [ 2, 'Hide labels','cmd.get_wizard().hide_Labels()'],
            [ 2, 'Done','cmd.get_wizard().btn_Done()'],                  
            ]

    #=======================================================================
    ''' Display the index of each atom of the ligand on the Pymol interface '''
    #=======================================================================
    def show_AtomsNumber(self):
        
        cmd.label(self.LigDisplay, "\"%d\" % ID")        
        cmd.refresh()

    #=======================================================================   
    ''' Display the name of each atom of the ligand on the Pymol interface '''
    #=======================================================================
    def show_AtomsName(self):
              
        cmd.label(self.LigDisplay, "\"%s\" % name")
        cmd.refresh()
    
    #=======================================================================   
    ''' Hide labels on the ligand '''
    #=======================================================================
    def hide_Labels(self):
              
        cmd.hide('labels', self.LigDisplay)
        cmd.refresh()

    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):
               
        self.Quit_Wizard()

    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):

        return self.ErrorStatus + \
                [ "Please click on a non-Hydrogen atom to define it as anchor." ]

    #=======================================================================
    ''' Pick an anchor atom then highlight it '''
    #=======================================================================
    def do_select(self, name):
        
        lt = cmd.index(name)
        
        for t in lt:
            if t[0] != self.LigDisplay:
                self.ErrorStatus = [ "You can only select an atom from the object " + self.LigDisplay + ". Try again." ]
            else:
                self.AnchorAtom = General_cmd.get_ID(t[1],t[0])
                self.RefreshDisplay()
            break
            
        cmd.deselect()
            