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

'''
@title: FlexAID - FlexSideChain.py

@summary: Permit to select the residu in the protein that will be used for the 
          Flexible Side Chain option.

@organization: Najmanovich Research Group
@creation date:  March 30, 2010
'''

from pymol.wizard import Wizard
from pymol import cmd

import pymol
import General_cmd

class flexSC(Wizard):

    FlexSCDisplay = 'FLEXIBLE_SIDE_CHAINS__'
    ResidueDisplay = 'HIGHLIGHT_RESIDUE__'

    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top, queue):
        
        #print "New instance of FlexSC Wizard"

        Wizard.__init__(self)
        self.top = top
        self.queue = queue
        self.FlexAID = self.top.top

        self.FlexAID.WizardError = False

        self.TargetName = self.FlexAID.IOFile.Target
        self.TargetFlex = self.top.Vars.TargetFlex

        self.ResidueName = '...'
        self.PanelResidue = '   Residue: '
        
        # Initial view
        self.View = cmd.get_view()
        self.State = cmd.get_state()
        self.auto_zoom = cmd.get("auto_zoom")
        
        # for Quit_Wizard
        self.ErrorCode = 1
        
        self.ErrorStatus = [ "Side-chains shown in white are already set as flexible.",
                             "Side-chains shown in orange are 'active' and are ready to be added." ]

        
    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):
        
        self.Quit_Wizard()
        
    #=======================================================================
    ''' Executes the first steps of the Wizard'''
    #=======================================================================    
    def Start(self):
        
        cmd.window('hide')
        self.queue.put(lambda: cmd.window('show'))
        cmd.refresh_wizard()
        # Display all Selected Flexible Bonds
        if self.show_SelectedSC():
            self.queue.put(lambda: self.FlexAID.DisplayMessage("  ERROR: Could not display selected flexible side-chains", 1))
            self.Quit_Wizard()
            return

        #self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        self.selection_mode = cmd.get("mouse_selection_mode")
        cmd.set("mouse_selection_mode", 1) # set selection mode to residue

        # Mask objects
        self.exc = [ self.TargetName ]
        General_cmd.mask_Objects(self.exc)
        cmd.zoom(self.TargetName)

        # remove any possible selection before selecting atoms
        cmd.deselect()
              
        self.ErrorCode += 1
        
        
    #=======================================================================
        ''' update_SelectedSC : updates the view of selected flexible sc. ''' 
    #=======================================================================       
    def show_SelectedSC(self):
        
        selString = ''

        try:
            cmd.delete(self.FlexSCDisplay)
            cmd.refresh()
        except:
            pass
        
        for sc in self.TargetFlex.listSideChain:
            
            resn = sc[0:3]
            resi = sc[3:3+len(sc)-4]
            chain = sc[len(sc)-1:len(sc)]

            if selString != '':
                selString += ' or '

            selString += ' ! n. C+O+N '
            selString += ' & resn ' + str(resn)
            selString += ' & resi ' + str(resi)
            
            if chain != '-':
                selString += ' & c. ' + str(chain)

            selString += ' & ' + self.TargetName
        
        
        if selString != '' and self.highlight_FlexibleSC(selString):
            return 1

        return 0

    #=======================================================================
    ''' Button Add Residue Clicked: Add the residue in the Flexible Side-Chains
                                   DDL in FlexAID '''
    #=======================================================================       
    def btn_AddResidue(self):

        self.TargetFlex.Add_SideChain(self.ResidueName)

        self.show_SelectedSC()
        self.btn_Reset()
        
    #=======================================================================
    ''' Button Add Residue Clicked: Add the residue in the Flexible Side-Chains
                                   DDL in FlexAID '''
    #=======================================================================       
    def btn_DelResidue(self):
                
        self.TargetFlex.Remove_SideChain(self.ResidueName)
        
        self.show_SelectedSC()
        self.btn_Reset()

    #=======================================================================
    ''' Button Clear Residues Clicked: Clears the side chains list '''
    #=======================================================================       
    def btn_ClearResidue(self):
                
        self.TargetFlex.Clear_SideChain()
        
        self.show_SelectedSC()
        self.btn_Reset()
        
    #=======================================================================
    ''' Button Reset Residues Clicked: Resets the highlighted '''
    #=======================================================================       
    def btn_Reset(self):
        
        self.ResidueName = '...'

        try:
            cmd.delete(self.ResidueDisplay)
            cmd.refresh()
        except:
            pass

        cmd.refresh_wizard()
        
    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):
     
        return self.ErrorStatus + \
                [ "Click on a side-chain of the protein object." ]

    #=======================================================================
    ''' Residue selection, then display the information related to it '''
    #=======================================================================
    def do_select(self, name):

        self.atom = self.get_Atom(name)

        if self.atom:

            resn = self.atom[1]
            resi = self.atom[2]
            chain = self.atom[3]

            if chain == '':
                chain = '-'

            self.ResidueName = resn + str(resi) + chain
            
            Error = self.top.Validate_EnterResidue(self.ResidueName)
            
            if Error != '':
                self.ErrorStatus = [ Error ]
                self.btn_Reset()
                
            else:
                self.ErrorStatus = [ ]
                self.highlight_Residue(name)
        
        cmd.deselect()
        cmd.refresh_wizard()

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
         [ 1, '* Flexible Side-chains *',''],
         [ 3, self.PanelResidue + self.ResidueName,''],
         [ 2, 'Add the side-chain','cmd.get_wizard().btn_AddResidue()'],
         [ 2, 'Delete the side-chain','cmd.get_wizard().btn_DelResidue()'],
         [ 2, 'Clear the side-chains','cmd.get_wizard().btn_ClearResidue()'],
         [ 2, 'Reset','cmd.get_wizard().btn_Reset()'],
         [ 2, 'Done','cmd.get_wizard().btn_Done()'],
         ]

    #=======================================================================
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self):
        
        try:
            
            #Delete the Residue objects
            cmd.delete(self.FlexSCDisplay)
            cmd.refresh()

            cmd.delete(self.ResidueDisplay)
            cmd.refresh()

            #cmd.delete(self.BackboneDisplay)
            #cmd.refresh()
            
        except:
            pass
            
        if self.ErrorCode != 1:
            General_cmd.unmask_Objects(self.exc)
            cmd.set('mouse_selection_mode', self.selection_mode)
            #cmd.config_mouse('three_button_editing', 1)

        if self.ErrorCode > 0:
            self.FlexAID.WizardError = True

        self.FlexAID.WizardResult = self.TargetFlex.Count_SideChain()

        self.queue.put(lambda: self.top.FlexSCRunning(False))
        self.FlexAID.ActiveWizard = None
        self.queue.put(lambda: self.FlexAID.root.lift())
        cmd.set_wizard()
        cmd.refresh()

        cmd.set_view(self.View)
        
     #=======================================================================   
    ''' gets atom information (coordinates and index)'''
    #=======================================================================    
    def get_Atom(self, name):

        info = []

        try:
            # test if the click is on the protein
            list = cmd.index(name + " & " + self.TargetName)

            if len(list) > 0:
                atoms = cmd.get_model(name, state=1)
                for at in atoms.atom:
                    info.extend([ at.index, at.resn, at.resi, at.chain, at.name,
                                  at.coord[0], at.coord[1], at.coord[2] ])
                    break

            else:
                self.ErrorStatus = [ "You can only select side-chains from the object " + self.TargetName + ". Try again." ]
                return info

            cmd.deselect()

        except:
            self.ErrorStatus = [ "An error occured while retrieving atom info. Try again." ]
            self.Quit_Wizard()
        
        return info

    #=======================================================================   
    ''' highlight_Residue: Highlight residue upon selection '''
    #=======================================================================    
    def highlight_Residue(self, name):

        try:            
            cmd.delete(self.ResidueDisplay)
            cmd.refresh()
        except:
            pass
            
        try:
            cmd.set("auto_zoom", 0)
            
            # Create new object from selection
            cmd.create(self.ResidueDisplay, name + ' & ! n. C+O+N', target_state=self.State)
            cmd.refresh()

            # Visual appearance
            cmd.hide('lines', self.ResidueDisplay)
            cmd.refresh()

            cmd.show('sticks', self.ResidueDisplay)
            cmd.refresh()

            cmd.color('orange', self.ResidueDisplay)
            cmd.refresh()

            cmd.mask(self.ResidueDisplay)

            # Toggle FlexSC obj to overlap ResidueDisplay
            cmd.disable(self.FlexSCDisplay)
            cmd.refresh()

            cmd.enable(self.FlexSCDisplay)
            cmd.refresh()
            
        except:
            self.queue.put(lambda: self.FlexAID.DisplayMessage("  ERROR: Could not highlight residue upon selection", 2))

        cmd.set("auto_zoom", self.auto_zoom)

    #=======================================================================   
    ''' highlight_FlexibleSC: Highlight flexible side-chains '''
    #=======================================================================    
    def highlight_FlexibleSC(self, selString):

        Error = 0
        
        try:
            cmd.delete(self.FlexSCDisplay)
            cmd.refresh()
        except:
            pass
            
        try:
            cmd.set("auto_zoom", 0)

            # Create new object from selection
            cmd.create(self.FlexSCDisplay, selString, target_state=self.State)
            cmd.refresh()
            
            # Visual appearance
            cmd.hide('lines', self.FlexSCDisplay)
            cmd.refresh()

            cmd.show('sticks', self.FlexSCDisplay)
            cmd.refresh()

            cmd.color('white', self.FlexSCDisplay)
            cmd.refresh()

            cmd.label(self.FlexSCDisplay + " & name CA", "resn+resi")

            cmd.mask(self.FlexSCDisplay)
            
        except:
            Error = 1
        
        cmd.set("auto_zoom", self.auto_zoom)

        return Error
