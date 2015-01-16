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
@title: FlexAID - AtomTypes.py

@summary: Permit to select to change the atom type from a wizard menu in pymol.

@organization: Najmanovich Research Group
@creation date:  April 28, 2011
'''

from pymol.wizard import Wizard
from pymol import cmd
from pymol.cgo import *
from pymol import util

import pymol
import General


class setType(Wizard):
    
    AtomDisplay = 'HIGHLIGHT_ATOM__'
    LigDisplay = 'ATOM_TYPES_LIGAND__'

    Translation = [1000,1000,1000]
    
    #=======================================================================
    ''' Initialization of the interface '''
    #=======================================================================
    def __init__(self, top, queue, UseOldTypes):
        
        #print "New instance of Wizard SetAtomType"

        Wizard.__init__(self)
        
        self.top = top
        self.FlexAID = self.top.top
        self.queue = queue
        # Save view
        self.View = cmd.get_view()
        self.State = cmd.get_state()
                
        #self.pdbPath = self.FlexAID.SimPath + 'tmp.pdb'
        self.pdbPath = self.FlexAID.IOFile.ProcessedLigandPath.get()
        
        self.ErrorCode = 1

        self.SelAtomName = '...'   
        self.SelAtomType = '...'

        self.Set_DDL_Colors(UseOldTypes)       

        smm = []
        smm.append([ 2, 'Type Selection', '' ])
        for a in self.types:
            smm.append([ 1, self.type_name[a], 'cmd.get_wizard().set_AtomType("' + a + '")'])
        self.menu['type']=smm
       
        self.atom = list()
        self.ID = ''

    #=======================================================================
    ''' Sets the number of colors and items depending on number of types '''
    #=======================================================================
    def Set_DDL_Colors(self, UseOldTypes):

        if UseOldTypes:

            self.types = [
                '1',
                '2',
                '3',
                '4',
                '5',
                '6',
                '7',
                '8',
                ]

            self.type_name = {
                '1':'\\099 1-Hydrophilic',
                '2':'\\009 2-Acceptor',
                '3':'\\900 3-Donor',
                '4':'\\090 4-Hydrophobic',
                '5':'\\990 5-Aromatic',
                '6':'\\555 6-Neutral',
                '7':'\\734 7-Neutral-Donor',
                '8':'\\559 8-Neutral-Acceptor',
                }

            self.type_color = {
                '1':'cyan',
                '2':'blue',
                '3':'red',
                '4':'green',
                '5':'yellow',
                '6':'grey',
                '7':'raspberry',
                '8':'slate',                          
                }
            
        else:

            self.types = [
                '1',
                '2',
                '3',
                '4',
                '5',
                '6',
                '7',
                '8',
                '9',
                '10'
                ]

            self.type_name = {
                '1':'\\099 1-Hydrophilic',
                '2':'\\009 2-Acceptor',
                '3':'\\900 3-Donor',
                '4':'\\090 4-Hydrophobic',
                '5':'\\990 5-Aromatic',
                '6':'\\555 6-Neutral',
                '7':'\\734 7-Electrorepulsor',
                '8':'\\559 8-Electroattractor',
                '9':'\\905 9-Positive',
                '10':'\\987 10-Negative',
                }

            self.type_color = {
                '1':'cyan',
                '2':'blue',
                '3':'red',
                '4':'green',
                '5':'yellow',
                '6':'grey',
                '7':'raspberry',
                '8':'slate',                          
                '9':'pink',
                '10':'wheat',
                }

    #=======================================================================
    ''' Executes the first steps of the Wizard '''
    #=======================================================================
    def Start(self):
        
        self.queue.put(lambda: self.FlexAID.root.withdraw())
        cmd.window('hide')
        self.queue.put(lambda: cmd.window('show'))
        cmd.refresh_wizard()

        # Display the ligand from the PDB file
        if self.DisplayLigand():
            self.FlexAID.DisplayMessage("Error while trying to display the ligand from PDB file", 1)
            self.Quit_Wizard(self.ErrorCode)
            return
                    
        # Color the Atoms based on their type
        if self.ColorByType():
            self.FlexAID.DisplayMessage("Error while coloring the ligand atoms by type", 1)
            self.Quit_Wizard(self.ErrorCode)
            return            

        # Mask everythign except modifyable object (ligand)
        self.exc = [self.LigDisplay]
        General.mask_Objects(self.exc)        

        self.selection_mode = cmd.get("mouse_selection_mode")
        cmd.set("mouse_selection_mode", 0) # set selection mode to atomic

        # remove any possible selection before selecting atoms
        cmd.deselect()

        self.ErrorCode += 1
        
    #=======================================================================
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self, rv):
        
        try:
            cmd.delete(self.LigDisplay)
            cmd.delete(self.AtomDisplay)         

            cmd.deselect()

            if rv != 1:
                General.unmask_Objects(self.exc)
                cmd.set('mouse_selection_mode', self.selection_mode)
                #cmd.config_mouse('three_button_editing', 1)
      
            if rv > 0:
                self.FlexAID.WizardError = True
        
            self.FlexAID.ActiveWizard = None

            cmd.set_wizard()
            cmd.set_view(self.View)
            cmd.refresh()

            self.top.SATRunning(False)

            self.queue.put(lambda: self.top.AnchorRunning(False))
            self.queue.put(lambda: self.FlexAID.root.deiconify())
            self.queue.put(lambda: self.FlexAID.root.update())

        except:
            pass
    
    #=======================================================================
    ''' Displays the ligand to be modified '''
    #=======================================================================    
    def DisplayLigand(self):
        
        try:
            cmd.load(self.pdbPath, self.LigDisplay, state=self.State)

            # Display the atoms spheres
            cmd.show('spheres', self.LigDisplay)
            cmd.alter(self.LigDisplay,'vdw=0.25')
            cmd.rebuild(self.LigDisplay)
        
            util.cbag(self.LigDisplay)
            cmd.translate(self.Translation,self.LigDisplay)
            cmd.zoom(self.LigDisplay)

        except:
            return 1

        return 0

    #=======================================================================
    ''' Reset all the variables '''
    #=======================================================================
    def reset(self):

        # Put back default types
        for atom in self.top.Vars.dictAtomTypes.keys():
            self.top.Vars.dictAtomTypes[atom][1] = self.top.Vars.dictAtomTypes[atom][0]

        self.ColorByType()
       
    #=======================================================================
    ''' Button Done selected '''
    #=======================================================================
    def btn_Done(self):
        
        Changed = 0
        Total = 0
        for atom in self.top.Vars.dictAtomTypes.keys():
            if self.top.Vars.dictAtomTypes[atom][1] != self.top.Vars.dictAtomTypes[atom][0]:
                Changed += 1
            Total += 1

        if Changed:
            if Changed == Total:
                self.top.SATStatus.set(' ALL atom types changed')
            else:
                self.top.SATStatus.set(' ' + str(Changed) + ' atom type(s) changed')

        self.Quit_Wizard(0)

    #=======================================================================
    ''' set_AtomType  '''
    #=======================================================================   
    def set_AtomType(self,type):

        if self.SelAtomName != '...':
            self.top.Vars.dictAtomTypes[self.ID][1] = type
                
            self.SelAtomName = '...'
            self.SelAtomType = '...'

        self.ColorByType()

        cmd.delete(self.AtomDisplay)
        cmd.refresh_wizard()
                       
    #=======================================================================
    ''' ColorByType: Change all atoms color based on their type '''
    #=======================================================================
    def ColorByType(self):
        
        try: 
            for atom in self.top.Vars.dictAtomTypes.keys():
                cmd.color(self.type_color[self.top.Vars.dictAtomTypes[atom][1]], self.LigDisplay + ' & ID ' + atom)
        except:
            return 1

        return 0

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
    ''' Highlight atom upon clicking '''
    #=======================================================================    
    def highlight_Atom(self, atom):

        try:
            cmd.pseudoatom(self.AtomDisplay, pos=atom[5:], vdw=0.30, color='white')
            cmd.hide('nonbonded', self.AtomDisplay)
            cmd.show('spheres', self.AtomDisplay)
            cmd.mask(self.AtomDisplay)
        except:
            self.FlexAID.DisplayMessage("Failed to highlight atom upon selecting atom", 1)
            return

    #=======================================================================
    ''' Pick an atom, then display the information related to it'''
    #=======================================================================
    def do_select(self, name):

        self.atom = self.get_Atom(name)
            
        if len(self.atom) > 0:
            if int(self.atom[2]) == self.FlexAID.IOFile.ResSeq.get() and \
                    self.atom[1] == 'LIG':

                cmd.delete(self.AtomDisplay)                    

                # highlight atom
                self.highlight_Atom(self.atom)

                self.ID = str(General.get_ID(self.atom[0], self.LigDisplay))
                self.SelAtomName = self.ID + ' (' + self.atom[4] + ')'
                self.SelAtomType = self.top.Vars.dictAtomTypes[self.ID][1]

                cmd.refresh_wizard()

            else:
                self.FlexAID.DisplayMessage("You can only select atoms in the object " + self.LigDisplay)

        else:
            self.FlexAID.DisplayMessage("No atom could be selected in " + self.LigDisplay)

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
                self.FlexAID.DisplayMessage("You must click in the object " + self.LigDisplay, 1)
                return info
                
            cmd.deselect()

        except:
            self.FlexAID.DisplayMessage("Error while retrieving atom info", 1)
            self.Quit_Wizard(self.ErrorCode)
        
        return info

    #=======================================================================
    ''' Panel displayed in the right side of the Pymol interface '''
    #=======================================================================
    def get_panel(self):

        return [
         [ 1, '* Set Atom Types *',''],
         [ 1, '\\999Atom Selected: ' + self.SelAtomName,''],
         [ 1, '\\999Type: ' + self.SelAtomType,''],
         [ 3, 'Change Atom Type','type'],
         [ 2, 'Show atom names','cmd.get_wizard().show_AtomsName()'],
         [ 2, 'Show atom IDs','cmd.get_wizard().show_AtomsNumber()'],
         [ 2, 'Hide labels','cmd.get_wizard().hide_Labels()'],
         [ 2, 'Reset','cmd.get_wizard().reset()'],  
         [ 2, 'Done','cmd.get_wizard().btn_Done()'],         
         ]

    #=======================================================================
    ''' Display a message in the interface '''
    #=======================================================================
    def get_prompt(self):

        return ['Click on the atom you wish to change the type']
