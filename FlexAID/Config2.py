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

import os
import functools
import General
import Vars
import Tabs

if __debug__:
    from pymol import cmd

    import General_cmd
    import Constraint
    import AtomTypes
    import FlexBonds


class Config2Vars(Vars.Vars):

    SATStatus = StringVar()
    FlexStatus = StringVar()
    UseReference = IntVar()
    CovDist = DoubleVar() 
    CovConsSelection = StringVar()
    IntTranslation = IntVar()
    IntRotation = IntVar()
    Anchor = IntVar()
    
    def __init__(self):
    
        self.dictCovConstraints = dict()
        self.dictAtomTypes = dict()
        self.dictFlexBonds = dict()
        self.dictNeighbours = dict()
    

class Config2(Tabs.Tab):

    DefaultCovDist = 1.5

    HIGHLIGHT = 'HIGHLIGHTCONS'
    COVCONS = 'COVCONS'

    
    def Def_Vars(self):

        self.SATStatus = self.Vars.SATStatus
        self.FlexStatus = self.Vars.FlexStatus
        self.UseReference = self.Vars.UseReference
        self.CovDist = self.Vars.CovDist
        self.CovConsSelection = self.Vars.CovConsSelection
        self.IntTranslation = self.Vars.IntTranslation
        self.IntRotation = self.Vars.IntRotation
        self.Anchor = self.Vars.Anchor

    def Init_Vars(self):
        
        self.dictCovConstraints.clear()

        self.UseReference.set(0)

        self.CovDist.set(0.0)
        self.CovConsSelection.set('')

        self.IntTranslation.set(1)
        self.IntRotation.set(1)

        self.FlexStatus.set('')
        self.SATStatus.set('')
        

    def Update_Vars(self):
    
        self.dictCovConstraints = self.Vars.dictCovConstraints
        self.dictAtomTypes = self.Vars.dictAtomTypes
        self.dictFlexBonds = self.Vars.dictFlexBonds
        self.dictNeighbours = self.Vars.dictNeighbours
        
    
    ''' ==================================================================================
    FUNCTION Btn_AddRemove_FlexBonds: Enables wizard to set flexible bond
    =================================================================================  '''    
    def Btn_AddRemove_FlexBonds(self):
        
        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:

            self.FlexStatus.set('')
            self.FlexBondsRunning(True)

            self.top.ActiveWizard = FlexBonds.flexbond(self)
            
            cmd.set_wizard(self.top.ActiveWizard)
            self.top.ActiveWizard.Start()
                
        else:
            self.top.DisplayMessage("A wizard is currently active", 2)

    ''' ==================================================================================
    FUNCTION FlexBondsRunning: Disables/enables controls related to Flexible lig. wizard
    =================================================================================  '''    
    def FlexBondsRunning(self,boolRun):

        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()
            
            #print "WizardResult=",self.top.WizardResult
            if self.top.WizardError or self.top.WizardResult == 0:
                self.FlexStatus.set('')
            
            elif self.top.WizardResult > 0:
                #print "IM HERE FUCK!"
                Status = ' (' + str(self.top.WizardResult) + ') flexible bond(s) set'
                self.FlexStatus.set(Status)

                    
    ''' ==================================================================================
    FUNCTION SAT_Toggle: Enables wizard to set atom types
    =================================================================================  '''    
    def Btn_Edit_AtomTypes(self):
        
        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:

            self.SATStatus.set('')

            self.SATRunning(True)
            self.top.ActiveWizard = AtomTypes.setType(self, self.top.IOFile.OldTypes.get())

            cmd.set_wizard(self.top.ActiveWizard)
            self.top.ActiveWizard.Start()
            
        else:
            self.top.DisplayMessage("A wizard is currently active", 2)

    ''' ==================================================================================
    FUNCTION SATRunning: Disables/enables controls related to Flexible lig. wizard
    =================================================================================  '''    
    def SATRunning(self,boolRun):

        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()

    ''' ==================================================================================
    FUNCTION Add_Constraint: Adds a new constraint type (starts wizard)
    =================================================================================  '''    
    def Add_Constraint(self, Type):

        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:
            
            self.ConsRunning(None, None, None, True)

            self.top.ActiveWizard = Constraint.constraint(self,Type)
            cmd.set_wizard(self.top.ActiveWizard)

            self.top.ActiveWizard.Start()
            
        else:
            self.top.DisplayMessage("A wizard is currently active", 2)


    ''' ==================================================================================
    FUNCTION Update_Constraints: Updates the list of constraints
    =================================================================================  '''    
    def Update_Constraints(self, dictCons, optCons, selection):

        optCons["menu"].delete(0, END)

        Cons = ''
        for key in sorted(dictCons): # odd sorting !
            optCons["menu"].add_command(label=key, command=lambda temp = key: optCons.setvar(optCons.cget("textvariable"), value = temp))
            Cons = key

        selection.set(Cons)


    ''' ==================================================================================
    FUNCTION ConsRunning: Toggles controls' state related to when Constraints are active
    =================================================================================  '''    
    def ConsRunning(self, dictCons, optCons, selection, boolRun):
        
        if boolRun:
            del self.StateList[:]

            General.saveState(self.fBonds,self.StateList)
            General.saveState(self.fSAT,self.StateList)

            General.setState(self.fBonds)
            General.setState(self.fSAT)
            
        else:
            self.Update_Constraints(dictCons, optCons, selection)

            General.backState(self.fBonds,self.StateList)
            General.backState(self.fSAT,self.StateList)
            

    ''' ==================================================================================
    FUNCTION Del_Constraint: Deletes a constraint type
    =================================================================================  '''    
    def Del_Constraint(self, dictCons, optCons, selection):

        key = selection.get()

        if key != '':
            
            try:
                cmd.delete(self.HIGHLIGHT)
                cmd.delete(dictCons[key][2])
            except:
                pass
            
            del dictCons[key]
            self.Update_Constraints(dictCons, optCons, selection)


    ''' ==================================================================================
    FUNCTION Clear_Constraints: Clear constraints of type i
    =================================================================================  '''    
    def Clear_Constraints(self, dictCons, optCons, selection):

        dictCons.clear()
        self.Update_Constraints(dictCons,optCons,selection)

    ''' ==================================================================================
    FUNCTION Frame: Generate the Configuration Options frame in the the middle 
                    frame section
    =================================================================================  '''    
    def Frame(self):
        
        self.fConfig2 = Frame(self.top.fMiddle)

        fC2Left = Frame(self.fConfig2)
        fC2Left.pack(side=LEFT, fill=Y)
        
        fC2Right = Frame(self.fConfig2)
        fC2Right.pack(side=LEFT, fill=BOTH, expand=True)

        #************************************************#
        #*         Ligand Translation/Rotation          *#
        #************************************************#

        self.fTransRot = Frame(fC2Left)
        self.fTransRot.pack(fill=X, side=TOP, padx=5, pady=5)
        fTransRotLine1 = Frame(self.fTransRot)
        fTransRotLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fTransRotLine2 = Frame(self.fTransRot)
        fTransRotLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fTransRotLine3 = Frame(self.fTransRot)
        fTransRotLine3.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fTransRotLine1, text='Ligand translation and rotation', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Checkbutton(fTransRotLine2, text='Translational degrees of freedom', variable=self.IntTranslation, font=self.top.font_Text).pack(fill=X, side=LEFT)
        Checkbutton(fTransRotLine3, text='Rotational degrees of freedom', variable=self.IntRotation, font=self.top.font_Text).pack(fill=X, side=LEFT)        

        #************************************************#
        #*             Ligand Flexible Bonds            *#
        #************************************************#

        self.fBonds = Frame(fC2Left)
        self.fBonds.pack(fill=X, side=TOP, padx=5, pady=5)
        fBondsLine1 = Frame(self.fBonds)
        fBondsLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fBondsLine2 = Frame(self.fBonds)
        fBondsLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fBondsLine3 = Frame(self.fBonds)
        fBondsLine3.pack(side=TOP, fill=X, padx=5, pady=2)
        
        Label(fBondsLine1, text='Ligand flexibility', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Button(fBondsLine2, text='Add/Remove flexible bonds', font=self.top.font_Text,command=self.Btn_AddRemove_FlexBonds, width=40).pack(side=LEFT, anchor=W)
        Entry(fBondsLine3, text='', state='disabled', textvariable=self.FlexStatus, font=self.top.font_Text, 
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT)

        #************************************************#
        #*             Set atom types                   *#
        #************************************************#

        self.fSAT = Frame(fC2Left)
        self.fSAT.pack(fill=X, side=TOP, padx=5, pady=5)
        fSATLine1 = Frame(self.fSAT)
        fSATLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fSATLine2 = Frame(self.fSAT)
        fSATLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fSATLine3 = Frame(self.fSAT)
        fSATLine3.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fSATLine1, text='Ligand atom typing', font=self.top.font_Title, state='disabled').pack(side=TOP,anchor=W)
        Button(fSATLine2, text='Edit atom types', font=self.top.font_Text,command=self.Btn_Edit_AtomTypes, width=40, state='disabled').pack(side=LEFT, anchor=W)
        Entry(fSATLine3, text='', state='disabled', textvariable=self.SATStatus, font=self.top.font_Text,
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT)

        #************************************************#
        #*                Reference pose (RMSD)         *#
        #************************************************#
        fRMSD = Frame(fC2Right)
        fRMSD.pack(fill=X, side=TOP, padx=5, pady=5)
        fRMSDLine1 = Frame(fRMSD)
        fRMSDLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fRMSDLine2 = Frame(fRMSD)
        fRMSDLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        
        Label(fRMSDLine1, text='RMSD structure', font=self.top.font_Title).pack(side=LEFT,anchor=W)        
        Checkbutton(fRMSDLine2, text='Ligand pose as reference', variable=self.UseReference, font=self.top.font_Text).pack(side=LEFT)

        #************************************************#
        #*               Add Constraints                *#
        #************************************************#
        #*                  COVALENT                    *#
        #************************************************#        
        
        fCovConstraint = Frame(fC2Right)#,bd=1,bg="red")
        fCovConstraint.pack(fill=X, side=TOP, padx=5, pady=5)
        fCovConstraintLine1 = Frame(fCovConstraint)
        fCovConstraintLine1.pack(fill=X, side=TOP, padx=5, pady=2)
        fCovConstraintLine2 = Frame(fCovConstraint)
        fCovConstraintLine2.pack(fill=X, side=TOP, padx=5, pady=2)
        fCovConstraintLine3 = Frame(fCovConstraint)
        fCovConstraintLine3.pack(fill=X, side=TOP, padx=5, pady=2)
        fCovConstraintLine4 = Frame(fCovConstraint)
        fCovConstraintLine4.pack(fill=X, side=TOP, padx=5, pady=2)
 
        #Label(fCovConstraintLine1, text = 'Covalent constraints', font=self.top.font_Title, justify=LEFT).pack(side=LEFT, anchor=W)
        Label(fCovConstraintLine1, text = 'Constraints', font=self.top.font_Title, justify=LEFT).pack(side=LEFT, anchor=W)

        Label(fCovConstraintLine2, text = 'Name:', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT, anchor=E)
        optionTuple = '',
        self.optCovCons = apply(OptionMenu, (fCovConstraintLine2, self.CovConsSelection) + optionTuple)
        self.optCovCons.config(bg=self.top.Color_White, font=self.top.font_Text)
        self.optCovCons.pack(side=LEFT, fill=X, expand=True)

        Label(fCovConstraintLine3, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT, anchor=E)
        self.btnAddCovConstraint = Button(fCovConstraintLine3, text = 'Add/Edit',font=self.top.font_Text,command=lambda: self.Add_Constraint('Covalent'))
        self.btnAddCovConstraint.pack(side=LEFT)
        self.btnDelCovConstraint = Button(fCovConstraintLine3, text = 'Delete', font=self.top.font_Text,
                                          command=lambda: self.Del_Constraint(self.dictCovConstraints,self.optCovCons,self.CovConsSelection))
        self.btnDelCovConstraint.pack(side=LEFT)
        self.btnClearCovConstraint = Button(fCovConstraintLine3, text = 'Clear', font=self.top.font_Text,
                                           command=lambda: self.Clear_Constraints(self.dictCovConstraints,self.optCovCons,self.CovConsSelection))
        self.btnClearCovConstraint.pack(side=LEFT)

        Label(fCovConstraintLine4, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT, anchor=E)
        Label(fCovConstraintLine4, text = 'Interaction distance (A):', font=self.top.font_Text).pack(side=LEFT, anchor=S)

        self.sclCovDist = Scale(fCovConstraintLine4, from_ = 0.0, to = 10.0, orient=HORIZONTAL, length=120, resolution=0.05, variable=self.CovDist)
        self.sclCovDist.pack(side=LEFT)
        self.sclCovDist.config(state='disabled')

        return self.fConfig2
        
    ''' ==================================================================================
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''        
    def Load_Message(self):
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < Ligand Cfg > Menu.', 2)
        self.DisplayMessage('  INFO:   Set different options for the Ligand.', 2)
        self.DisplayMessage('          1) Include flexibility in the ligand, if necessary.', 2)
        self.DisplayMessage('          2) Alter definition of atom types, if necessary.', 2)
        self.DisplayMessage('          3) Add contraints of optimization, if necessary.', 2)

