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

if __debug__:
    from pymol import cmd

    import General_cmd
    import Constraint
    import AtomTypes
    import FlexBonds

class Config2:

    def __init__(self,top,PyMOL):
        
        #print "New instance of Config2"
        
        self.PyMOL = PyMOL
        
        self.top = top
        self.Tab = self.top.Btn_Config2
        self.FrameName = 'Config2'

        self.Def_Vars()
        self.Init_Vars()

        self.Frame()
        self.Trace()

    def Def_Vars(self):

        self.dictCovConstraints = dict()
        self.dictIntConstraints = dict()

        self.SetAtmTypeOpt = StringVar()
        self.FlexBondsPar = StringVar()

        self.SATStatus = StringVar()
        self.FlexStatus = StringVar()
        self.UseReference = IntVar()

        self.CovDist = DoubleVar() 
        self.CovConsSelection = StringVar()
        self.IntFactor = DoubleVar()
        self.IntConsSelection = StringVar()
        self.IntForced = IntVar()

        self.IntTranslation = IntVar()
        self.IntRotation = IntVar()
        
        self.dictAtomTypes = dict()
        self.dictFlexBonds = dict()
        self.dictNeighbours = dict()

        self.Translation = list()
        self.Validator = list()
        self.StateList = list()

    def Init_Vars(self):
        
        self.dictCovConstraints.clear()
        self.dictIntConstraints.clear()

        self.SetAtmTypeOpt.set('DEFAULT')
        self.FlexBondsPar.set('RIGID')
        self.UseReference.set(0)

        self.CovDist.set(0.0)
        self.CovConsSelection.set('')
        self.IntFactor.set(0.0)
        self.IntConsSelection.set('')
        self.IntForced.set(0)

        self.IntTranslation.set(1)
        self.IntRotation.set(1)

        self.FlexStatus.set('')
        self.SATStatus.set('')

        self.Translation = [1000.000, 1000.000, 1000.000]
        self.StateList = []

        self.DefaultCovDist = 1.5
        self.DefaultIntFactor = 5.0

        self.HIGHLIGHT = 'HIGHLIGHTCONS'
        self.COVCONS = 'COVCONS'
        self.INTCONS = 'INTCONS'

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        self.FlexTrace = self.FlexBondsPar.trace('w',self.Flex_Toggle)
        self.SATTrace = self.SetAtmTypeOpt.trace('w',self.SAT_Toggle)

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''    
    def Del_Trace(self):

        self.FlexBondsPar.trace_vdelete('w',self.FlexTrace)
        self.SetAtmTypeOpt.trace_vdelete('w',self.SATTrace)

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
        
        self.fConfig2.pack_forget()
        #self.fConfig2.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    def Validator_Fail(self):

        return

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
        
        self.LoadMessage()

        self.fConfig2.pack(fill=BOTH, expand=True)

    ''' ==================================================================================
    FUNCTION Flex_Toggle: Enables wizard to set flexible bond
    =================================================================================  '''    
    def Flex_Toggle(self,*args):
        
        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:

            if self.FlexBondsPar.get() == 'CUSTOM':
                self.FlexStatus.set('')
                self.FlexBondsRunning(True)

                self.top.ActiveWizard = FlexBonds.flexbond(self)
                
                cmd.set_wizard(self.top.ActiveWizard)
                self.top.ActiveWizard.Start()
                
            elif self.FlexBondsPar.get() == 'RIGID':
                # Clear the data in the text box
                self.FlexStatus.set('')
                
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

            if self.top.WizardError or self.top.WizardResult == 0:
                self.FlexBondsPar.set('RIGID')
                self.FlexStatus.set('')
            elif self.top.WizardResult > 0:
                self.FlexStatus.set(str(self.top.WizardResult) + ' flexible bond(s) set')

    ''' ==================================================================================
    FUNCTION SAT_Toggle: Enables wizard to set atom types
    =================================================================================  '''    
    def SAT_Toggle(self,*args):
        
        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:

            if self.SetAtmTypeOpt.get() == 'CUSTOM':
                self.SATStatus.set('')

                self.SATRunning(True)
                self.top.ActiveWizard = AtomTypes.setType(self, self.top.IOFile.OldTypes.get())

                cmd.set_wizard(self.top.ActiveWizard)
                self.top.ActiveWizard.Start()

            elif self.SetAtmTypeOpt.get() == 'DEFAULT':
                # Clear the data in the text box
                self.SATStatus.set('')
            
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
            
            print "CONSRUNNING"
            self.ConsRunning(None, None, None, True)

            print "CONSTRUCTOR"
            self.top.ActiveWizard = Constraint.constraint(self,Type)
            cmd.set_wizard(self.top.ActiveWizard)

            print "START"
            self.top.ActiveWizard.Start()
            
        else:
            self.top.DisplayMessage("A wizard is currently active", 2)


    ''' ==================================================================================
    FUNCTION Update_Constraints: Updates the list of constraints
    =================================================================================  '''    
    def Update_Constraints(self, dictCons, optCons, selection):

        optCons["menu"].delete(0, END)

        Cons = ''
        for key in sorted(dictCons):
            optCons["menu"].add_command(label=key, command=lambda temp = key: optCons.setvar(optCons.cget("textvariable"), value = temp))
            Cons = key

        selection.set(Cons)


    ''' ==================================================================================
    FUNCTION ConsRunning: Toggles controls' state related to when Constraints are active
    =================================================================================  '''    
    def ConsRunning(self, dictCons, optCons, selection, boolRun):
        
        if boolRun:
            self.StateList = []

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
    FUNCTION Disable_Frame: Disables all controls on main frame
    =================================================================================  '''    
    def Disable_Frame(self):

        self.StateList = []
        General.saveState(self.fConfig2, self.StateList)
        General.setState(self.fConfig2)

    ''' ==================================================================================
    FUNCTION Enable_Frame: Enables all controls on main frame
    =================================================================================  '''    
    def Enable_Frame(self):

        General.backState(self.fConfig2, self.StateList)

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
        fTransRotLine1.pack(side=TOP, fill=X)
        fTransRotLine2 = Frame(self.fTransRot)
        fTransRotLine2.pack(side=TOP, fill=X)
        fTransRotLine3 = Frame(self.fTransRot)
        fTransRotLine3.pack(side=TOP, fill=X)

        Label(fTransRotLine1, text='Ligand translation and rotation', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Checkbutton(fTransRotLine2, text='Translational degrees of freedom', width=25, variable=self.IntTranslation, font=self.top.font_Text).pack(fill=X, anchor=W)
        Checkbutton(fTransRotLine3, text='Rotational degrees of freedom', width=25, variable=self.IntRotation, font=self.top.font_Text).pack(fill=X, anchor=W)        

        #************************************************#
        #*             Ligand Flexible Bonds            *#
        #************************************************#

        self.fBonds = Frame(fC2Left)
        self.fBonds.pack(fill=X, side=TOP, padx=5, pady=5)
        fBondsLine1 = Frame(self.fBonds)
        fBondsLine1.pack(side=TOP, fill=X)
        fBondsLine2 = Frame(self.fBonds)
        fBondsLine2.pack(side=TOP, fill=X)
        fBondsLine3 = Frame(self.fBonds)
        fBondsLine3.pack(side=TOP, fill=X)
        
        self.fBonds = Frame(fC2Left)
        self.fBonds.pack(fill=X, side=TOP, padx=5, pady=5)
        fBondsLine1 = Frame(self.fBonds)
        fBondsLine1.pack(side=TOP, fill=X)
        fBondsLine2 = Frame(self.fBonds)
        fBondsLine2.pack(side=TOP, fill=X)
        fBondsLine3 = Frame(self.fBonds)
        fBondsLine3.pack(side=TOP, fill=X)

        Label(fBondsLine1, text='Ligand flexibility', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Radiobutton(fBondsLine2, text='Rigid bonds', width=15, value='RIGID', variable=self.FlexBondsPar, font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Radiobutton(fBondsLine3, text='Customize:', width=15, value='CUSTOM', variable=self.FlexBondsPar, font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Entry(fBondsLine3, text='', state='disabled', textvariable=self.FlexStatus, font=self.top.font_Text, disabledforeground=self.top.Color_Black, disabledbackground=self.top.Color_White).pack(side=LEFT)

        #************************************************#
        #*             Set atom types                   *#
        #************************************************#

        self.fSAT = Frame(fC2Left)
        self.fSAT.pack(fill=X, side=TOP, padx=5, pady=5)
        fSATLine1 = Frame(self.fSAT)
        fSATLine1.pack(side=TOP, fill=X)
        fSATLine2 = Frame(self.fSAT)
        fSATLine2.pack(side=TOP, fill=X)
        fSATLine3 = Frame(self.fSAT)
        fSATLine3.pack(side=TOP, fill=X)

        Label(fSATLine1, text='Ligand atom types', font=self.top.font_Title, state='disabled').pack(side=TOP,anchor=W)
        Radiobutton(fSATLine2, text='Default types', width=15, value='DEFAULT', variable=self.SetAtmTypeOpt, font=self.top.font_Text, state='disabled', justify=LEFT).pack(side=LEFT, anchor=W)
        Radiobutton(fSATLine3, text='Customize:', width=15, value='CUSTOM', variable=self.SetAtmTypeOpt, font=self.top.font_Text, state='disabled').pack(side=LEFT, anchor=W)
        Entry(fSATLine3, text='', state='disabled',textvariable=self.SATStatus, font=self.top.font_Text, disabledforeground=self.top.Color_Black, disabledbackground=self.top.Color_White).pack(side=LEFT)

        #************************************************#
        #*                Reference pose (RMSD)         *#
        #************************************************#
        fRMSD = Frame(fC2Left)
        fRMSD.pack(fill=X, side=TOP, padx=5, pady=5)
        fRMSDLine1 = Frame(fRMSD)
        fRMSDLine1.pack(side=TOP, fill=X)
        fRMSDLine2 = Frame(fRMSD)
        fRMSDLine2.pack(side=TOP, fill=X)
        
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
        fCovConstraintLine1.pack(fill=X, side=TOP)
        fCovConstraintLine2 = Frame(fCovConstraint)
        fCovConstraintLine2.pack(fill=X, side=TOP)
        fCovConstraintLine3 = Frame(fCovConstraint)
        fCovConstraintLine3.pack(fill=X, side=TOP)
        fCovConstraintLine4 = Frame(fCovConstraint)
        fCovConstraintLine4.pack(fill=X, side=TOP)
 
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

        """
        #************************************************#
        #*                  INTERACTION                 *#
        #************************************************#        
        
        fIntConstraint = Frame(fC2Right)#,bd=1,bg="orange")
        fIntConstraint.pack(fill=X, side=TOP, padx=5, pady=5)
        fIntConstraintLine1 = Frame(fIntConstraint)
        fIntConstraintLine1.pack(fill=X, side=TOP)
        fIntConstraintLine2 = Frame(fIntConstraint)
        fIntConstraintLine2.pack(fill=X, side=TOP)
        fIntConstraintLine3 = Frame(fIntConstraint)
        fIntConstraintLine3.pack(fill=X, side=TOP)
        fIntConstraintLine4 = Frame(fIntConstraint)
        fIntConstraintLine4.pack(fill=X, side=TOP)
        fIntConstraintLine5 = Frame(fIntConstraint)
        fIntConstraintLine5.pack(fill=X, side=TOP)
 
        Label(fIntConstraintLine1, text = 'Interaction Constraints', font=self.top.font_Title, justify=LEFT).pack(side=LEFT, anchor=W)

        Label(fIntConstraintLine2, text = 'Name:', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT, anchor=E)
        optionTuple = '',
        self.optIntCons = apply(OptionMenu, (fIntConstraintLine2, self.IntConsSelection) + optionTuple)
        self.optIntCons.config(bg=self.top.Color_White, font=self.top.font_Text)
        self.optIntCons.pack(side=LEFT, fill=X, expand=True)

        Label(fIntConstraintLine3, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT)
        self.btnAddIntConstraint = Button(fIntConstraintLine3, text = 'Add/Edit',font=self.top.font_Text,command=lambda: self.Add_Constraint('Interaction'))
        self.btnAddIntConstraint.pack(side=LEFT)
        self.btnDelIntConstraint = Button(fIntConstraintLine3, text = 'Delete', font=self.top.font_Text,
                                          command=lambda: self.Del_Constraint(self.dictIntConstraints,self.optIntCons,self.IntConsSelection))
        self.btnDelIntConstraint.pack(side=LEFT)
        self.btnClearIntConstraint = Button(fIntConstraintLine3, text = 'Clear',font=self.top.font_Text,
                                           command=lambda: self.Clear_Constraints(self.dictIntConstraints,self.optIntCons,self.IntConsSelection))
        self.btnClearIntConstraint.pack(side=LEFT)

        Label(fIntConstraintLine4, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT)
        Label(fIntConstraintLine4, text = 'Multiplier:', width=10, font=self.top.font_Text).pack(side=LEFT, anchor=S)
        self.sclIntFactor = Scale(fIntConstraintLine4, from_ = -10.0, to = 10.0, orient=HORIZONTAL, length=120, resolution=1.0, variable=self.IntFactor)
        self.sclIntFactor.pack(side=LEFT)
        self.sclIntFactor.config(state='disabled')

        #Label(fIntConstraintLine5, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT)
        #Label(fIntConstraintLine5, text = '', width=10, font=self.top.font_Text, justify=RIGHT).pack(side=LEFT)
        #self.chkIntForced = Checkbutton(fIntConstraintLine5, text = 'Force', font=self.top.font_Text, variable=self.IntForced)
        #self.chkIntForced.pack(side=LEFT)
        """

    ''' ==================================================================================
    FUNCTION MenuLoadMessage: Display the message based on the menu selected
    ================================================================================== '''        
    def LoadMessage(self):
        
        self.top.DisplayMessage('', 0)
        self.top.DisplayMessage('  FlexAID < Ligand Cfg > Menu.', 2)
        self.top.DisplayMessage('  INFO:   Set different options for the Ligand.', 2)
        self.top.DisplayMessage('          1) Include flexibility in the ligand, if necessary.', 2)
        self.top.DisplayMessage('          2) Alter definition of atom types, if necessary.', 2)
        self.top.DisplayMessage('          3) Add contraints of optimization, if necessary.', 2)

