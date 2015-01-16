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

from Tkinter import *

import os
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

    IntTranslation = IntVar()
    IntRotation = IntVar()
    FlexStatus = StringVar()
    ConsStatus = StringVar()
    SATStatus = StringVar()
    UseReference = IntVar()
    ConsDist = DoubleVar() 
    ActiveCons = StringVar()
    
    def __init__(self):
    
        self.dictConstraints = dict()
    
class Config2(Tabs.Tab):
    
    MIN_DIST_CONSTRAINT = 1.00
    MAX_DIST_CONSTRAINT = 7.00
    RESOLUTION_CONSTRAINT = 0.05
    
    def Def_Vars(self):

        self.SATStatus = self.Vars.SATStatus
        self.FlexStatus = self.Vars.FlexStatus
        self.UseReference = self.Vars.UseReference
        self.ConsDist = self.Vars.ConsDist
        self.ConsStatus = self.Vars.ConsStatus
        self.IntTranslation = self.Vars.IntTranslation
        self.IntRotation = self.Vars.IntRotation
        self.ActiveCons = self.Vars.ActiveCons
        
    def Init_Vars(self):
        
        self.ActiveCons.set('')      
        self.UseReference.set(0)
        self.ConsDist.set(0.25)
        self.IntTranslation.set(1)
        self.IntRotation.set(1)
        self.SATStatus.set('')
        self.ConsStatus.set('No constraint(s) set')
        
        self.Vars.dictConstraints.clear()
        self.ResetFlexBonds()
        self.ResetAtomTypes()
    
    def Trace(self):
    
        try:
            self.ActiveConsTrace = self.ActiveCons.trace('w', self.ActiveCons_Toggle)
            self.ConsDistTrace = self.ConsDist.trace('w', self.ConsDist_Toggle)
        except:
            pass
    
    def Del_Trace(self):
    
        try:
            self.ActiveCons.trace_vdelete('w',self.ActiveConsTrace)
            self.ConsDist.trace_vdelete('w',self.ConsDistTrace)
        except:
            pass
    
    ''' ==================================================================================
    FUNCTION Btn_AddRemove_FlexBonds: Enables wizard to set flexible bond
    =================================================================================  '''    
    def Btn_AddRemove_FlexBonds(self):
        
        if not self.PyMOL:
            return

        self.FlexStatus.set('')
        self.FlexBondsRunning(True)
        
        self.top.ActiveWizard = FlexBonds.flexbond(self, self.queue, self.top.IOFile.ResSeq.get(),self.top.IOFile.ProcessedLigandPath.get())
        
        cmd.set_wizard(self.top.ActiveWizard)
        cmd.refresh()

        self.top.ActiveWizard.Start()
            
    ''' ==================================================================================
    FUNCTION ResetAtomTypes: Resets the atom types of the ligand to defaults
    =================================================================================  '''    
    def ResetAtomTypes(self):
    
        return

    ''' ==================================================================================
    FUNCTION ResetFlexBonds: Resets the selected flexible bonds of the ligand
    =================================================================================  '''    
    def ResetFlexBonds(self):

        for index in self.top.IOFile.Vars.dictFlexBonds.keys():
            self.top.IOFile.Vars.dictFlexBonds[index][0] = 0
        self.Update_FlexStatus()

    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):
                
        #Show the list of selection/objects in case the user already worked in PyMOL
        self.Update_FlexStatus()        

    ''' ==================================================================================
    FUNCTION FlexBondsRunning: Disables/enables controls related to Flexible lig. wizard
    =================================================================================  '''    
    def FlexBondsRunning(self,boolRun):

        if boolRun:
            self.Start_Update()
            self.Disable_Frame()
            
        else:
            self.End_Update()
            self.Enable_Frame()
            
            if self.top.WizardError or self.top.WizardResult == 0:
                Status = 'No flexible bond(s) set'
            
            elif self.top.WizardResult > 0:
                Status = ' (' + str(self.top.WizardResult) + ') flexible bond(s) set'

            self.FlexStatus.set(Status)
    
    ''' ==================================================================================
    FUNCTION SAT_Toggle: Enables wizard to set atom types
    =================================================================================  '''    
    def Btn_Edit_AtomTypes(self):
        
        if not self.PyMOL:
            return
            
        self.SATStatus.set('')

        self.SATRunning(True)
        self.top.ActiveWizard = AtomTypes.setType(self, self.queue, self.top.IOFile.OldTypes.get())

        cmd.set_wizard(self.top.ActiveWizard)
        cmd.refresh()

        self.top.ActiveWizard.Start()
            
    ''' ==================================================================================
    FUNCTION SATRunning: Disables/enables controls related to Flexible lig. wizard
    =================================================================================  '''    
    def SATRunning(self,boolRun):

        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()
    
    ''' ==================================================================================
    FUNCTION AddEditDel_Constraint: Opens the Constraint Wizard.
    =================================================================================  '''    
    def AddEditDel_Constraint(self):

        if not self.PyMOL:
            return
        
        self.ConsRunning(True)
        
        self.top.ActiveWizard = Constraint.constraint(self, self.queue)
        cmd.set_wizard(self.top.ActiveWizard)
        cmd.refresh()

        self.top.ActiveWizard.Start()

    ''' ==================================================================================
    FUNCTION ConsRunning: Toggles controls' state related to when Constraints are active
    =================================================================================  '''    
    def ConsRunning(self, boolRun):
        
        if boolRun:
            self.Start_Update()
            
            if self.ActiveCons.get():
                self.Disable_Frame(self.lblInteraction, self.sclConsDist)
            else:
                self.Disable_Frame(self.lblInteraction)
        else:
            self.Enable_Frame()

            if self.top.WizardResult == 0:
                Status = 'No constraint(s) set'
            else:
                Status = ' (' + str(self.top.WizardResult) + ') constraint(s) set'

            self.ConsStatus.set(Status)
            
            self.End_Update()
                        
    ''' ==================================================================================
    FUNCTION Update_FlexStatus: Update self.FlexStatus StringVar
    =================================================================================  '''                
    def Update_FlexStatus(self):

        nflexbonds = 0
        for k in self.top.IOFile.Vars.dictFlexBonds.keys():
            if self.top.IOFile.Vars.dictFlexBonds[k][0]:
                nflexbonds += 1

        if nflexbonds:
            status = '(' + str(nflexbonds) + ') flexible bond(s) set'
            #status = '(' + 'all' + ') flexible bond(s) set'
        else:
            status = 'No flexible bond(s) set'

        self.FlexStatus.set(status)

    ''' ==================================================================================
    FUNCTION Frame: Generate the Configuration Options frame in the the middle 
                    frame section
    =================================================================================  '''    
    def Frame(self):
        
        self.fConfig2 = Frame(self.top.fMiddle)

        fC2Left = Frame(self.fConfig2)
        fC2Left.pack(side=LEFT, fill=BOTH, expand=True)
        
        fC2Right = Frame(self.fConfig2)
        fC2Right.pack(side=LEFT, fill=BOTH, expand=True)

        #************************************************#
        #*         Ligand Translation/Rotation          *#
        #************************************************#

        fDOF = Frame(fC2Left, border=1, relief=RAISED)
        fDOF.pack(fill=BOTH, expand=True, padx=5, pady=5)

        Label(fDOF, text='Degrees of freedom', font=self.top.font_Title_H).pack(side=TOP, fill=X, expand=True, pady=2)

        self.fTransRot = Frame(fDOF)
        self.fTransRot.pack(fill=X, side=TOP, padx=5, pady=5)
        fTransRotLine1 = Frame(self.fTransRot)
        fTransRotLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fTransRotLine2 = Frame(self.fTransRot)
        fTransRotLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fTransRotLine3 = Frame(self.fTransRot)
        fTransRotLine3.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fTransRotLine1, text='Ligand translation and rotation', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Checkbutton(fTransRotLine2, text=' Include translational', variable=self.IntTranslation, font=self.top.font_Text).pack(fill=X, side=LEFT)
        Checkbutton(fTransRotLine3, text=' Include rotational', variable=self.IntRotation, font=self.top.font_Text).pack(fill=X, side=LEFT)        

        #************************************************#
        #*             Ligand Flexible Bonds            *#
        #************************************************#

        self.fBonds = Frame(fDOF)
        self.fBonds.pack(fill=X, side=TOP, padx=5, pady=5)
        fBondsLine1 = Frame(self.fBonds)
        fBondsLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fBondsLine2 = Frame(self.fBonds)
        fBondsLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fBondsLine3 = Frame(self.fBonds)
        fBondsLine3.pack(side=TOP, fill=X, padx=5, pady=2)
        self.Update_FlexStatus()
        Label(fBondsLine1, text='Ligand flexibility', font=self.top.font_Title).pack(side=TOP, anchor=W)
        Button(fBondsLine2, text='Add/Delete flexible bonds', font=self.top.font_Text,command=self.Btn_AddRemove_FlexBonds, width=40).pack(side=LEFT, anchor=NW,expand=True, fill=X)
        Entry(fBondsLine3, text='', state='disabled', textvariable=self.FlexStatus, font=self.top.font_Text, 
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT, fill=X, expand=True, anchor=NW)

        #************************************************#
        #*             Set atom types                   *#
        #************************************************#

        self.fSAT = Frame(fC2Left)
        #self.fSAT.pack(fill=X, side=TOP, padx=5, pady=5)
        fSATLine1 = Frame(self.fSAT)
        fSATLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fSATLine2 = Frame(self.fSAT)
        fSATLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fSATLine3 = Frame(self.fSAT)
        fSATLine3.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fSATLine1, text='Ligand atom typing', font=self.top.font_Title, state='disabled').pack(side=TOP,anchor=W)
        Button(fSATLine2, text='Edit atom types', font=self.top.font_Text,command=self.Btn_Edit_AtomTypes, width=40, state='disabled').pack(side=LEFT, anchor=NW,expand=True)
        Entry(fSATLine3, text='', state='disabled', textvariable=self.SATStatus, font=self.top.font_Text,
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT,anchor=NW,expand=True)

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
        Checkbutton(fRMSDLine2, text=' Ligand pose as reference', variable=self.UseReference, font=self.top.font_Text).pack(side=LEFT)

        #************************************************#
        #*               Add Constraints                *#
        #************************************************#
        
        fConstraint = Frame(fC2Right)#,bd=1,bg="red")
        fConstraint.pack(fill=X, side=TOP, padx=5, pady=5)
        fConstraintLine1 = Frame(fConstraint)
        fConstraintLine1.pack(fill=X, side=TOP, padx=5, pady=2)
        fConstraintLine2 = Frame(fConstraint)
        fConstraintLine2.pack(fill=X, side=TOP, padx=5, pady=2)
        fConstraintLine3 = Frame(fConstraint)
        fConstraintLine3.pack(fill=X, side=TOP, padx=5, pady=2)
        fConstraintLine4 = Frame(fConstraint)
        fConstraintLine4.pack(fill=X, side=TOP, padx=5, pady=2)
 
        Label(fConstraintLine1, text = 'Constraints', font=self.top.font_Title, justify=LEFT).pack(side=LEFT, anchor=W)

        self.btnAddConstraint = Button(fConstraintLine2, text = 'Add/Edit/Delete constraints',font=self.top.font_Text,command=self.AddEditDel_Constraint)
        self.btnAddConstraint.pack(side=LEFT, fill=X, expand=True, anchor=NE)
        
        Entry(fConstraintLine3, text='', state='disabled', textvariable=self.ConsStatus, font=self.top.font_Text, 
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT, fill=X, expand=True, anchor=NE)

        self.lblInteraction = Label(fConstraintLine4, text = 'Interaction distance (A):', font=self.top.font_Text)
        self.lblInteraction.pack(side=LEFT, anchor=SW)

        self.sclConsDist = Scale(fConstraintLine4, from_ = self.MIN_DIST_CONSTRAINT, to = self.MAX_DIST_CONSTRAINT,
                                 orient=HORIZONTAL, length=120, resolution=self.RESOLUTION_CONSTRAINT, variable=self.ConsDist, showvalue=1)
        self.sclConsDist.pack(side=LEFT, fill=X, expand=True)
        self.sclConsDist.config(state='disabled')
        

        return self.fConfig2
    
    ''' ==================================================================================
    FUNCTION parse_cons: parses the constraint from the dictConstraints
    ================================================================================== '''        
    def parse_cons(self, constraint):

        l = []
        
        l.append(constraint[1:constraint.find(' ')])

        st = constraint.find(' ') + 1
        rnc = constraint[st:]

        l.append(rnc[0:3])
        l.append(rnc[3:len(rnc)-1])
        l.append(rnc[len(rnc)-1:len(rnc)])

        return l
    
    ''' ==================================================================================
    FUNCTION ActiveCons_Toggle: The active constraint is changed
    ================================================================================== '''        
    def ActiveCons_Toggle(self, *args):
        
        if self.top.WizardRunning():
            
            if self.ActiveCons.get():
                self.sclConsDist.config(state='normal')
                self.ConsDist.set(self.Vars.dictConstraints[self.ActiveCons.get()][5])
            else:
                self.sclConsDist.config(state='disabled')
            
            self.top.ActiveWizard.refresh_display()
            
    ''' ==================================================================================
    FUNCTION ConsDist_Toggle: The active constraint's interaction distance is changed
    ================================================================================== '''        
    def ConsDist_Toggle(self, *args):

        if self.top.WizardRunning():
            if self.ActiveCons.get() in self.Vars.dictConstraints.keys():
                self.Vars.dictConstraints[self.ActiveCons.get()][5] = self.ConsDist.get()
            self.top.ActiveWizard.refresh_distance()
            
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

