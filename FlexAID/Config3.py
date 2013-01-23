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

if __debug__:
    from pymol import cmd

class Config3:

    def __init__(self,top,PyMOL):

        #print "New instance of Config3"
        self.PyMOL = PyMOL

        self.top = top
        self.AtomTypes = self.top.IOFile.AtomTypes
        self.Tab = self.top.Btn_Config3
        
        self.FrameName = 'Config3'

        self.Def_Vars()
        self.Init_Vars()

        self.Frame()
        self.Trace()

    def Def_Vars(self):

        self.StateList = list()
        self.CompFct = StringVar()

        self.UseDEE = IntVar()
        self.IncludeHET = IntVar()
        self.ExcludeHOH = IntVar()

        self.DEE_Clash_Threshold = StringVar()
        self.Permeability = StringVar()

        self.SolventType = StringVar()
        self.SolventTypeIndex = IntVar()        
        self.SolventTerm = StringVar()

        self.DeltaAngle = StringVar()
        self.DeltaDihedral = StringVar()
        self.DeltaDihedralFlex = StringVar()

        self.ValidDEE = list()
        self.ValidSolventTerm = list()

        self.ValidDeltaAngle = list()
        self.ValidDeltaDihedral = list()
        self.ValidDeltaDihedralFlex = list()
        self.Validator = list()


    def Init_Vars(self):
        
        self.StateList = []

        self.CompFct.set('VCT')

        self.UseDEE.set(0)
        self.IncludeHET.set(1)
        self.ExcludeHOH.set(1)

        self.DEE_Clash_Threshold.set('0.50')
        self.Permeability.set('0.1')
        
        if self.AtomTypes.get() == 'Sobolev':
            self.SolventType.set('< No type >')
            self.SolventTypeIndex.set(0)
        elif self.AtomTypes.get() == 'Gaudreault':
            self.SolventType.set('< Type-based >')
            self.SolventTypeIndex.set(13)
        elif self.AtomTypes.get() == 'Sybyl':
            self.SolventType.set('< No type >')
            self.SolventTypeIndex.set(0)

        self.SolventTerm.set('-5.0')

        self.DeltaAngle.set('2.5')
        self.DeltaDihedral.set('2.5')
        self.DeltaDihedralFlex.set('10.0')

        self.ValidDEE = [True, 1, 0, None]
        self.ValidPermeability = [True, 1, 0, None]
        self.ValidSolventTerm = [True, 1, 0, None]

        self.ValidDeltaAngle = [True, 1, 0, None]
        self.ValidDeltaDihedral = [True, 1, 0, None]
        self.ValidDeltaDihedralFlex = [True, 1, 0, None]

        self.Validator = [ self.ValidDEE, self.ValidPermeability,
                           self.ValidDeltaAngle, self.ValidDeltaDihedral, 
                           self.ValidDeltaDihedralFlex, self.ValidSolventTerm ]


    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        try:
            self.IncludeHETTrace = self.IncludeHET.trace('w',self.IncludeHET_Toggle)
            self.AtomTypesTrace = self.AtomTypes.trace('w',self.AtomTypes_Toggle)
            self.SolventTypeTrace = self.SolventType.trace('w',self.SolventType_Toggle)
            self.DEETrace = self.UseDEE.trace('w',self.DEE_Toggle)
        except:
            pass
        
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''
    def Del_Trace(self):

        try:
            self.IncludeHET.trace_vdelete('w',self.IncludeHETTrace)
            self.AtomTypes.trace_vdelete('w',self.AtomTypesTrace)
            self.SolventType.trace_vdelete('w',self.SolventTypeTrace)
            self.UseDEE.trace_vdelete('w',self.DEETrace)
        except:
            pass        

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''
    def Kill_Frame(self):
        
        self.fConfig3.pack_forget()
        #self.fConfig3.destroy()

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

        self.fConfig3.pack(fill=BOTH, expand=True)

    ''' ==================================================================================
    FUNCTION Disable_Frame: Disables all controls on main frame
    =================================================================================  '''    
    def Disable_Frame(self):

        self.StateList = []
        General.saveState(self.fConfig3, self.StateList)
        General.setState(self.fConfig3)

    ''' ==================================================================================
    FUNCTION Enable_Frame: Enables all controls on main frame
    =================================================================================  '''    
    def Enable_Frame(self):

        General.backState(self.fConfig3, self.StateList)

    ''' ==================================================================================
    FUNCTION IncludeHET_Toggle: Toggle the controls related to Including HET Groups
    =================================================================================  '''    
    def IncludeHET_Toggle(self, *args):
        
        if self.IncludeHET.get():
            self.chkHOH.config(state='normal')
        else:
            self.chkHOH.config(state='disabled')

    ''' ==================================================================================
    FUNCTION AtomTypes_Toggle: Toggle the controls related to Atom Types
    =================================================================================  '''    
    def AtomTypes_Toggle(self, *args):
        
        if self.AtomTypes.get() != 'Gaudreault':
            self.optSolventType.config(state='normal')
            self.SolventType.set('< No type >')
        else:
            self.optSolventType.config(state='disabled')
            self.SolventType.set('< Type-based >')

        # Need for processing the ligand again if atom typing is changed
        self.top.IOFile.fProcessLigand = False

        # Reset atom type definition
        self.top.Config2.SetAtmTypeOpt.set('DEFAULT')
        self.top.Config2.SATStatus.set('')

    ''' ==================================================================================
    FUNCTION SolventType_Toggle: Toggle the controls related to Solvent Types
    =================================================================================  '''    
    def SolventType_Toggle(self, *args):
        
        if self.SolventType.get() == '< No type >':
            self.entSolventTerm.config(state='normal')
            self.SolventTypeIndex.set(0)

        elif self.SolventType.get() == '< Type-based >':
            self.entSolventTerm.config(state='disabled')
            if self.AtomTypes.get() == 'Gaudreault':
                self.SolventTypeIndex.set(13)
            elif self.AtomTypes.get() == 'Sybyl':
                self.SolventTypeIndex.set(27)
            
    ''' ==================================================================================
    FUNCTION DEE_Toggle: Disables/enables the edit box according to checkstate
    =================================================================================  '''    
    def DEE_Toggle(self, *args):
        
        if self.UseDEE.get() == 1:
            self.entDEE.config(state='normal')
        else:
            self.entDEE.config(state='disabled')

    ''' ==================================================================================
    FUNCTION Frame: Generate the Configuration Options frame in the the middle 
                    frame section
    =================================================================================  '''    
    def Frame(self):
        
        self.fConfig3 = Frame(self.top.fMiddle)

        fC3Left = Frame(self.fConfig3)
        fC3Left.pack(side=LEFT, fill=BOTH, expand=True)
        
        fC3Right = Frame(self.fConfig3)
        fC3Right.pack(side=RIGHT, fill=BOTH, expand=True)

        #==================================================================================
        # Include/Exclude HET Groups
        #==================================================================================
        fHET = Frame(fC3Left)
        fHET.pack(fill=X, side=TOP, padx=5, pady=5)
        fHETLine1 = Frame(fHET)
        fHETLine1.pack(fill=X, padx=5, pady=2)
        fHETLine2 = Frame(fHET)
        fHETLine2.pack(fill=X, padx=5, pady=2)
        fHETLine3 = Frame(fHET)
        fHETLine3.pack(fill=X, padx=5, pady=2)

        Label(fHETLine1, text='HET groups', font=self.top.font_Title).pack(side=LEFT, anchor=W)
        Checkbutton(fHETLine2, text='Include bound molecules', variable=self.IncludeHET, font=self.top.font_Text).pack(side=LEFT)
        
        self.chkHOH = Checkbutton(fHETLine3, text='Exclude water molecules', variable=self.ExcludeHOH, font=self.top.font_Text)
        self.chkHOH.pack(side=LEFT)

        #==================================================================================
        # Permeability
        #==================================================================================
        fPermea = Frame(fC3Left)#, bd=1, relief=SUNKEN)
        fPermea.pack(fill=X, side=TOP, padx=5, pady=5)
        fPermeaLine1 = Frame(fPermea)
        fPermeaLine1.pack(fill=X, padx=5, pady=2)
        fPermeaLine2 = Frame(fPermea)
        fPermeaLine2.pack(fill=X, padx=5, pady=2)

        Label(fPermeaLine1, text='Permeability', font=self.top.font_Title).pack(side=LEFT, anchor=W)
        Label(fPermeaLine2, text='Van der Waals ratio: ', font=self.top.font_Text).pack(side=LEFT)
        entPermea = Entry(fPermeaLine2, textvariable=self.Permeability, font=self.top.font_Text, justify=CENTER, width=4)
        entPermea.pack(side=RIGHT)
        args_list = [entPermea, self.Permeability, 0.00, 1.00, 2, self.ValidPermeability,'VDW Permeability','float']
        entPermea.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidPermeability[3] = entPermea


        #==================================================================================
        # Delta (variations of distances/angles)
        #==================================================================================
        fSearchSpace = Frame(fC3Right, border=2, relief=RAISED)
        fSearchSpace.pack(fill=BOTH, expand=True, padx=5, pady=5)

        Label(fSearchSpace, text='Search space', font=self.top.font_Title).pack(side=TOP, fill=X, expand=True)
        
        fDelta = Frame(fSearchSpace)#, bd=1, relief=SUNKEN)
        fDelta.pack(fill=X, side=TOP, padx=5, pady=5)
        fDeltaLine1 = Frame(fDelta)
        fDeltaLine1.pack(fill=X, side=TOP, padx=5, pady=2)
        fDeltaLine2 = Frame(fDelta)
        fDeltaLine2.pack(fill=X, side=TOP, padx=5, pady=2)
        fDeltaLine3 = Frame(fDelta)
        fDeltaLine3.pack(fill=X, side=TOP, padx=5, pady=2)
        fDeltaLine4 = Frame(fDelta)
        fDeltaLine4.pack(fill=X, side=TOP, padx=5, pady=2)
        
        Label(fDeltaLine1, text='DELTA parameters', font=self.top.font_Title).pack(side=LEFT, anchor=W)

        # Angles
        Label(fDeltaLine2, text='Angle from reference: ', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        entDAng = Entry(fDeltaLine2, textvariable=self.DeltaAngle, font=self.top.font_Text, justify=CENTER, width=5)
        entDAng.pack(side=RIGHT, anchor=W)
        args_list = [entDAng, self.DeltaAngle, 1.0, 10.0, 2, self.ValidDeltaAngle,'DELTA Angle','float']
        entDAng.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidDeltaAngle[3] = entDAng

        # Dihedrals
        Label(fDeltaLine3, text='Dihedrals from reference: ', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        entDDih = Entry(fDeltaLine3, textvariable=self.DeltaDihedral, font=self.top.font_Text, justify=CENTER, width=5)
        entDDih.pack(side=RIGHT, anchor=W)
        args_list = [entDDih, self.DeltaDihedral, 1.0, 10.0, 2, self.ValidDeltaDihedral,'DELTA Dihedral','float']
        entDDih.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidDeltaDihedral[3] = entDDih

        # Dihedrals (Flex Bonds)
        Label(fDeltaLine4, text='Dihedrals of flexible bonds: ', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        entDDihFlex = Entry(fDeltaLine4, textvariable=self.DeltaDihedralFlex, font=self.top.font_Text, width=5, justify=CENTER)
        entDDihFlex.pack(side=RIGHT, anchor=W)
        args_list = [entDDihFlex, self.DeltaDihedralFlex, 1.0, 30.0, 2, self.ValidDeltaDihedralFlex,'DELTA Dihedral Flexible Bonds','float']
        entDDihFlex.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidDeltaDihedralFlex[3] = entDDihFlex

        #==================================================================================
        # Side-chain optimization (DEE)
        #==================================================================================
        fDEE = Frame(fSearchSpace)#, bd=1, relief=SUNKEN)
        fDEE.pack(fill=X, side=TOP, padx=5, pady=5)
        fDEELine1 = Frame(fDelta)
        fDEELine1.pack(fill=X, side=TOP, padx=5, pady=2)
        fDEELine2 = Frame(fDelta)
        fDEELine2.pack(fill=X, side=TOP, padx=5, pady=2)
        fDEELine3 = Frame(fDelta)
        fDEELine3.pack(fill=X, side=TOP, padx=5, pady=2)
        fDEELine4 = Frame(fDelta)
        fDEELine4.pack(fill=X, side=TOP, padx=5, pady=2)

        Label(fDEELine1, text='Side-chain optimization', font=self.top.font_Title).pack(side=TOP, anchor=W)        
        Checkbutton(fDEELine2, text='Dead-end elimination theorem', variable=self.UseDEE, font=self.top.font_Text).pack(side=LEFT)
        
        Label(fDEELine3, text='Clashing Threshold:', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        self.entDEE = Entry(fDEELine3, width=5, font=self.top.font_Text, textvariable=self.DEE_Clash_Threshold, state='disabled', justify=CENTER)
        self.entDEE.pack(side=RIGHT)
        args_list = [self.entDEE, self.DEE_Clash_Threshold, 0.01, 1.00, 2, self.ValidDEE,'DEE Clash Threshold','float']
        self.entDEE.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidDEE[3] = self.entDEE

        #==================================================================================
        # Implicit Solvent Type
        #==================================================================================

        fSolvent = Frame(fC3Left)#, bd=1, relief=SUNKEN)
        fSolvent.pack(fill=X, side=TOP, padx=5, pady=5)
        fSolventLine1 = Frame(fSolvent)
        fSolventLine1.pack(fill=X, side=TOP, padx=5, pady=2)
        fSolventLine2 = Frame(fSolvent)
        fSolventLine2.pack(fill=X, side=TOP, padx=5, pady=2)
        fSolventLine3 = Frame(fSolvent)
        fSolventLine3.pack(fill=X, side=TOP, padx=5, pady=2)

        Label(fSolventLine1, text='Implicit solvent properties', font=self.top.font_Title).pack(side=LEFT)
        Label(fSolventLine2, text='Solvent type:', font=self.top.font_Text).pack(side=LEFT)

        optionTuple = '< No type >', '< Type-based >',
        self.optSolventType = apply(OptionMenu, (fSolventLine2, self.SolventType) + optionTuple)
        self.optSolventType.config(bg=self.top.Color_White, font=self.top.font_Text, width=10)
        self.optSolventType.pack(side=RIGHT, fill=X, expand=True)

        Label(fSolventLine3, text='Solvent term:', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        self.entSolventTerm = Entry(fSolventLine3, textvariable=self.SolventTerm, font=self.top.font_Text, width=4, justify=CENTER)
        self.entSolventTerm.pack(side=RIGHT)
        args_list = [self.entSolventTerm, self.SolventTerm, -10.0, 10.0, 1, self.ValidSolventTerm,'Solvent Term','float']
        self.entSolventTerm.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidSolventTerm[3] = self.entSolventTerm

    ''' ==================================================================================
    FUNCTION MenuLoadMessage: Display the message based on the menu selected
    ================================================================================== '''        
    def LoadMessage(self):
        
        self.top.DisplayMessage('', 0)
        self.top.DisplayMessage('  FlexAID < Scoring Cfg > Menu.', 2)
        self.top.DisplayMessage('  INFO:   Set different options of the scoring function.', 2)


