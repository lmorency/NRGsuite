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
import General
import Vars
import Tabs 

if __debug__:
    from pymol import cmd

class Config3Vars(Vars.Vars):

    CompFct = StringVar()
    UseDEE = IntVar()
    IncludeHET = IntVar()
    ExcludeHOH = IntVar()
    DEE_Clash_Threshold = StringVar()
    Permeability = StringVar()
    SolventType = StringVar()
    SolventTypeIndex = IntVar()        
    SolventTerm = StringVar()
    DeltaAngle = StringVar()
    DeltaDihedral = StringVar()
    DeltaDihedralFlex = StringVar()


class Config3(Tabs.Tab):

    def Def_Vars(self):

        self.CompFct = self.Vars.CompFct
        self.UseDEE = self.Vars.UseDEE
        self.IncludeHET = self.Vars.IncludeHET
        self.ExcludeHOH = self.Vars.ExcludeHOH
        self.DEE_Clash_Threshold = self.Vars.DEE_Clash_Threshold
        self.Permeability = self.Vars.Permeability
        self.SolventType = self.Vars.SolventType
        self.SolventTypeIndex = self.Vars.SolventTypeIndex
        self.SolventTerm = self.Vars.SolventTerm
        self.DeltaAngle = self.Vars.DeltaAngle
        self.DeltaDihedral = self.Vars.DeltaDihedral
        self.DeltaDihedralFlex = self.Vars.DeltaDihedralFlex

    def Init_Vars(self):
        
        self.CompFct.set('VCT')

        self.UseDEE.set(0)
        self.IncludeHET.set(0)
        self.ExcludeHOH.set(1)

        self.DEE_Clash_Threshold.set('0.50')
        self.Permeability.set('0.1')
        
        if self.top.IOFile.AtomTypes.get() == 'Sobolev':
            self.SolventType.set('< No type >')
            self.SolventTypeIndex.set(0)
        elif self.top.IOFile.AtomTypes.get() == 'Gaudreault':
            self.SolventType.set('< Type-based >')
            self.SolventTypeIndex.set(13)
        elif self.top.IOFile.AtomTypes.get() == 'Sybyl':
            self.SolventType.set('< No type >')
            self.SolventTypeIndex.set(0)

        self.SolventTerm.set('-5.0')

        self.DeltaAngle.set('2.5')
        self.DeltaDihedral.set('2.5')
        self.DeltaDihedralFlex.set('10.0')


    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        try:
            self.IncludeHETTrace = self.IncludeHET.trace('w',self.IncludeHET_Toggle)
            self.top.IOFile.AtomTypesTrace = self.top.IOFile.AtomTypes.trace('w',self.AtomTypes_Toggle)
            self.SolventTypeTrace = self.SolventType.trace('w',self.SolventType_Toggle)
            self.DEETrace = self.UseDEE.trace('w',self.DEE_Toggle)

            #self.PermeabilityTrace = self.Permeability.trace('w', lambda args=args_list: self.Validate_Field(self.Permeability_args_list))
            #self.DeltaAngleTrace = self.DeltaAngle.trace('w', lambda args=args_list: self.Validate_Field(args_list))            
            #self.DeltaDihedralTrace = self.DeltaDihedral.trace('w', lambda args=args_list: self.Validate_Field(args_list))
        except:
            pass
        
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''
    def Del_Trace(self):

        try:
            self.IncludeHET.trace_vdelete('w',self.IncludeHETTrace)
            self.top.IOFile.AtomTypes.trace_vdelete('w',self.top.IOFile.AtomTypesTrace)
            self.SolventType.trace_vdelete('w',self.SolventTypeTrace)
            self.UseDEE.trace_vdelete('w',self.DEETrace)
        except:
            pass        

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
        
        if self.top.IOFile.AtomTypes.get() != 'Gaudreault':
            self.optSolventType.config(state='normal')
            self.SolventType.set('< No type >')
        else:
            self.optSolventType.config(state='disabled')
            self.SolventType.set('< Type-based >')

        # Need for processing the ligand again if atom typing is changed
        self.top.IOFile.fProcessLigand = False

        # Reset atom type definition
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
            if self.top.IOFile.AtomTypes.get() == 'Gaudreault':
                self.SolventTypeIndex.set(13)
            elif self.top.IOFile.AtomTypes.get() == 'Sybyl':
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
        args_list = [entPermea, self.Permeability, 0.00, 1.00, 2, 'VDW Permeability','float']
        #entPermea.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidPermeability = [1, 0, entPermea]
        #self.PermeabilityTrace = self.Permeability.trace('w', lambda args=args_list: self.Validate_Field(args_list))

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
        args_list = [entDAng, self.DeltaAngle, 1.0, 10.0, 2, 'DELTA Angle','float']
        #entDAng.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidDeltaAngle = [1, 0, entDAng]
        #self.DeltaAngleTrace = self.DeltaAngle.trace('w', lambda args=args_list: self.Validate_Field(args_list))

        # Dihedrals
        Label(fDeltaLine3, text='Dihedrals from reference: ', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        entDDih = Entry(fDeltaLine3, textvariable=self.DeltaDihedral, font=self.top.font_Text, justify=CENTER, width=5)
        entDDih.pack(side=RIGHT, anchor=W)
        args_list = [entDDih, self.DeltaDihedral, 1.0, 10.0, 2, 'DELTA Dihedral','float']
        #entDDih.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidDeltaDihedral = [1, 0, entDDih]
        #self.DeltaDihedralTrace = self.DeltaDihedral.trace('w', lambda args=args_list: self.Validate_Field(args_list))

        # Dihedrals (Flex Bonds)
        Label(fDeltaLine4, text='Dihedrals of flexible bonds: ', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        entDDihFlex = Entry(fDeltaLine4, textvariable=self.DeltaDihedralFlex, font=self.top.font_Text, width=5, justify=CENTER)
        entDDihFlex.pack(side=RIGHT, anchor=W)
        args_list = [entDDihFlex, self.DeltaDihedralFlex, 1.0, 30.0, 2, 'DELTA Dihedral Flexible Bonds','float']
        #entDDihFlex.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidDeltaDihedralFlex = [1, 0, entDDihFlex]
        #self.DeltaDihedralFlexTrace = self.DeltaDihedralFlex.trace('w', lambda args=args_list: self.Validate_Field(args_list))

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
        args_list = [self.entDEE, self.DEE_Clash_Threshold, 0.01, 1.00, 2,'DEE Clash Threshold','float']
        #self.entDEE.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidDEE = [1, 0, self.entDEE]
        #self.DEE_Clash_ThresholdTrace = self.DEE_Clash_Threshold.trace('w', lambda args=args_list: self.Validate_Field(args_list))

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
        args_list = [self.entSolventTerm, self.SolventTerm, -10.0, 10.0, 1, 'Solvent Term','float']
        #self.entSolventTerm.config(validate='key', validatecommand=lambda args=args_list: self.Validate_Field(args_list))
        self.ValidSolventTerm = [1, 0, self.entSolventTerm]
        #self.SolventTermTrace = self.SolventTerm.trace('w', lambda args=args_list: self.Validate_Field(args_list))


        self.Validator = [ self.ValidDEE, self.ValidPermeability,
                           self.ValidDeltaAngle, self.ValidDeltaDihedral, 
                           self.ValidDeltaDihedralFlex, self.ValidSolventTerm ]


        return self.fConfig3
        
    ''' ==================================================================================
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''        
    def Load_Message(self):
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < Scoring Cfg > Menu.', 2)
        self.DisplayMessage('  INFO:   Set different options of the scoring function.', 2)


