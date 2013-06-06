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

import math
import os
import shutil

import Vars
import Tabs
import tkFileDialog
import tkMessageBox
import General
import Smiles
import Constants
import ProcessLigand

if __debug__:
    from pymol import cmd
    
    import Anchor
    import General_cmd


'''
@title: FlexAID - IOFile tab - Interface

@summary: This is the IOFile tab interface of FlexAID application

@organization: Najmanovich Research Group
@creation date:  Nov. 24 2011
'''


class IOFileVars(Vars.Vars):
    
    Complex = StringVar()
    #TargetMD5 = StringVar()
    ProcessedTargetPath = StringVar()
    TargetPath = StringVar()
    TargetName = StringVar()
    #LigandMD5 = StringVar()
    ProcessedLigandPath = StringVar()
    ProcessedLigandINPPath = StringVar()
    ProcessedLigandICPath = StringVar()
    LigandPath = StringVar()
    LigandName = StringVar()
    AtomTypes = StringVar()
    SmilesString = StringVar()
    SmilesName = StringVar()
    Gen3D = IntVar()
    Anchor = IntVar()
    
    def __init__(self):
        
        self.dictAtomTypes = dict()
        self.dictNeighbours = dict()
        self.dictFlexBonds = dict()
    
    
class IOFile(Tabs.Tab):
    
    ATOM_INDEX = 90000
    RESIDUE_NUMBER = 9999
    
    LigandSupportedFormats = [ ('All supported formats', ('*.pdb','*.mol','*.mol2','*.sdf','*.smi')),
                               ('PDB File','*.pdb'),
                               ('MOL File','*.mol'),
                               ('MOL2 File','*.mol2'),
                               ('SDF File','*.sdf'),
                               ('SMI File','*.smi') ]

    TargetSupportedFormats = [ ('PDB File','*.pdb') ]
    
    SmilesLigand = 'SMILES_LIGAND__'
    ReferenceLigand = 'REFERENCE_LIGAND__'
    Ligand = 'LIGAND__'
    Target = 'TARGET__'
    ExtractObject = 'EXTRACTED_OBJECT__'
    
    def Def_Vars(self):
        
        self.defaultOption = StringVar()
        self.FetchPDB = StringVar()
        #self.ReferencePath = StringVar()
        self.ResSeq = IntVar()
                
        # vars class objects
        self.Complex = self.Vars.Complex

        self.ProcessedTargetPath = self.Vars.ProcessedTargetPath
        #self.TargetMD5 = self.Vars.TargetMD5
        self.TargetPath = self.Vars.TargetPath
        self.TargetName = self.Vars.TargetName

        self.ProcessedLigandPath = self.Vars.ProcessedLigandPath
        self.ProcessedLigandINPPath = self.Vars.ProcessedLigandINPPath
        self.ProcessedLigandICPath = self.Vars.ProcessedLigandICPath
        #self.LigandMD5 = self.Vars.LigandMD5
        self.LigandPath = self.Vars.LigandPath
        self.LigandName = self.Vars.LigandName

        self.AtomTypes = self.Vars.AtomTypes
        self.SmilesString = self.Vars.SmilesString
        self.SmilesName = self.Vars.SmilesName

        self.Gen3D = self.Vars.Gen3D
        self.Anchor = self.Vars.Anchor

    def Init_Vars(self):

        self.ProcessedTargetPath.set('')
        #self.TargetMD5.set('')
        self.TargetPath.set('')
        self.TargetName.set('')

        self.ProcessedLigandPath.set('')
        self.ProcessedLigandINPPath.set('')
        self.ProcessedLigandICPath.set('')
        #self.LigandMD5.set('')
        self.LigandPath.set('')
        self.LigandName.set('')

        self.SmilesString.set('')
        self.SmilesName.set('')

        #self.ReferencePath.set('')
        self.defaultOption.set('')
        self.FetchPDB.set('')
        self.AtomTypes.set('Sybyl')
        self.Complex.set('')
        
        self.ResSeq.set(0)
        self.Gen3D.set(0)
        self.Anchor.set(-1)
        
        # flags
        self.ProcessOnly = False
        
        self.Vars.dictAtomTypes.clear()
        self.Vars.dictNeighbours.clear()
        self.Vars.dictFlexBonds.clear()
    
    """
    ''' ==================================================================================
    FUNCTION Check_Integrity: Checks if the referenced files are the same upon loading    
    =================================================================================  '''
    def Check_Integrity(self):
        
        if (self.TargetName.get() != '' and General.hashfile(self.TargetPath.get()) != self.TargetMD5.get()) or \
           (self.LigandName.get() != '' and General.hashfile(self.LigandPath.get()) != self.LigandMD5.get()):
            return 1
            
        return 0
    """
    
    ''' ==================================================================================
    FUNCTION Load_Session: Actions related to when a new session is loaded
    =================================================================================  '''    
    def Load_Session(self):
        
        self.Btn_DisplayObject_Clicked('Ligand')
        self.Btn_DisplayObject_Clicked('Target')
    
        if self.ProcessedLigandPath.get():
            self.Load_ProcConvLigand(self.ProcessedLigandPath.get(), self.Ligand, True)

        if self.ProcessedTargetPath.get():
            self.Load_ProcConvLigand(self.ProcessedTargetPath.get(), self.Target, True)

    ''' ==================================================================================
    FUNCTION Before_Kill_Frame: Actions related before killing the frame
    =================================================================================  '''    
    def Before_Kill_Frame(self):
               
        # Process ligand
        if not self.ProcessedLigandPath.get():
            self.top.ProcessError = False
            
            self.ProcessLigand( True, self.LigandPath.get(), self.ATOM_INDEX, self.AtomTypes.get(), 
                                self.Anchor.get(), False, self.ProcessOnly, self.Gen3D.get(), False )
            
            if self.top.ProcessError:
                return False
            else:
                self.ProcessedLigandPath.set(os.path.join(self.top.FlexAIDTempProject_Dir,
                                                          os.path.split(os.path.splitext(self.LigandPath.get())[0])[1]) + '_ref.pdb')
                self.ProcessedLigandINPPath.set(os.path.join(self.top.FlexAIDTempProject_Dir,
                                                             os.path.split(os.path.splitext(self.LigandPath.get())[0])[1]) + '.inp')
                self.ProcessedLigandICPath.set(os.path.join(self.top.FlexAIDTempProject_Dir,
                                                            os.path.split(os.path.splitext(self.LigandPath.get())[0])[1]) + '.ic')
                
                self.ResSeq.set(self.RESIDUE_NUMBER)
                
                #ReferencePath = os.path.join(self.top.FlexAIDTempProject_Dir, 'LIGAND_ref.pdb')
                #self.ReferencePath.set(ReferencePath)
        
                # Store content of ligand input files
                if self.store_InpFile():
                    return False
                else:
                    if self.Load_ProcConvLigand(self.ProcessedLigandPath.get(), self.Ligand, False):
                        return False
                        
        if not self.ProcessedTargetPath.get():
            self.top.ProcessError = False
            
            self.ProcessTarget( True, self.TargetPath.get(), self.ATOM_INDEX, self.AtomTypes.get(),
                                self.Anchor.get(), False, False, False, True )
            
            if self.top.ProcessError:
                return False
            else:
                self.ProcessedTargetPath.set(os.path.join(self.top.FlexAIDTempProject_Dir,
                                                          os.path.split(os.path.splitext(self.TargetPath.get())[0])[1]) + '.inp.pdb')
                
        return True
    
    ''' ==================================================================================
    FUNCTION ProcessLigand: Processes ligand PDB file using Process_Ligand
    ================================================================================== '''    
    def ProcessLigand(self, boolRun, MoleculeFile, StartAtomIndex, AtomTypes,
                      Anchor, ConvertOnly, ProcessOnly, Gen3D, Target):
        
        if boolRun:
            self.Disable_Frame()
            
            self.top.ProcessRunning = True
            self.DisplayMessage("  Processing the ligand...", 1)
            
            p = ProcessLigand.ProcLig(self, MoleculeFile, StartAtomIndex, AtomTypes, Anchor,
                                      ConvertOnly, ProcessOnly, Gen3D, Target, self.ProcessLigand)
            p.join()
            
        else:
            self.Enable_Frame()
        
    ''' ==================================================================================
    FUNCTION ProcessTarget: Processes target PDB file
    ================================================================================== '''
    def ProcessTarget(self, boolRun, MoleculeFile, StartAtomIndex, AtomTypes, 
                      Anchor, ConvertOnly, ProcessOnly, Gen3D, Target):
        
        if boolRun:
            self.Disable_Frame()
            
            self.top.ProcessRunning = True
            self.DisplayMessage("  Processing the target...", 1)
            
            p = ProcessLigand.ProcLig(self, MoleculeFile, StartAtomIndex, AtomTypes, Anchor,
                                      ConvertOnly, ProcessOnly, Gen3D, Target, self.ProcessTarget)
            p.join()
            
        else:
            self.Enable_Frame()
    
    ''' ==================================================================================
    FUNCTION Convert_Smiles: Converts the smiles string to a viewable molecule format
    ================================================================================== '''
    def Convert_Smiles(self):

        self.ProcessLigand( True, self.LigandPath.get(), 0, '', 0, True, False, 1, False)
        
        if self.top.ProcessError:
            return 1
        
        return 0
        
    ''' ==================================================================================
    FUNCTION Convert_Smiles: Write the Smiles to a .smi file
    ================================================================================== '''    
    def Write_Smiles(self, Filebase):
        
        LigandPath = os.path.join(self.top.FlexAIDTempProject_Dir, Filebase + '.smi')
        
        try:
            handle = open(LigandPath, 'w')
            handle.write(self.SmilesString.get() + '\n')
            handle.close()
        except IOError:
            self.DisplayMessage("  ERROR: Could not write Smiles string to a file.", 2)
            return 1
            
        self.LigandPath.set(LigandPath)
        
        return 0
    
    ''' ==================================================================================
    FUNCTION SmilesRunning: Disables all controls when smiles windows is opened
    ================================================================================== '''    
    def SmilesRunning(self, boolRun, Convert):

        if boolRun:
            self.Disable_Frame()
            
        else:
            # When a user inputs a SMILES string, need to first write to a .smi file
            
            self.Enable_Frame()

            if Convert and self.SmilesString.get():
                
                Filebase = self.SmilesLigand
                if self.SmilesName.get():
                    Filebase = self.SmilesName.get()
                    
                if not self.Write_Smiles(Filebase) and not self.Convert_Smiles() and \
                   not self.Move_TempLigand(Filebase):
                    
                    LigandFile = os.path.join(self.top.FlexAIDTempProject_Dir, Filebase + '.mol2')
                    
                    try:
                        cmd.delete(Filebase)
                        cmd.refresh()
                    except:
                        pass
                    
                    self.LigandPath.set(LigandFile)
                    
                    if not self.Load_ProcConvLigand(LigandFile, Filebase, True) and \
                       not self.Validate_ObjectSelection(Filebase, 'Ligand', 1):
                           
                        self.LigandName.set(Filebase)
                        
                        if self.ForceSaveObject(LigandFile, Filebase, 'Ligand'):
                            self.LigandPath.set('')
                            self.LigandName.set('')
                        
                        #self.ProcessOnly = True
                                                
                    else:
                        self.LigandPath.set('')
                        self.LigandName.set('')
                else:
                    self.LigandPath.set('')
                    self.LigandName.set('')
    
                self.ProcessedLigandPath.set('')
                self.ProcessedLigandINPPath.set('')
                self.ProcessedLigandICPath.set('')

            self.top.ChildWindow = None
            
    ''' ==================================================================================
                         Reset tabs when the target/ligand is modified
    ================================================================================== '''             
    def Target_Toggle(self, *args):

        # Reset binding-site and target flexibility
        self.top.Config1.Init_Vars()

        # Clear constraint because there might be ligand-target constraints
        self.top.Config2.Vars.dictConstraints.clear()

        self.top.Simulate.Init_Vars()

        self.ValidateLigTarget()
        
    def Ligand_Toggle(self, *args):
        
        self.ValidateLigTarget()

    ''' ==================================================================================
                         ENABLE / DISABLE - Buttons
    ================================================================================== '''             
    def ValidateLigTarget(self):
        
        if self.TargetName.get() and self.LigandName.get():
            self.top.Go_Step2()
        else:
            self.top.Go_Step1()
        
        self.Complex.set(self.TargetName.get() + '-' + self.LigandName.get())
        
    ''' ==================================================================================
    FUNCTION Reset_Ligand: Rests all parameters variables to processing of the ligand
    ================================================================================== '''    
    def Reset_Ligand(self):
        
        self.Gen3D.set(0)
        self.Anchor.set(-1)
        
        self.Vars.dictAtomTypes.clear()
        self.Vars.dictFlexBonds.clear()
        self.Vars.dictNeighbours.clear()
        
        self.top.Config2.Init_Vars()
        
        self.top.Simulate.Init_Vars()
        
        self.ProcessOnly = False
        
        self.ProcessedLigandPath.set('')
        self.ProcessedLigandINPPath.set('')
        self.ProcessedLigandICPath.set('')

    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):
                
        #Show the list of selection/objects in case the user already worked in PyMOL
        self.Btn_RefreshOptMenu_Clicked()
        
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        try:
            self.TargetNameTrace = self.TargetName.trace('w',self.Target_Toggle)
            self.LigandNameTrace = self.LigandName.trace('w',self.Ligand_Toggle)
            #self.AtomTypesTrace = self.AtomTypes.trace('w',self.AtomTypes_Toggle)
            #self.Gen3DTrace = self.Gen3D.trace('w',self.Gen3D_Toggle)
        except:
            pass
            
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        try:
            self.TargetName.trace_vdelete('w',self.TargetNameTrace)
            self.LigandName.trace_vdelete('w',self.LigandNameTrace)
            #self.AtomTypes.trace_vdelete('w',self.AtomTypesTrace)
            #self.Gen3D.trace_vdelete('w',self.Gen3DTrace)
        except:
            pass
    
    ''' ==================================================================================
    FUNCTION AtomTypes_Toggle: Toggle the controls related to Atom Types
    =================================================================================  '''    
    #def AtomTypes_Toggle(self, *args):
    def AtomTypes_Toggle(self):
        
        if self.AtomTypes.get() != 'Gaudreault':
            self.top.Config3.optSolventType.config(state='normal')
            self.top.Config3.SolventType.set('< No type >')
        else:
            self.top.Config3.optSolventType.config(state='disabled')
            self.top.Config3.SolventType.set('< Type-based >')

        # Reset atom type definition
        self.Vars.dictAtomTypes.clear()

        # Need for processing the ligand again if atom typing is changed
        self.ProcessedTargetPath.set('')

        self.ProcessedLigandPath.set('')
        self.ProcessedLigandINPPath.set('')
        self.ProcessedLigandICPath.set('')

    ''' ==================================================================================
    FUNCTION Gen3D_Toggle: Toggles the generation of a 3D molecule
    =================================================================================  '''    
    #def Gen3D_Toggle(self, *args):
    def Gen3D_Toggle(self):
        
        self.ProcessedLigandPath.set('')
        self.ProcessedLigandINPPath.set('')
        self.ProcessedLigandICPath.set('')
        
    ''' ==================================================================================
    FUNCTION Frame: Generate the Input / Output Files frame in the the middle 
                    frame section    
    ==================================================================================  '''  
    def Frame(self):
        
        self.fIOFile = Frame(self.top.fMiddle)

        #==================================================================================
        #                                 Object/selections
        #==================================================================================
        fPDB_options2 = Frame(self.fIOFile)#, border=1, relief=SUNKEN)
        fPDB_options2.pack(fill=X, side=TOP, padx=5, pady=5)

        fPDB_options2Line1 = Frame(fPDB_options2)#, border=1, relief=SUNKEN)
        fPDB_options2Line1.pack(side=TOP, fill=X)
        fPDB_options2Line2 = Frame(fPDB_options2)#, border=1, relief=SUNKEN)
        fPDB_options2Line2.pack(side=TOP, fill=X)
        fPDB_options2Line3 = Frame(fPDB_options2)#, border=1, relief=SUNKEN)
        fPDB_options2Line3.pack(side=TOP, fill=X)

        # Header Select
        Label(fPDB_options2Line1, text='Select an object or selection', font=self.font_Title).pack(side=LEFT)

        # List of selectionsSave
        Label(fPDB_options2Line2, text='PyMOL objects/selections:', width=25, justify=RIGHT, font=self.font_Text).pack(side=LEFT, anchor=E)

        optionTuple = ('',)
        self.optionMenuWidget = apply(OptionMenu, (fPDB_options2Line2, self.defaultOption) + optionTuple)
        self.optionMenuWidget.config(bg=self.Color_White, width=15, font=self.font_Text)
        self.optionMenuWidget.pack(side=LEFT)
        
        # Refresh the list with the selections in Pymol        
        Button(fPDB_options2Line2, text='Refresh', command=self.Btn_RefreshOptMenu_Clicked, font=self.font_Text).pack(side=LEFT)

        # List of selections
        #Button(fPDB_options2Line2, text='LIGAND', font=self.font_Text, relief=RIDGE, command=self.Btn_SetLigand_Clicked).pack(side=RIGHT)
        #Button(fPDB_options2Line2, text='TARGET', font=self.font_Text, relief=RIDGE, command=self.Btn_SetTarget_Clicked).pack(side=RIGHT)

        #Label(fPDB_options2Line2, text='Set as...', justify=RIGHT, font=self.font_Text).pack(side=RIGHT, anchor=E)

        # List of selections
        Button(fPDB_options2Line2, text='Save as target', command=lambda objtype='Target': self.Btn_SaveObject_Clicked(objtype), 
                                                          font=self.font_Text).pack(side=LEFT, padx=10)
        
        Button(fPDB_options2Line2, text='Save as ligand', command=lambda objtype='Ligand': self.Btn_SaveObject_Clicked(objtype),
                                                          font=self.font_Text).pack(side=LEFT)
        
        Button(fPDB_options2Line2, text='Extract as ligand', command=self.Btn_ExtractLigand_Clicked, font=self.font_Text).pack(side=LEFT)

        #Label(fPDB_options2Line2, text='Save as...', justify=RIGHT, font=self.font_Text).pack(side=RIGHT, anchor=E)
        

               
        #==================================================================================
        #                                SET TARGET
        #==================================================================================                
        
        fPDBTarget = Frame(self.fIOFile, border=1, relief=RAISED, width=500, height=70)
        fPDBTarget.pack(side=TOP, pady=10, padx=10)#, fill=X, expand=True)
        fPDBTarget.pack_propagate(0)

        fPDBTargetLine1 = Frame(fPDBTarget)
        fPDBTargetLine1.pack(side=TOP, fill=X, padx=3, pady=3)

        fPDBTargetLine2 = Frame(fPDBTarget)
        fPDBTargetLine2.pack(side=TOP, fill=X, padx=3, pady=3)
        
        # First line
        Label(fPDBTargetLine1, width=20, text='THE TARGET', font=self.font_Title).pack(side=LEFT)
        Button(fPDBTargetLine1, text='Load', command=lambda objtype='Target': self.Btn_LoadObject_Clicked(objtype), 
                                              font=self.font_Text).pack(side=LEFT)
        Button(fPDBTargetLine1, text='Display', command=lambda objtype='Target': self.Btn_DisplayObject_Clicked(objtype),
                                                 font=self.font_Text).pack(side=LEFT)
        Button(fPDBTargetLine1, text='Reset', command=self.Btn_ResetTarget_Clicked, font=self.font_Text).pack(side=LEFT)

        # Second line
        Label(fPDBTargetLine2, width=20, text='', font=self.font_Title).pack(side=LEFT)
        EntTarget = Entry(fPDBTargetLine2, textvariable=self.TargetName, disabledbackground=self.Color_White, 
                            disabledforeground=self.Color_Black, font=self.font_Text, justify=CENTER, width=20)
        EntTarget.pack(side=LEFT, fill=X)
        EntTarget.config(state='disabled')
        #Checkbutton(fPDBTargetLine2, variable=self.TargetRNA, width=10, text='RNA', font=self.font_Text, justify=LEFT).pack(side=LEFT)
        
        #==================================================================================
        #                               SET LIGAND
        #==================================================================================   

        fPDBligand = Frame(self.fIOFile, border=1, relief=RAISED, width=500, height=70)
        fPDBligand.pack(side=TOP, pady=10, padx=10)
        fPDBligand.pack_propagate(0)

        fPDBligandLine1 = Frame(fPDBligand)
        fPDBligandLine1.pack(side=TOP, fill=X, padx=3, pady=3)

        fPDBligandLine2 = Frame(fPDBligand)
        fPDBligandLine2.pack(side=TOP, fill=X, padx=3, pady=3)

        # First line
        Label(fPDBligandLine1, width=20, text='THE LIGAND', font=self.font_Title).pack(side=LEFT)
        Button(fPDBligandLine1, text='Load', command=lambda objtype='Ligand': self.Btn_LoadObject_Clicked(objtype), 
                                             font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Display', command=lambda objtype='Ligand': self.Btn_DisplayObject_Clicked(objtype),
                                                font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Reset', command=self.Btn_ResetLigand_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Input', command=self.Btn_Input_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Anchor', command=self.Btn_Anchor_Clicked, font=self.font_Text).pack(side=LEFT)
        
        # Second line
        Label(fPDBligandLine2, width=20, text='', font=self.font_Title).pack(side=LEFT)
        EntLigand = Entry(fPDBligandLine2, disabledbackground=self.Color_White, disabledforeground=self.Color_Black, 
                            textvariable=self.LigandName, font=self.font_Text, justify=CENTER, width=20)
        EntLigand.pack(side=LEFT, fill=X)
        EntLigand.config(state='disabled')
        Checkbutton(fPDBligandLine2, variable=self.Gen3D, text=' Generate 3D conformation',
            font=self.font_Text, justify=RIGHT, command=self.Gen3D_Toggle).pack(side=LEFT)

        #==================================================================================
        #                                 Processing of molecules
        #==================================================================================

        fProcessing = Frame(self.fIOFile)
        #fProcessing.pack(side=TOP, fill=X, pady=5, padx=5)

        fProcessingLine1 = Frame(fProcessing)#, border=1, relief=SUNKEN)
        fProcessingLine1.pack(side=TOP, fill=X)
        fProcessingLine2 = Frame(fProcessing)#, border=1, relief=SUNKEN)
        fProcessingLine2.pack(side=TOP, fill=X)
        fProcessingLine3 = Frame(fProcessing)#, border=1, relief=SUNKEN)
        fProcessingLine3.pack(side=TOP, fill=X)
        
        Label(fProcessingLine1, text='Processing of the ligand', font=self.font_Title).pack(side=LEFT)

        Label(fProcessingLine2, text='Atom typing:', justify=RIGHT, font=self.font_Text).pack(side=LEFT, anchor=E, padx=10)

        Radiobutton(fProcessingLine2, text='Sobolev', variable=self.AtomTypes,
                    value="Sobolev", font=self.font_Text).pack(side=LEFT)
        Radiobutton(fProcessingLine2, text='Gaudreault', variable=self.AtomTypes, 
                    value="Gaudreault", font=self.font_Text).pack(side=LEFT, padx=10)
        Radiobutton(fProcessingLine2, text='Sybyl', variable=self.AtomTypes, value="Sybyl", font=self.font_Text).pack(side=LEFT)
        
        return self.fIOFile
    
    ''' ==================================================================================
    FUNCTIONS Reset Ligand and Target textbox fields
    ================================================================================== '''
    def Btn_ResetTarget_Clicked(self):
    
        self.TargetPath.set('')
        self.TargetName.set('')

    def Btn_ResetLigand_Clicked(self):

        self.LigandPath.set('')
        self.LigandName.set('')

    ''' ==================================================================================
    FUNCTION Btn_Anchor_Clicked: Selects the anchor atom of the ligand
    ==================================================================================  '''    
    def Btn_Anchor_Clicked(self):

        if not self.PyMOL:
            return

        if self.LigandPath.get():
                            
            self.top.ActiveWizard = Anchor.anchor(self, self.LigandPath.get(), self.Anchor.get())
            
            cmd.set_wizard(self.top.ActiveWizard)
            self.top.ActiveWizard.Start()
    
    ''' ==================================================================================
    FUNCTION Btn_Input_Clicked: Selects the anchor atom of the ligand
    ==================================================================================  '''    
    def Btn_Input_Clicked(self):
        
        if self.PyMOL:
            
            self.SmilesRunning(True, False)
            
            self.top.ChildWindow = Smiles.Smiles(self, self.SmilesString, self.SmilesName)
        
    ''' ==================================================================================
    FUNCTION Set_Object_Variables: Sets class level variables
    ==================================================================================  '''    
    def Set_Object_Variables(self, objtype):

        if objtype == 'Ligand':
            self.savepath = self.top.FlexAIDLigandProject_Dir

            self.VarPath = self.LigandPath
            self.VarName = self.LigandName
            self.VarProc = self.ProcessedLigandPath
            #self.VarMD5 = self.LigandMD5
            
        elif objtype == 'Target':
            self.savepath = self.top.TargetProject_Dir        

            self.VarPath = self.TargetPath
            self.VarName = self.TargetName
            self.VarProc = self.ProcessedTargetPath
            #self.VarMD5 = self.TargetMD5

    ''' ==================================================================================
    FUNCTION Validate_ObjectSelection: Validates whether the obj exists on the current state
                                       and has at least N atoms
    ==================================================================================  '''    
    def Validate_ObjectSelection(self, sele, objtype, state):
        
        try:
            n = cmd.count_atoms(sele + ' & ! n. H*' , state=state)
            
            if n < 5:
                self.DisplayMessage("  ERROR for object/selection '" + sele + "': The object must have at least (5) heavy atoms)", 1)
                return 1
            
            elif objtype == 'Ligand' and n > Constants.MAX_LIGAND_ATOMS:
                self.DisplayMessage("  ERROR for object/selection '" + sele + "': The ligand must have a maximum of (" + 
                                    str(Constants.MAX_LIGAND_ATOMS) + ") atoms", 1)
                return 1
            
        except:
            self.DisplayMessage("  ERROR: object/selection '" + sele + "' does not exist on current state", 1)
            return 1
        
        return 0
        
    ''' ==================================================================================
    FUNCTION ForceSaveObject: Forces saving of the object in the corresponding folder
    ==================================================================================  '''    
    def ForceSaveObject(self, file, objname, objtype):
        
        basepath, filename = os.path.split(file)
        filebasename, fileextname = os.path.splitext(filename)
        
        self.Set_Object_Variables(objtype)
        
        if self.savepath != basepath:
        
            newfile = os.path.join(self.savepath,filebasename + '.pdb')
            
            if os.path.isfile(newfile):
                answer = tkMessageBox.askquestion("Question", 
                                                  message=  "An object with that name already exists in your '" + \
                                                            objtype + "' folder. Would you like to overwrite it?",
                                                  icon='warning')
                if str(answer) == 'no':
                    return 2
            
            try:
                cmd.save(newfile, objname, state=1)
            except:
                self.DisplayMessage("  ERROR: The file could not be saved in your '" + objtype + "' folder.", 2)
                return 1
            
            self.VarPath.set(newfile)
            
        return 0
    
    ''' ==================================================================================
    FUNCTION Btn_SaveObjectClicked: Save Ligand and Target objects
                                    Object is reloaded when renamed
    ==================================================================================  '''    
    def Btn_SaveObject_Clicked(self, objtype):

        # Get the Drop Down List Selection Name
        ddlSelection = self.defaultOption.get()

        state = cmd.get_state()

        self.Set_Object_Variables(objtype)
        
        if ddlSelection == '' or self.Validate_ObjectSelection(ddlSelection, objtype, state):
            return
        
        Path = tkFileDialog.asksaveasfilename(filetypes=[('PDB File','*.pdb')],
                                              initialdir=self.savepath, title='Save the PDB File',
                                              initialfile=ddlSelection, defaultextension='.pdb')
        
        if len(Path) > 0:

            Path = os.path.normpath(Path)

            if General.validate_String(Path, '.pdb', True, False, True):
                self.DisplayMessage("  ERROR: Could not save the file because you entered an invalid name.", 2)
                return

            if self.top.ValidateSaveProject(Path, objtype):
                self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                return
            
            try:
                cmd.save(Path, ddlSelection, state)
                Name = os.path.basename(os.path.splitext(Path)[0])

                cmd.load(Path, Name, state=1)
                cmd.refresh()
                
                if ddlSelection != Name:
                    # as if the object was renamed, delete the object
                    cmd.delete(ddlSelection)
                
            except:
                self.DisplayMessage("  ERROR: An error occured while saving the object.", 1)
                return
            
            self.VarPath.set(os.path.normpath(Path))
            self.VarName.set(Name)
            self.VarProc.set('')
            self.top.SaveSessionFile = ''
            #self.VarMD5.set(General.hashfile(self.VarPath.get()))
            
            if objtype == 'Ligand':
                self.Reset_Ligand()
            
            self.DisplayMessage('  Successfully saved and loaded the object:  ' + self.VarName.get() + "'", 0)
        
    ''' ==================================================================================
    FUNCTIONS Get_SupportedFormats: Returns the supported formats depending on the object type
    ================================================================================== '''        
    def Get_SupportedFormats(self, objtype):
    
        if objtype == 'Ligand':
            return self.LigandSupportedFormats
        else:
            return self.TargetSupportedFormats

    ''' ==================================================================================
    FUNCTIONS Load Ligand and Target and display the filename in the textbox
    ================================================================================== '''        
    def Btn_LoadObject_Clicked(self, objtype):
        
        self.Set_Object_Variables(objtype)
        
        SupportedFormats = self.Get_SupportedFormats(objtype)
        
        Path = tkFileDialog.askopenfilename(filetypes=SupportedFormats,
                                            initialdir=self.savepath, title='Select a file to Load')

        if len(Path) > 0:
        
            Path = os.path.normpath(Path)
            
            if General.validate_String(Path, '', True, False, True):
                self.DisplayMessage("  ERROR: Could not load the file because it has an invalid name.", 2)
                return

            if Path == self.VarPath.get():
                self.DisplayMessage("  Loading skipped. File is the same as the one already loaded.", 2)
                return
            
            try:
                Name = os.path.basename(os.path.splitext(Path)[0])

                cmd.load(Path, Name, state=1)
                cmd.refresh()
                
                if self.Validate_ObjectSelection(Name, objtype, 1):
                    return
                                        
            except:
                self.DisplayMessage("  ERROR: An error occured while loading the file.", 1)
                return
            
            self.VarPath.set(Path)
            self.VarName.set(Name)
            self.VarProc.set('')
            self.top.SaveSessionFile = ''

            if objtype == 'Ligand':
                self.Reset_Ligand()

            if self.ForceSaveObject(Path,Name,objtype):
                self.VarPath.set('')
                self.VarName.set('')
            
            #self.VarMD5.set(General.hashfile(self.VarPath.get()))

            self.DisplayMessage("  Successfully loaded the object: '" + self.VarName.get() + "'", 0)
        
    ''' ==================================================================================
    FUNCTIONS Btn_ExtractLigand_Clicked: Extracts the selection out of the current object
    ================================================================================== '''        
    def Btn_ExtractLigand_Clicked(self):
                
        state = cmd.get_state()

        # Get the Drop Down List Selection Name
        ddlSelection = self.defaultOption.get()
        
        if ddlSelection == '' or self.Validate_ObjectSelection(ddlSelection, 'Ligand', state):
            return
        
        LigandPath = tkFileDialog.asksaveasfilename(filetypes=[('PDB File','*.pdb')],
                                                    initialdir=self.top.FlexAIDLigandProject_Dir, title='Save the PDB File', 
                                                    initialfile=ddlSelection, defaultextension='.pdb')
        
        if len(LigandPath) > 0:
            
            LigandPath = os.path.normpath(LigandPath)

            if General.validate_String(LigandPath, '.pdb', True, False, True):
                self.DisplayMessage("  ERROR: Could not save the file because you entered an invalid name.", 2)
                return

            if self.top.ValidateSaveProject(LigandPath, 'Ligand'):
                self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                return

            try:
                cmd.save(LigandPath, ddlSelection, state)
                LigandName = os.path.basename(os.path.splitext(LigandPath)[0])

                cmd.extract(self.ExtractObject, ddlSelection)
                cmd.set_name(self.ExtractObject, LigandName)
                
                if ddlSelection != LigandName:
                    # as if the object was renamed, delete the object
                    cmd.delete(ddlSelection)
                
            except:
                self.DisplayMessage("  ERROR: An error occured while extracting the ligand object.", 1)
                return
            
            self.LigandPath.set(os.path.normpath(LigandPath))
            self.LigandName.set(LigandName)
            #self.LigandMD5.set(General.hashfile(self.LigandPath.get()))

            self.Reset_Ligand()
            
            self.DisplayMessage('  Successfully extracted the ligand:  ' + self.LigandName.get() + "'", 0)
    
    def Load_ProcConvLigand(self, LigandFile, ObjectName, Zoom):
        
        Error = 0
        
        try:
            cmd.delete(ObjectName)
            cmd.refresh()
        except:
            pass
        
        auto_zoom = cmd.get("auto_zoom")
        
        if self.PyMOL:
            try:
                cmd.set("auto_zoom", 0)                
                cmd.load(LigandFile, ObjectName, state=1)
                cmd.refresh()
                
                if Zoom:
                    cmd.zoom(ObjectName)
                    cmd.refresh()
            except:
                self.DisplayMessage('  ERROR: Could not load the ligand file in PyMOL', 1)
                Error = 1
        
        cmd.set("auto_zoom", auto_zoom)
        
        return Error

    ''' ==================================================================================
    FUNCTION Btn_RefreshOptMenu_Clicked: Refresh the selections list in the application
                                         with the selections in Pymol 
    ==================================================================================  '''                
    def Btn_RefreshOptMenu_Clicked(self): 
        
        if not self.PyMOL:
            return
            
        General_cmd.Refresh_DDL(self.optionMenuWidget, self.defaultOption, [], None)
   
    ''' ==================================================================================
    FUNCTIONS Move_TempLigand(self):
    ================================================================================== '''        
    def Move_TempLigand(self, Filebase):

        LigandFile = os.path.join(self.top.FlexAIDTempProject_Dir, Filebase + '.mol2')
        TempLigandFile = LigandFile + '.tmp'
        
        try:
            shutil.move(TempLigandFile, LigandFile)
        except IOError:
            self.DisplayMessage("  ERROR: Could not move the temporary ligand file.", 2)
            return 1
        except shutil.Error:
            pass
        
        return 0

    ''' ==================================================================================
    FUNCTIONS Display Ligand and Target in pymol
    ================================================================================== '''        
    def Btn_DisplayObject_Clicked(self, objtype):

        self.Set_Object_Variables(objtype)

        if self.VarName.get() != '':
            try:
                if not General_cmd.object_Exists(self.VarName.get()):
                    cmd.load(self.VarPath.get(), state=1)
                    cmd.refresh()
                else:                    
                    cmd.zoom(self.VarName.get())    
                    cmd.refresh()
            except:
                self.VarName.set('')
                self.VarPath.set('')
                self.VarProc.set('')
                self.top.SaveSessionFile = ''
                
                self.DisplayMessage("  ERROR: An error occured while displaying the object.", 1)
    
    #=======================================================================
    ''' Store inp file information (flexible bonds, atom types)  '''
    #=======================================================================   
    def store_InpFile(self):
            
        inpInfo = dict()
        flexInfo = dict()

        #Read the inp file and get the flexible bonds
        try:
            file = open(self.ProcessedLigandINPPath.get(), 'r')
            inpLines = file.readlines()
            file.close()

            for Line in inpLines:

                #HETTYP  902 1  N1  m   909  910  903    0
                if Line.startswith('HETTYP'):
                    
                    ATOM = Line[6:11].strip()

                    list = []
                    list.append(Line[11:13].strip()) # Type
                    list.append(Line[21:26].strip()) # Neighbour 1
                    list.append(Line[26:31].strip()) #           2
                    list.append(Line[31:36].strip()) #           3
                    inpInfo[ATOM] = list

                #FLEDIH  4   916  917   
                elif Line.startswith('FLEDIH'):

                    INDEX = Line[7:9].strip()

                    list = []
                    for i in range(0,len(Line[10:])/5):
                        list.append(Line[(10+i*5):(10+5+i*5)].strip())
                    flexInfo[INDEX] = list
                
        except:
            self.DisplayMessage('  ERROR: Could not retrieve ligand input file', 1)
            return 1
        
        if not len(self.Vars.dictNeighbours): 
            self.store_Neighbours(inpInfo)
            
        if not len(self.Vars.dictAtomTypes):
            self.store_AtomTypes(inpInfo)
        
        if not len(self.Vars.dictFlexBonds):
            self.store_FlexBonds(flexInfo)
        
        return 0

    #=======================================================================
    ''' Store Neighbours Dictionary'''
    #=======================================================================   
    def store_Neighbours(self, inpInfo):
        
        for atom in inpInfo.keys():
            self.Vars.dictNeighbours[atom] = inpInfo[atom][1:]
                    
    #=======================================================================
    ''' Store Flexible Bonds Dictionary'''
    #=======================================================================   
    def store_FlexBonds(self, flexInfo):
                
        for index in flexInfo.keys():
            
            ''' [ Selected as flexible,
                  Forced as flexible,
                  Number of atoms defining the bond,
                  Atom list defining bond ] '''
            
            dictList = [ 0, 0, len(flexInfo[index]) ]
            for i in range(0, len(flexInfo[index])):
                dictList.append(flexInfo[index][i])
                
            self.Vars.dictFlexBonds[index] = dictList

    #=======================================================================
    ''' Store Atom Types Dictionary'''
    #=======================================================================   
    def store_AtomTypes(self, inpInfo):
        
        for atom in inpInfo.keys():
            self.Vars.dictAtomTypes[atom] = [inpInfo[atom][0], inpInfo[atom][0]]
    
    #=======================================================================
    ''' Enables/disables controls related to AnchorWizard  '''
    #=======================================================================   
    def AnchorRunning(self, boolRun):
        
        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()
            
            if self.Anchor.get() != self.top.WizardResult:
                self.ProcessedLigandPath.set('')
                self.ProcessedLigandINPPath.set('')
                self.ProcessedLigandICPath.set('')
                
                self.Vars.dictAtomTypes.clear()
                self.Vars.dictNeighbours.clear()
                self.Vars.dictFlexBonds.clear()
            
                self.top.Config2.Init_Vars()
            
            self.Anchor.set(self.top.WizardResult)

    # Welcome menu message
    def Load_Message(self):

        self.DisplayMessage('' ,0)
        self.DisplayMessage('  FlexAID < Input Files > Menu', 2)
        self.DisplayMessage('  INFO:   Select a < TARGET > and a < LIGAND > by:', 2)
        self.DisplayMessage('          1) Saving a PyMOL object/selection to your project directory', 2)
        self.DisplayMessage('          2) Loading an existing object file from your project directory', 2)
