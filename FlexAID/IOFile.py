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
import hashlib
import Vars
import Tabs
import tkFileDialog
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
    
    ProtPath = StringVar()
    ProtName = StringVar()
    LigandPath = StringVar()
    LigandName = StringVar()
    AtomTypes = StringVar()

    def __init__(self):
        
        self.LigandPathMD5 = ''
    
    
class IOFile(Tabs.Tab):

    ProcessError = False
    
    # flags
    fProcessLigand = False
    fStoreInfo = False
    fLoadProcessed = False
    
    SupportedFormats = [ ('PDB File','*.pdb'),
                         ('MOL File','*.mol'),
                         ('MOL2 File','*.mol2'),
                         ('SDF File','*.sdf'),
                         ('SMI File','*.smi') ]
    
    def Def_Vars(self):
        
        self.defaultOption = StringVar()
        self.FetchPDB = StringVar()
        self.ReferencePath = StringVar()
        self.ResSeq = IntVar()
        
        # vars class objects
        self.ProtPath = self.Vars.ProtPath
        self.ProtName = self.Vars.ProtName
        self.LigandPath = self.Vars.LigandPath
        self.LigandName = self.Vars.LigandName
        self.AtomTypes = self.Vars.AtomTypes
        

    def Init_Vars(self):

        #print("Init Vars for IOFile")
        self.ProtPath.set('')
        self.ProtName.set('')
        self.LigandPath.set('')
        self.LigandName.set('')

        self.ReferencePath.set('')
        self.defaultOption.set('')
        self.FetchPDB.set('')
        self.AtomTypes.set('Sybyl')
        
        #self.TargetRNA.set(0)
        self.ResSeq.set(0)
        #self.top.Config2.Vars.Anchor.set(-1)

        self.Vars.LigandPathMD5 = ''
    
    
    ''' ==================================================================================
    FUNCTION Load_Session: Actions related to when a new session is loaded
    =================================================================================  '''    
    def Load_Session(self):

        self.Btn_DisplayLigand_Clicked()
        self.Btn_DisplayProtein_Clicked()

    ''' ==================================================================================
    FUNCTION Before_Kill_Frame: Actions related before killing the frame
    =================================================================================  '''    
    def Before_Kill_Frame(self):
        
        self.ProcessError = False
        
        # Process ligand
        if not self.fProcessLigand:
            
            AtomIndex = General.store_Residues(self.top.Config1.listResidues, self.ProtPath.get(), 0)
            if AtomIndex == -1:
                return False

            self.ProcessLigand(True, AtomIndex + 1)

            if self.ProcessError:
                return False

        # Store content of ligand input files
        if not self.fStoreInfo:
            if self.store_InpFile():
                return False

        # Loads the PDB of the processed ligand
        if not self.fLoadProcessed:
            if self.Load_ProcessedLigand():
                return False

        return True
    
    ''' ==================================================================================
    FUNCTION ProcessLigand: Processes ligand PDB file using lig_extractor
    ================================================================================== '''    
    def ProcessLigand(self, boolRun, StartAtomIndex):

        if boolRun:
            self.Disable_Frame()
            
            self.top.ProcessRunning = True
            p = ProcessLigand.ProcLig(self, StartAtomIndex, self.AtomTypes.get(), self.top.Config2.Anchor.get())

        else:
            self.Enable_Frame()

    ''' ==================================================================================
    FUNCTION SmilesRunning: Disables all controls when smiles windows is opened
    ================================================================================== '''    
    def SmilesRunning(self, boolRun, SmilesString):

        if boolRun:
            self.Disable_Frame()
        else:
            print "SmilesString is", SmilesString
            
            self.Enable_Frame()

    ''' ==================================================================================
                         ENABLE / DISABLE - Buttons
    ================================================================================== '''             
    def ValidateLigProt(self,*args):

        if self.ProtName.get() and self.LigandName.get():
            self.top.Go_Step2()
        else:
            self.top.Go_Step1()
        
        self.fProcessLigand = False
        self.top.Reset_Step2()

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
            self.ProtNameTrace = self.ProtName.trace('w',self.ValidateLigProt)
            self.LigandNameTrace = self.LigandName.trace('w',self.ValidateLigProt)
        except:
            pass

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        try:
            self.ProtName.trace_vdelete('w',self.ProtNameTrace)
            self.LigandName.trace_vdelete('w',self.LigandNameTrace)
        except:
            pass
            
    ''' ==================================================================================
    FUNCTION Frame: Generate the Input / Output Files frame in the the middle 
                    frame section    
    ==================================================================================  '''  
    def Frame(self):
        
        self.fIOFile = Frame(self.top.fMiddle)

        #==================================================================================
        #                              PDB Options
        #==================================================================================
        fPDB_options = Frame(self.fIOFile)#, border=1, relief=SUNKEN)
        #fPDB_options.pack(fill=X, side=TOP, padx=5, pady=5)

        fPDB_optionsLine1 = Frame(fPDB_options)#, border=1, relief=SUNKEN)
        fPDB_optionsLine1.pack(side=TOP, fill=X)
        fPDB_optionsLine2 = Frame(fPDB_options)#, border=1, relief=SUNKEN)
        fPDB_optionsLine2.pack(side=TOP, fill=X)

        # Header Get PDB
        Label(fPDB_optionsLine1, text='Retrieve a molecule', font=self.font_Title).pack(side=LEFT)
        
        # Get a PDB File from a file on your harddrive
        Button(fPDB_optionsLine2, text='Open file', command=self.Btn_OpenPDB_Clicked, font=self.font_Text).pack(side=LEFT, padx=5)
        
        # Download a PDB File from the internet
        Button(fPDB_optionsLine2, text='Download', command=self.Btn_DownloadPDB_Clicked, font=self.font_Text, width=10).pack(side=RIGHT, padx=5)
        entFetchPDB = Entry(fPDB_optionsLine2, textvariable=self.FetchPDB, width=10, background='white', font=self.font_Text, justify=CENTER)#,
        #                    validate="key", validatecommand=lambda v=self.FetchPDB: len(v.get()) < 5)
        #entFetchPDB.pack(side=RIGHT)
        Label(fPDB_optionsLine2, text='Enter the PDB code:', font=self.font_Text, justify=CENTER).pack(side=RIGHT, padx=5)
        

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

        # List of selections
        Label(fPDB_options2Line2, text='PyMOL objects/selections:', width=25, justify=RIGHT, font=self.font_Text).pack(side=LEFT, anchor=E)

        optionTuple = ('',)
        self.optionMenuWidget = apply(OptionMenu, (fPDB_options2Line2, self.defaultOption) + optionTuple)
        self.optionMenuWidget.config(bg=self.Color_White, width=15, font=self.font_Text)
        self.optionMenuWidget.pack(side=LEFT)
        
        # Refresh the list with the selections in Pymol        
        Button(fPDB_options2Line2, text='Refresh', command=self.Btn_RefreshOptMenu_Clicked, font=self.font_Text).pack(side=LEFT)

        # List of selections
        #Button(fPDB_options2Line2, text='LIGAND', font=self.font_Text, relief=RIDGE, command=self.Btn_SetLigand_Clicked).pack(side=RIGHT)
        #Button(fPDB_options2Line2, text='TARGET', font=self.font_Text, relief=RIDGE, command=self.Btn_SetProt_Clicked).pack(side=RIGHT)

        #Label(fPDB_options2Line2, text='Set as...', justify=RIGHT, font=self.font_Text).pack(side=RIGHT, anchor=E)

        # List of selections
        Button(fPDB_options2Line2, text='LIGAND', command=self.Btn_SaveLigand_Clicked, font=self.font_Text, relief=RIDGE).pack(side=RIGHT)
        Button(fPDB_options2Line2, text='TARGET', command=self.Btn_SaveProt_Clicked, font=self.font_Text, relief=RIDGE).pack(side=RIGHT)

        Label(fPDB_options2Line2, text='Save as...', justify=RIGHT, font=self.font_Text).pack(side=RIGHT, anchor=E)
        

               
        #==================================================================================
        #                                SET TARGET
        #==================================================================================                
        fPDBsep = Frame(self.fIOFile, border=1, relief=RAISED, width=500, height=3)

        fPDBprotein = Frame(self.fIOFile, border=1, relief=RAISED, width=500, height=70)
        fPDBprotein.pack(side=TOP, pady=10)
        fPDBprotein.pack_propagate(0)

        fPDBproteinLine1 = Frame(fPDBprotein)
        fPDBproteinLine1.pack(side=TOP, fill=X, padx=3, pady=3)

        fPDBproteinLine2 = Frame(fPDBprotein)
        fPDBproteinLine2.pack(side=TOP, fill=X, padx=3, pady=3)

        # First line
        Label(fPDBproteinLine1, width=30, text='THE TARGET', font=self.font_Title).pack(side=LEFT)
        Button(fPDBproteinLine1, text='Load', command=self.Btn_LoadProt_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBproteinLine1, text='Display', command=self.Btn_DisplayProtein_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBproteinLine1, text='Reset', command=self.Btn_ResetProt_Clicked, font=self.font_Text).pack(side=LEFT)

        # Second line
        Label(fPDBproteinLine2, width=30, text='', font=self.font_Title).pack(side=LEFT)
        EntProtein = Entry(fPDBproteinLine2, textvariable=self.ProtName, disabledbackground=self.Color_White, 
                            disabledforeground=self.Color_Black, font=self.font_Text, justify=CENTER)
        EntProtein.pack(side=LEFT, fill=X)
        EntProtein.config(state='disabled')
        #Checkbutton(fPDBproteinLine2, variable=self.TargetRNA, width=10, text='RNA', font=self.font_Text, justify=LEFT).pack(side=LEFT)
        
        #==================================================================================
        #                               SET LIGAND
        #==================================================================================   

        fPDBligand = Frame(self.fIOFile, border=1, relief=RAISED, width=500, height=70)
        fPDBligand.pack(side=TOP, pady=10)
        fPDBligand.pack_propagate(0)

        fPDBligandLine1 = Frame(fPDBligand)
        fPDBligandLine1.pack(side=TOP, fill=X, padx=3, pady=3)

        fPDBligandLine2 = Frame(fPDBligand)
        fPDBligandLine2.pack(side=TOP, fill=X, padx=3, pady=3)

        # First line
        Label(fPDBligandLine1, width=30, text='THE LIGAND', font=self.font_Title).pack(side=LEFT)
        Button(fPDBligandLine1, text='Load', command=self.Btn_LoadLigand_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Display', command=self.Btn_DisplayLigand_Clicked,font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Reset', command=self.Btn_ResetLigand_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Input', command=self.Btn_Input_Clicked, font=self.font_Text).pack(side=LEFT)
        Button(fPDBligandLine1, text='Anchor', command=self.Btn_Anchor_Clicked, font=self.font_Text).pack(side=LEFT)
        
        # Second line
        Label(fPDBligandLine2, width=30, text='', font=self.font_Title).pack(side=LEFT)
        EntLigand = Entry(fPDBligandLine2, disabledbackground=self.Color_White, disabledforeground=self.Color_Black, 
                            textvariable=self.LigandName, font=self.font_Text, justify=CENTER)
        EntLigand.pack(side=LEFT, fill=X)
        EntLigand.config(state='disabled')
        Label(fPDBligandLine2, width=10, text='', font=self.font_Text).pack(side=LEFT)

        #==================================================================================
        #                                 Processing of molecules
        #==================================================================================

        fProcessing = Frame(self.fIOFile)
        fProcessing.pack(side=TOP, fill=X, pady=5, padx=5)

        fProcessingLine1 = Frame(fProcessing)#, border=1, relief=SUNKEN)
        fProcessingLine1.pack(side=TOP, fill=X)
        fProcessingLine2 = Frame(fProcessing)#, border=1, relief=SUNKEN)
        fProcessingLine2.pack(side=TOP, fill=X)
        
        Label(fProcessingLine1, text='Processing of molecules', font=self.font_Title).pack(side=LEFT)

        Label(fProcessingLine2, text='Atom typing:', width=30, justify=RIGHT, font=self.font_Text).pack(side=LEFT, anchor=E)

        Radiobutton(fProcessingLine2, text='Sobolev', variable=self.AtomTypes, value="Sobolev", font=self.font_Text).pack(side=LEFT)
        Radiobutton(fProcessingLine2, text='Gaudreault', variable=self.AtomTypes, value="Gaudreault", font=self.font_Text).pack(side=LEFT, padx=10)
        Radiobutton(fProcessingLine2, text='Sybyl', variable=self.AtomTypes, value="Sybyl", font=self.font_Text).pack(side=LEFT)

        return self.fIOFile
        
    ''' ==================================================================================
    FUNCTIONS Reset Ligand and Protein textbox fields
    ================================================================================== '''
    def Btn_ResetProt_Clicked(self):
    
        self.ProtPath.set('')
        self.ProtName.set('')

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
            
            if self.top.ActiveWizard is None:
                
                self.top.ActiveWizard = Anchor.anchor(self, self.LigandPath.get(), self.top.Config2.Anchor.get())
                
                cmd.set_wizard(self.top.ActiveWizard)
                self.top.ActiveWizard.Start()
                
            else:
                self.DisplayMessage("A wizard is currently active", 2)

    ''' ==================================================================================
    FUNCTION Btn_Input_Clicked: Selects the anchor atom of the ligand
    ==================================================================================  '''    
    def Btn_Input_Clicked(self):

        #if not self.PyMOL:
        #    return

        Smi = Smiles.Smiles(self)
        
    ''' ==================================================================================
    FUNCTION Btn_Save : Save Ligand and Protein objects, object is reloaded automatically
    ==================================================================================  '''    
    def Btn_SaveLigand_Clicked(self):
        
        if not self.PyMOL:
            return

        # Get the Drop Down List Selection Name
        ddlSelection = self.defaultOption.get()

        if ddlSelection == '':
            return

        try:
            state = cmd.get_state()
            n = cmd.count_atoms(ddlSelection, state=state)

            if n < 3:
                self.DisplayMessage("  ERROR for object/selection '" + ddlSelection + "': The ligand must have at least (3) atoms)", 1)
                return

            elif n > Constants.MAX_LIGAND_ATOMS:
                self.DisplayMessage("  ERROR for object/selection '" + ddlSelection + "': The ligand must have a maximum of (" + 
                                    Constants.MAX_LIGAND_ATOMS + ') atoms', 1)
                return

        except:
            self.DisplayMessage("  ERROR: object/selection '" + ddlSelection + "' does not exist on current state", 1)
            return

        LigandPath = tkFileDialog.asksaveasfilename(initialdir=self.top.FlexAIDLigandProject_Dir, title='Save the PDB File', 
                                                    initialfile=ddlSelection, filetypes=[('PDB File','*.pdb')])
        
        if len(LigandPath) > 0:
            
            self.LigandPath.set(os.path.normpath(LigandPath))
            
            if self.LigandPath.get().find('.pdb') == -1:
                self.LigandPath.set(self.LigandPath.get() + '.pdb')

            cmd.save(self.LigandPath.get(), ddlSelection, state, 'pdb') # Save the Selection
            self.LigandName.set(os.path.basename(os.path.splitext(self.LigandPath.get())[0]))
            cmd.load(self.LigandPath.get(), state=1)

            self.DisplayMessage('  Successfully saved and loaded the ligand:  ' + self.LigandName.get() + "'", 0)
    
            

    def Btn_SaveProt_Clicked(self):
        
        if not self.PyMOL:
            return

        # Get the Drop Down List Selection Name
        ddlSelection = self.defaultOption.get()

        if ddlSelection == '':
            return

        try:
            state = cmd.get_state()
            n = cmd.count_atoms(ddlSelection, state=state)

            if n < 3:
                self.DisplayMessage("  ERROR for object/selection '" + ddlSelection + "': The protein must have at least (3) atoms)", 1)
                return
        except:
            self.DisplayMessage(  "The object/selection '" + ddlSelection + "' does not exist on current state", 1)
            return
                    
        ProtPath = tkFileDialog.asksaveasfilename(initialdir=self.top.TargetProject_Dir, title='Save the PDB File', initialfile=ddlSelection, filetypes=[('PDB File','*.pdb')])
        
        if len(ProtPath) > 0:

            self.ProtPath.set(os.path.normpath(ProtPath))
            
            if self.ProtPath.get().find('.pdb') == -1:
                self.ProtPath.set(self.ProtPath.get() + '.pdb')

            cmd.save(self.ProtPath.get(), ddlSelection, state, 'pdb') # Save the Selection
            self.ProtName.set(os.path.basename(os.path.splitext(self.ProtPath.get())[0]))
            cmd.load(self.ProtPath.get(), state=1)
            
            self.DisplayMessage('  Successfully saved and loaded the target: ' + self.ProtName.get(), 0)                                    
                
            

    ''' ==================================================================================
    FUNCTIONS Load Ligand and Protein and display the filename in the textbox
    ================================================================================== '''        
    def Btn_LoadLigand_Clicked(self):        
        
        LigandPath = tkFileDialog.askopenfilename(  filetypes=self.SupportedFormats,
                                                    initialdir=self.top.FlexAIDLigandProject_Dir, title='Select a Ligand File to Load')

        if len(LigandPath) > 0:
            LigandPath = os.path.normpath(LigandPath)

            if LigandPath == self.LigandPath.get():
                return
            
            self.LigandPath.set(LigandPath)

            try:

                Name = os.path.basename(os.path.splitext(self.LigandPath.get())[0])

                if self.PyMOL:
                    cmd.load(self.LigandPath.get(), Name, state=1)
                    n = cmd.count_atoms(Name, state=1)
                else:
                    # For testing purposes only
                    n = 50

                if n < 3:
                    self.DisplayMessage("  ERROR for object '" + Name + "': The ligand must have at least (3) atoms)", 1)
                    return

                elif n > Constants.MAX_LIGAND_ATOMS:
                    self.DisplayMessage("  ERROR for object '" + Name + "': The ligand must have a maximum of (" + 
                                        Constants.MAX_LIGAND_ATOMS + ') atoms', 1)
                    return

            except:
                self.DisplayMessage("  ERROR for file '" + LigandPath + "': Could not load the ligand file", 1)
                return
                
            self.LigandName.set(Name)
            self.DisplayMessage("  Successfully loaded the ligand: '" + self.LigandName.get() + "'", 0)


    def Load_ProcessedLigand(self):
        
        try:
            if self.PyMOL:
                cmd.load(self.ReferencePath.get(), self.LigandName.get(), state=1)
        except:
            self.DisplayMessage('  ERROR: Could not load the PDB of the processed ligand', 1)
            return 1
        
        self.fLoadProcessed = True

        return 0

           
    def Btn_LoadProt_Clicked(self):
        
        ProtPath = tkFileDialog.askopenfilename(filetypes=self.SupportedFormats, 
                                                initialdir=self.top.TargetProject_Dir, title='Select a Target File to Load')
        
        if len(ProtPath) > 0:
            ProtPath = os.path.normpath(ProtPath)

            if ProtPath == self.ProtPath.get():
                return
            
            self.ProtPath.set(ProtPath)
            
            try:
                Name = os.path.basename(os.path.splitext(self.ProtPath.get())[0])
                if self.PyMOL:
                    cmd.load(self.ProtPath.get(), Name, state=1)

            except:
                self.DisplayMessage("  ERROR for object '" + ProtPath + "': Could not load the target file", 1)
                return

            self.ProtName.set(Name)
            self.DisplayMessage("  Successfully loaded the target: '" + self.ProtName.get() + "'", 0)

    ''' ==================================================================================
    FUNCTION openPDB: Import PDB file
    ================================================================================== '''
    def Btn_OpenPDB_Clicked(self):
        
        FilePath = tkFileDialog.askopenfilename(filetypes=[('PDB File','*.pdb')],
                                                initialdir=self.top.Project_Dir, 
                                                title='Select a PDB File to Import')
        
        if len(FilePath) > 0:

            FilePath = os.path.normpath(FilePath)
            
            if self.PyMOL:
                cmd.load(FilePath, state=1)

                
    ''' ==================================================================================
    FUNCTION Btn_DownloadPDB_Clicked: Download a PDB from the internet and display the
                                      result in Pymol 
    ==================================================================================  '''    
    def Btn_DownloadPDB_Clicked(self):       

        PdbCode = self.FetchPDB.get()
        
        try:            
            if self.PyMOL:
                cmd.fetch(PdbCode, async=0)

        except:
            self.DisplayMessage('You entered an invalid pdb code.', 1)
            
        self.FetchPDB.set('')        

    ''' ==================================================================================
    FUNCTION Btn_RefreshOptMenu_Clicked: Refresh the selections list in the application
                                         with the selections in Pymol 
    ==================================================================================  '''                
    def Btn_RefreshOptMenu_Clicked(self): 
        
        if not self.PyMOL:
            return
            
        General_cmd.Refresh_DDL(self.optionMenuWidget, self.defaultOption, [], None)
   
    ''' ==================================================================================
    FUNCTIONS Display Ligand and Protein in pymol
    ================================================================================== '''        
    def Btn_DisplayLigand_Clicked(self):
        
        if not self.PyMOL:
            return

        if self.LigandName.get() != '':

            if not General_cmd.object_Exists(self.LigandName.get()):
                cmd.load(self.LigandPath.get(), state=1)                      # Load the pdb file in Pymol                     
            else:
                cmd.center(self.LigandName.get())
                cmd.zoom(self.LigandName.get())    
    
        
    def Btn_DisplayProtein_Clicked(self):
        
        if not self.PyMOL:
            return

        if self.ProtName.get() != '':

            if not General_cmd.object_Exists(self.ProtName.get()):
                cmd.load(self.ProtPath.get(), state=1)                        # Load the pdb file in Pymol                     
            else:
                cmd.center(self.ProtName.get())
                cmd.zoom(self.ProtName.get())        
            


    #=======================================================================
    ''' Store inp file information (flexible bonds, atom types)  '''
    #=======================================================================   
    def store_InpFile(self):
                
        inpFilePath = os.path.join(self.top.FlexAIDSimulationProject_Dir,'LIG.inp')

        inpInfo = dict()
        flexInfo = dict()

        #Read the inp file and get the flexible bonds
        try:
            file = open(inpFilePath, 'r')
            inpFile = file.readlines()
            file.close()

            nbLines = len(inpFile)      
            for line in range(0, nbLines):

                #HETTYP  902 1  N1  m   909  910  903    0
                if inpFile[line].startswith('HETTYP'):
                    
                    ATOM = inpFile[line][7:11].strip()

                    list = []
                    list.append(inpFile[line][11:13].strip()) # Type
                    list.append(inpFile[line][22:26].strip()) # Neighbour 1
                    list.append(inpFile[line][27:31].strip()) #           2
                    list.append(inpFile[line][32:36].strip()) #           3
                    inpInfo[ATOM] = list

                #FLEDIH  4   916  917   
                elif inpFile[line].startswith('FLEDIH'):

                    INDEX = inpFile[line][7:9].strip()

                    list = []
                    for i in range(0,len(inpFile[line][10:])/5):
                        list.append(inpFile[line][(10+i*5):(10+5+i*5)].strip())
                    flexInfo[INDEX] = list
                
        except:
            self.DisplayMessage('  ERROR: Could not retrieve ligand input file', 1)
            return 1

        if not self.Check_LigandPathMD5(inpFilePath):
            self.store_Neighbours(inpInfo)
            self.store_AtomTypes(inpInfo)
            self.store_FlexBonds(flexInfo)
        
        self.Vars.LigandPathMD5 = self._LigandPathMD5
        self.fStoreInfo = True

        return 0

    #=======================================================================
    ''' Store Neighbours Dictionary'''
    #=======================================================================   
    def store_Neighbours(self, inpInfo):
        
        self.top.Config2.Vars.dictNeighbours.clear()
        
        for atom in inpInfo.keys():
            self.top.Config2.Vars.dictNeighbours[atom] = inpInfo[atom][1:]
                    
    #=======================================================================
    ''' Store Flexible Bonds Dictionary'''
    #=======================================================================   
    def store_FlexBonds(self, flexInfo):
        
        self.top.Config2.Vars.dictFlexBonds.clear()
        
        for index in flexInfo.keys():
            
            ''' [ Selected as flexible,
                  Forced as flexible,
                  Number of atoms defining the bond,
                  Atom list defining bond ] '''

            dictList = [ 0, 0, len(flexInfo[index]) ]
            for i in range(0, len(flexInfo[index])):
                dictList.append(flexInfo[index][i])
                
            self.top.Config2.Vars.dictFlexBonds[index] = dictList

    #=======================================================================
    ''' Store Atom Types Dictionary'''
    #=======================================================================   
    def store_AtomTypes(self, inpInfo):
        
        self.top.Config2.Vars.dictAtomTypes.clear()
        
        for atom in inpInfo.keys():
            self.top.Config2.Vars.dictAtomTypes[atom] = [inpInfo[atom][0], inpInfo[atom][0]]

    
    #=======================================================================
    ''' Compares the file content of the ligand when the session was saved '''
    #=======================================================================   
    def Check_LigandPathMD5(self, fname):
        
        self._LigandPathMD5 = self.hashfile(open(fname, 'r'), hashlib.md5())
        
        return self.Vars.LigandPathMD5 == self._LigandPathMD5

    #=======================================================================
    ''' Returns the Signature of a file contents '''
    #=======================================================================   
    def hashfile(self, afile, hasher, blocksize=65536):
    
        buf = afile.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = afile.read(blocksize)
        
        return hasher.digest()

    #=======================================================================
    ''' Enables/disables controls related to AnchorWizard  '''
    #=======================================================================   
    def AnchorRunning(self, boolRun):
        
        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()

            if self.top.Config2.Anchor.get() != self.top.WizardResult:
                self.fProcessLigand = False

            self.top.Config2.Anchor.set(self.top.WizardResult)


    # Welcome menu message
    def Load_Message(self):

        self.DisplayMessage('' ,0)
        self.DisplayMessage('  FlexAID < Input Files > Menu', 2)
        self.DisplayMessage('  INFO:   Select a < TARGET > and a < LIGAND > by:', 2)
        self.DisplayMessage('          1) Saving a PyMOL object/selection to your project directory', 2)
        self.DisplayMessage('          2) Loading an existing object file from your project directory', 2)
