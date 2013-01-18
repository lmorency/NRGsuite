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

'''
@title: GetCleft - Interface - Default tab

@summary: This is the interface of GetCleft application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Oct. 19, 2010
'''

from Tkinter import *

from subprocess import Popen, PIPE
from time import time

import os
import tkFileDialog
import General
import CleftObj
import BindingSite

import functools
import threading
import Color
import pickle

if __debug__:
    from pymol import cmd
    import General_cmd

#=========================================================================================
'''                        ---   STARTING GETCLEFT  ---                                '''
#=========================================================================================
class RunThread(threading.Thread):

    def __init__(self, top, selection, cmdline):
                 
        threading.Thread.__init__(self)

        self.top = top

        self.Selection = selection
        self.GetCleft = self.top.top

        self.cmdline = cmdline

        self.StateList = list()
        self.top.GetCleftRunning(True)

        # Start the thread
        self.start()


    def run(self): 
       
        self.timeBegin = time()

        self.top.ProcessRunning = True
        if self.GetCleft.OSid == 'WIN':
            self.GetCleft.Run = Popen(self.cmdline, shell=False, stderr=PIPE)
        else:
            self.GetCleft.Run = Popen(self.cmdline, shell=True, stderr=PIPE)
        self.GetCleft.Run.wait()

        
        if self.GetCleft.Run.returncode != 0:

            self.top.DisplayMessage("  ERROR: An error occured while executing GetCleft:", 1)
            self.top.DisplayMessage(self.GetCleft.Run.stderr.read(), 1)
            
        else:

            duration = str(time() - self.timeBegin)
            duration = '%.2f' % float(duration)
            
            self.top.DisplayMessage('  Process completed in ' + str(duration) + ' seconds ...', 2)

        self.GetCleft.Run = None

        # Store clefts and show them
        self.GetCleft.Manage.delete_Temp()
        #self.GetCleft.Manage.rename_Temp()
        self.GetCleft.Manage.store_Temp(self.top.LastdefaultOption)
            
        nCleft = self.top.TempBindingSite.Count_Cleft()
        if nCleft > 0:
            self.top.DisplayMessage("  Stored (" + str(nCleft) + 
                                    ") cleft objects from object/selection '" +
                                    self.Selection + "'", 0)

            self.top.Display_Temp()
            self.GetCleft.Go_Step2()
        else:
            self.top.DisplayMessage("  No clefts found for object/selection '" +
                                    self.Selection + "'", 0)
            
        self.top.GetCleftRunning(False)
        
#=========================================================================================
'''                        ---  GETCLEFT's DEFAULT FRAME  ---                          '''
#=========================================================================================     
class Default:

    def __init__(self,top, PyMOL):

        #print "New instance of Default"

        self.PyMOL = PyMOL

        self.top = top
        self.FrameName = 'Default'
        self.Tab = self.top.Btn_Config

        self.font_Text = self.top.font_Text
        self.font_Title = self.top.font_Title

        self.Color_Black = self.top.Color_Black
        self.Color_Blue = self.top.Color_Blue
        self.Color_White = self.top.Color_White

        self.DisplayMessage = self.top.DisplayMessage

        self.Def_Vars()
        self.Init_Vars()
        
        self.Frame()
        self.Trace()


    def Def_Vars(self):

        self.TempBindingSite = BindingSite.BindingSite()

        self.NbCleft = StringVar()
        self.MinRadius = DoubleVar()
        self.MaxRadius = DoubleVar()
        self.MessageValue = StringVar()

        self.FetchPDB = StringVar()
        self.ResiduValue = StringVar()
        self.OutputFileValue = StringVar()
                        
        self.defaultOption = StringVar()
        self.listResidues = list()

        self.ColorList = list()
        self.ColorRGB = list()

        self.ValidNbCleft = list()
        self.ValidOutput = list()
        self.Validator = list()


    def Init_Vars(self):

        self.FetchPDB.set('')
        self.NbCleft.set('10')
        self.MinRadius.set('1.50')
        self.MaxRadius.set('4.00')
        self.OutputFileValue.set('')
        self.ResiduValue.set('')
        self.MessageValue.set('')
        self.defaultOption.set('')

        self.listResidues = []
        self.PartitionColor = 'partition'
        self.LastdefaultOption = ''
        
        self.ValidNbCleft = [ True, 1, 0, None ]
        self.ValidMinRadius = [ True, 1, 0, None ]
        self.ValidMaxRadius = [ True, 1, 0, None ]
        self.ValidOutput = [ True, 1, 0, None ]

        self.Validator = [ self.ValidNbCleft, self.ValidMinRadius, self.ValidMaxRadius, self.ValidOutput ]
        
    
    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
                
        self.fDefault.pack_forget()
        #self.fDefault.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        #self.intChkAtomsTrace = self.intChkAtoms.trace('w',self.Check_Atoms)
        #self.intChkCleftsTrace = self.intChkClefts.trace('w',self.Check_Spheres)
        #self.RadioCLFTrace = self.RadioCLF.trace('w',self.Toggle_CleftView)

        self.defaultOptionTrace = self.defaultOption.trace('w',self.Toggle_defaultOption)

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        #self.intChkAtoms.trace_vdelete('w',self.intChkAtomsTrace)
        #self.intChkClefts.trace_vdelete('w',self.intChkCleftsTrace)
        #self.RadioCLF.trace_vdelete('w',self.RadioCLFTrace)

        self.defaultOption.trace_vdelete('w',self.defaultOptionTrace)

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
        
        self.fDefault.pack(fill=BOTH, expand=True)
        self.LoadMessage()
                
        self.Btn_RefreshOptMenu_Clicked()
        
        if self.LastdefaultOption != '' and General_cmd.object_Exists(self.LastdefaultOption):
            self.defaultOption.set(self.LastdefaultOption)
        

    ''' ==================================================================================
    FUNCTION Frame: Displays the Default Frame
    ==================================================================================  '''          
    def Frame(self):
                
        self.fDefault = Frame(self.top.fMiddle, relief=RIDGE)

        #==================================================================================
        #                           SELECTION OF STRUCTURE
        #==================================================================================                
        fStructure = Frame(self.fDefault)
        fStructure.pack(side=TOP, padx=5, pady=5, fill=X, expand=True)
        fStructureLine1 = Frame(fStructure)
        fStructureLine1.pack(side=TOP, fill=X, expand=True)
        fStructureLine2 = Frame(fStructure)
        fStructureLine2.pack(side=TOP, fill=X, expand=True)

        Label(fStructureLine1, text='Retrieve a structure', font=self.top.font_Title).pack(side=LEFT)

        # Get a PDB File from a file on your harddrive
        Button(fStructureLine2, text='Open target', command=self.Btn_OpenPDB_Clicked, font=self.font_Text).pack(side=LEFT, padx=5)
        
        # Download a PDB File from the internet
        Button(fStructureLine2, text='Download', command=self.Btn_DownloadPDB_Clicked, font=self.font_Text, width=10).pack(side=RIGHT, padx=5)
        Entry(fStructureLine2, textvariable=self.FetchPDB, width=10, background='white', font=self.font_Text, justify=CENTER).pack(side=RIGHT)
        Label(fStructureLine2, text='Enter the PDB code:', font=self.font_Text, justify=CENTER).pack(side=RIGHT, padx=5)

        fSelection = Frame(self.fDefault)
        fSelection.pack(side=TOP, padx=5, pady=5, fill=X, expand=True)
        fSelectionLine1 = Frame(fSelection)
        fSelectionLine1.pack(side=TOP, fill=X, expand=True)
        fSelectionLine2 = Frame(fSelection)
        fSelectionLine2.pack(side=TOP, fill=X, expand=True)

        Label(fSelectionLine1, text='Select a structure', font=self.top.font_Title).pack(side=LEFT)

        # List of selections
        Label(fSelectionLine2, text='PyMOL objects/selections:', font=self.font_Text, justify=LEFT).pack(side=LEFT, padx=5)
        
        # Refresh the list with the selections in Pymol
        Button(fSelectionLine2, text='Refresh', command=self.Btn_RefreshOptMenu_Clicked, font=self.font_Text, width=10).pack(side=RIGHT, padx=5)

        optionTuple = ('',)
        self.optionMenuWidget = apply(OptionMenu, (fSelectionLine2, self.defaultOption) + optionTuple)
        self.optionMenuWidget.config(bg=self.Color_White, font=self.font_Text, width=15)
        self.optionMenuWidget.pack(side=RIGHT)
                
        #==================================================================================
        '''                     --- Basic GetCleft Options ---                          '''
        #==================================================================================
        #-----------------------------------------------------------------------------------------        
        fBasic = Frame(self.fDefault, border=0, relief=RAISED)
        fBasic.pack(fill=X, expand=True, padx=5, pady=5, ipady=3)
        fBasicLine1 = Frame(fBasic)
        fBasicLine1.pack(fill=X, side=TOP)
        fBasicLine2 = Frame(fBasic)
        fBasicLine2.pack(fill=X, side=TOP)
        fBasicLine3 = Frame(fBasic)
        fBasicLine3.pack(fill=X, side=TOP)
        fBasicLine4 = Frame(fBasic)
        fBasicLine4.pack(fill=X, side=TOP)
        fBasicLine5 = Frame(fBasic)
        fBasicLine5.pack(fill=X, side=TOP)

        Label(fBasicLine1, text='Parametrization', font=self.top.font_Title).pack(side=LEFT)

        Label(fBasicLine2, text='Insert spheres radii:', font=self.top.font_Text).pack(side=LEFT)

        EntryMaxRadius = Entry(fBasicLine2, width=5, textvariable=self.MaxRadius, background='white', justify=CENTER, font=self.top.font_Text)
        EntryMaxRadius.pack(side=RIGHT)
        Label(fBasicLine2, text='Max:', font=self.top.font_Text).pack(side=RIGHT, padx=2)

        EntryMinRadius = Entry(fBasicLine2, width=5, textvariable=self.MinRadius, background='white', justify=CENTER, font=self.top.font_Text)
        EntryMinRadius.pack(side=RIGHT)
        Label(fBasicLine2, text='Min:', font=self.top.font_Text).pack(side=RIGHT, padx=2)

        args_list = [EntryMaxRadius, self.MaxRadius, EntryMinRadius, 10.00, 2, self.ValidMaxRadius, 'Maximum radius', 'float']
        EntryMaxRadius.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidMaxRadius[3] = EntryMaxRadius

        args_list = [EntryMinRadius, self.MinRadius, 0.25, EntryMaxRadius, 2, self.ValidMinRadius, 'Minimum radius', 'float']
        EntryMinRadius.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidMinRadius[3] = EntryMinRadius
        
        Label(fBasicLine3, text='Residue in contact (e.g. ALA13A):', font=self.top.font_Text).pack(side=LEFT)
        self.EntryResidu = Entry(fBasicLine3,textvariable=self.ResiduValue, background='white', justify=CENTER, font=self.top.font_Text)
        self.EntryResidu.pack(side=RIGHT)

        Label(fBasicLine4, text='Number of clefts to show:', font=self.top.font_Text).pack(side=LEFT)
        EntryNbCleft = Entry(fBasicLine4, width=8, textvariable=self.NbCleft, background='white', justify=CENTER, font=self.top.font_Text)
        EntryNbCleft.pack(side=RIGHT)
        args_list = [EntryNbCleft, self.NbCleft, 1, 100, -1, self.ValidNbCleft, 'Number of clefts', 'int']
        EntryNbCleft.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidNbCleft[3] = EntryNbCleft

        Label(fBasicLine5, text='Output cleft basename:', font=self.top.font_Text).pack(side=LEFT)
        EntryOutput = Entry(fBasicLine5, textvariable=self.OutputFileValue, background='white', font=self.top.font_Text, justify=CENTER)
        EntryOutput.pack(side=RIGHT)
        args_list = [EntryOutput, self.OutputFileValue, -1, -1, -1, self.ValidOutput, 'Output filename', 'str']
        EntryOutput.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidOutput[3] = EntryOutput

        #Button(fBasicLine5, text='Advanced parameters', font=self.top.font_Text, width=20, relief=RIDGE, command=self.Btn_AdvancedOptions).pack(side=RIGHT)

        #==================================================================================
        '''                           --- BUTTONS AREA ---                              '''
        #==================================================================================
        fButtons = Frame(self.fDefault)
        fButtons.pack(fill=X, expand=True, padx=5, pady=5)
                                                  
        # Starts the GetCleft application
        self.Btn_StartGetCleft = Button(fButtons, text='Start', command=self.Btn_StartGetCleft_Clicked, font=self.top.font_Text)
        self.Btn_StartGetCleft.pack(side=LEFT)
        self.Btn_StartGetCleft.config(state='disabled')

        # Clear PyMOL elements
        Btn_Clear = Button(fButtons, text='Clear', command=self.Btn_Clear_Clicked, font=self.top.font_Text)
        Btn_Clear.pack(side=LEFT)

        #==================================================================================
        '''                         --- COLOR CHART AREA ---                            '''
        #==================================================================================
        #-----------------------------------------------------------------------------------------        
        fChart = Frame(self.fDefault)
        fChart.pack(fill=X, expand=True, padx=5, side=BOTTOM)
        fChartLine1 = Frame(fChart)
        fChartLine1.pack(fill=X, expand=True)
        fChartLine2 = Frame(fChart)
        fChartLine2.pack(fill=X, expand=True)
        fChartLine3 = Frame(fChart)
        fChartLine3.pack(fill=X, expand=True)

        #Label(fChartLine1, text='Clefts color chart', font=self.top.font_Title).pack(side=TOP)
        #Label(fChartLine2, text='Index', width=10, font=self.top.font_Text).pack(side=LEFT)        

        fSim_Prog = Frame(fChartLine2, border=1, relief=SUNKEN, width=400, height=25)
        fSim_Prog.pack(side=TOP)#, anchor=W) 
        self.ChartBar = Canvas(fSim_Prog, bg=self.top.Color_Grey, width=400, height=25, highlightthickness=0, relief='flat', bd=0)
        self.ChartBar.pack(fill=BOTH, anchor=W)

        #==================================================================================
        '''        --- The DISPLAY options for the CLEFT and SPHERES ---                '''
        #==================================================================================
        fDisplay = Frame(self.fDefault)
        fDisplay.pack(side=TOP, fill=X, padx=5, pady=5)
        fDisplayLine1 = Frame(fDisplay)
        fDisplayLine1.pack(side=TOP, fill=X)
        fDisplayLine2 = Frame(fDisplay)
        fDisplayLine2.pack(side=TOP, fill=X)
        fDisplayLine3 = Frame(fDisplay)
        fDisplayLine3.pack(side=TOP, fill=X)
        
        #Label(fDisplayLine1, text='Display options', font=self.top.font_Title).pack(side=LEFT)
        #Checkbutton(fDisplayLine2, text='Atoms', width=20, variable=self.intChkAtoms, font=self.top.font_Text, justify=LEFT).pack(side=LEFT) 
        #Radiobutton(fDisplayLine2, text='Spheres', variable=self.RadioCLF, value='sphere', font=self.top.font_Text).pack(side=LEFT)
        #Radiobutton(fDisplayLine2, text='Surface', variable=self.RadioCLF, value='surface', font=self.top.font_Text).pack(side=LEFT)
        #Checkbutton(fDisplayLine3, text='Clefts', width=20, variable=self.intChkClefts, font=self.top.font_Text, justify=LEFT).pack(side=LEFT)
       
    ''' ==================================================================================
    FUNCTION Btn_Clear_Clicked: Clear the temp dir of clefts and delete associated clefts 
    ==================================================================================  '''
    def Btn_Clear_Clicked(self):
        
        # Delete Cleft/Sphere objects in PyMOL
        for Cleft in iter(self.TempBindingSite.listClefts):
            if self.PyMOL:
                try:
                    cmd.delete(Cleft.CleftName)
                except:
                    pass
        
        self.TempBindingSite.Clear_Cleft()

        # Clean temporary files
        self.top.Manage.Clean()

        self.top.Go_Step1()
        self.ChartBar.delete('all')
        
    ''' ==================================================================================
    FUNCTION Btn_StartGetCleft_Clicked: Run GetCleft and display the result in Pymol 
    ==================================================================================  '''
    def Btn_StartGetCleft_Clicked(self):
        
        if not self.top.ProcessRunning:
            
            # Trigger lost_focus event for validation
            self.top.fMiddle.focus_set()
            self.top.fMiddle.update_idletasks()

            rv = self.top.Validate_Entries(self.Validator)
            if rv > 0:
                self.DisplayMessage("  Cannot run GetCleft: Not all fields are validated", 2)
                return
            
            try:
                TmpFile = os.path.join(self.top.GetCleftProject_Dir,'tmp.pdb')

                if self.PyMOL:
                    if cmd.count_atoms(self.defaultOption.get()) == 0:
                        self.DisplayMessage("  ERROR: No atoms found for object/selection '" + self.defaultOption.get() + "'", 2)
                        return

                    cmd.save(TmpFile, self.defaultOption.get())

            except:
                self.DisplayMessage("  ERROR: Could not save the object/selection '" + self.defaultOption.get() + "'", 2)
                return

            if self.ResiduValue.get() != '':
                
                # with HETATM groups
                if General.store_Residues(self.listResidues, TmpFile, 1) == -1:
                    self.DisplayMessage("  ERROR: Could not retrieve list of residues for object/selection '" + self.defaultOption.get() + "'", 2)
                    return            
                
                if self.listResidues.count(self.ResiduValue.get()) == 0:
                    self.DisplayMessage("  ERROR: The residue entered could not be found in the object/selection '" + self.defaultOption.get() + "'", 2)
                    self.EntryResidu.config(bg=self.top.Color_Red)
                    return

                self.EntryResidu.config(bg=self.top.Color_White)
            
            self.DisplayMessage("  Analyzing clefts for object/selection '" + self.defaultOption.get() + "'...", 0)
            
            #TmpFile = '/Users/francisgaudreault/1stp.pdb'
            Command_Line = '"' + self.top.GetCleftExecutable + '" -p "' + TmpFile + '"' + self.Get_Arguments()
            
            # Clear temporary clefts
            self.Btn_Clear_Clicked()
            
            # Run GetCleft
            StartRun = RunThread(self, self.defaultOption.get(), Command_Line)
            

    ''' ========================================================
                 Display all temporary clefts
        ========================================================'''
    def Display_Temp(self):

        self.SetColorList()
        self.DisplayColorChart()

        View = cmd.get_view()
        
        self.Load_Clefts()
        self.Show_Clefts(self.ColorList)

        # reset view
        cmd.set_view(View)
            
    ''' ========================================================
                 Gets all arguments for the cmdline
        ========================================================'''
    def Get_Arguments(self):

        Args = ''
        
        # Centralized on a residue
        if self.ResiduValue.get() != '':
            Args += ' -a ' + self.ResiduValue.get() + '-'
            
        # Output location
        OutputPath = self.top.GetCleftTempProject_Dir
        if self.OutputFileValue.get() != '':
            OutputPath = os.path.join(OutputPath,self.OutputFileValue.get())
        else:
            OutputPath = os.path.join(OutputPath,self.defaultOption.get())

        Args += ' -o "' + OutputPath + '"'

        # Number of clefts maximum
        Args += ' -t ' + self.NbCleft.get()          
            
        # Size of spheres
        Args += ' -l ' + str(self.MinRadius.get())
        Args += ' -u ' + str(self.MaxRadius.get())
        
        """    
        # Misc.
        Args += ' -k ' + self.top.AdvOptions.Entry_K.get()

        # Booleans
        if self.top.AdvOptions.Entry_R_Value.get():
            Args += ' -r'
                
        if self.top.AdvOptions.Entry_h_Value.get():
            Args += ' -h'
                
        if self.top.AdvOptions.Entry_H_Value.get():
            Args += ' -H'
                
        if self.top.AdvOptions.Entry_C.get() != 'ALL':
            Args += ' -c ' + self.Entry_C.get()
        
        """
        #Args += ' -ca'
        
        Args += ' -s'

        self.LastdefaultOption = self.defaultOption.get()

        return Args

    ''' ========================================================
                  Toggles the state of the Start button
        ========================================================'''
    def Toggle_defaultOption(self, *args):

        if self.defaultOption.get() != '':
            self.Btn_StartGetCleft.config(state='normal')
        else:
            self.Btn_StartGetCleft.config(state='disabled')

    ''' ==================================================================================
    FUNCTION Btn_RefreshOptMenu_Clicked: Refresh the selections list in the application
                                         with the selections in Pymol 
    ==================================================================================  '''                
    def Btn_RefreshOptMenu_Clicked(self):
        
        if self.PyMOL:
            exc = []
            General_cmd.Refresh_DDL(self.optionMenuWidget, self.defaultOption, exc, None)

    """
    ''' ==================================================================================
    FUNCTION Btn_AdvancedOptions: Opens the Advanced menu of GetCleft
    ==================================================================================  '''                
    def Btn_AdvancedOptions(self):

        self.top.SetActiveFrame(self.top.AdvOptions)
       
    """
    
    ''' ========================================================
                  Welcome message upon frame built
    ========================================================='''
    def LoadMessage(self):

        self.DisplayMessage('', 0)
        self.DisplayMessage('  Opened the default menu... ',0)        


    # Disable/enables the whole frames
    def GetCleftRunning(self, boolRun):
        
        if boolRun:
            self.StateList = []
            General.saveState(self.fDefault, self.StateList)
            General.setState(self.fDefault)

        else:
            General.backState(self.fDefault, self.StateList)

    ''' ==================================================================================
    FUNCTION DisplayColorChart: Display the color chart in the GetCleft application 
    ================================================================================== '''   
    def DisplayColorChart(self):
        
        nClefts = self.TempBindingSite.Count_Cleft()
        
        self.ChartBar.delete('all')
        
        if nClefts > Color.NBCOLOR:
            nClefts = Color.NBCOLOR
        
        CellWidth = int(self.ChartBar.cget('width')) / nClefts
        
        LeftCoord = 0
        RightCoord = 0

        for i in range(0, nClefts):
            
            if i == (nClefts-1):
                RightCoord = self.ChartBar.cget('width')
            else:    
                RightCoord = LeftCoord + CellWidth
                
            self.ChartBar.create_rectangle(LeftCoord, 0,
                                           RightCoord, self.ChartBar.cget('height'),
                                           fill=self.ColorRGB[i])
            LeftCoord = RightCoord

    ''' ==================================================================================
    FUNCTION Get_CleftPath: Retrieves the default path of the clefts
    ==================================================================================  '''        
    def Get_CleftPath(self):
                
        if self.LastdefaultOption != '':
            TARGETNAME = self.LastdefaultOption.upper()
        else:
            TARGETNAME = self.defaultOption.get().upper()

        CleftPath = os.path.join(self.top.CleftProject_Dir,TARGETNAME)
        
        return CleftPath

    ''' ==================================================================================
    FUNCTION Get_BindingSitePath: Retrieves the default path of the bindingsite
    ==================================================================================  '''        
    def Get_BindingSitePath(self):
        
        if self.LastdefaultOption != '':
            TARGETNAME = self.LastdefaultOption.upper()
        else:
            TARGETNAME = self.defaultOption.get().upper()

        BindingSitePath = os.path.join(self.top.BindingSiteProject_Dir,TARGETNAME)
        
        return BindingSitePath

    ''' ==================================================================================
    FUNCTION Btn_Load_Clefts: Asks for user to load clefts
    ==================================================================================  '''        
    def Btn_Load_Clefts(self):

        CleftPath = self.Get_CleftPath()
        
        if not os.path.isdir(CleftPath):
            self.DisplayMessage("  Could not find a Cleft folder for your target:", 2)
            self.DisplayMessage("  The default Cleft folder is used.", 2)
        
            CleftPath = self.top.CleftProject_Dir

        
        LoadFiles = tkFileDialog.askopenfilename(filetypes=[('Cleft file','*.nrgclf')],
                                                 initialdir=CleftPath, title='Select cleft file(s) to load',
                                                 multiple=1)
        
        if len(LoadFiles) > 0:
        
            TempBindingSite = BindingSite.BindingSite()
            
            for LoadFile in iter(LoadFiles):

                LoadFile = os.path.normpath(LoadFile)
                
                try:
                    in_ = open(LoadFile, 'r')
                    Cleft = pickle.load(in_)
                    in_.close()
                    
                    TempBindingSite.Add_Cleft(Cleft)
                except:
                    self.top.DisplayMessage("  ERROR: Could not the cleft '" + LoadFile + "'", 1)
                    pass
                
            if TempBindingSite.Count_Cleft() > 0:
                self.Btn_Clear_Clicked()
                self.TempBindingSite = TempBindingSite
                
                self.Display_Temp()
                self.top.Go_Step2()

    ''' ==================================================================================
    FUNCTION Btn_Save_Clefts: Asks for user to save clefts
    ==================================================================================  '''        
    def Btn_Save_Clefts(self):

        if not self.PyMOL:
            return
            
        if self.TempBindingSite.Count_Cleft() > 0:
                        
            self.Update_TempBindingSite()
            
            self.top.Manage.save_Temp()
            
            for Cleft in self.TempBindingSite.listClefts:
            
                CleftPath = os.path.join(self.top.CleftProject_Dir,Cleft.UTarget)
                
                if not os.path.isdir(CleftPath):
                    os.makedirs(CleftPath)

                CleftSaveFile = os.path.join(CleftPath,os.path.basename(os.path.splitext(Cleft.CleftFile)[0]) + '.nrgclf')

                try:
                    out = open(CleftSaveFile, 'w')
                    pickle.dump(Cleft, out)
                    out.close()
                    
                    self.top.DisplayMessage("  Successfully saved '" + CleftSaveFile + "'", 0)
                except:
                    self.top.DisplayMessage("  ERROR: Could not save binding-site configuration", 1)
        
        else:
            self.top.DisplayMessage("  No clefts to save as 'clefts'", 2)

    ''' ==================================================================================
    FUNCTION Btn_Load_BindingSite: Asks for user to load binding-site
    ==================================================================================  '''        
    def Btn_Load_BindingSite(self):

        BindingSitePath = self.Get_BindingSitePath()

        if not os.path.isdir(BindingSitePath):
            self.DisplayMessage("  Could not find a BindingSite folder for your target:", 2)
            self.DisplayMessage("  The default BindingSite folder is used.", 2)
            
            BindingSitePath = self.top.BindingSiteProject_Dir
        
        LoadPath = tkFileDialog.askopenfilename(filetypes=[('NRG BindingSite','*.nrgbs')],
                                                initialdir=BindingSitePath, title='Select a NRG BindingSite file to load')
        
        if len(LoadPath) > 0:
            LoadPath = os.path.normpath(LoadPath)
            
            try:
                in_ = open(LoadPath, 'r')
                TempBindingSite = pickle.load(in_)
                in_.close()
                
                if TempBindingSite.Count_Cleft() > 0:
                    self.Btn_Clear_Clicked()
                    self.TempBindingSite = TempBindingSite

                    self.Display_Temp()
                    self.top.Go_Step2()
                    
            except:
                self.DisplayMessage("  ERROR: Could not find read the BindingSite", 2)


    ''' ==================================================================================
    FUNCTION Btn_Save_BindingSite: Asks for user to save binding-site
    ==================================================================================  '''        
    def Btn_Save_BindingSite(self):

        BindingSitePath = self.Get_BindingSitePath()

        if self.TempBindingSite.Count_Cleft() > 0:
            
            if not os.path.isdir(BindingSitePath):
                os.makedirs(BindingSitePath)

            SaveFile = tkFileDialog.asksaveasfilename(initialdir=BindingSitePath,
                                                      title='Save the BindingSite file', initialfile='def_cleft_bindingsite',
                                                      filetypes=[('NRG BindingSite','*.nrgbs')])
            
            if len(SaveFile) > 0:
                SaveFile = os.path.norm(SaveFile)
                
                if SaveFile.find('.nrgbs') == -1:
                    SaveFile = SaveFile + '.nrgbs'

                self.Update_TempBindingSite()
                self.TempBindingSite.Type = 2
                
                # Save all cleft 
                self.top.Manage.save_Temp()

                try:
                    out = open(SaveFile, 'w')
                    pickle.dump(self.TempBindingSite, out)
                    out.close()
                    
                    self.top.DisplayMessage("  Successfully saved '" + SaveFile + "'", 0)
                except:
                    self.top.DisplayMessage("  ERROR: Could not save binding-site configuration", 1)

        else:
            self.top.DisplayMessage("  No clefts to save as 'binding-site'", 2)
    
    ''' ==================================================================================
    FUNCTION Load_Clefts: Loads the list of temp clefts
    ==================================================================================  '''        
    def Load_Clefts(self):

        if not self.PyMOL:
            return
            
        for Cleft in iter(self.TempBindingSite.listClefts):
                cmd.load(Cleft.CleftFile, Cleft.CleftName, state=1)
                
                if Cleft.Partition and Cleft.PartitionParent != None and \
                        General_cmd.object_Exists(Cleft.PartitionParent.CleftName):
                    
                    General_cmd.Oscillate(Cleft.PartitionParent.CleftName, 0.0)

    ''' ==================================================================================
    FUNCTION Update_TempBindingSite: Update cleft object with only those undeleted in PyMOL
    ==================================================================================  '''                 
    def Update_TempBindingSite(self):
        
        if not self.PyMOL:
            return

        for Cleft in iter(self.TempBindingSite.listClefts):
            if not General_cmd.object_Exists(Cleft.CleftName):
                self.TempBindingSite.Remove_Cleft(Cleft)

    ''' ==================================================================================
    FUNCTION Btn_OpenPDB_Clicked: Open a PDB file and display it in PyMOL 
    ==================================================================================  '''                 
    def Btn_OpenPDB_Clicked(self):
        
        FilePath = tkFileDialog.askopenfilename(filetypes=[('PDB File','*.pdb'),('MOL File','*.mol'),('MOL2 File','*.mol2'),('SDF File','*.sdf')], initialdir=self.top.TargetProject_Dir, title='Select a PDB File to Load')
        
        if len(FilePath) > 0:
            FilePath = os.path.normpath(FilePath)
            
            Name = os.path.basename(os.path.splitext(FilePath)[0])
        
            if self.PyMOL:
                # Load the pdb file in Pymol
                cmd.load(FilePath, state=1)
                self.defaultOption.set(Name)
                
    ''' ==================================================================================
    FUNCTION SetColorList: Reset the color lists
    ==================================================================================  '''        
    def SetColorList(self):
        
        nColor = self.TempBindingSite.Count_Cleft()
        
        if nColor > Color.NBCOLOR:
            nColor = Color.NBCOLOR
            
        self.ColorList = Color.GetHeatColorList(nColor, False)
        self.ColorRGB = Color.GetHeatColorList(nColor, True)
                
    ''' ==========================================================
    SetPartitionColor: sets a new color to the partition based on its parent
    ==========================================================='''           
    def SetPartitionColor(self, Sel):
        
        if not self.PyMOL:
            return

        try:
            mycolors = {'colors': []}
            cmd.iterate( Sel, 'colors.append(color)', space=mycolors)

            for color in mycolors['colors']:
                lRGB = list( cmd.get_color_tuple(color) )
                break

            for i in range(0,3):
                if lRGB[i] <= 0.80:
                    lRGB[i] += 0.20
                else:
                    lRGB[i] -= 0.20
         
            cmd.set_color(self.PartitionColor, lRGB)
            
        except:
            return 1

        return 0
    
    ''' ==================================================================================
    FUNCTION Show_Clefts: Shows the 'Cleft' object a cavity
    ================================================================================== '''   
    def Show_Clefts(self, ColorList):
        
        i = 0
        for Cleft in iter(self.TempBindingSite.listClefts):
            if self.PyMOL:
                try:
                    cmd.hide('everything', Cleft.CleftName)
                    cmd.show('surface', Cleft.CleftName)
                    cmd.color(ColorList[i], Cleft.CleftName)
                    cmd.refresh()
            
                    i += 1
                except:
                    self.DisplayMessage("  ERROR: Could not display sphere object '" + key + "'", 2)
                    continue
                    

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
            self.DisplayMessage('  You entered an invalid PDB code.', 1)
            
        self.FetchPDB.set('')
