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

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
    import tkFileDialog
else:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog

from pymol import util
from subprocess import Popen

import os
import time
import re
import pickle
import shutil
import glob

import MultiList
import General
import Color
import ManageFiles
import Result
import Vars
import Tabs 

if __debug__:
    from pymol import *
    from pymol import cmd

    import pymol
    import Simulation


class SimulateVars(Vars.Vars):
    
    SimLigDisplay = StringVar()
    SimCartoonDisplay = IntVar()
    SimLinesDisplay = IntVar()

    def __init__(self):
        
        return

class Simulate(Tabs.Tab):
    
    # 100 msec
    INTERVAL = 0.10
    
    # 1 minute timeout
    TIMEOUT = INTERVAL * 300
    
    SimStatus = StringVar()
    ProgBarText = StringVar()
    
    def Def_Vars(self):

        self.ResultsName = StringVar()
        
        # vars class objects
        self.SimLigDisplay = self.Vars.SimLigDisplay
        self.SimCartoonDisplay = self.Vars.SimCartoonDisplay
        self.SimLinesDisplay = self.Vars.SimLinesDisplay
        
        self.ResultsContainer = Result.ResultsContainer()
        self.Manage = ManageFiles.Manage(self)
        
    def Init_Vars(self):

        if self.Condition_Update():
            return
        
        self.ResultsName.set('')
        self.ProgBarText.set('... / ...')
        self.SimLigDisplay.set('sticks')
        self.hasConstraints = bool(self.top.Config2.Vars.dictConstraints)
        self.SimCartoonDisplay.set(1)
        self.SimLinesDisplay.set(0)
        
        self.BarWidth = 0
        self.BarHeight = 0
        
        self.dictSimData = {}
        
        self.ResultsContainer.Clear()

        self.Paused = False
        
    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):

        self.IdleStatus()
        self.Reset_Buttons()
        
        self.ConfigMD5 = self.Manage.Hash_CONFIG()
        
        self.ResultsName_Toggle()
        
        self.BarWidth = self.ProgressBar.cget("width")
        self.BarHeight = self.ProgressBar.cget("height")

        nConstraints = bool(self.top.Config2.Vars.dictConstraints)
        if self.hasConstraints != nConstraints:
            self.hasConstraints = nConstraints
            self.Table.Undraw()
            if self.hasConstraints:
                self.Table = MultiList.Table(self.fTable, 6,
                                           [ 'Color', 'TOP', 'CF', 'Apparent CF', 'Fitness', 'Last RMSD' ],
                                           [ 40, 40, 167, 167, 120, 97 ],
                                           [ 0, 6, 6, 6, 6, 6 ],
                                           [ False, True, True, True, True, True ],
                                           self.top.font_Text,
                                           self.top.Color_Blue)
            else:
                self.Table = MultiList.Table(self.fTable, 5,
                                           [ 'Color', 'TOP', 'CF', 'Fitness', 'Last RMSD' ],
                                           [ 40, 40, 167, 120, 97 ],
                                           [ 0, 6, 6, 6, 6, 6 ],
                                           [ False, True, True, True, True ],
                                           self.top.font_Text,
                                           self.top.Color_Blue)
            # self.fTable.pack(side=BOTTOM, fill=BOTH, expand=True)
            # self.fTable.pack_propagate(0)
            self.Table.Draw()
        
    ''' =============================================================================== 
    FUNCTION Frame: Generate the CSimulation frame in the the middle frame 
                    section.
    ===============================================================================  '''      
    def Frame(self):
        
        self.fSimulate = Frame(self.top.fMiddle)

        #==================================================================================
        '''                         --- SIMULATION TABLE ---                            '''
        #==================================================================================
        
        self.fRes = Frame(self.fSimulate)#, relief=RAISED, border=1)
        self.fRes.pack(fill=BOTH, side=BOTTOM, padx=5, pady=10)

        Label(self.fRes, text='Simulation results', font=self.top.font_Title_H).pack(side=TOP, fill=X, pady=3)
        
        fConfRes = Frame(self.fRes)
        fConfRes.pack(side=TOP, fill=X, expand=True, padx=5, pady=5)
        
        Label(fConfRes, text='Pre-generated results:', width=30, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfRes, text='Load', command=self.Btn_Load_Results_Clicked, width=10, font=self.top.font_Text).pack(side=LEFT)
        #Button(fConfRes, text='Save', command=self.Btn_Save_Results_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Entry(fConfRes, textvariable=self.ResultsName, font=self.top.font_Text, state='disabled', disabledbackground=self.top.Color_White,
                        disabledforeground=self.top.Color_Black, justify=CENTER).pack(side=LEFT, fill=X, expand=True)

        self.BtnViewReport = Button(self.fRes, text='View report', command=self.Btn_View_Report, font=self.top.font_Text)
        self.BtnViewReport.pack(side=TOP, fill=X, padx=3)
        self.BtnViewReport.config(state='disabled')
        
        fNaviRes = Frame(self.fRes)
        #fNaviRes.pack(side=TOP, fill=X, expand=True, padx=5, pady=5)
        
        Label(fNaviRes, text='Results navigation', width=30, font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show parent', font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show next child', font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show previous child', font=self.top.font_Text).pack(side=LEFT)
        
        self.fTable = Frame(self.fSimulate, height=150, pady=5, padx=5)
        self.fTable.pack(side=BOTTOM, fill=BOTH, expand=True)
        self.fTable.pack_propagate(0)
        
        if self.hasConstraints:
            self.Table = MultiList.Table(self.fTable, 6,
                                       [ 'Color', 'TOP', 'CF', 'Apparent CF', 'Fitness', 'Last RMSD' ],
                                       [ 40, 40, 167, 167, 120, 97 ],
                                       [ 0, 6, 6, 6, 6, 6 ],
                                       [ False, True, True, True, True, True ],
                                       self.top.font_Text,
                                       self.top.Color_Blue)
        else:
            self.Table = MultiList.Table(self.fTable, 5,
                                       [ 'Color', 'TOP', 'CF', 'Fitness', 'Last RMSD' ],
                                       [ 40, 40, 167, 120, 97 ],
                                       [ 0, 6, 6, 6, 6, 6 ],
                                       [ False, True, True, True, True ],
                                       self.top.font_Text,
                                       self.top.Color_Blue)
        self.Table.Draw()
        
        #fSep = Frame(self.fSimulate, height=3, relief=RAISED).pack(side=BOTTOM, fill=X, expand=True, pady=5)
        
        #==================================================================================
        '''                         --- LEFT AND RIGHT ---                              '''
        #==================================================================================       

        fSimulate_Left = Frame(self.fSimulate)
        fSimulate_Left.pack(side=LEFT, fill=X, expand=True)
        
        fSimulate_Right = Frame(self.fSimulate)
        fSimulate_Right.pack(side=RIGHT, fill=X, expand=True)
 
        #==================================================================================
        '''                         --- PLAYER BUTTONS ---                             '''
        #==================================================================================       
        
        fSim_Player = Frame(fSimulate_Left)
        fSim_Player.pack(fill=X, expand=True, side=TOP, padx=5, pady=5)
        fSim_PlayerLine1 = Frame(fSim_Player)
        fSim_PlayerLine1.pack(side=TOP, fill=X)
        fSim_PlayerLine2 = Frame(fSim_Player)
        fSim_PlayerLine2.pack(side=TOP, fill=X)
        
        Label(fSim_PlayerLine1, text='Simulation controls', font=self.top.font_Title).pack(side=LEFT, anchor=W)        
        self.Btn_Start = Button(fSim_PlayerLine2, text='Start', command=self.Btn_StartSim, font=self.top.font_Text, state='normal')
        self.Btn_Start.pack(side=LEFT)
        self.Btn_Continue = Button(fSim_PlayerLine2, text='Continue', command=lambda bContinue=True: self.Btn_StartSim(bContinue),
                                   font=self.top.font_Text, state='disabled')
        self.Btn_Continue.pack(side=LEFT)
        self.Btn_PauseResume = Button(fSim_PlayerLine2, text='Pause', command=self.Btn_PauseResumeSim, font=self.top.font_Text, state='disabled')
        self.Btn_PauseResume.pack(side=LEFT)
        self.Btn_Stop = Button(fSim_PlayerLine2, text='Stop', command=self.Btn_StopSim, font=self.top.font_Text, state='disabled')
        self.Btn_Stop.pack(side=LEFT)
        self.Btn_Abort = Button(fSim_PlayerLine2, text='Abort', command=self.Btn_AbortSim, font=self.top.font_Text, state='disabled')
        self.Btn_Abort.pack(side=LEFT)
                
        #==================================================================================
        '''                           --- PROGRESSION BAR---                            '''
        #==================================================================================
        fSim_Progress = Frame(fSimulate_Left)
        fSim_Progress.pack(side=TOP, expand=True, padx=5, pady=5, anchor=W)
        fSim_ProgressLine1 = Frame(fSim_Progress)
        fSim_ProgressLine1.pack(side=TOP, fill=X)
        fSim_ProgressLine2 = Frame(fSim_Progress)
        fSim_ProgressLine2.pack(side=TOP, fill=X)
        fSim_ProgressLine3 = Frame(fSim_Progress)
        fSim_ProgressLine3.pack(side=TOP, fill=X)
        
        self.lblSimStatus = Label(fSim_ProgressLine1, bg=self.Color_Black, textvariable=self.SimStatus, font=self.top.font_Title, fg=self.top.Color_Blue)
        self.lblSimStatus.pack(side=LEFT, fill=X, expand=True, anchor=W)

        self.ProgressBar = Canvas(fSim_ProgressLine2, bg=self.top.Color_White, width=325, height=25, relief=RAISED, bd=1)
        self.ProgressBar.pack(side=LEFT, fill=BOTH, expand=True, anchor=W)
        
        self.RectPB = self.ProgressBar.create_rectangle(0, 0, 0, self.ProgressBar.winfo_reqheight(), fill=self.top.Color_Blue, width=0)
        self.TextPB = self.ProgressBar.create_text(self.ProgressBar.winfo_reqwidth()/2, self.ProgressBar.winfo_reqheight()/2, text='', fill='black')
        
        Label(fSim_ProgressLine3, text='Genetic algorithm progress:', font=self.top.font_Text).pack(side=LEFT)
        Label(fSim_ProgressLine3, textvariable=self.ProgBarText, font=self.top.font_Text).pack(side=RIGHT)
                                        
        #==================================================================================
        '''                           --- DISPLAY OPTIONS ---                           '''
        #==================================================================================
        fSim_Display = Frame(fSimulate_Right, padx=5, pady=5)
        fSim_Display.pack(fill=BOTH,expand=True)
        fSim_DisplayLine1 = Frame(fSim_Display)
        fSim_DisplayLine1.pack(side=TOP, fill=X)
        fSim_DisplayLine2 = Frame(fSim_Display)
        fSim_DisplayLine2.pack(side=TOP, fill=X)
        fSim_DisplayLine3 = Frame(fSim_Display)
        fSim_DisplayLine3.pack(side=TOP, fill=X)
        
        Label(fSim_DisplayLine1, text='Display options of TOP objects', font=self.top.font_Title).pack(side=LEFT)
        
        Label(fSim_DisplayLine2, text='Ligand:', font=self.top.font_Text).pack(side=LEFT)
        
        Radiobutton(fSim_DisplayLine2, text='Sticks', variable=self.SimLigDisplay, value='sticks', 
                    font=self.top.font_Text).pack(side=RIGHT)
                    
        Radiobutton(fSim_DisplayLine2, text='Spheres', variable=self.SimLigDisplay, value='spheres',
                    font=self.top.font_Text).pack(side=RIGHT)
                        
        Label(fSim_DisplayLine3, text='Target:', font=self.top.font_Text).pack(side=LEFT)
        
        Checkbutton(fSim_DisplayLine3, text=' Cartoon', variable=self.SimCartoonDisplay, 
                    command=lambda var=self.SimCartoonDisplay, display='cartoon': self.Modify_Display(var, display),
                    font=self.top.font_Text).pack(side=RIGHT)
                    
        Checkbutton(fSim_DisplayLine3, text=' Lines', variable=self.SimLinesDisplay,
                    command=lambda var=self.SimLinesDisplay, display='lines': self.Modify_Display(var, display),
                    font=self.top.font_Text).pack(side=RIGHT)
        
        return self.fSimulate
        
    ''' =============================================================================== 
    FUNCTION Btn_View_Report: Views the final report of a simulation
    ===============================================================================  '''     
    def Btn_View_Report(self):

        if self.ResultsContainer.Report:

            if self.top.OSid == 'WIN':
                apps = ['notepad']
            elif self.top.OSid == 'MAC':
                apps = ['open']
            elif self.top.OSid == 'LINUX':
                apps = ['gedit','kate','kedit']

            for app in apps:
                try:
                    Popen([app, self.ResultsContainer.Report])
                    return
                except OSError:
                    pass
                except:
                    pass

            self.top.DisplayMessage("  ERROR: No text editor found for your operating system", 2)
        
    ''' =============================================================================== 
    FUNCTION Btn_ContinueSim: Continue a simulation from an existing one
    ===============================================================================  '''     
    def Btn_ContinueSim(self):

        return

    ''' =============================================================================== 
    FUNCTION Btn_StartSim: Start the simulation (If requirements are meet...)
    ===============================================================================  '''     
    def Btn_StartSim(self, bContinue=False):
        
        self.Manage.Reference_Folders()

        if not self.Manage.Clean():
            self.DisplayMessage('  ERROR: Could not clean files before the simulation', 1)
            return
                
        if not self.Manage.Create_Folders():
            self.DisplayMessage('  ERROR: Could not create complex folder', 1)
            return
        
        #if not self.Manage.Move_Files():
        #    self.DisplayMessage('  ERROR: Could not move input files to RESULT folder', 1)
        #    return
        
        if not self.Manage.Executable_Exists():
            self.DisplayMessage('  ERROR: Could not find FlexAID executable', 1)
            return

        self.DisplayMessage('   Creating input files...', 2)
        
        self.DisplayMessage('   CONFIG.inp...', 2)
        self.Manage.Create_CONFIG()
        
        self.DisplayMessage('   ga_inp.dat...', 2)
        self.Manage.Create_ga_inp(bContinue)

        self.DisplayMessage('   Writing report...', 2)
        self.Manage.Write_Report(bContinue)
        
        self.DisplayMessage('   Saving and modifying input files...', 2)
        #self.Manage.Modify_Input()
        self.Manage.CreateTempPDB()
        self.Manage.Get_CoordRef()
        self.Manage.Get_VarAtoms()
        self.Manage.Get_DisAngDih()
        self.Manage.Get_RecAtom()
        
        self.DisplayMessage('   RESULT file(s) will be saved in: ', 0)
        self.DisplayMessage(self.Manage.FlexAIDRunSimulationProject_Dir, 0)

        # Start the simulation
        self.Btn_Start.config(state='disabled')
        self.Btn_Continue.config(state='disabled')
        self.Btn_PauseResume.config(state='normal')
        self.Btn_Stop.config(state='normal')
        self.Btn_Abort.config(state='normal')

        self.ColorList = Color.GetHeatColorList(int(self.top.GAParam.NbTopChrom.get()), True)
        self.PymolColorList = Color.GetHeatColorList(int(self.top.GAParam.NbTopChrom.get()), False)

        self.Init_Table()

        # START FLEXAID AS THREAD
        commandline =   '"%s" "%s" "%s" "%s"' % (   self.top.FlexAIDExecutable,
                                                    self.Manage.CONFIG,
                                                    self.Manage.ga_inp,
                                                    os.path.join(self.Manage.FlexAIDRunSimulationProject_Dir,'RESULT') )
        self.Results = False
        self.Paused = False
        
        if bContinue:
            ParentResultsContainer = self.ResultsContainer
        
        self.ResultsContainer = Result.ResultsContainer()
        self.ResultsContainer.ConfigMD5 = self.ConfigMD5

        if bContinue:
            self.ResultsContainer.ParentResult = ParentResultsContainer
        
        self.ResultsContainer.Report = self.Manage.Report
        
        # return codes description
        #       -1: Parse thread has not started yet
        #        0: Parsing the logfile
        #        1: NRGsuite update error
        #        2: FlexAID error parsed
        #       10: Done parsing
        self.top.ParseState = -1
        
        # return codes description
        #      -1: FlexAID Simulation thread not started yet
        #       0: FlexAID is running or done when FlexAID.Run is None
        # [1-100[: FlexAID return codes error
        #     100: IOError exception when running Simulate thread
        #     200: OSError exception when running Simulate thread
        #     300: Other exception when running Simulate thread
        self.top.SimulateState = -1

        self.Start_Update()
        
        # START PARSING AS THREAD
        self.DisplayMessage('  Starting parsing thread.', 2)
        self.Parse = Simulation.Parse(self, self.queue)
        
        while self.top.ParseState < 0:
            time.sleep(self.INTERVAL)

        # START SIMULATION
        self.DisplayMessage('  Starting executable thread.', 2)
        self.Start = Simulation.Start(self, commandline)
        
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        try:
            self.ResultsNameTrace = self.ResultsName.trace('w',self.ResultsName_Toggle)
            self.SimLigDisplayTrace = self.SimLigDisplay.trace('w',self.Modify_LigDisplay)
        except:
            pass    

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        try:
            self.ResultsName.trace_vdelete('w',self.ResultsNameTrace)
            self.SimLigDisplay.trace_vdelete('w',self.SimLigDisplayTrace)
        except:
            pass

    ''' =============================================================================== 
    FUNCTION IdleSim: Simulation was not started yet
    ===============================================================================  '''     
    def IdleStatus(self):

        self.lblSimStatus.config(fg='white')
        self.SimStatus.set('Idle.')

    ''' =============================================================================== 
    FUNCTION InitStatus: Initializes necessary variables from FlexAID parsing thread
    ===============================================================================  '''   
    def InitStatus(self):

        self.lblSimStatus.config(fg='cyan')
        self.SimStatus.set('Initializing...')

    ''' ===============================================================================
    FUNCTION RunStatus: Signal given to parse genetic algorithm
    ===============================================================================  '''
    def RunStatus(self):

        self.lblSimStatus.config(fg='blue')
        self.SimStatus.set('Running...')

    ''' ===============================================================================
    FUNCTION PauseStatus: Signal given to pause FlexAID
    ===============================================================================  '''
    def PauseStatus(self):

        self.lblSimStatus.config(fg='purple')
        self.SimStatus.set('Paused.')

    ''' =============================================================================== 
    FUNCTION StopStatus: Signal given to stop FlexAID
    ===============================================================================  '''     
    def StopStatus(self):

        self.lblSimStatus.config(fg='yellow')
        self.SimStatus.set('Stopped.')

    ''' =============================================================================== 
    FUNCTION AbortStatus: Signal given to pause FlexAID
    ===============================================================================  '''     
    def AbortStatus(self):

        self.lblSimStatus.config(fg='orange')
        self.SimStatus.set('Aborting...')

    ''' =============================================================================== 
    FUNCTION ClusterStatus: GA has ended - clustering all individuals
    ===============================================================================  '''     
    def ClusterStatus(self):

        self.lblSimStatus.config(fg='blue violet')
        self.SimStatus.set('Clustering...')

    ''' =============================================================================== 
    FUNCTION SuccessStatus: GA has ended - clustering all individuals
    ===============================================================================  '''     
    def SuccessStatus(self):

        self.lblSimStatus.config(fg='chartreuse')
        
        if self.Results:
            self.SimStatus.set('Simulation ended successfully.')
        else:
            self.SimStatus.set('Simulation aborted successfully.')
    
    ''' =============================================================================== 
    FUNCTION ErrorStatus: An error occured.
    ===============================================================================  '''     
    def ErrorStatus(self, ErrorMsg):

        self.lblSimStatus.config(fg='red')
        self.SimStatus.set(ErrorMsg)
    
    ''' ==================================================================================
    FUNCTION Init_Table: Initialisation the dictionary that contain the energy 
                                                 and fitness for each solution.
    ==================================================================================  '''                
    def Init_Table(self):
                
        self.dictSimData.clear()
        
        # Empty table list
        self.Table.Clear()
        
        for key in range(1, len(self.ColorList) + 1):
            if self.hasConstraints:
                self.Table.Add( [ '', key, 0.000, 0.000, 0.000, 0.000 ],
                            [ self.ColorList[key-1], None, None, None, None ] )
            else:
                self.Table.Add( [ '', key, 0.000, 0.000, 0.000 ],
                            [ self.ColorList[key-1], None, None, None, None ] )

            self.dictSimData[key] = [ 0.0, 0.0, 0.0, 'N/A' ]
    
    ''' ==================================================================================
    FUNCTION update_DataResults: Updates the dictionary with the values of the results
    ==================================================================================  '''                
    def update_DataResults(self):
        
        NRes = 0
        self.dictSimData.clear()
        
        for Result in self.ResultsContainer.Results:
            
            self.dictSimData[Result.ResultID] = [ Result.CF, Result.CFapp, 'N/A', Result.RMSD ]
            
            if Result.ResultID != 'REF':
                NRes = NRes + 1
            
        return NRes
    
    ''' ==================================================================================
    FUNCTION update_DataList: Update the displayed Data List informations.
    ==================================================================================  '''                
    def update_DataList(self):

        self.Table.Clear()
        
        i = 0
        for key in sorted(self.dictSimData.keys()):
            if key == 'REF':
                if self.hasConstraints:
                    self.Table.Add( [ '', key, self.dictSimData[key][0], self.dictSimData[key][1], 
                                  self.dictSimData[key][2], '0.000' ],
                                [ self.top.Color_White, None, None, None, None, None ] )
                else:
                    self.Table.Add( [ '', key, self.dictSimData[key][0], 
                                  self.dictSimData[key][2], '0.000' ],
                                [ self.top.Color_White, None, None, None, None, None ] )
            else:
                if self.hasConstraints:
                    self.Table.Add( [ '', key, self.dictSimData[key][0], self.dictSimData[key][1],
                                      self.dictSimData[key][2], self.dictSimData[key][3] ],
                                    [ self.ColorList[i], None, None, None, None, None ] )
                else:
                    self.Table.Add( [ '', key, self.dictSimData[key][0],
                                      self.dictSimData[key][2], self.dictSimData[key][3] ],
                                    [ self.ColorList[i], None, None, None, None, None ] )
                i += 1
    
    ''' =============================================================================== 
    FUNCTION Reset_Buttons(self): resets button states back to defaults
    ===============================================================================  '''        
    def Reset_Buttons(self):
        
        self.Btn_Start.config(state='normal')
        self.Btn_Continue.config(state='disabled')
        self.Btn_PauseResume.config(state='disabled')
        self.Btn_Stop.config(state='disabled')
        self.Btn_Abort.config(state='disabled')
    
    ''' =============================================================================== 
    FUNCTION Clean_Update: Remove the .update file is present
    ===============================================================================  '''        
    def Clean_Update(self):
    
        if os.path.isfile(self.Manage.UPDATE):
            if self.Parse.Remove_UPDATE():
                self.Parse.Error = True
                self.Parse.ErrorMsg = '*NRGsuite ERROR: Could not successfully remove .update file'

    ''' =============================================================================== 
    FUNCTION Btn_PauseResumeSim: Pauses/Resumes the simulation   
    ===============================================================================  '''        
    def Btn_PauseResumeSim(self):
                
        if self.SimStatus.get() == 'Running...':
        
            try:
                #Create the .pause file            
                pause_file = open(self.Manage.PAUSE, 'w')
                pause_file.close()
                
            except OSError:
                self.DisplayMessage('  ERROR: An error occured while trying to pause the simulation.', 0)
                return
                
            self.Btn_PauseResume.config(text='Resume')
            self.Btn_Stop.config(state='disabled')
            self.Btn_Abort.config(state='disabled')

            self.PauseStatus()
            self.Paused = True

        elif self.SimStatus.get() == 'Paused.':
        
            if os.path.isfile(self.Manage.PAUSE):
                try:
                    os.remove(self.Manage.PAUSE)
                except OSError:
                    self.DisplayMessage('  ERROR: An error occured while trying to resume the simulation.', 0)
                    pass

            self.Btn_PauseResume.config(text='Pause')
            self.Btn_Stop.config(state='normal')
            self.Btn_Abort.config(state='normal')

            self.RunStatus()
            self.Paused = False

    ''' =============================================================================== 
    FUNCTION Btn_AbortSim: Abort the simulation 
    ===============================================================================  '''    
    def Btn_AbortSim(self):
        
        if self.SimStatus.get() == 'Running...':
            
            try:
                #Create the .abort file
                abort_file = open(self.Manage.ABORT, 'w')
                abort_file.close()
            except OSError:
                self.DisplayMessage('  ERROR: An error occured while trying to abort the simulation.', 0)
                return
            try:  
                self.Btn_PauseResume.config(state='disabled')
                self.Btn_Stop.config(state='disabled')
                self.Btn_Abort.config(state='disabled')

                self.Clean_Update()
                self.AbortStatus()
                
                self.Parse.ParseFile = self.Manage.LOGFILE
                self.Paused = False
            except:
                self.DisplayMessage('  ERROR: An errror occured while reseting parameters after .abort.',0)
    
    ''' =============================================================================== 
    FUNCTION Btn_StopSim: Stop the simulation 
    ===============================================================================  '''    
    def Btn_StopSim(self):
        
        if self.SimStatus.get() == 'Running...':
 
            try:
                #Create the .stop file
                stop_file = open(self.Manage.STOP, 'w')
                stop_file.close()
                
            except OSError:
                self.DisplayMessage('  ERROR: An error occured while trying to stop the simulation.', 0)
                return


            try:
                self.Btn_PauseResume.config(state='disabled')
                self.Btn_Stop.config(state='disabled')
                self.Btn_Abort.config(state='disabled')

                self.Clean_Update()
                self.StopStatus()
                self.Parse.ParseFile = self.Manage.LOGFILE
                self.Paused = False
            except:
                self.DisplayMessage('  ERROR: An error occured while trying to reset paraments after .stop.',0)
                return

    ''' ==================================================================================
    FUNCTION: Loads all result files
    ==================================================================================  '''               
    def Load_Results(self):
        
        self.Manage.Load_ResultFiles()
        
    ''' ==================================================================================
    FUNCTION: Displays nicely the result complex
    ==================================================================================  '''               
    def Nice_Display(self, Result, ResultName):
        
        cmd.hide('everything', ResultName)
        cmd.refresh()
        
        cmd.show('cartoon', ResultName)
        cmd.refresh()
                
        for opt in Result.Optimizable:
            if opt.rnc[:3] == 'LIG':
                cmd.show('lines', 'byres(resn LIG around 5.0) & ' + ResultName)
            
            res = opt.rnc[:3].replace('-','')
            num = opt.rnc[3:len(opt.rnc)-1].replace('-','')
            chn = opt.rnc[-1:].replace('-','')
            
            sele =  'resn ' + res + ' & resi ' + num + ' & chain \''  + chn + '\' & ' + ResultName
            cmd.show('sticks', sele)
            cmd.color('white', sele)
        
        util.cnc(ResultName)
        cmd.refresh()
    
    ''' ==================================================================================
    FUNCTION: Displays the h-bonds
    ==================================================================================  '''               
    def Highlight_HBonds(self, Result, ResultName, ResultHBondsName):
        
        cmd.distance(ResultHBondsName, 
                     'resn LIG & ' + ResultName,
                     ResultName + ' & !resn LIG',
                     3.5, # distance in A
                     2)   # mode (2=polar atoms)

    ''' ==================================================================================
    FUNCTION: Shows the result in the PyMOL viewer
    ==================================================================================  '''               
    def Show_Results(self):

        i = 0

        for key in sorted(self.dictSimData.keys()):
            
            Result = self.ResultsContainer.Get_ResultID(key)
            if Result is not None:
                try:
                    ResultName = 'RESULT_' + str(Result.ResultID) + '__'
                    ResultHBondsName = 'RESULT_' + str(Result.ResultID) + '_H_BONDS__'

                    cmd.load(Result.ResultFile, ResultName, state=1)
                    cmd.refresh()

                    cmd.color(self.PymolColorList[i], ResultName)
                    util.cnc(ResultName)
                    cmd.refresh()

                    self.Nice_Display(Result, ResultName)
                    self.Highlight_HBonds(Result, ResultName, ResultHBondsName)
                    
                except:
                    continue

                i += 1
                    
        #self.Modify_LigDisplay()
        #self.Modify_Display(self.SimCartoonDisplay, 'cartoon')
        #self.Modify_Display(self.SimLinesDisplay, 'lines')
    
    ''' ==================================================================================
    FUNCTION Modify_LigDisplay: Modifies how the ligand is visualized in the TOP*/RESULT* objects
    ==================================================================================  '''
    def Modify_LigDisplay(self, *args):
        
        display = self.SimLigDisplay.get()
        
        try:
            cmd.hide('everything', 'TOP_*__ & resn LIG')
            cmd.refresh()

            #cmd.hide('everything', 'RESULT_*__ & resn LIG')
            #cmd.refresh()
            
            cmd.show(display, 'TOP_*__ & resn LIG')
            cmd.refresh()

            #cmd.show(display, 'RESULT_*__ & resn LIG')
            #cmd.refresh()

        except:
            pass
                        
    ''' ==================================================================================
    FUNCTION Modify_Display: Modifies how the target is visualized in the TOP*/RESULT* objects
    ==================================================================================  '''
    def Modify_Display(self, var, display):
    
        if var.get():
        
            try:
                cmd.show(display, 'TOP_*__ & ! resn LIG')
                cmd.refresh()

                #cmd.show(display, 'RESULT_*__ & ! resn LIG')
                #cmd.refresh()
                
            except:
                pass
            
        else:   
            try:
                cmd.hide(display, 'TOP_*__ & ! resn LIG')
                cmd.refresh()

                #cmd.hide(display, 'RESULT_*__ & ! resn LIG')
                #cmd.refresh()
            except:
                pass                
    
    '''
    @summary: SUBROUTINE progressBarHandler: Update the progression bar in the interface                  
    '''
    def progressBarHandler(self, Generation, Total):

        try:
            percentage = float(Generation)/float(Total)

            if percentage > 1.0:
                percentage = 1.0

            elif percentage < 0.0:
                percentage = 0.0

            # Change the Generation(s) done display
            self.ProgBarText.set(str(Generation).rjust(4, ' ') + "/" + str(Total))

            text = str(int(round(100.0 * percentage))) + " %"
            barValue = int(float(self.BarWidth) * percentage)

            self.ProgressBar.coords(self.RectPB, 0, 0, barValue, self.BarHeight)
            self.ProgressBar.itemconfigure(self.TextPB, text=text)
            
        except:
            pass

        return
    
    ''' ==================================================================================
    FUNCTION Condition_Update: Tests tab specific conditions to trigger stopping
    ==================================================================================  '''               
    def Condition_Update(self):
        
        if self.top.ParseState != 10:
            return True
        else:
            return False
            
    ''' ==================================================================================
    FUNCTION After_Update: Executes tasks when done updating Tkinter
    ==================================================================================  '''               
    def After_Update(self):
        
        self.Reset_Buttons()

        # Results were generated.
        if self.Results:
            #print "Results were generated!"
            self.Load_Results()

            NRes = self.Manage.NUMBER_RESULTS
            #if self.top.Config2.UseReference.get():
            #    NRes = NRes + 1

            self.Process_ResultsContainer()

            Results_Dir = os.path.join(self.top.FlexAIDResultsProject_Dir, self.top.IOFile.Complex.get().upper())
            if not os.path.isdir(Results_Dir):
                os.makedirs(Results_Dir)

            Results_File = os.path.join(Results_Dir,self.Manage.Now + '.nrgfr')

            self.Save_Results(Results_File)
            self.ResultsName.set(os.path.splitext(os.path.split(Results_File)[1])[0])

        else:
            self.ResultsContainer.Clear()
            try:
                pass
                # shutil.rmtree(self.Manage.FlexAIDRunSimulationProject_Dir, True)
            except shutil.Error:
                pass

            self.ResultsName.set('')
                
    ''' ==================================================================================
    FUNCTION Process_ResultsContainer: updates the data and show the results
    ==================================================================================  '''
    def Process_ResultsContainer(self):
        
        NRes = self.update_DataResults()
        self.ColorList = Color.GetHeatColorList(NRes, True)
        self.PymolColorList = Color.GetHeatColorList(NRes, False)
        
        self.Show_Results()
        self.update_DataList()

    ''' ==================================================================================
    FUNCTION Btn_Load_Results_Clicked: Loads a previously saved results
    ==================================================================================  '''        
    def Btn_Load_Results_Clicked(self):

        if self.top.ValidateProcessRunning() or self.top.ValidateWizardRunning() or \
            self.top.ValidateWindowRunning():
            return

        Results_Dir = os.path.join(self.top.FlexAIDResultsProject_Dir, self.top.IOFile.Complex.get().upper())
        if not os.path.isdir(Results_Dir):
            Results_Dir = self.top.FlexAIDResultsProject_Dir

        LoadFile = tkFileDialog.askopenfilename(initialdir=Results_Dir,
                                                filetypes=[('NRG FlexAID Results','*.nrgfr')],
                                                title='Load a Results file')
        # LoadFile = self.top.root.master.splitlist(LoadFile)
        if len(LoadFile) > 0:
                                
            LoadFile = os.path.normpath(LoadFile)
            
            try:
                in_ = open(LoadFile, 'rb')
                TmpResultsContainer = pickle.load(in_)
                in_.close()
                
                if TmpResultsContainer is not None and len(TmpResultsContainer.Results):
                
                    self.ResultsContainer = TmpResultsContainer
                    self.Process_ResultsContainer()
                    
                    self.ResultsName.set(os.path.basename(os.path.splitext(LoadFile)[0]))
                else:
                    self.DisplayMessage("  ERROR: The Results file is empty or has unknown format.", 2)

                self.DisplayMessage("  Successfully loaded the Results file '" + os.path.basename(os.path.splitext(LoadFile)[0]) + "'", 2)
                
            except:
                self.DisplayMessage("  ERROR: Could not read the Results file.", 2)
    
    ''' ==================================================================================
    FUNCTION Btn_Save_Results_Clicked: Saves the current results
    ==================================================================================  '''        
    def Btn_Save_Results_Clicked(self):

        if self.top.ValidateProcessRunning() or self.top.ValidateWizardRunning() or \
           self.top.ValidateWindowRunning():
            return
        
        if self.ResultsContainer is None or not len(self.ResultsContainer.Results):
            self.DisplayMessage("  No result file(s) to save.", 2)
            return
        
        Results_Dir = os.path.join(self.top.FlexAIDResultsProject_Dir, self.top.IOFile.Complex.get().upper())
        if not os.path.isdir(Results_Dir):
            os.makedirs(Results_Dir)
        
        SaveFile = tkFileDialog.asksaveasfilename(filetypes=[('NRG FlexAID Results','*.nrgfr')],
                                                  initialdir=Results_Dir,
                                                  title='Save the Results file', initialfile='default_results',
                                                  defaultextension='.nrgfr')
        
        if len(SaveFile) > 0:

            SaveFile = os.path.normpath(SaveFile)
            
            if General.validate_String(SaveFile, '.nrgfr', True, True, False):
                self.DisplayMessage("  ERROR: Could not save the file because you entered an invalid name.", 2)
                return
                
            if self.top.ValidateSaveProject(SaveFile, 'Results'):
                self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                return

            self.Save_Results(SaveFile)
            
    ''' ==================================================================================
    FUNCTION ResultsName_Toggle: Toggles the continue button
    ==================================================================================  '''        
    def Save_Results(self, SaveFile):
                            
        try:
            out = open(SaveFile, 'wb')
            pickle.dump(self.ResultsContainer, out)
            out.close()
            
            #self.ResultsName.set(os.path.basename(os.path.splitext(SaveFile)[0]))
            
            self.DisplayMessage("  Successfully saved '" + SaveFile + "'", 2)
            
        except:
            self.DisplayMessage("  ERROR: Could not save Results file.", 2)

    ''' ==================================================================================
    FUNCTION ResultsName_Toggle: Toggles the continue button
    ==================================================================================  '''        
    def ResultsName_Toggle(self, *args):

        self.Btn_Continue.config(state='disabled')
        self.BtnViewReport.config(state='disabled')
        
        if self.ResultsName.get():
            self.BtnViewReport.config(state='normal')
            if self.ConfigMD5 == self.ResultsContainer.ConfigMD5:
                self.Btn_Continue.config(state='normal')

    ''' ==================================================================================
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''        
    def Load_Message(self):
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < Simulate > Menu.', 2)
        self.DisplayMessage('  INFO:   Press the < Start > button to begin the simulation.', 2)
        self.DisplayMessage('           After a simulation is started, you can press the < Pause > or the < Stop > button.', 2)
