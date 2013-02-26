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
import re
import pickle
import shutil
import glob
import Queue
import tkFileDialog

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
    
    ResultsName = StringVar()
    SimDefDisplay = StringVar()
    SimCartoonDisplay = IntVar()

    def __init__(self):
        
        self.ResultsContainer = Result.ResultsContainer()


class Simulate(Tabs.Tab):
    
    def Def_Vars(self):
        
        self.SimStatus = StringVar()
        self.ProgBarText = StringVar()

        # vars class objects
        self.ResultsName = self.Vars.ResultsName
        self.SimDefDisplay = self.Vars.SimDefDisplay
        self.SimCartoonDisplay = self.Vars.SimCartoonDisplay

        self.Manage = ManageFiles.Manage(self)
    
    def Init_Vars(self):
        
        self.ResultsName.set('')
        self.ProgBarText.set('... / ...')
        self.SimDefDisplay.set('sticks')
        self.SimCartoonDisplay.set(0)

        self.BarWidth = 0
        self.BarHeight = 0

        self.dictSimData = {}
        
        self.Vars.ResultsContainer.Clear()
        
        self.ProcessParsing = False
        
    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):

        self.IdleStatus()
        self.Reset_Buttons()
    
        self.ConfigMD5 = self.Manage.Hash_CONFIG()
        
        self.ResultsName_Toggle()
        
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
        Button(fConfRes, text='Load', command=self.Btn_Load_Results_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfRes, text='Save', command=self.Btn_Save_Results_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Entry(fConfRes, textvariable=self.ResultsName, font=self.top.font_Text, state='disabled', disabledbackground=self.top.Color_White,
                        disabledforeground=self.top.Color_Black, justify=CENTER).pack(side=LEFT, fill=X, expand=True)

        fNaviRes = Frame(self.fRes)
        #fNaviRes.pack(side=TOP, fill=X, expand=True, padx=5, pady=5)

        Label(fNaviRes, text='Results navigation', width=30, font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show parent', font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show next child', font=self.top.font_Text).pack(side=LEFT)
        Button(fNaviRes, text='Show previous child', font=self.top.font_Text).pack(side=LEFT)


        fTable = Frame(self.fSimulate, height=150, pady=5, padx=5)
        fTable.pack(side=BOTTOM, fill=BOTH, expand=True)
        fTable.pack_propagate(0)

        self.Table = MultiList.Table(fTable, 5,
                                   [ 'Color', 'TOP', 'CF', 'Fitness', 'Last RMSD' ],
                                   [ 65, 65, 167, 167, 167 ],
                                   [ 0, 6, 6, 6, 6 ],
                                   [ False, True, True, True, True ],
                                   self.top.font_Text,
                                   self.top.Color_Blue)
        self.Table.Draw()
        

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
       
        self.BarWidth = self.ProgressBar.cget("width")
        self.BarHeight = self.ProgressBar.cget("height")
        
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
        
        Label(fSim_DisplayLine1, text='Display options of TOP and RESULT objects', font=self.top.font_Title).pack(side=LEFT)
        
        Label(fSim_DisplayLine2, text='Ligand:', font=self.top.font_Text).pack(side=LEFT)
        Radiobutton(fSim_DisplayLine2, text='Sticks', variable=self.SimDefDisplay, value='sticks', command=self.Click_RadioSIM, font=self.top.font_Text).pack(side=RIGHT)
        Radiobutton(fSim_DisplayLine2, text='Spheres', variable=self.SimDefDisplay, value='spheres', command=self.Click_RadioSIM, font=self.top.font_Text).pack(side=RIGHT)
                        
        Label(fSim_DisplayLine3, text='Target:', font=self.top.font_Text).pack(side=LEFT)
        #Checkbutton(fSim_DisplayLine3, text='Mesh', variable=self.SimMeshDisplay, command=self.Check_MeshSIM,font=self.top.font_Text).pack(side=RIGHT)
        Checkbutton(fSim_DisplayLine3, text=' Cartoon', variable=self.SimCartoonDisplay, command=self.Check_CartoonSIM, font=self.top.font_Text).pack(side=RIGHT)
                
        return self.fSimulate
    
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

        if not self.Manage.Move_Files():
            self.DisplayMessage('  ERROR: Could not move input files to RESULT folder', 1)
            return

        if not self.Manage.Executable_Exists():
            self.DisplayMessage('  ERROR: Could not find FlexAID executable', 1)
            return

        self.DisplayMessage('   Creating input files...', 2)
        
        self.DisplayMessage('   CONFIG.inp...', 2)
        self.Manage.Create_CONFIG()
        
        self.DisplayMessage('   ga_inp.dat...', 2)
        self.Manage.Create_ga_inp(bContinue)
        
        self.DisplayMessage('   Saving and modifying input files...', 2)
        #self.Manage.Modify_Input()
        self.Manage.CreateTempPDB()
        self.Manage.Get_CoordRef()
        self.Manage.Get_VarAtoms()
        self.Manage.Get_DisAngDih()
        self.Manage.Get_RecAtom()

        #self.DisplayMessage('   RESULT file(s) will be saved in: ', 0)
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
        
        self.Results = True
        
        if bContinue:
            ParentResultsContainer = self.Vars.ResultsContainer
        
        self.Vars.ResultsContainer = Result.ResultsContainer()
        self.Vars.ResultsContainer.ConfigMD5 = self.ConfigMD5
        
        if bContinue:
            self.Vars.ResultsContainer.ParentResult = ParentResultsContainer
        
        # START SIMULATION
        self.DisplayMessage('  Starting executable thread.', 2)
        self.Start = Simulation.Start(self, commandline)
        
        # Stacking up tasks related to updating Tkinter
        self.queue = Queue.Queue()
        
        # START PARSING AS THREAD
        self.DisplayMessage('  Starting parsing thread.', 2)
        self.Parse = Simulation.Parse(self, self.queue)
        
        self.Update_Tkinter()
        
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        try:
            self.ResultsNameTrace = self.ResultsName.trace('w',self.ResultsName_Toggle)
            self.SimDefDisplayTrace = self.SimDefDisplay.trace('w',self.Click_RadioSIM)
            self.SimCartoonDisplayTrace = self.SimCartoonDisplay.trace('w',Check_CartoonSIM)
        except:
            pass    

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        try:
            self.ResultsName.trace_vdelete('w',self.ResultsNameTrace)
            self.SimDefDisplay.trace_vdelete('w',self.SimDefDisplayTrace)
            self.SimCartoonDisplay.trace_vdelete('w',self.SimCartoonDisplayTrace)
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
                
        # Empty table list
        self.Table.Clear()
        
        for key in range(1, len(self.ColorList) + 1):
            self.Table.Add( [ '', key, 0.000, 0.000, 0.000 ], 
                            [ self.ColorList[key-1], None, None, None, None ] )
                            
            self.dictSimData[key] = [ 0.0, 0.0, 0.0 ]
    
    ''' ==================================================================================
    FUNCTION update_DataResults: Updates the dictionary with the values of the results
    ==================================================================================  '''                
    def update_DataResults(self):
        
        NRes = 0
        self.dictSimData.clear()
        
        for Result in self.Vars.ResultsContainer.Results:
            
            self.dictSimData[Result.ResultID] = [ Result.CF, 'N/A', Result.RMSD ]
            
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
                self.Table.Add( [ '', key, self.dictSimData[key][0], self.dictSimData[key][1], '0.000' ],
                                [ self.top.Color_White, None, None, None, None ] )
            else:
                self.Table.Add( [ '', key, self.dictSimData[key][0], self.dictSimData[key][1], self.dictSimData[key][2] ],
                                [ self.ColorList[i], None, None, None, None ] )
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

            self.Btn_PauseResume.config(state='disabled')
            self.Btn_Stop.config(state='disabled')
            self.Btn_Abort.config(state='disabled')

            self.Clean_Update()
            self.AbortStatus()
            self.Parse.ParseFile = self.Manage.LOGFILE
            self.Results = False
    
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

            self.Btn_PauseResume.config(state='disabled')
            self.Btn_Stop.config(state='disabled')
            self.Btn_Abort.config(state='disabled')

            self.Clean_Update()
            self.StopStatus()
            self.Parse.ParseFile = self.Manage.LOGFILE
    
    ''' ==================================================================================
    FUNCTION: Loads all result files
    ==================================================================================  '''               
    def Load_Results(self):
        
        self.Manage.Load_ResultFiles()

    ''' ==================================================================================
    FUNCTION: Shows the result in the PyMOL viewer
    ==================================================================================  '''               
    def Show_Results(self):
        
        i = 0
        
        for key in sorted(self.dictSimData.keys()):
            
            Result = self.Vars.ResultsContainer.Get_ResultID(key)
            if Result is not None:
                try:
                    ResultName = 'RESULT_' + str(Result.ResultID) + '__'
                    
                    cmd.load(Result.ResultFile, ResultName, state=1)
                    cmd.refresh()
                    
                    cmd.color(self.PymolColorList[i], ResultName)
                    cmd.refresh()
                    
                except:
                    continue
                
                i += 1
        
        self.Refresh_LigDisplay()
        self.Refresh_CartoonDisplay()
    
    ''' ==================================================================================
    FUNCTION Click_RadioSIM: Change the way the ligand is displayed during a Simulation
    ==================================================================================  '''
    def Click_RadioSIM(self, *args):
        
        self.Refresh_LigDisplay()
    
    ''' ==================================================================================
    FUNCTION Refresh_LigDisplay: Refreshes the visual appearance of the ligand in TOP_* objects
    ==================================================================================  '''
    def Refresh_LigDisplay(self):
    
        if self.SimDefDisplay.get() == 'spheres':
            try:
                cmd.show('spheres', 'TOP_*__ & resn LIG')
                cmd.refresh()

                cmd.hide('sticks', 'TOP_*__ & resn LIG')
                cmd.refresh()

                cmd.show('spheres', 'RESULT_*__ & resn LIG')
                cmd.refresh()

                cmd.hide('sticks', 'RESULT_*__ & resn LIG')
                cmd.refresh()

            except:
                pass
            
        elif self.SimDefDisplay.get() == 'sticks':
            try:
                cmd.hide('spheres', 'TOP_*__ & resn LIG')
                cmd.refresh()

                cmd.show('sticks', 'TOP_*__ & resn LIG')
                cmd.refresh()

                cmd.hide('spheres', 'RESULT_*__ & resn LIG')
                cmd.refresh()

                cmd.show('sticks', 'RESULT_*__ & resn LIG')
                cmd.refresh()
                
            except:
                pass
    
    ''' ==================================================================================
    FUNCTION Check_CartoonSIM: Change the protein display adding or removing the cartoon
    ==================================================================================  '''
    def Check_CartoonSIM(self, *args):

        self.Refresh_CartoonDisplay()
        
    ''' ==================================================================================
    FUNCTION Refresh_CartoonDisplay: Refreshes the visual appearance of the protein
    ==================================================================================  '''
    def Refresh_CartoonDisplay(self):

        # Display the Cartoon
        if self.SimCartoonDisplay.get():
            try:
                cmd.show('cartoon', 'TOP_*__ & ! resn LIG')
                cmd.refresh()

                cmd.show('cartoon', 'RESULT_*__ & ! resn LIG')
                cmd.refresh()
                
            except:
                pass
            
        else:   
            # Remove the Cartoon
            try:
                cmd.hide('cartoon', 'TOP_*__ & ! resn LIG')
                cmd.refresh()

                cmd.hide('cartoon', 'RESULT_*__ & ! resn LIG')
                cmd.refresh()
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
    FUNCTION Update_Tkinter: update the tkinter interface (tasks queued from the worker)
    ==================================================================================  '''               
    def Update_Tkinter(self):
        
        # Check every 100 ms if there is something new in the queue.
        while self.queue.qsize():
            try:
                func = self.queue.get()
                func()
            except Queue.Empty:
                pass
        
        if self.top.ProcessRunning or self.ProcessParsing:
            self.top.root.after(self.top.TKINTER_UPDATE_INTERVAL, self.Update_Tkinter)
            
        else:
            # Results were generated.
            if self.Results:
                self.Load_Results()
                
                NRes = self.Manage.NUMBER_RESULTS
                #if self.top.Config2.UseReference.get():
                #    NRes = NRes + 1
                
                self.Process_ResultsContainer()
            else:
                self.Vars.ResultsContainer = None
                shutil.rmtree(self.Manage.FlexAIDRunSimulationProject_Dir, True)
            
            self.Reset_Buttons()
    
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

        if len(LoadFile) > 0:
                                
            LoadFile = os.path.normpath(LoadFile)
            
            try:
                in_ = open(LoadFile, 'r')
                TmpResultsContainer = pickle.load(in_)
                in_.close()
                
                if TmpResultsContainer is not None and len(TmpResultsContainer.Results):
                
                    self.Vars.ResultsContainer = TmpResultsContainer
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
        
        if self.Vars.ResultsContainer is None or not len(self.Vars.ResultsContainer.Results):
            self.DisplayMessage("  No result file(s) to save.", 2)
            return
        
        Results_Dir = os.path.join(self.top.FlexAIDResultsProject_Dir, self.top.IOFile.Complex.get().upper())
        if not os.path.isdir(Results_Dir):
            os.makedirs(Results_Dir)
        
        SaveFile = tkFileDialog.asksaveasfilename(initialdir=Results_Dir,
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
    
            try:
                out = open(SaveFile, 'w')
                pickle.dump(self.Vars.ResultsContainer, out)
                out.close()

                self.ResultsName.set(os.path.basename(os.path.splitext(SaveFile)[0]))

                self.DisplayMessage("  Successfully saved '" + SaveFile + "'", 2)
            except:
                self.DisplayMessage("  ERROR: Could not save Results file.", 2)
            
    ''' ==================================================================================
    FUNCTION ResultsName_Toggle: Toggles the continue button
    ==================================================================================  '''        
    def ResultsName_Toggle(self, *args):

        self.Btn_Continue.config(state='disabled')
        
        if self.ResultsName.get():
            if self.ConfigMD5 == self.Vars.ResultsContainer.ConfigMD5:
                self.Btn_Continue.config(state='normal')

    ''' ==================================================================================
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''        
    def Load_Message(self):
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < Simulate > Menu.', 2)
        self.DisplayMessage('  INFO:   Press the < Start > button to begin the simulation.', 2)
        self.DisplayMessage('           After a simulation is started, you can press the < Pause > or the < Stop > button.', 2)
