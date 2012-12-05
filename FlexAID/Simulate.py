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
import tkTable
import General
import Color
import ManageFiles

if __debug__:
    from pymol import *
    from pymol import cmd

    import pymol
    import Simulation


class Simulate:

    def __init__(self,top,PyMOL):
	        
		self.PyMOL = PyMOL
		self.top = top
		
		self.Tab = self.top.Btn_Simulate
		self.FrameName = 'Simulate'
		
		self.Def_Vars()
		self.Init_Vars()
		
		self.Frame()
		self.Trace()

    def Def_Vars(self):
        
        self.SimStatus = StringVar()
        self.ProgBarText = StringVar()
        self.SimDefDisplay = StringVar()
        self.SimMeshDisplay = IntVar()
        self.SimLinesDisplay = IntVar()
        self.SimCartoonDisplay = IntVar()
        self.dictSimData = dict()
        self.Validator = list()


    def Init_Vars(self):

        self.currentClick = None
        self.ProcessError = False

        self.SimStatus.set('Not running.')
        self.ProgBarText.set('... / ...')
        self.SimDefDisplay.set('sticks')
        self.SimMeshDisplay.set(0)
        self.SimLinesDisplay.set(1)
        self.SimCartoonDisplay.set(0)

        self.BarWidth = 0
        self.BarHeight = 0

        self.dictSimData = {}

        self.Manage = None

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
        
        self.fSimulate.pack_forget()
        #self.fSimulate.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    def Validator_Fail(self):
        
        pass

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):

        self.LoadMessage()

        self.fSimulate.pack(fill=BOTH, expand=True)

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        return

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''    
    def Del_Trace(self):

        return
    ''' =============================================================================== 
    FUNCTION Frame: Generate the CSimulation frame in the the middle frame 
                    section.   
    ===============================================================================  '''      
    def Frame(self):
        
        self.fSimulate = Frame(self.top.fMiddle)

        #==================================================================================
        '''                         --- SIMULATION TABLE ---                            '''
        #==================================================================================
        
        fTable = Frame(self.fSimulate, height=200, pady=5, padx=5)
        fTable.pack(side=BOTTOM, fill=BOTH, expand=True)
        fTable.pack_propagate(0)
        
        self.Table = tkTable.Table(fTable, 5,
                                   [ 'Color', 'TOP', 'Energy', 'Fitness', 'RMSD' ],
                                   [ 65, 65, 167, 167, 167 ],
                                   [ 0, 6, 15, 15, 15 ],
                                   [ False, True, True, True, True ],
                                   self.top.font_Text,
                                   self.top.Color_Blue)
        self.Table.Draw()
        
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
        self.Btn_Start = Button(fSim_PlayerLine2, text='Start', command=self.Btn_StartSim, font=self.top.font_Text)
        self.Btn_Start.pack(side=LEFT)
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
        
        Label(fSim_ProgressLine1, text='Status:', font=self.top.font_Title).pack(side=LEFT) 
        Label(fSim_ProgressLine1, textvariable=self.SimStatus, font=self.top.font_Title, fg=self.top.Color_Blue).pack(side=LEFT) 

        self.ProgressBar = Canvas(fSim_ProgressLine2, bg='white', width=325, height=25, relief=RAISED, bd=1)
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
        
        Label(fSim_DisplayLine1, text='PyMOL display options', font=self.top.font_Title).pack(side=LEFT)
        
        Label(fSim_DisplayLine2, text='Ligand:', font=self.top.font_Text).pack(side=LEFT)
        Radiobutton(fSim_DisplayLine2, text='Sticks', variable=self.SimDefDisplay, value='sticks', command=self.Click_RadioSIM, font=self.top.font_Text).pack(side=RIGHT)
        Radiobutton(fSim_DisplayLine2, text='Spheres', variable=self.SimDefDisplay, value='spheres', command=self.Click_RadioSIM, font=self.top.font_Text).pack(side=RIGHT)
                        
        Label(fSim_DisplayLine3, text='Target:', font=self.top.font_Text).pack(side=LEFT)
        Checkbutton(fSim_DisplayLine3, text='Mesh', variable=self.SimMeshDisplay, command=self.Check_MeshSIM,font=self.top.font_Text).pack(side=RIGHT)
        Checkbutton(fSim_DisplayLine3, text=' Cartoon', variable=self.SimCartoonDisplay, command=self.Check_CartoonSIM, font=self.top.font_Text).pack(side=RIGHT)
                
                       
    ''' =============================================================================== 
    FUNCTION Btn_StartSim: Start the simulation (If requirements are meet...)   
    ===============================================================================  '''     
    def Btn_StartSim(self):

        self.Manage = ManageFiles.Manage(self)

        if not self.Manage.Clean():
            self.top.DisplayMessage('   Fatal error: Could not clean files before the simulation',1)
            self.top.DisplayMessage('   Please contact the developers of the NRGsuite',1)
            return
        
        if not self.Manage.Create_Folders():
            self.top.DisplayMessage('   Fatal error: Could not create RESULT folder',1)
            self.top.DisplayMessage('   Please contact the developers of the NRGsuite',1)
            return

        if not self.Manage.Move_Files():
            self.top.DisplayMessage('   Fatal error: Could not move input files to RESULT folder',1)
            self.top.DisplayMessage('   Please contact the developers of the NRGsuite',1)
            return

        if not self.Manage.Executable_Exists():
            self.top.DisplayMessage('   Fatal error: Could not find FlexAID executable',1)
            self.top.DisplayMessage('   Please contact the developers of the NRGsuite',1)
            return

                
        self.top.DisplayMessage('   Will create input files...', 2)
        self.Manage.Create_CONFIG()

        self.top.DisplayMessage('   CONFIG.inp file created.', 2)
        self.Manage.Create_ga_inp()
        self.top.DisplayMessage('   ga_inp.dat file created.', 2)
        self.Manage.Modify_Input()

        self.top.DisplayMessage('   RESULT file(s) will be saved in: ', 0)
        self.top.DisplayMessage(self.Manage.FlexAIDRunSimulationProject_Dir, 0)

        # Start the simulation
        self.Btn_Start.config(state='disabled') 
        self.Btn_PauseResume.config(state='normal')
        self.Btn_Stop.config(state='normal')
        self.Btn_Abort.config(state='normal')
        self.SimStatus.set('Running...')

        # START FLEXAID AS THREAD
        commandline =   '"%s" "%s" "%s" "%s"' % (   self.top.FlexAIDExecutable,
                                                    self.Manage.CONFIG,
                                                    self.Manage.ga_inp,
                                                    self.Manage.FlexAIDRunSimulationProject_Dir )

        # START PARSING AS THREAD
        self.top.DisplayMessage('   Starting parsing thread.', 2)
        Parse = Simulation.Parse(self)
        
        # START SIMULATION
        self.top.DisplayMessage('   Starting executable thread.', 2)
        Start = Simulation.Start(self, commandline)
        
        self.top.DisplayMessage("For better performance you can disable object SPHERE_AREA__ or GRID_AREA__", 0)
    
    
    ''' ==================================================================================
    FUNCTION Init_Table: Initialisation the dictionary that contain the energy 
                                                 and fitness for each solution.
    ==================================================================================  '''                
    def Init_Table(self):
        
        nbTopChrom = int(self.top.GAParam.NbTopChrom.get())
        self.ColorList = Color.GetHeatColorList(nbTopChrom, True)

        # Empty table list
        self.Table.Clear()

        for key in range(1, nbTopChrom + 1):
            self.Table.Add( [ '', key, 0.000, 0.000, 0.000 ], 
                            [ self.ColorList[key-1], None, None, None, None ] )
            
        #self.Table.Delete('1', 'TOP')
        #self.Table.Set('5', 'TOP', 1.000, 'Energy')

    ''' ==================================================================================
    FUNCTION update_DataList: Update the displayed Data List informations.
    ==================================================================================  '''                
    def update_DataList(self):

        self.Table.Clear()

        i = 0
        for key in iter(self.dictSimData):
            self.Table.Add( [ '', key, self.dictSimData[key][0], self.dictSimData[key][1], self.dictSimData[key][2] ], 
                            [ self.ColorList[i], None, None, None, None ] )
            i += 1
                
    ''' =============================================================================== 
    FUNCTION Btn_PauseResumeSim: Pauses/Resumes the simulation   
    ===============================================================================  '''        
    def Btn_PauseResumeSim(self):
                
        if self.SimStatus.get() == 'Running...':
        
          self.Btn_PauseResume.config(text='Resume...')
          self.Btn_Stop.config(state='disabled')
          self.Btn_Abort.config(state='disabled')
          self.SimStatus.set('Paused.')
          
          #Create the .pause file            
          pause_file = open(self.Manage.PAUSE, 'w')
          pause_file.close()
          
          self.top.DisplayMessage('   the simulation...', 0)
          
        elif self.SimStatus.get() == 'Paused.':
        
          self.Btn_PauseResume.config(text='Pause...')
          self.Btn_Stop.config(state='normal')
          self.Btn_Abort.config(state='normal')
          self.SimStatus.set('Running...')
        
          if os.path.isfile(self.Manage.PAUSE):
              try:
                  os.remove(self.Manage.PAUSE)
              except OSError:
                  time.sleep(0.1)
                  os.remove(self.Manage.PAUSE)
          
          self.top.DisplayMessage('  Resuming the simulation...', 0)
        
    ''' =============================================================================== 
    FUNCTION Btn_AbortSim: Abort the simulation 
    ===============================================================================  '''    
    def Btn_AbortSim(self):
        
        if self.SimStatus.get() == 'Running.':
 
            self.Btn_PauseResume.config(state='disabled')
            self.Btn_Stop.config(state='disabled')
            self.Btn_Abort.config(state='disabled')
            
            self.SimStatus.set('Not running: Aborted.')
            
            #Create the .stop file
            abort_file = open(self.Manage.ABORT, 'w')
            abort_file.close()
            
            self.top.DisplayMessage('  Abort the simulation.', 0)
            

    ''' =============================================================================== 
    FUNCTION Btn_StopSim: Stop the simulation 
    ===============================================================================  '''    
    def Btn_StopSim(self):
        
        if self.SimStatus.get() == 'Running...':
 
            self.Btn_PauseResume.config(state='disabled')
            self.Btn_Stop.config(state='disabled')
            self.Btn_Abort.config(state='disabled')
            
            self.SimStatus.set('Not running: Stopped.')
            
            #Create the .stop
            stop_file = open(self.Manage.STOP, 'w')
            stop_file.close()
            
            self.top.DisplayMessage('  Stop the simulation.', 0)

               
        
    ''' ==================================================================================
    FUNCTION Click_RadioSIM: Change the way the ligand is displayed in Pymol
                             during a Simulation
    ==================================================================================  '''               
    def Click_RadioSIM(self):
        
        if self.SimDefDisplay.get() == 'spheres':
            try:
                cmd.show('spheres', 'TOP_*__ & resn LIG')
                cmd.hide('sticks', 'TOP_*__ & resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
        
        elif self.SimDefDisplay.get() == 'sticks':
            try:
                cmd.hide('spheres', 'TOP_*__ & resn LIG')
                cmd.show('sticks', 'TOP_*__ & resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
        

    ''' ==================================================================================
    FUNCTION Check_MeshSIM: Change the protein display adding or removing the mesh.
    ==================================================================================  '''
    def Check_MeshSIM(self):

        # Display the MESH
        if self.SimMeshDisplay.get():
            try:
                cmd.show('mesh', 'TOP_*__ & ! resn LIG')
                #cmd.set('mesh_color', 'gray30')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
            
        else:   
            # Remove the MESH
            try:
                cmd.hide('mesh', 'TOP_*__ & ! resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
       
    ''' ==================================================================================
    FUNCTION Check_LinesSIM: Change the protein display adding or removing the lines
    ==================================================================================  '''
    def Check_LinesSIM(self):

        # Display the Lines
        if self.SimLinesDisplay.get():
            try:
                cmd.show('lines', 'TOP_*__ & ! resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
            
        else:   
            # Remove the Lines
            try:
                cmd.hide('lines', 'TOP_*__ & ! resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
         
    ''' ==================================================================================
    FUNCTION Check_CartoonSIM: Change the protein display adding or removing the cartoon
    ==================================================================================  '''
    def Check_CartoonSIM(self):

        # Display the Cartoon
        if self.SimCartoonDisplay.get():
            try:
                cmd.show('cartoon', 'TOP_*__ & ! resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
            
        else:   
            # Remove the Cartoon
            try:
                cmd.hide('cartoon', 'TOP_*__ & ! resn LIG')
            except:
                self.top.DisplayMessage("Could not find object to modify", 1)
                
    

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
            #self.ProgressBar.update_idletasks()

        except:
            pass

        return

    ''' ==================================================================================
    FUNCTION LoadMessage: Display the message based on the menu selected
    ================================================================================== '''        
    def LoadMessage(self):
        
        self.top.DisplayMessage('', 0)
        self.top.DisplayMessage('  FlexAID < Simulate > Menu.', 2)
        self.top.DisplayMessage('  INFO:   Press the < Start > button to begin the simulation.', 2)
        self.top.DisplayMessage('           After a simulation is started, you can press the < Pause > or the < Stop > button.', 2)
