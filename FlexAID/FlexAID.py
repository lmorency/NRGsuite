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
@title: FlexAID - Interface

@summary: This is the interface of FlexAID application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Aug. 25, 2010
'''

from Tkinter import *

import os, sys
import pickle
import time
import tkFont
import tkFileDialog

import Prefs
import Color
import General

import Base
import IOFile
import Config1
import Config2
import Config3
import GAParam
import Simulate

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayFlexAID(Base.Base):
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, top, ActiveWizard, Project_Dir, Install_Dir, AlreadyRunning_Dir, OSid, PyMOL):
        
        #print("New instance of FlexAID")
        self.PyMOL = PyMOL

        self.Name = 'FlexAID'
        self.RunFile = '.frun'

        self.WINDOWWIDTH = 700
        self.WINDOWHEIGHT = 600

        #================================================================================== 
        ''' USER SAVE PATH DIRECTORY  '''
        #================================================================================== 
        self.AlreadyRunning_Dir = AlreadyRunning_Dir

        self.Install_Dir = Install_Dir
        self.FlexAIDInstall_Dir = os.path.join(self.Install_Dir,'FlexAID')

        self.OSid = OSid
        if self.OSid == 'WIN':
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID.exe')
        else:
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID')

        self.Project_Dir = Project_Dir
        self.FlexAIDProject_Dir = os.path.join(self.Project_Dir,'FlexAID')
        self.GetCleftProject_Dir = os.path.join(self.Project_Dir,'GetCleft')
        self.GetCleftSaveProject_Dir = os.path.join(self.GetCleftProject_Dir,'Save')

        self.CleftProject_Dir = os.path.join(self.Project_Dir,'Cleft')
        self.TargetProject_Dir = os.path.join(self.Project_Dir,'Target')

        self.FlexAIDLigandProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Ligand')
        self.FlexAIDSimulationProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Simulation')
        self.FlexAIDSessionProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Session')
        self.FlexAIDResultsProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Results')
        self.FlexAIDBindingSiteProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Binding_Site')
        self.FlexAIDTargetFlexProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Target_Flexibility')

        self.LOGFILE = os.path.join(self.FlexAIDProject_Dir,'logfile.txt')
        #self.LOGCOMMANDS = os.path.join(self.UserFlexAID,'logcmd.pml')
        
        self.ValidateFolders()      

        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''
        #================================================================================== 
        self.OSid = OSid

        self.BindingSiteDisplay = 'BINDING_SITE_AREA'
        self.SphereDisplay = 'SPHERE_AREA'
                
        self.AppIsRunning()
        
        #self.MsgLineCounter = 2
        
        self.Color_Green = '#CCFFCC'
        self.Color_Grey = '#EDEDED'
        self.Color_Blue = '#6699FF'
        self.Color_Red = '#FF9999'
        self.Color_White = '#FFFFFF'
        self.Color_Black = 'black'                

        #print("FlexAID: initializing window")
        # Initialize the window
        self.top = top        
        self.top.title(self.Name)

        #print("FlexAID: center window")
        General.CenterWindow(self.top,self.WINDOWWIDTH,self.WINDOWHEIGHT)


        #self.top.geometry()   # Interface DIMENSIONS
        #self.top.maxsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        #self.top.minsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        self.top.protocol('WM_DELETE_WINDOW', self.Quit)

        #================================================================================== 
        #                 SET the default fonts of the interface
        #==================================================================================
        FontType = Prefs.GetFontType()
        FontSize = Prefs.GetFontSize()
        
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)        
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_I = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        self.font_Text_U = tkFont.Font(family=FontType,size=FontSize, underline=True)       
        self.Img_Msg = list()        
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================
                
        self.frame = Frame(top)
        self.frame.pack(expand=True)
        self.curCursor = self.frame['cursor']
        
        #print "Main frame generated"

        self.ActiveFrame = None
        self.ActiveWizard = ActiveWizard
        self.WizardError = False
        self.WizardResult = 0

        self.ProcessRunning = False
        self.Run = None

        self.Frame_Main()
        
        # Build class objects of each tab
        self.IOFile = IOFile.IOFile(self, self.PyMOL, self.Btn_IOFiles, 'IOFile', IOFile.IOFileVars())
        self.Config1 = Config1.Config1(self, self.PyMOL, self.Btn_Config1, 'Config1', Config1.Config1Vars())
        self.Config2 = Config2.Config2(self, self.PyMOL, self.Btn_Config2, 'Config2', Config2.Config2Vars())
        self.Config3 = Config3.Config3(self, self.PyMOL, self.Btn_Config3, 'Config3', Config3.Config3Vars())
        self.GAParam = GAParam.GAParam(self, self.PyMOL, self.Btn_GAParam, 'GAParam', GAParam.GAParamVars())
        self.Simulate = Simulate.Simulate(self, self.PyMOL, self.Btn_Simulate, 'Simulate', Simulate.SimulateVars())

        self.MakeMenuBar()

        self.listTabs = [self.IOFile, self.Config1, self.Config2, self.Config3, self.GAParam, self.Simulate]
        self.listBtnTabs = [self.Btn_IOFiles, self.Btn_Config1, self.Btn_Config2, self.Btn_Config3, self.Btn_GAParam, self.Btn_Simulate]

        #print "Created instances of all tab classes"

        # default Tab (IOFile)
        self.Btn_IOFiles_Clicked()

        # By default hide advanced tabs
        self.bAdvancedView = False
        self.Btn_Toggle_AdvView()


    #=====================================================================================
    '''                     --- FRAME DISPLAY SETTINGS ---                             '''
    #=====================================================================================
    
    ''' ==================================================================================
    FUNCTION Btn_Toggle_AdvView: Hides/shows the advanced tabs (scoring + ga)
    ==================================================================================  '''  
    def Btn_Toggle_AdvView(self):
    
        if not self.bAdvancedView:
            self.Btn_Config3.pack_forget()
            self.Btn_GAParam.pack_forget()
            self.bAdvancedView = True
            
            if self.ActiveFrame == self.Config3 or self.ActiveFrame == self.GAParam:
                self.Btn_IOFiles_Clicked()

            self.menubar.entryconfig(3,label="Show")
            
        else:
            self.Btn_Simulate.pack_forget()
            self.Btn_Config3.pack(side=LEFT, fill=BOTH, expand=True)
            self.Btn_GAParam.pack(side=LEFT, fill=BOTH, expand=True)
            self.Btn_Simulate.pack(side=LEFT, fill=BOTH, expand=True)
            self.bAdvancedView = False
            
            self.menubar.entryconfig(3,label="Hide")
            
    
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''  
    def Frame_Main(self):
   
        #==================================================================================
        '''                  --- TOP MENU OF THE INTERFACE ---                          '''
        #==================================================================================

        self.fTop = Frame(self.frame, relief=RIDGE, border=4, height=50)
        self.fTop.pack(fill=BOTH, expand=True)#, padx=10, pady=10, ipady=10, ipadx=10, side=TOP)
        self.fTop.pack_propagate(0)

        self.Btn_IOFiles = Button(self.fTop, text='Input Files', bg=self.Color_White, command=self.Btn_IOFiles_Clicked, font=self.font_Text)
        self.Btn_IOFiles.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_IOFiles.config(state='normal')

        self.Btn_Config1 = Button(self.fTop, text='Target Cfg', bg=self.Color_Grey, command=self.Btn_Config1_Clicked, font=self.font_Text)
        self.Btn_Config1.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config1.config(state='disabled')

        self.Btn_Config2 = Button(self.fTop, text='Ligand Cfg', bg=self.Color_Grey, command=self.Btn_Config2_Clicked, font=self.font_Text)
        self.Btn_Config2.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config2.config(state='disabled')

        self.Btn_Config3 = Button(self.fTop, text='Scoring Cfg', bg=self.Color_Grey, command=self.Btn_Config3_Clicked, font=self.font_Text)
        self.Btn_Config3.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config3.config(state='disabled')
        
        self.Btn_GAParam = Button(self.fTop, text='GA Param', bg=self.Color_Grey, command=self.Btn_GAParam_Clicked, font=self.font_Text)
        self.Btn_GAParam.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_GAParam.config(state='disabled')    

        self.Btn_Simulate = Button(self.fTop, text='Simulate', bg=self.Color_Grey, command=self.Btn_Simulate_Clicked, font=self.font_Text)
        self.Btn_Simulate.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Simulate.config(state='disabled') 
                
        #==================================================================================
        '''                MIDDLE DISPLAY SECTION OF THE INTERFACE                      '''
        #==================================================================================
           
        self.fMiddle = Frame(self.frame)#, height=355)
                
        self.fMiddle.pack(fill=X, expand=True, padx=10, side=TOP)
        #self.fMiddle.config(takefocus=1)
        #self.fMiddle.bind('<FocusIn>', lambda e: None)
        #self.fMiddle.pack_propagate(0)
        
        #==================================================================================
        '''                 BOTTOM DISPLAY SECTION OF THE INTERFACE                     '''
        #==================================================================================
 
        fBottom = Frame(self.frame, border=1, relief=SUNKEN)
        fBottom.pack(fill=BOTH, expand=True, padx=10, pady=10, ipadx=10, ipady=10, side=TOP)

        fBottomRight = Frame(fBottom)
        fBottomRight.pack(fill=Y, side=RIGHT)

        Btn_Default = Button(fBottomRight, text='Default', width=15, command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Default.pack(side=TOP, fill=X)

        #Btn_SaveDefault = Button(fBottomRight, text='Save as default', command=self.Btn_SaveDefault_Clicked, font=self.font_Text)
        #Btn_SaveDefault.pack(side=TOP, fill=X)

        #Btn_Restore = Button(fBottomRight, text='Restore', command=self.Btn_Restore_Clicked, font=self.font_Text)
        #Btn_Restore.pack(side=TOP, fill=X)

        Btn_Quit = Button(fBottomRight, text='Close', width=15, command=self.Quit, font=self.font_Text)
        Btn_Quit.pack(side=BOTTOM, fill=X)

        fBottomLeft = Frame(fBottom)
        fBottomLeft.pack(side=LEFT, fill=Y, ipadx=20, ipady=20)

        scrollBar = Scrollbar(fBottomLeft)
        scrollBar.pack(side=RIGHT, fill=Y)

        self.TextMessage = Text(fBottomLeft, border=1, background=self.Color_Grey, font=self.font_Text)
        self.TextMessage.pack(side=RIGHT, fill=BOTH, expand=True)

        scrollBar.config(command=self.TextMessage.yview)
        self.TextMessage.config(state='disabled', yscrollcommand=scrollBar.set)                                       

        #self.Btn_Dummy = Button(fBottomRight, text='Close', font=self.font_Text)
        #self.Btn_Dummy.bind('<Button-1>', lambda event, arg=self.ActiveFrame: self.SwitchTab(event,arg))
    
    ''' ==================================================================================
    FUNCTION MakeMenuBar: Builds the menu on the upper left corner    
    ==================================================================================  '''        
    def MakeMenuBar(self):
        
        self.menubar = Menu(self.top)
        
        loadmenu = Menu(self.menubar, tearoff=0)
        loadmenu.add_command(label="Load Session", command=self.Btn_Load_Session)
        loadmenu.add_separator()
        loadmenu.add_command(label="Load Results", command=self.Btn_Load_Results)
        self.menubar.add_cascade(label="Load", menu=loadmenu)

        savemenu = Menu(self.menubar, tearoff=0)        
        savemenu.add_command(label="Save Session", command=self.Btn_Save_Session)
        savemenu.add_separator()
        savemenu.add_command(label="Save Results", command=self.Btn_Save_Results)
        self.menubar.add_cascade(label="Save", menu=savemenu)

        viewmenu = Menu(self.menubar, tearoff=0)
        viewmenu.add_command(label="Advanced view", command=self.Btn_Toggle_AdvView)
        self.menubar.add_cascade(label="View", menu=viewmenu)
        
        self.top.config(menu=self.menubar)
        
    ''' ==================================================================================
    FUNCTION Btn_Load_Session: Loads a previously saved session
    ==================================================================================  '''        
    def Btn_Load_Session(self):

        if self.Run is None and not self.ProcessRunning:
            if self.ActiveWizard is None:

                LoadFile = tkFileDialog.askopenfilename(initialdir=self.FlexAIDSessionProject_Dir,
                                                        filetypes=[('NRG FlexAID Session','*.nrgfs')],
                                                        title='Select the Session to load')

                if len(LoadFile) > 0:
                    
                    self.Btn_IOFiles_Clicked()
                    
                    LoadFile = os.path.normpath(LoadFile)
                    
                    try:
                        in_ = open(LoadFile, 'r')
                        for Tab in self.listTabs:
                            try:
                                Tab.Vars = pickle.load(in_)
                                Tab.Vars.refresh()
                                Tab.Load_Session()
                            except:
                                pass
                        in_.close()
                    except:
                        self.DisplayMessage("  ERROR: Could not properly load the session", 2)
                        self.DisplayMessage("  Unexpected error: " + str(sys.exc_info()), 2)

            else:
                self.DisplayMessage("  Cannot save session while a wizard is active", 2)            
        else:
            self.DisplayMessage("  Cannot save session while a process is active", 2)
        
    ''' ==================================================================================
    FUNCTION Btn_Load_Results: Loads a previously saved results
    ==================================================================================  '''        
    def Btn_Load_Results(self):

        if self.Run is None and not self.ProcessRunning:

            if self.ActiveWizard is None:

                LoadFile = tkFileDialog.askopenfilename(initialdir=self.FlexAIDResultsProject_Dir,
                                                        filetypes=[('NRG FlexAID Results','*.nrgfr')],
                                                        title='Select the Results to load')

                if len(LoadFile) > 0:
                                        
                    LoadFile = os.path.normpath(LoadFile)
                    
            else:
                self.DisplayMessage("  Cannot save session while a wizard is active", 2)            
        else:
            self.DisplayMessage("  Cannot save session while a process is active", 2)
        
    ''' ==================================================================================
    FUNCTION Btn_Save_Results: Saves the current results
    ==================================================================================  '''        
    def Btn_Save_Results(self):

        if self.Run is None and not self.ProcessRunning:
        
            if self.ActiveWizard is None:
                
                SaveFile = tkFileDialog.asksaveasfilename(initialdir=self.FlexAIDResultsProject_Dir,
                                          title='Save the Results file', initialfile='default_results',
                                          filetypes=[('NRG FlexAID Results','*.nrgfr')])
            
                if len(SaveFile) > 0:

                    SaveFile = os.path.normpath(SaveFile)
            
            else:
                self.DisplayMessage("  Cannot save results while a wizard is active", 2)
        else:
            self.DisplayMessage("  Cannot save results while a process is active", 2)
            

    ''' ==================================================================================
    FUNCTION Btn_Save_Session: Saves the current session
    ==================================================================================  '''        
    def Btn_Save_Session(self):

        if self.Run is None and not self.ProcessRunning:
        
            if self.ActiveWizard is None:
                
                SaveFile = tkFileDialog.asksaveasfilename(initialdir=self.FlexAIDSessionProject_Dir,
                                          title='Save the Session file', initialfile='default_session',
                                          filetypes=[('NRG FlexAID Session','*.nrgfs')])
            
                if len(SaveFile) > 0:

                    SaveFile = os.path.normpath(SaveFile)
                    
                    if SaveFile.find('.nrgfs') == -1:
                        SaveFile = SaveFile + '.nrgfs'

                    try:
                        out = open(SaveFile, 'w')
                        for Tab in self.listTabs:
                            try:
                                pickle.dump(Tab.Vars, out)
                            except:
                                pass
                        out.close()
                    except:
                        self.DisplayMessage("  ERROR: Could not properly save the session:", 2)
                        self.DisplayMessage("  Unexpected error: " + str(sys.exc_info()), 2)
            else:
                self.DisplayMessage("  Cannot save session while a wizard is active", 2)
        else:
            self.DisplayMessage("  Cannot save session while a process is active", 2)
    
    ''' ==================================================================================
    FUNCTION Btn_*_Clicked: Display the Tab options menu
    ================================================================================== '''    
    def Btn_IOFiles_Clicked(self):
    
        self.SetActiveFrame(self.IOFile)        

    def Btn_Config1_Clicked(self):
        
        self.SetActiveFrame(self.Config1)

    def Btn_Config2_Clicked(self):
        
        self.SetActiveFrame(self.Config2)

    def Btn_Config3_Clicked(self):
        
        self.SetActiveFrame(self.Config3)

    def Btn_GAParam_Clicked(self):
        
        self.SetActiveFrame(self.GAParam)

    def Btn_Simulate_Clicked(self):
        
        self.SetActiveFrame(self.Simulate)

    ''' ==================================================================================
    FUNCTION Go_Step1: Enables/Disables buttons for step 1
    ================================================================================== '''    
    def Go_Step1(self):

        #print "Setting Tab buttons to Step 1"        
        for Btn in self.listBtnTabs:
            Btn.config(state='disabled',bg=self.Color_Grey)

        self.Btn_IOFiles.config(state='normal')

    ''' ==================================================================================
    FUNCTION Go_Step2: Enables/Disables buttons for step 2
    ================================================================================== '''    
    def Go_Step2(self):
        
        #print "Setting Tab buttons to Step 2"        
        for Btn in self.listBtnTabs:
            Btn.config(state='normal',bg=self.Color_White)
        
    ''' ==================================================================================
    FUNCTION Reset_Step2: Reset ALL the parameters
    ==================================================================================  '''    
    def Reset_Step2(self):
        
        for Tab in self.listTabs:
            if Tab != self.IOFile:
                Tab.Init_Vars()

    #=====================================================================================
    '''                            --- BUTTONS EVENT ---                               '''
    #=====================================================================================     
        
    ''' ==================================================================================
    FUNCTION ValidateResiduValue: Validate the residue entered 
    ================================================================================== '''
    def ValidateResiduValue(self, event):        
        
        term = self.ResiduValue.get().upper()
        self.ResiduValue.set(term)
        
        NbChar = len(term)
        self.ValidResn = False
        
        if NbChar > 0:
            # Be sure the number of characters do not exceed 10
            if NbChar > 10:
                term = term[0:10]
                self.ResiduValue.set(term)
        
            Notfound = True            
            
            self.NoResn = 0 
            for elem in self.listResidu:                
                if elem.startswith(term):                    
                    if (NbChar == len(elem)) or (elem.startswith(term + '-')):
                        self.ValidResn = True
                    Notfound = False
                    break
                self.NoResn += 1
                
            if Notfound:
                self.EntryResidu.config(bg=self.Color_Red)
                self.Btn_AddResidu.config(state='disabled')
                
            else:
                if self.ValidResn:
                    self.EntryResidu.config(bg=self.Color_Green)
                    self.Btn_AddResidu.config(state='normal')
                else:
                    self.EntryResidu.config(bg=self.Color_White)
                    self.Btn_AddResidu.config(state='disabled')
        else:
            self.EntryResidu.config(bg=self.Color_White)
            
    ''' ==================================================================================
    FUNCTION ValidateFolders: Be sure the folders Exists 
    ==================================================================================  '''    
    def ValidateFolders(self):
       
        if not os.path.isdir(self.FlexAIDProject_Dir):
            os.makedirs(self.FlexAIDProject_Dir)
            
        if not os.path.isdir(self.FlexAIDLigandProject_Dir):
            os.makedirs(self.FlexAIDLigandProject_Dir)
            
        if not os.path.isdir(self.TargetProject_Dir):
            os.makedirs(self.TargetProject_Dir)
            
        if not os.path.isdir(self.CleftProject_Dir):
            os.makedirs(self.CleftProject_Dir)

        if not os.path.isdir(self.FlexAIDSimulationProject_Dir):
            os.makedirs(self.FlexAIDSimulationProject_Dir)
            
        if not os.path.isdir(self.FlexAIDSessionProject_Dir):
            os.makedirs(self.FlexAIDSessionProject_Dir)

        if not os.path.isdir(self.FlexAIDResultsProject_Dir):
            os.makedirs(self.FlexAIDResultsProject_Dir)

        if not os.path.isdir(self.FlexAIDBindingSiteProject_Dir):
            os.makedirs(self.FlexAIDBindingSiteProject_Dir)
            
        if not os.path.isdir(self.FlexAIDTargetFlexProject_Dir):
            os.makedirs(self.FlexAIDTargetFlexProject_Dir)

