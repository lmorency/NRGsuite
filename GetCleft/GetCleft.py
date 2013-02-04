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
@title: GetCleft - Interface

@summary: This is the interface of GetCleft application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Oct. 19, 2010
'''

from Tkinter import *

import os, sys
import pickle
import time
import tkFont
import tkFileDialog
import tkMessageBox

import Prefs
import Color
import General

import Base
import CleftObj
import ManageFiles2
import Default
import CropCleft
import Volume

if __debug__:
	from pymol import cmd
	from pymol.cgo import *
	from pymol.vfont import plain

	import General_cmd

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayGetCleft(Base.Base):
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, top, ActiveWizard, Project_Dir, Install_Dir, AlreadyRunning_Dir, OSid, PyMOL):
        
        self.PyMOL = PyMOL
        if self.PyMOL:
            from pymol import cmd
            import General_cmd

        self.Name = 'GetCleft'
        self.RunFile = '.grun'
        self.OSid = OSid
        
        self.WINDOWWIDTH = 500
        self.WINDOWHEIGHT = 550

        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''
        self.AlreadyRunning_Dir = AlreadyRunning_Dir
        
        self.Install_Dir = Install_Dir
        self.Project_Dir = Project_Dir

        self.GetCleftInstall_Dir = os.path.join(self.Install_Dir,'GetCleft')

        if self.OSid == 'WIN':
            self.GetCleftExecutable = os.path.join(self.GetCleftInstall_Dir,'WRK','GetCleft.exe')
            self.VolumeExecutable = os.path.join(self.GetCleftInstall_Dir,'WRK','volume_calc.exe')
        else:
            self.GetCleftExecutable = os.path.join(self.GetCleftInstall_Dir,'WRK','GetCleft')
            self.VolumeExecutable = os.path.join(self.GetCleftInstall_Dir,'WRK','volume_calc')

        self.Project_Dir = Project_Dir
        self.GetCleftProject_Dir = os.path.join(self.Project_Dir,'GetCleft')
        self.TargetProject_Dir = os.path.join(self.Project_Dir,'Target')
        self.CleftProject_Dir = os.path.join(self.Project_Dir,'Cleft')

        self.GetCleftSaveProject_Dir = os.path.join(self.GetCleftProject_Dir,'Save')
        self.GetCleftTempProject_Dir = os.path.join(self.GetCleftProject_Dir,'Temp')

        #self.GetCleftTempCleftsProject_Dir = os.path.join(self.GetCleftTempProject_Dir,'Clefts')
        #self.GetCleftTempAtomsProject_Dir = os.path.join(self.GetCleftTempProject_Dir,'Atoms')

        # Create the folders if the arent exist
        self.ValidateFolders()
               
        self.AppIsRunning()

        #self.MsgLineCounter = 2
        
        self.Color_Green = '#CCFFCC'
        self.Color_Grey = '#EDEDED'
        self.Color_Blue = '#6699FF'
        self.Color_Red = '#FF9999'
        self.Color_White = '#FFFFFF'
        self.Color_Black = 'black'

        self.top = top
        self.top.title(self.Name)

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
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================        
        self.fMain = Frame(self.top)  # Main frame
        self.fMain.pack(expand=True)

        self.ActiveFrame = None
        self.ActiveWizard = ActiveWizard
        self.WizardError = False
        self.WizardResult = 0

        self.ProcessRunning = False
        self.Run = None

        # flag for saving unsaved clefts
        self.CopySession = True
        
        # flag for saving unsaved calculated cleft volumes
        self.EditSession = True

        self.Frame_Main()

        # Build class objects of each tab
        self.Default = Default.Default(self, self.PyMOL, self.Btn_Config, 'Default', None)
        self.Manage = ManageFiles2.Manage(self)
        self.Crop = CropCleft.CropCleft(self, self.PyMOL, self.Btn_CropCleft, 'Partition', None)
        self.Volume = Volume.EstimateVolume(self, self.PyMOL, self.Btn_Volume, 'Volume', None)
        
        self.MakeMenuBar()

        self.listBtnTabs = [ self.Btn_Config, self.Btn_Volume, self.Btn_CropCleft ]
        self.listTabs = [ self.Default, self.Crop, self.Volume ]

        # Default view
        self.Btn_Config_Clicked()

        # Remove all temporary clefts in Temporary Dir
        self.Manage.Clean()

    ''' ==================================================================================
    FUNCTION MakeMenuBar: Builds the menu on the upper left corner    
    ==================================================================================  '''        
    def MakeMenuBar(self):
        
        menubar = Menu(self.top)
        
        loadmenu = Menu(menubar, tearoff=0)
        loadmenu.add_command(label="Load Clefts", command=self.Default.Btn_Load_Clefts)
        menubar.add_cascade(label="Load", menu=loadmenu)

        savemenu = Menu(menubar, tearoff=0)        
        savemenu.add_command(label="Save Clefts", command=self.Default.Btn_Save_Clefts)
        menubar.add_cascade(label="Save", menu=savemenu)
        
        self.top.config(menu=menubar)
        
        
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP MENU OF THE INTERFACE ---                          '''
        #==================================================================================        
        fTop = Frame(self.fMain, relief=RIDGE, border=4, height=50)
        fTop.pack(fill=BOTH, expand=True)#, padx=10, pady=10, ipady=10, ipadx=10, side=TOP)
        fTop.pack_propagate(0)
        
        self.Btn_Config = Button(fTop, text='Generate', bg=self.Color_White, command=self.Btn_Config_Clicked, font=self.font_Text)
        self.Btn_Config.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config.config(state='normal')

        self.Btn_CropCleft = Button(fTop, text='Partition', bg=self.Color_Grey, command=self.Btn_CropCleft_Clicked, font=self.font_Text)
        self.Btn_CropCleft.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_CropCleft.config(state='disabled')

        self.Btn_Volume = Button(fTop, text='Volume', bg=self.Color_Grey, command=self.Btn_Volume_Clicked, font=self.font_Text)
        self.Btn_Volume.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Volume.config(state='disabled')

        #==================================================================================
        '''                  --- MIDDLE FRAME OF THE INTERFACE ---                      '''
        #==================================================================================
        
        self.fMiddle = Frame(self.fMain, relief=RIDGE)
        self.fMiddle.pack(fill=X, expand=True, padx=10, pady=10)
        
        #==================================================================================
        '''                 BOTTOM DISPLAY SECTION OF THE INTERFACE                     '''
        #==================================================================================
 
        fBottom = Frame(self.fMain, border=1, relief=SUNKEN)
        fBottom.pack(fill=BOTH, expand=True, padx=10, pady=5, ipadx=10, ipady=5, side=TOP)

        fBottomRight = Frame(fBottom)
        fBottomRight.pack(fill=Y, side=RIGHT)

        Btn_Default = Button(fBottomRight, text='Default', command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Default.pack(side=TOP, fill=X)

        #Btn_SaveDefault = Button(fBottomRight, text='Save as default', command=self.Btn_SaveDefault_Clicked, font=self.font_Text)
        #Btn_SaveDefault.pack(side=TOP, fill=X)

        #Btn_Restore = Button(fBottomRight, text='Restore', command=self.Btn_Restore_Clicked, font=self.font_Text)
        #Btn_Restore.pack(side=TOP, fill=X)

        Btn_Quit = Button(fBottomRight, text='Close', command=self.Quit, font=self.font_Text)
        Btn_Quit.pack(side=BOTTOM, fill=X)

        fBottomLeft = Frame(fBottom)
        fBottomLeft.pack(side=LEFT, fill=Y)
        
        #self.MessageBox = MessageBox.scrollableContainer(fBottomLeft, bd=2, bg="black")
        #self.MessageBox.pack(fill=BOTH, expand=True)

        scrollBar = Scrollbar(fBottomLeft)
        scrollBar.pack(side=RIGHT, fill=Y)

        self.TextMessage = Text(fBottomLeft, border=1, background=self.Color_Grey, font=self.font_Text)
        self.TextMessage.pack(side=RIGHT, fill=BOTH, expand=True)

        scrollBar.config(command=self.TextMessage.yview)
        self.TextMessage.config(state='disabled', yscrollcommand=scrollBar.set)
    
    ''' ==================================================================================
    FUNCTION Btn_Config_Clicked: Default configuration view
    ==================================================================================  '''                 
    def Btn_Config_Clicked(self):
        
        self.SetActiveFrame(self.Default)

    ''' ==================================================================================
    FUNCTION Btn_Volume_Clicked: Estimates the volume of the loaded Clefts
    ==================================================================================  '''                 
    def Btn_Volume_Clicked(self):
        
        self.SetActiveFrame(self.Volume)
        
    ''' ==================================================================================
    FUNCTION Btn_CropCleft_Clicked: Opens the crop cleft menu 
    ==================================================================================  '''                 
    def Btn_CropCleft_Clicked(self):
        
        self.SetActiveFrame(self.Crop)

    ''' ==================================================================================
    FUNCTION Before_Quit: Execute tasks before exitting the application
    ==================================================================================  '''
    def Before_Quit(self):
    
        if not self.CopySession:
            if tkMessageBox.askquestion("Question", message="One or more cleft(s) are unsaved. Would you like to save them before leaving?",
                                        icon='warning') == 'yes':
                self.Default.Btn_Save_Clefts()
    
    ''' ==================================================================================
    FUNCTION Go_Step1: Enables/Disables buttons for step 1
    ================================================================================== '''    
    def Go_Step1(self):

        #print "Setting Tab buttons to Step 1"        
        for Btn in self.listBtnTabs:
            Btn.config(state='disabled', bg=self.Color_Grey)

        self.Btn_Config.config(state='normal')
        self.Btn_Config.config(bg=self.Color_Blue)

    ''' ==================================================================================
    FUNCTION Go_Step2: Enables/Disables buttons for step 2
    ================================================================================== '''    
    def Go_Step2(self):
        
        #print "Setting Tab buttons to Step 2"
        for Btn in self.listBtnTabs:
            Btn.config(state='normal',bg=self.Color_White)

        self.Btn_Config.config(bg=self.Color_Blue)

    ''' ==================================================================================
    FUNCTION ValidateFolders: Be sure the folders Exists 
    ==================================================================================  '''    
    def ValidateFolders(self):
        
        if not os.path.isdir(self.GetCleftProject_Dir):
            os.makedirs(self.GetCleftProject_Dir)
            
        if not os.path.isdir(self.TargetProject_Dir):
            os.makedirs(self.TargetProject_Dir)

        if not os.path.isdir(self.CleftProject_Dir):
            os.makedirs(self.CleftProject_Dir)
        
        if not os.path.isdir(self.GetCleftSaveProject_Dir):
            os.makedirs(self.GetCleftSaveProject_Dir)

        if not os.path.isdir(self.GetCleftTempProject_Dir):
            os.makedirs(self.GetCleftTempProject_Dir)
    