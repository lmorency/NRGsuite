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
@title: NewProject - Interface

@summary: This is the interface of the New Project application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec 8, 2010
'''

from Tkinter import *
from pymol import cmd

import General
import tkFileDialog
import tkFont
import re
import os
import time
import Prefs

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayNewProject:
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, parent, inThread, NewProjectPath, HomeDir, RunDir_Path, OSid):      
        
        self.WINDOWWIDTH = 400
        self.WINDOWHEIGHT = 320

        self.OSid = OSid
        self.RunDir_Path = RunDir_Path
        self.HomeDir = HomeDir
        self.ValidateBaseDir()       
        
        self.inThread = inThread       
            
        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''        
        self.NewProject_Path = NewProjectPath
        
        self.NewProjectIsRunning()
               
        #================================================================================== 
        #                 SET the default fonts of the interface
        #================================================================================== 
        FontType = Prefs.GetFontType()
        FontSize = Prefs.GetFontSize()
       
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)        
        self.font_Title_H = tkFont.Font(family=FontType,size=FontSize + 2, weight=tkFont.BOLD)        
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_H = tkFont.Font(family=FontType,size=FontSize + 2)
        self.font_Text_U = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        
        self.Color_Gray = '#EDEDED'
        self.Color_White = '#FFFFFF'
        self.Color_Red = '#FF9999'
        self.Color_Green = '#CCFFCC'
        
        #================================================================================== 
        #                        Initialize the window
        #================================================================================== 
        top = self.top = Toplevel(parent)
        top.title('NRGsuite - Project')
        #top.geometry('')   # Interface DIMENSIONS
        General.CenterWindow(self.top,self.WINDOWWIDTH,self.WINDOWHEIGHT)
        top.maxsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        top.minsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        self.top.protocol('WM_DELETE_WINDOW', self.Btn_Cancel_Clicked)
        
        
        self.NameSpacer = '  '
        self.ActualDirPath = self.HomeDir
        
        self.ProjDirPath = StringVar()
        self.ProjDirPath.set(os.path.join(self.NameSpacer,self.ActualDirPath)) 
        
        self.useDefault = IntVar()
        self.useDefault.set(1)
        
        self.ProjName = StringVar()
        self.ProjName.set('')
        
        self.ClearMsg = False
        self.BtnCreateStatus = 'disabled'       
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================
        self.frame = Frame(top)  # Main frame
        self.frame.pack(expand='true')

        self.Frame_Main() 

    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        self.fTop = Frame(self.frame, relief=RIDGE, border=0, width=400, height=260)        

        fTitle = Frame(self.fTop, relief=RIDGE, border=1, width=400, height=40, bg=self.Color_White)
        
        Title = Label(fTitle, text='Create a new NRG suite Project.', bg=self.Color_White)
        Title.pack(padx=10, side=LEFT, anchor=W)
        Title['font'] = self.font_Title
        
        fTitle.pack(fill=X, expand=True)
        fTitle.pack_propagate(0)
        
        fProjName = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=50)
        
        Name = Label(fProjName, text='Project name:')
        Name.pack(padx=10, side=LEFT, anchor=W)
        Name['font'] = self.font_Text
        
        self.NameVal = Entry(fProjName, width=28, textvariable=self.ProjName, background='white')
        self.NameVal.pack(side=LEFT, anchor=W)
        self.NameVal['font'] = self.font_Text
        self.NameVal.focus_set()
        
        fProjName.pack(fill=X, expand=True, padx=10)
        fProjName.pack_propagate(0)
        
        fProjContent = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=30)
        
        Content = Label(fProjContent, text='Project contents:')
        Content.pack(padx=10, side=LEFT, anchor=S)
        Content['font'] = self.font_Text
        
        fProjContent.pack(fill=X, expand=True, padx=10)
        fProjContent.pack_propagate(0)
        
        fDefault = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=30)
        
        Default = Checkbutton(fDefault, text='  Use default', variable=self.useDefault, command=self.Default_Checked)
        Default.pack(padx=25, side=LEFT, anchor=W)
        Default['font'] = self.font_Text_U
        
        fDefault.pack(fill=X, expand=True, padx=10)
        fDefault.pack_propagate(0)
        
        fDirectory = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=30)
        
        Directory = Label(fDirectory, text='Create project directory under:')
        Directory.pack(padx=10, side=LEFT, anchor=W)
        Directory['font'] = self.font_Text
        
        fDirectory.pack(fill=X, expand=True, padx=10)
        fDirectory.pack_propagate(0)
        
        fDirPath = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=30)
        
        PadDir = Label(fDirPath, text='')
        PadDir.pack(padx=13, side=LEFT, anchor=W)       
        
        self.DirPath = Entry(fDirPath, width=38, textvariable=self.ProjDirPath, background='white', \
                             disabledforeground='black')
        self.DirPath.pack(side=LEFT, anchor=W)
        self.DirPath['font'] = self.font_Text_U
        self.DirPath.config(state='disabled')
        
        self.Btn_Browse = Button(fDirPath, text='Browse', width=7, command=self.Btn_Browse_Clicked)
        self.Btn_Browse.pack(side=RIGHT, anchor=E)
        self.Btn_Browse['font'] = self.font_Text
        self.Btn_Browse.config(state='disabled')

        fDirPath.pack(fill=X, expand=True, padx=10)
        fDirPath.pack_propagate(0)
        
        fButtons = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=50)        
    
        Btn_Cancel = Button(fButtons, text='Cancel', width=7, command=self.Btn_Cancel_Clicked)
        Btn_Cancel.pack(side=RIGHT, anchor=SE, pady=3)
        Btn_Cancel['font'] = self.font_Text
        
        self.Btn_Create = Button(fButtons, text='Create', width=7, command=self.Btn_Create_Clicked)
        self.Btn_Create.pack(side=RIGHT, anchor=SE, padx=4, pady=3)
        self.Btn_Create['font'] = self.font_Text
        self.Btn_Create.config(state='normal')
        
        fButtons.pack(fill=X, expand=True, padx=10)
        fButtons.pack_propagate(0)
        
        self.fTop.pack(fill=X, expand=True)
        self.fTop.pack_propagate(0) 
        
        #==================================================================================
        '''                  --- BOTTOM FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        fBottom = Frame(self.frame, relief=RIDGE, border=0, width=400, height=60)
        
        # Messages Box
        fMsg = Frame(fBottom, border=1, width=380, height=50, relief=SUNKEN)        

        self.TextMessage = Text(fMsg, width=380, height=5, background=self.Color_Gray)
        self.TextMessage.pack(side=LEFT, anchor=S)
        self.TextMessage['font'] = self.font_Text
        self.TextMessage.config(state='disabled')

        fMsg.pack(fill=X, expand=True, side=LEFT, anchor=S, padx=10, pady= 5)
        fMsg.pack_propagate(0)

        fBottom.pack(fill=X, expand=True, side=BOTTOM, anchor=S)
        fBottom.pack_propagate(0)
        
    ''' ==================================================================================
    FUNCTION NewProjectIsRunning: Update or Create the Running File to BLOCK multiple GUI 
    ==================================================================================  '''       
    def NewProjectIsRunning(self):
        
        #Create the .run file
        RunPath = os.path.join(self.RunDir_Path,'.nrun')
        RunFile = open(RunPath, 'w')
        RunFile.write(str(os.getpid()))
        RunFile.close()
        
    ''' ==================================================================================
    FUNCTION Delete_RunFile: Delete the run file.
    ==================================================================================  '''    
    def Delete_RunFile(self):
        
        #Delete the .run file
        RunPath = os.path.join(self.RunDir_Path,'.nrun')
        if (os.path.isfile(RunPath)):
            try:
                os.remove(RunPath)
            except OSError:
                time.sleep(0.1)
                os.remove(RunPath)
        
    ''' ==================================================================================
    FUNCTION Btn_Cancel_Clicked: Cancel the project creation then quit the application.
    ==================================================================================  '''        
    def Btn_Cancel_Clicked(self):
        
        print('   CANCELLED the creation of a New Project...')
        self.inThread.EnableMenu = False
        self.inThread.StopThread = True
        self.inThread.UserPath = ''
        
        self.Delete_RunFile()
        
        self.top.destroy()        
        
    ''' ==================================================================================
    FUNCTION Btn_Create_Clicked: Create the new project then quit the application.
    ==================================================================================  '''        
    def Btn_Create_Clicked(self):
        
        name = self.ProjName.get().strip()
        if name != '':
        
            words = name.split()
            for word in words:
                if not re.match('^[A-Za-z0-9_\-\.]*$', word):
                    print('  ERROR: Invalid name for project.')
                    return
            
            userPath = os.path.join(self.ActualDirPath,name)
            
            if not os.path.isdir(userPath):
                try:
                    os.makedirs(userPath)
                except OSError:
                    print('  ERROR: Could not create the project: \'' + name + '\' under that directory.')
                    return 
                    
                self.Update_ProjectFile(userPath, name)
                
                self.inThread.EnableMenu = True
                self.inThread.StopThread = True
                self.inThread.UserPath = userPath
                self.inThread.ProjectName = name       
                
                self.Delete_RunFile()
                self.top.destroy()

                print('  Successfully created the project: \'' + name + '\'')
    
            else:
                print('  The project: \'' + name + '\' already exists.' )
                
                        
    ''' ==================================================================================
    FUNCTION Btn_Browse_Clicked: Browse to specify the directory to install the project.
    ==================================================================================  '''        
    def Btn_Browse_Clicked(self):

        src = tkFileDialog.askdirectory(initialdir=self.HomeDir, title='Select a directory')

        if len(src) > 0:
            if os.path.isdir(src):
                self.ActualDirPath = 'src'
                self.ProjDirPath.set(os.path.join(self.NameSpacer,src))
        
    ''' ==================================================================================
    FUNCTION Step3_Checked: The check box in the Step3 in Check or NOT.
    ==================================================================================  '''        
    def Default_Checked(self):
        
        if self.useDefault.get() == 0:
            self.Btn_Browse.config(state='normal')
            self.DirPath.config(state='normal')
            self.DirPath['font'] = self.font_Text
        else:
            self.Btn_Browse.config(state='disabled')
            self.DirPath.config(state='disabled')
            self.DirPath['font'] = self.font_Text_U
            self.ActualDirPath = self.HomeDir
            self.ProjDirPath.set(os.path.join(self.NameSpacer,self.ActualDirPath))
            
    ''' ==================================================================================
    FUNCTION Update_ProjectFile: Add the created project to the Project File list.
    ==================================================================================  '''            
    def Update_ProjectFile(self, userPath, ProjName):
        
        MainFilePath = os.path.join(self.HomeDir,'.proj')
        UserFilePath = os.path.join(userPath,'.proj')
        Date = self.Get_Date()
        Seconds = time.time()
        
        LineToAdd = 'Name: ' + ProjName + ' Seconds: ' + str(Seconds) + ' LastUsed: ' + Date + ' Creation: ' + Date + ' Path: ' + userPath + '\n'
        
        NotFound = True
        
        if os.path.isfile(MainFilePath):
            # Add the project path to the .proj file
            file = open(MainFilePath, 'r')
            ProjFile = file.readlines()
            file.close()
            
            for line in ProjFile:
                if LineToAdd == line:
                    NotFound = False
                    break

        if NotFound:
            # Add the project path to the .proj file
            ProjFileMain = open(MainFilePath, 'a')
            ProjFileMain.write(LineToAdd)
            ProjFileMain.close()
        
        # Add the project path to the .proj file
        ProjFileUser = open(UserFilePath, 'w')
        ProjFileUser.write(LineToAdd)
        ProjFileUser.close()
        
    ''' ==================================================================================
    FUNCTION Get_Date: Return the actual date (mm/dd/yyyy - hh:mm:ss).
    ==================================================================================  '''    
    def Get_Date(self):
            
        timeNow = time.localtime()

        Date = ''
        
        month = int(timeNow.tm_mon)
        day = int(timeNow.tm_mday)
        hour = int(timeNow.tm_hour)
        minute = int(timeNow.tm_min)
        second = int(timeNow.tm_sec)
        
        if month < 10:
            Date += '0'
            
        Date += str(month) + '/'
        
        if day < 10:
            Date += '0'
            
        Date += str(day) + '/' + str(timeNow.tm_year) + ' - '

        if hour < 10:
            Date += '0'
            
        Date += str(hour) + ':'
        
        if minute < 10:
            Date += '0'
            
        Date += str(minute) + ':'
            
        if second < 10:
            Date += '0'
            
        Date += str(second) 
               
        return Date
                
    ''' ==================================================================================
    FUNCTION ValidateBaseDir: Create if didnt exist the based path directory.
    ==================================================================================  '''            
    def ValidateBaseDir(self):
        
        self.HomeDir = os.path.join(self.HomeDir,'Documents','NRGsuite')
        
        if not os.path.isdir(self.HomeDir):
            os.makedirs(self.HomeDir)
            
        # Be sure the Default directory exist
        Default = os.path.join(self.HomeDir,'Default')
        if not os.path.isdir(Default):
            os.makedirs(Default)       # ADD the Directory       
            self.Update_ProjectFile(Default, 'Default')     # ADD the .proj Files         
                                
    ''' ==================================================================================
    FUNCTION DisplayMessage: Display the message  
    ==================================================================================  '''    
    def DisplayMessage(self, msg, priority):
        
        NBCHARMAX = 45
        self.TextMessage.config(state='normal')
 
        lineNo = '1'
        insNo = lineNo + '.0'      

        if msg != '':
            # Verify if the message need to be split
            NbChar = len(msg)

            #self.TextMessage.config(font='red')
            self.TextMessage.insert(INSERT, '\n' + msg)
            
            if priority == 1:
                self.TextMessage.tag_add('warn', '0.0', '2.' + str(NbChar))
                self.TextMessage.tag_config('warn', foreground='red')
            elif priority == 2:
                self.TextMessage.tag_add('notice', '0.0', '2.' + str(NbChar))
                self.TextMessage.tag_config('notice', foreground='blue')   
                
            self.TextMessage.yview(INSERT)

        else:
            # Skip a line...
            self.TextMessage.delete(1.0,END)   
        
        self.TextMessage.config(state='disabled')
            


