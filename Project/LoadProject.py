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
@title: LoadProject - Interface

@summary: This is the interface of the Load Project application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec 10, 2010
'''

from Tkinter import *
from pymol import cmd

import General
import MultiList
import tkFileDialog
import tkFont
import time
import os
import Prefs

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayLoadProject:
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, parent, inThread, LoadProjectPath, HomeDir, RunDir_Path, OSid):      
        
        self.WINDOWHEIGHT = 350
        self.WINDOWWIDTH = 460

        self.OSid = OSid
        self.RunDir_Path = RunDir_Path
        self.HomeDir = HomeDir
        self.ValidateBaseDir()

        self.inThread = inThread       
            
        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''        
        self.LoadProject_Path = LoadProjectPath
        
        self.LoadProjectIsRunning()
               
        #================================================================================== 
        #                 SET the default fonts of the interfaceself.current = None

        #==================================================================================
        FontType = Prefs.GetFontType()
        FontSize = Prefs.GetFontSize()
       
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)        
        self.font_Title_H = tkFont.Font(family=FontType,size=FontSize + 2, weight=tkFont.BOLD)        
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_H = tkFont.Font(family=FontType,size=FontSize + 2)
        self.font_Text_I = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        
        self.Color_Gray = '#EDEDED'
        self.Color_White = '#FFFFFF'
        self.Color_Red = '#FF9999'
        self.Color_Green = '#CCFFCC'
        self.Color_Blue = '#4477CC'
        
        #================================================================================== 
        #                        Initialize the window
        #================================================================================== 
        top = self.top = Toplevel(parent)
        top.title('NRG suite - Project')

        General.CenterWindow(self.top,self.WINDOWWIDTH,self.WINDOWHEIGHT)
        #top.geometry('460x350')   # Interface DIMENSIONS
        
        top.maxsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        top.minsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        self.top.protocol('WM_DELETE_WINDOW', self.Btn_Cancel_Clicked)
        
        self.NameSpacer = '  '
        self.ActualDirPath = ''
        self.ActualFileName = ''
        self.ActualKey = ''
       
        self.useDefault = IntVar()
        self.useDefault.set(1)
        
        self.ProjectDir = StringVar()
        self.ProjectDir.set('')
        
        self.dictProjectList = {}
        self.LoadProjectList()
        
        self.ClearMsg = False
        self.BtnCreateStatus = 'disabled'
        self.current = None
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================
        self.frame = Frame(top)  # Main frame
        self.frame.pack(expand='true')

        self.Frame_Main() 
        self.Fill_Table()

        self.Trace()

    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        fTitle = Frame(self.frame, relief=RIDGE, border=1, width=460, height=40, bg=self.Color_White)
        fTitle.pack(fill=X, expand=True)
        fTitle.pack_propagate(0)        
        
        Title = Label(fTitle, text='Load an existing NRG suite project.', bg=self.Color_White, font=self.font_Title)
        Title.pack(padx=10, side=LEFT, anchor=W)        
        
        self.fTop = Frame(self.frame, relief=RIDGE, border=0, width=460, height=240)        
        self.fTop.pack(fill=X, expand=True, padx=10)
        self.fTop.pack_propagate(0) 
        
        fLblList = Frame(self.fTop, relief=RIDGE, border=0, width=460, height=30)
        
        LblList = Label(fLblList, text='Project(s) list:', font=self.font_Text)
        LblList.pack(side=LEFT, anchor=S, pady=2)
        
        fLblList.pack(fill=X, expand=True)
        fLblList.pack_propagate(0)
        
        #==================================================================================
        '''                       --- LISTBOX FRAME SECTION ---                         '''
        #==================================================================================

        fProjList = Frame(self.fTop, relief=RIDGE, border=1, width=460, height=110)
        fProjList.pack(fill=X, expand=True)
        fProjList.pack_propagate(0)

        self.Table = MultiList.Table(fProjList, 3,
                                   [ 'Name', 'Last used', 'Creation date' ],
                                   [ 140, 140, 140 ],
                                   [ 1, 1, 1 ],
                                   [ True, True, True ],
                                   self.font_Text,
                                   self.Color_Blue)
        self.Table.Draw()

        # References to StringVar'
        self.ProjectName = self.Table.Columns['Name']['StringVar']
        self.LastUsed = self.Table.Columns['Last used']['StringVar']
        self.CreationDate = self.Table.Columns['Creation date']['StringVar']

        #==================================================================================
        '''              --- NAME & PATH SELECTED - FRAME SECTION ---                   '''
        #==================================================================================        
        fProjName = Frame(self.fTop, relief=RIDGE, border=0, width=460, height=35)
        
        fProjName_L = Frame(fProjName, relief=RIDGE, border=0, width=100, height=35)
        fProjName_R = Frame(fProjName, relief=RIDGE, border=0, width=360, height=35)
                
        Name = Label(fProjName_L, text='Name:', font=self.font_Text)
        Name.pack(side=RIGHT, anchor=S, pady=4, padx=8)
        
        fProjName_L.pack(fill=X, expand=True, side=LEFT)
        fProjName_L.pack_propagate(0)
        
        self.NameVal = Entry(fProjName_R, width=20, textvariable=self.ProjectName, background='white', \
                             disabledbackground='white', disabledforeground='black', font=self.font_Text)
        self.NameVal.pack(side=LEFT, anchor=S, padx=4, pady=4)
        self.NameVal.config(state='disabled')
        
        fProjName_R.pack(fill=X, expand=True, side=LEFT)
        fProjName_R.pack_propagate(0)
                
        fProjName.pack(fill=X, expand=True)
        fProjName.pack_propagate(0)
        
        fProjPath = Frame(self.fTop, relief=RIDGE, border=0, width=460, height=35)
        
        fProjPath_L = Frame(fProjPath, relief=RIDGE, border=0, width=100, height=35)
        fProjPath_R = Frame(fProjPath, relief=RIDGE, border=0, width=360, height=35)
        
        Path = Label(fProjPath_L, text='Path:', font=self.font_Text)
        Path.pack(side=RIGHT, anchor=N, padx=8)
        
        fProjPath_L.pack(fill=X, expand=True, side=LEFT)
        fProjPath_L.pack_propagate(0)     
        
        self.PathVal = Entry(fProjPath_R, textvariable=self.ProjectDir, background='white',
                             disabledbackground='white', disabledforeground='black')
        self.PathVal.pack(side=LEFT, anchor=N, padx=4, fill=X, expand=True)
        self.PathVal['font'] = self.font_Text
        self.PathVal.config(state='disabled')
        
        fProjPath_R.pack(fill=X, expand=True, side=LEFT)
        fProjPath_R.pack_propagate(0)
                
        fProjPath.pack(fill=X, expand=True)
        fProjPath.pack_propagate(0)
        
        fButtons = Frame(self.fTop, relief=RIDGE, border=0, width=460, height=40)
        
        #self.Btn_Browse = Button(fButtons, text='Browse', width=7, command=self.Btn_Browse_Clicked)
        #self.Btn_Browse.pack(side=LEFT, anchor=SW, pady=3)
        #self.Btn_Browse['font'] = self.font_Text
    
        Btn_Cancel = Button(fButtons, text='Cancel', width=7, command=self.Btn_Cancel_Clicked, font=self.font_Text)
        Btn_Cancel.pack(side=RIGHT, anchor=SE, pady=3)
        
        self.Btn_Load = Button(fButtons, text='Load', width=7, command=self.Btn_Load_Clicked, font=self.font_Text)
        self.Btn_Load.pack(side=RIGHT, anchor=SE, padx=4, pady=3)
        self.Btn_Load.config(state='disabled')
        
        fButtons.pack(fill=X, expand=True)
        fButtons.pack_propagate(0)
                
        #==================================================================================
        '''                  --- BOTTOM FRAME OF THE INTERFACE ---                      '''
        #==================================================================================
        
        fBottom = Frame(self.frame, relief=RIDGE, border=0, width=460, height=60)
        
        # Messages Box
        fMsg = Frame(fBottom, border=1, width=440, height=50, relief=SUNKEN)        

        self.TextMessage = Text(fMsg, width=440, height=5, background=self.Color_Gray, font=self.font_Text)
        self.TextMessage.pack(side=LEFT, anchor=S)
        self.TextMessage.config(state='disabled')

        fMsg.pack(fill=X, expand=True, side=LEFT, anchor=S, padx=10, pady= 5)
        fMsg.pack_propagate(0)

        fBottom.pack(fill=X, expand=True, side=BOTTOM, anchor=S)
        fBottom.pack_propagate(0)        
        

    ''' ==================================================================================
    FUNCTION LoadProjectIsRunning: Update or Create the Running File to BLOCK multiple GUI 
    ==================================================================================  '''       
    def LoadProjectIsRunning(self):
        
        #Create the .run file
        RunPath = os.path.join(self.RunDir_Path,'.lrun')
        RunFile = open(RunPath, 'w')
        RunFile.write(str(os.getpid()))
        RunFile.close()

        #Load up NRG built-in functions
        #cmd.do("run 
        
        
    ''' ==================================================================================
    FUNCTION Delete_RunFile: Delete the run file.
    ==================================================================================  '''    
    def Delete_RunFile(self):
        
        #Delete the .run file
        RunPath = os.path.join(self.RunDir_Path,'.lrun')
        if (os.path.isfile(RunPath)):
            try:
                os.remove(RunPath)
            except OSError:
                time.sleep(0.1)
                os.remove(RunPath)
                
                
    ''' ==================================================================================
    FUNCTION Fill_Table: Insert the possible project in the list boxes.
    ==================================================================================  '''                
    def Fill_Table(self):
        
        Keys = self.dictProjectList.keys()
        Keys.sort(reverse=True)

        for key in Keys:
            self.Table.Add( [ self.dictProjectList[key][0], self.dictProjectList[key][1], self.dictProjectList[key][2] ], 
                            [ None, None, None ] )
        
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):
        
        self.ProjectNameTrace = self.ProjectName.trace('w',self.ProjectName_Toggle)
        self.ProjectDirTrace = self.ProjectDir.trace('w',self.ProjectDir_Toggle)

    ''' ==================================================================================
    FUNCTION ProjectDir_Toggle: Callback function when the path is changed
    =================================================================================  '''    
    def ProjectDir_Toggle(self, *args):

        if self.ProjectDir.get() != '':
            # Loadable project
            self.Btn_Load.config(state='normal')
        else:
            self.Btn_Load.config(state='disabled')

    ''' ==================================================================================
    FUNCTION ProjectName_Toggle: Callback function when the user clicks on a project
    =================================================================================  '''    
    def ProjectName_Toggle(self, *args):

        self.ActualKey = ''
        for k in self.dictProjectList.keys():
            if self.dictProjectList[k][0] == self.ProjectName.get():
                self.ActualKey = k
                break

        if self.ActualKey != '':
            self.ActualDirPath = self.dictProjectList[self.ActualKey][3]
            self.ActualFileName = self.dictProjectList[self.ActualKey][0]
            self.ProjectDir.set(str(' ' + self.ActualDirPath))

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
            
            MainFilePath = os.path.join(self.HomeDir,'.proj')
            UserFilePath = os.path.join(Default,'.proj')
            Date = self.Get_Date()
            Seconds = time.time()

            LineToAdd = 'Name: Default Seconds: ' + str(Seconds) + ' LastUsed: ' + \
                        Date + ' Creation: ' + Date + ' Path: ' + Default + '\n'

            # Add the project path to the .proj file
            ProjFileMain = open(MainFilePath, 'w')
            ProjFileMain.write(LineToAdd)
            ProjFileMain.close()

            # Add the project path to the .proj file
            ProjFileUser = open(UserFilePath, 'w')
            ProjFileUser.write(LineToAdd)
            ProjFileUser.close()
            

    ''' ==================================================================================
    FUNCTION Update_ProjectFile: Add the created project to the Project File list.
    ==================================================================================  '''            
    def Update_ProjectFile(self):
       
        MainFilePath = os.path.join(self.HomeDir,'.proj')
        UserFilePath = os.path.join(self.ActualDirPath,'.proj')
        Date = self.Get_Date()
        Seconds = time.time()
        
        # Add the project path to the .proj file
        file = open(MainFilePath, 'r')
        ProjFile = file.readlines()
        file.close()
        
        # Update the file        
        LineToAdd = 'Name: ' + self.ActualFileName + ' Seconds: ' + str(Seconds) + ' LastUsed: ' + Date + \
                    ' Creation: ' + self.dictProjectList[self.ActualKey][2] + ' Path: ' + self.ActualDirPath + '\n'
        StartLine = 'Name: ' + self.ActualFileName
        
        WriteFile = open(MainFilePath, 'w')
        NotFound = True
        
        for line in ProjFile:
            if NotFound and line.startswith(StartLine):
                WriteFile.write(LineToAdd)
                NotFound = False
            else:
                WriteFile.write(line)
        
        WriteFile.close()  

        # Add the project path to the .proj file
        ProjFileUser = open(UserFilePath, 'w')
        ProjFileUser.write(LineToAdd)
        ProjFileUser.close()

        
    ''' ==================================================================================
    FUNCTION Btn_Cancel_Clicked: Cancel the project creation then quit the application.
    ==================================================================================  '''        
    def Btn_Cancel_Clicked(self):
        
        print('   CANCELLED the loading of a Project...')
        self.inThread.EnableMenu = False
        self.inThread.StopThread = True
        self.inThread.UserPath = ''
        
        self.Delete_RunFile()
        
        self.top.destroy()        
        
    
    ''' ==================================================================================
    FUNCTION Btn_Load_Clicked: Load a project then quit the application.
    ==================================================================================  '''        
    def Btn_Load_Clicked(self):
        
        print('   Loaded the Project ' + self.ActualFileName)

        self.Update_ProjectFile()
        
        self.inThread.EnableMenu = True
        self.inThread.StopThread = True
        self.inThread.UserPath = self.ActualDirPath
        self.inThread.ProjectName = self.ActualFileName
        
        self.Delete_RunFile()
        
        self.top.destroy()
        

    ''' ==================================================================================
    FUNCTION Btn_Browse_Clicked: Browse to specify the directory to install the project.
    ==================================================================================  '''        
    def Btn_Browse_Clicked(self):

        src = tkFileDialog.askdirectory(initialdir=self.HomeDir, title='Select a directory')

        if len(src) > 0:
            if os.path.isdir(src):
                self.ActualDirPath = src             
                self.ProjectDir.set(self.NameSpacer + src)
                
                # Get Project Name   
        

    ''' ==================================================================================
    FUNCTION LoadProjectList: Load the projects list from the .proj file.
    ==================================================================================  '''   
    def LoadProjectList(self):
            
        MainFilePath = os.path.join(self.HomeDir,'.proj')
        
        if os.path.isfile(MainFilePath):
        
            # Add the project path to the .proj file
            file = open(MainFilePath, 'r')
            ProjFile = file.readlines()
            file.close()
            
            for line in ProjFile:

                if line.startswith('Name:'):
                    
                    endSec = line.rfind('Seconds:') - 1
                    startSec = endSec + 10
                    
                    endUsed = line.rfind('LastUsed:') - 1
                    startUsed = endUsed + 11
                    
                    endCreate = line.rfind('Creation:') - 1
                    startCreate = endCreate + 11
                    
                    endPath = line.rfind('Path:') - 1
                    startPath = endPath + 7                    
                    
                    Name = line[6:endSec]
                    Seconds = line[startSec:endUsed]                  
                    LastUsed = line[startUsed:endCreate]
                    Creation = line[startCreate:endPath]
                    Path = line[startPath:].replace('\n', '')
                    
                    # Validate the Path existence
                    if os.path.isdir(Path):
                        self.dictProjectList[Seconds] = [ Name, LastUsed, Creation, Path ]
                               
               
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
            
#===================================================================================
# REMOVE THE COMMENTS OF THIS SECTION TO BE ABLE TO SEE THIS INTERFACE IN CONSOLE
#===================================================================================
# root = Tk()
# 
# root.path_FlexAID = '/home/eugene/FlexAID/'
# Button(root, text='Hello!').pack()
# root.update()
# 
# d = displayCfgFile(root)
# 
# root.wait_window(d.top)
#===================================================================================


