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

'''
@title: LoadProject - Interface

@summary: This is the interface of the Load Project application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec 10, 2010
'''

from Tkinter import *

import os
import time
import tkFileDialog

import Base
import NRGsuite
import General
import MultiList

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayLoadProject(Base.Base):
    
    ''' ==================================================================================
    FUNCTION After_Init: Initialization of extra variables for the class
    ================================================================================== '''
    def After_Init(self):
        
        self.LoadProjectList()
        self.Fill_Table()    
        self.Trace()

    ''' ==================================================================================
    FUNCTION Def_Vars: Define extra variables for the class
    ==================================================================================  '''  
    def Def_Vars(self):
            
        self.ProjectPath = StringVar()        
        self.dictProjectList = {}
        
    ''' ==================================================================================
    FUNCTION Init_Vars: Initialize the extra variables
    ==================================================================================  '''  
    def Init_Vars(self):

        self.NameSpacer = '  '
        self.ActualDirPath = ''
        self.ActualFileName = ''
        self.ActualKey = ''

        self.ProjectPath.set('')
        self.dictProjectList.clear()

    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        fTitle = Frame(self.fMain, relief=RIDGE, border=1, width=self.WINDOWWIDTH, height=40, bg=self.Color_White)
        fTitle.pack(fill=X, expand=True)
        fTitle.pack_propagate(0)
        
        Title = Label(fTitle, text='Load an existing NRG suite project.', bg=self.Color_White, font=self.font_Title)
        Title.pack(padx=10, side=LEFT, anchor=W)        
        
        self.fTop = Frame(self.fMain, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=240)
        self.fTop.pack(fill=X, expand=True, padx=10)
        self.fTop.pack_propagate(0)
        
        fLblList = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=30)
        
        LblList = Label(fLblList, text='Project(s) list:', font=self.font_Text)
        LblList.pack(side=LEFT, anchor=S, pady=2)
        
        fLblList.pack(fill=X, expand=True)
        fLblList.pack_propagate(0)
        
        #==================================================================================
        '''                       --- LISTBOX FRAME SECTION ---                         '''
        #==================================================================================

        fProjList = Frame(self.fTop, relief=RIDGE, border=1, width=self.WINDOWWIDTH, height=110)
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
        fProjName = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=35)
        
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
        
        fProjPath = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=35)
        
        fProjPath_L = Frame(fProjPath, relief=RIDGE, border=0, width=100, height=35)
        fProjPath_R = Frame(fProjPath, relief=RIDGE, border=0, width=360, height=35)
        
        Path = Label(fProjPath_L, text='Path:', font=self.font_Text)
        Path.pack(side=RIGHT, anchor=N, padx=8)
        
        fProjPath_L.pack(fill=X, expand=True, side=LEFT)
        fProjPath_L.pack_propagate(0)     
        
        self.PathVal = Entry(fProjPath_R, textvariable=self.ProjectPath, background='white',
                             disabledbackground='white', disabledforeground='black')
        self.PathVal.pack(side=LEFT, anchor=N, padx=4, fill=X, expand=True)
        self.PathVal['font'] = self.font_Text
        self.PathVal.config(state='disabled')
        
        fProjPath_R.pack(fill=X, expand=True, side=LEFT)
        fProjPath_R.pack_propagate(0)
                
        fProjPath.pack(fill=X, expand=True)
        fProjPath.pack_propagate(0)
        
        fButtons = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=40)
        
        #self.Btn_Browse = Button(fButtons, text='Browse', width=7, command=self.Btn_Browse_Clicked)
        #self.Btn_Browse.pack(side=LEFT, anchor=SW, pady=3)
        #self.Btn_Browse['font'] = self.font_Text
    
        Btn_Cancel = Button(fButtons, text='Cancel', width=7, command=self.Btn_Cancel_Clicked, font=self.font_Text)
        Btn_Cancel.pack(side=RIGHT, anchor=SE, pady=3)
        
        self.Btn_Load = Button(fButtons, text='Load', width=7, command=self.Btn_Load_Clicked, font=self.font_Text)
        self.Btn_Load.pack(side=RIGHT, anchor=SE, padx=4, pady=3)
        self.Btn_Load.config(state='disabled')
        

        for col in self.Table.Columns.keys():
            self.Table.Columns[col]['List'].bind('<Double-Button-1>',self.Btn_Load_Clicked)
            self.Table.Columns[col]['List'].bind('<Return>',self.Btn_Load_Clicked)

        fButtons.pack(fill=X, expand=True)
        fButtons.pack_propagate(0)
                
        #==================================================================================
        '''                  --- BOTTOM FRAME OF THE INTERFACE ---                      '''
        #==================================================================================
        
        fBottom = Frame(self.fMain, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=60)
        
        # Messages Box
        fMsg = Frame(fBottom, border=1, width=440, height=50, relief=SUNKEN)        

        self.TextMessage = Text(fMsg, width=440, height=5, background=self.Color_Grey, font=self.font_Text)
        self.TextMessage.pack(side=LEFT, anchor=S)
        self.TextMessage.config(state='disabled')

        fMsg.pack(fill=X, expand=True, side=LEFT, anchor=S, padx=10, pady= 5)
        fMsg.pack_propagate(0)

        fBottom.pack(fill=X, expand=True, side=BOTTOM, anchor=S)
        fBottom.pack_propagate(0)        
    
    ''' ==================================================================================
    FUNCTION Fill_Table: Insert the possible project in the list boxes.
    ==================================================================================  '''                
    def Fill_Table(self):

        for key in reversed(sorted(self.dictProjectList.keys())):
        
            self.Table.Add( [ self.dictProjectList[key][0], self.dictProjectList[key][1], self.dictProjectList[key][2] ], 
                            [ None, None, None ] )
        
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):
        
        try:
            self.ProjectNameTrace = self.ProjectName.trace('w',self.ProjectName_Toggle)
            self.ProjectPathTrace = self.ProjectPath.trace('w',self.ProjectPath_Toggle)
        except:
            pass

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes the trace of some variables
    ================================================================================== '''    
    def Del_Trace(self):
    
        try:
            self.ProjectName.trace_vdelete('w',self.ProjectNameTrace)
            self.ProjectPath.trace_vdelete('w',self.ProjectPathTrace)
        except:
            pass
        
    ''' ==================================================================================
    FUNCTION ProjectPath_Toggle: Callback function when the path is changed
    =================================================================================  '''    
    def ProjectPath_Toggle(self, *args):

        if self.ProjectPath.get() != '':
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
            self.ProjectPath.set(str(' ' + self.ActualDirPath))
                
    ''' ==================================================================================
    FUNCTION Update_ProjectFile: Add the created project to the Project File list.
    ==================================================================================  '''            
    def Update_ProjectFile(self):
       
        MainFilePath = os.path.join(self.NRGsuite_Dir,'.proj')
        UserFilePath = os.path.join(self.ActualDirPath,'.proj')
        Date = General.Get_Date()
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
                
        self.Quit()
        
    ''' ==================================================================================
    FUNCTION Btn_Load_Clicked: Load a project then quit the application.
    ==================================================================================  '''        
    def Btn_Load_Clicked(self, *args):
        
        self.Update_ProjectFile()
        
        self.top.Project_Dir = self.ActualDirPath
        self.top.ProjectName = self.ActualFileName
        
        self.Quit()

        self.After_Quit()

        print('  Successfully loaded the project \'' + self.ActualFileName + '\'')
            
    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        
        NRGsuite.EnableDisableMenu(self.top, [ 'disabled', 'disabled', 'normal', 'normal', 'normal', 'normal', 'normal' ] )
    
    ''' ==================================================================================
    FUNCTION Btn_Browse_Clicked: Browse to specify the directory to install the project.
    ==================================================================================  '''        
    def Btn_Browse_Clicked(self):

        src = tkFileDialog.askdirectory(initialdir=self.NRGsuite_Dir, title='Select a directory')

        if len(src) > 0:
            if os.path.isdir(src):
                self.ActualDirPath = src             
                self.ProjectPath.set(self.NameSpacer + src)

    ''' ==================================================================================
    FUNCTION LoadProjectList: Load the projects list from the .proj file.
    ==================================================================================  '''   
    def LoadProjectList(self):
            
        MainFilePath = os.path.join(self.NRGsuite_Dir,'.proj')
        
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
                               
