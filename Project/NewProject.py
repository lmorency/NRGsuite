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
@title: NewProject - Interface

@summary: This is the interface of the New Project application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec 8, 2010
'''

from Tkinter import *

import os
import time
import tkFileDialog

import Base
import NRGsuite
import General

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayNewProject(Base.Base):
    
    ''' ==================================================================================
    FUNCTION Def_Vars: Define extra variables for the class
    ==================================================================================  '''  
    def Def_Vars(self):
            
        self.ProjDirPath = StringVar()
        self.ProjName = StringVar()
        self.useDefault = IntVar()
    
    ''' ==================================================================================
    FUNCTION Init_Vars: Initialize the extra variables
    ==================================================================================  '''  
    def Init_Vars(self):
        
        self.NameSpacer = '  '
        self.ActualDirPath = self.NRGsuite_Dir

        self.ProjName.set('')
        self.ProjDirPath.set(os.path.join(self.NameSpacer,self.ActualDirPath))         
        
        self.useDefault.set(1)
    
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        self.fTop = Frame(self.fMain, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=260)        

        fTitle = Frame(self.fTop, relief=RIDGE, border=1, width=self.WINDOWWIDTH, height=40, bg=self.Color_White)
        
        Title = Label(fTitle, text='Create a new NRG suite Project.', bg=self.Color_White)
        Title.pack(padx=10, side=LEFT, anchor=W)
        Title['font'] = self.font_Title
        
        fTitle.pack(fill=X, expand=True)
        fTitle.pack_propagate(0)
        
        fProjName = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=50)
        
        Name = Label(fProjName, text='Project name:')
        Name.pack(padx=10, side=LEFT, anchor=W)
        Name['font'] = self.font_Text
        
        self.NameVal = Entry(fProjName, width=28, textvariable=self.ProjName, background='white')
        self.NameVal.pack(side=LEFT, anchor=W)
        self.NameVal['font'] = self.font_Text
        self.NameVal.focus_set()
        
        fProjName.pack(fill=X, expand=True, padx=10)
        fProjName.pack_propagate(0)
        
        fProjContent = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=30)
        
        Content = Label(fProjContent, text='Project contents:')
        Content.pack(padx=10, side=LEFT, anchor=S)
        Content['font'] = self.font_Text
        
        fProjContent.pack(fill=X, expand=True, padx=10)
        fProjContent.pack_propagate(0)
        
        fDefault = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=30)
        
        Default = Checkbutton(fDefault, text='  Use default', variable=self.useDefault, command=self.Default_Checked)
        Default.pack(padx=25, side=LEFT, anchor=W)
        Default['font'] = self.font_Text_U
        
        fDefault.pack(fill=X, expand=True, padx=10)
        fDefault.pack_propagate(0)
        
        fDirectory = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=30)
        
        Directory = Label(fDirectory, text='Create project directory under:')
        Directory.pack(padx=10, side=LEFT, anchor=W)
        Directory['font'] = self.font_Text
        
        fDirectory.pack(fill=X, expand=True, padx=10)
        fDirectory.pack_propagate(0)
        
        fDirPath = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=30)
        
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
        
        fButtons = Frame(self.fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=50)        
    
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
        
        fBottom = Frame(self.fMain, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=60)
        
        # Messages Box
        fMsg = Frame(fBottom, border=1, width=380, height=50, relief=SUNKEN)        

        self.TextMessage = Text(fMsg, width=380, height=5, background=self.Color_Grey)
        self.TextMessage.pack(side=LEFT, anchor=S)
        self.TextMessage['font'] = self.font_Text
        self.TextMessage.config(state='disabled')

        fMsg.pack(fill=X, expand=True, side=LEFT, anchor=S, padx=10, pady= 5)
        fMsg.pack_propagate(0)

        fBottom.pack(fill=X, expand=True, side=BOTTOM, anchor=S)
        fBottom.pack_propagate(0)
    
    ''' ==================================================================================
    FUNCTION Btn_Cancel_Clicked: Cancel the project creation then quit the application.
    ==================================================================================  '''        
    def Btn_Cancel_Clicked(self):
            
        self.top.Project_Dir = ''
        
        self.Quit()
        
        print('  The creation of a New Project was cancelled.')

    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        
        NRGsuite.EnableDisableMenu(self.top, [ 'disabled', 'disabled', 'normal', 'normal', 'normal', 'normal', 'normal' ] )

    ''' ==================================================================================
    FUNCTION Btn_Create_Clicked: Create the new project then quit the application.
    ==================================================================================  '''        
    def Btn_Create_Clicked(self):
        
        name = self.ProjName.get().strip()
        if name != '':
            
            if General.validate_String(name, '', False, True, False):
                print('  ERROR: The name \'' + name + '\' is invalid. Try again.')
                return
            
            Project_Dir = os.path.join(self.ActualDirPath,name)
            
            if not os.path.isdir(Project_Dir):
                try:
                    os.makedirs(Project_Dir)
                except OSError:
                    print('  ERROR: Could not create the project: \'' + name + '\' under that directory.')
                    return 
                    
                self.Update_ProjectFile(Project_Dir, name)
                
                self.top.Project_Dir = Project_Dir
                self.top.ProjectName = name       
                
                self.Quit()
                self.After_Quit()
                
                print('  Successfully created the project: \'' + name + '\'')
                print('  The project was loaded by default.')
    
            else:
                print('  The project: \'' + name + '\' already exists.')
    
    ''' ==================================================================================
    FUNCTION Btn_Browse_Clicked: Browse to specify the directory to install the project.
    ==================================================================================  '''        
    def Btn_Browse_Clicked(self):

        src = tkFileDialog.askdirectory(initialdir=self.NRGsuite_Dir, title='Select a directory')

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
            self.DirPath.config(state='normal', font=self.font_Text)
        else:
            self.Btn_Browse.config(state='disabled')
            self.DirPath.config(state='disabled', font=self.font_Text)

            self.ActualDirPath = self.self.NRGsuite_Dir
            self.ProjDirPath.set(os.path.join(self.NameSpacer,self.ActualDirPath))
            
    ''' ==================================================================================
    FUNCTION Update_ProjectFile: Add the created project to the Project File list.
    ==================================================================================  '''            
    def Update_ProjectFile(self, Project_Dir, ProjName):
        
        MainFilePath = os.path.join(self.NRGsuite_Dir,'.proj')
        UserFilePath = os.path.join(Project_Dir,'.proj')
        Date = General.Get_Date()
        Seconds = time.time()
        
        LineToAdd = 'Name: ' + ProjName + ' Seconds: ' + str(Seconds) + ' LastUsed: ' + Date + ' Creation: ' + Date + ' Path: ' + Project_Dir + '\n'
        
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
        
