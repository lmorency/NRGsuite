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

from __future__ import print_function

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
    import tkFont
else:
    from tkinter import *
    import tkinter.font as tkFont

from subprocess import Popen, PIPE
from pymol import cmd

import ctypes
import os

import Prefs
import General

class Base(object):

    ActiveWizard = None
    
    WizardError = False
    WizardResult = 0

    # in milliseconds
    TKINTER_UPDATE_INTERVAL = 100
    
    WEBSITE = 'http://bcb.med.usherbrooke.ca/'
    VERSION = '2.48l'
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, root, top, menuindex, Project_Dir, Install_Dir, NRGsuite_Dir, OSid, PyMOL, Name, WINDOWWIDTH, WINDOWHEIGHT, RootPrefs):
        
        self.Name = Name
        self.OSid = OSid
        self.PyMOL = PyMOL
        self.pymol_major_version = int(cmd.get_version()[0][0])
        
        self.NRGsuite_Dir = NRGsuite_Dir
        self.Install_Dir = Install_Dir
        self.Project_Dir = Project_Dir
        
        self.Folders = [ self.NRGsuite_Dir ]
        if self.Project_Dir:
            self.Folders.append(self.Project_Dir)
        
        self.Set_Folders()
        self.Validate_Folders()
        
        self.Prefs = RootPrefs

        self.Color_Green = '#CCFFCC'
        self.Color_Grey = '#EDEDED'
        self.Color_Gray = '#EDEDED'
        self.Color_Blue = '#6699FF'
        self.Color_Red = '#FF9999'
        self.Color_White = '#FFFFFF'
        self.Color_Black = '#000000'
        
        self.root = root
        self.root.title(self.Name)

        self.top = top

        self.menuindex = menuindex
        if self.menuindex != -1:
            self.Disable_MenuItem()
            
        self.WINDOWWIDTH = WINDOWWIDTH
        self.WINDOWHEIGHT = WINDOWHEIGHT
        
        #self.root.geometry()   # Interface DIMENSIONS
        #self.root.maxsize(WINDOWWIDTH,WINDOWHEIGHT)
        self.root.minsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        self.root.protocol('WM_DELETE_WINDOW', self.Quit)
        
        General.CenterWindow(self.root,self.WINDOWWIDTH,self.WINDOWHEIGHT)
        
        #================================================================================== 
        #                 SET the default fonts of the interface
        #==================================================================================
        
        FontType = self.Prefs.GetFontType()
        FontSize = self.Prefs.GetFontSize()
        
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)
        self.font_Title_H = tkFont.Font(family=FontType,size=FontSize + 1, weight=tkFont.BOLD)
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_H = tkFont.Font(family=FontType,size=FontSize + 1)
        self.font_Text_I = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        self.font_Text_U = tkFont.Font(family=FontType,size=FontSize, underline=True)
        
        self.ChildWindow = None
        self.ActiveFrame = None
        self.Run = None
        
        self.ProcessError = False
        self.ProcessRunning = False
        self.ProcessParsing = False
        self.ProcessDone = False
        self.ParseState = 10
        self.SimulateState = -1

        self.fMain = Frame(self.root)
        self.fMain.pack(expand=True)
        
        self.Def_Vars()
        self.Init_Vars()
        
        self.Frame_Main()
                
        self.Build_Tabs()
        self.MakeMenuBar()
        
        self.After_Init()
        
        self.Clean()
                
    ''' ==================================================================================
    FUNCTION Toggle_MenuItem: Toggle the state of the menu item associated to the GUI
    ==================================================================================  '''  
    def Toggle_MenuItem(self, state):

        if self.pymol_major_version == 2:
            k = 0
            for mi in self.top.menuBar._menudict['NRGsuite'].actions():
                print(mi,k,self.menuindex)
                if k == self.menuindex:
                    mi.setEnabled(state)
                    break
                k += 1
        else:
            str_state = 'disabled' if state == False else 'normal'
            self.top.menuBar.component('NRGsuite-menu').entryconfig(self.menuindex, state=str_state)
    
    ''' ==================================================================================
    FUNCTION Disable_MenuItem: Disable the menu item associated to the GUI
    ==================================================================================  '''  
    def Disable_MenuItem(self):

        self.Toggle_MenuItem(False)
        
    ''' ==================================================================================
    FUNCTION Enable_MenuItem: Enable the menu item associated to the GUI
    ==================================================================================  '''  
    def Enable_MenuItem(self):

        self.Toggle_MenuItem(True)
        
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''  
    def Frame_Main(self):

        return
    
    ''' ==================================================================================
    FUNCTION After_Init: Generates actions after initializing variables
    ==================================================================================  '''  
    def After_Init(self):

        return    
    
    ''' ==================================================================================
    FUNCTION Clean: Cleans the temporary folders of the application
    ==================================================================================  '''  
    def Clean(self):

        return    

    ''' ==================================================================================
    FUNCTION Def_Vars: Define extra variables for the class
    ==================================================================================  '''  
    def Def_Vars(self):
    
        return
        
    ''' ==================================================================================
    FUNCTION Init_Vars: Initialize the extra variables
    ==================================================================================  '''  
    def Init_Vars(self):
    
        return

    ''' ==================================================================================
    FUNCTION Build_Tabs: Builds the tab' classes of the main frame
    ==================================================================================  '''  
    def Build_Tabs(self):

        return

    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Default_Clicked(self):
        
        if self.ValidateProcessRunning() or self.ValidateWizardRunning() or \
            self.ValidateWindowRunning():
            return
            
        self.ActiveFrame.Init_Vars()

    ''' ==================================================================================
    FUNCTION SetActiveFrame: Switch up tabs in the uppper menu
    ================================================================================== '''    
    def SetActiveFrame(self, Frame):

        if self.ValidateProcessRunning() or self.ValidateWizardRunning() or \
            self.ValidateWindowRunning():
            return
        
        if self.ActiveFrame != Frame:

            if not self.ActiveFrame is None:
                
                rv = self.ActiveFrame.Validate_Fields()
                if rv > 0:
                    if rv == 1:
                        self.DisplayMessage("Cannot switch tab: Not all fields are validated", 2)
                    elif rv == 2:
                        self.ActiveFrame.Validator_Fail()
                    return

                if not self.ActiveFrame.Before_Kill_Frame() or not self.ActiveFrame.Kill_Frame():
                    self.DisplayMessage("Cannot switch tab: Could not kill the frame", 2)
                    return

                #print "Killed Frame " + self.ActiveFrame.FrameName
                self.ActiveFrame.Tab.config(bg=self.Color_White)

            self.ActiveFrame = Frame
            
            self.ActiveFrame.Show()
            self.ActiveFrame.After_Show()
            self.ActiveFrame.Tab.config(bg=self.Color_Blue)

    ''' ==================================================================================
    FUNCTION Restore: Restore the original default configuration
    ================================================================================== '''    
    def Restore(self):
        
        return
    
    ''' ==================================================================================
    FUNCTION SaveDefault: Saves the current configuration as default
    ================================================================================== '''    
    def SaveDefault(self):
        
        return

    ''' ==================================================================================
    FUNCTION DisplayMessage: Display the message  
    ==================================================================================  '''    
    def DisplayMessage(self, msg, priority):
        
        # Prepend text at end of control text
        self.TextMessage.mark_set(INSERT, END)

        self.TextMessage.config(state='normal')
         
        #self.TextMessage.config(font='red')
        self.TextMessage.insert(INSERT, '\n' + msg)

        if priority == 1:
            #self.TextMessage.tag_add('warn', lineNo + '.0', lineNo + '.' + str(NbChar))
            self.TextMessage.tag_config('warn', foreground='red')
        elif priority == 2:
            #self.TextMessage.tag_add('notice', lineNo + '.0', lineNo + '.' + str(NbChar))
            self.TextMessage.tag_config('notice', foreground='blue')

        self.TextMessage.yview(INSERT)        
        self.TextMessage.config(state='disabled')
                
    ''' ==================================================================================
    FUNCTION Before_Quit: Execute tasks before exitting the application
    ==================================================================================  '''
    def Before_Quit(self):
    
        return

    ''' ==================================================================================
    FUNCTION Kill: Kills a process by PID
    ==================================================================================  '''
    def Kill(self, pid):
    
         print("   The following process will be killed", self.Run)

         if self.OSid == 'WIN':
            kernel32 = ctypes.windll.kernel32
            handle = kernel32.OpenProcess(1, 0, self.Run.pid)
            return (0 != kernel32.TerminateProcess(handle, 0))

         else:
            killcmd = 'kill ' + str(self.Run.pid)
            pkill = Popen(killcmd, shell=True, stdout=PIPE, stderr=PIPE)
            (out, err) = pkill.communicate()
            if err:
                print("   An error occured while killing the following process", self.Run)

    ''' ==================================================================================
    FUNCTION Quit: Exit the application
    ==================================================================================  '''
    def Quit(self):
        
        if self.Before_Quit():
            return
        
        # Cannot quit while process is running
        if self.ProcessRunning and self.Run is not None:
            #self.Kill(self.Run.pid)
            self.Kill(self.Run.pid)

        # Close any Wizard interface in Pymol if started
        if self.ActiveWizard is not None:
            self.ActiveWizard.btn_Done()
        
        # Close child windows if any
        if self.ChildWindow is not None:
            self.ChildWindow.Quit()

        self.Del_Trace()
        
        if self.root is not None:
            self.root.destroy()
            self.root = None

        if self.menuindex != -1:
            print("Enable !")
            self.Enable_MenuItem()
                
        print('  Closed ' + self.Name)

    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        
        return
        
    ''' ==================================================================================
    FUNCTION Trace: Adds tracer to Tk variables
    ================================================================================== '''    
    def Trace(self):
    
        return

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes the trace of some variables
    ================================================================================== '''    
    def Del_Trace(self):
    
        return

    ''' ==================================================================================
    FUNCTION MakeMenuBar: Builds the menu on the upper left corner    
    ==================================================================================  '''        
    def MakeMenuBar(self):
    
        return

    ''' ==================================================================================
    FUNCTION ValidateWizardRunning: Validates if a wizard is active
    ================================================================================== '''    
    def ValidateWizardRunning(self):

        if self.ActiveWizard is not None:
            self.DisplayMessage("  Cannot execute task because a wizard is currently running.", 2)
            return 1

        return 0

    ''' ==================================================================================
    FUNCTION WizardRunning: Checks if a wizard is running
    ================================================================================== '''    
    def WizardRunning(self):

        if self.ActiveWizard is not None:
            return 1

        return 0
        
    ''' ==================================================================================
    FUNCTION ValidateProcessRunning: Validates if a process is currently running
    ================================================================================== '''    
    def ValidateProcessRunning(self):

        if self.ProcessRunning or self.Run is not None:
            self.DisplayMessage("  Cannot execute task because a process is currently running.", 2)
            return 1

        return 0
        
    ''' ==================================================================================
    FUNCTION ValidateWindowRunning: Validates if a child window is opened
    ================================================================================== '''    
    def ValidateWindowRunning(self):

        if self.ChildWindow is not None:
            self.DisplayMessage("  Cannot execute task because a child window is currently opened.", 2)
            return 1

        return 0
    
    ''' ==================================================================================
    FUNCTION Set_Folders: Builds the list of folders that will be built
    ==================================================================================  '''
    def Set_Folders(self):
        
        return

    ''' ==================================================================================
    FUNCTION Validate_Folders: Make all the folders necessary for running the project
    ==================================================================================  '''
    def Validate_Folders(self):

        for Folder in self.Folders:
            if not os.path.isdir(Folder):
                os.makedirs(Folder)
