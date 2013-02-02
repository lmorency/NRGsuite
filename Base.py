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

import os, sys
import pickle
import time
import tkFont
import tkFileDialog

import Prefs
import Color
import General

if __debug__:
	from pymol import cmd
	from pymol.cgo import *
	from pymol.vfont import plain

	import General_cmd


class Base:

    def __init__(self):
    
        return
        
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''  
    def Frame_Main(self):

            return
    
    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Default_Clicked(self):
        
        if self.ActiveWizard != None:
            self.DisplayMessage("Cannot reset values while a Wizard is active", 2)
            return

        if self.ProcessRunning is True:
            self.DisplayMessage("Cannot reset values while a Process is running", 2)
            return
            
        self.ActiveFrame.Init_Vars()

    ''' ==================================================================================
    FUNCTION SetActiveFrame: Switch up tabs in the uppper menu
    ================================================================================== '''    
    def SetActiveFrame(self, Frame):

        if not self.ActiveWizard is None:
            self.DisplayMessage("Cannot switch tab: A wizard is currently running...", 2)
            return

        if self.ProcessRunning:
            self.DisplayMessage("Cannot switch tab: A process is currently running...", 2)
            return

        if self.ActiveFrame != Frame:

            if not self.ActiveFrame is None:

                # Trigger lost_focus event for validation
                self.fMiddle.focus_set()
                self.fMiddle.update_idletasks()

                rv = self.ActiveFrame.Validate_Entries(self.ActiveFrame.Validator)
                if rv > 0:
                    if rv == 1:
                        self.DisplayMessage("Cannot switch tab: Not all fields are validated", 2)
                    elif rv == 2:
                        self.ActiveFrame.Validator_Fail()
                    return

                if not self.ActiveFrame.Before_Kill_Frame() or not self.ActiveFrame.Kill_Frame():
                    self.DisplayMessage("Cannot switch tab: Not all fields are validated", 2)
                    return

                self.fMiddle.update_idletasks()

                #print "Killed Frame " + self.ActiveFrame.FrameName
                self.ActiveFrame.Tab.config(bg=self.Color_White)

            self.ActiveFrame = Frame
            
            self.ActiveFrame.Show()
            self.ActiveFrame.After_Show()
            self.ActiveFrame.Tab.config(bg=self.Color_Blue)

            self.fMiddle.update_idletasks()

        return

    ''' ==================================================================================
    FUNCTION AppIsRunning: Update or Create the Running File to BLOCK multiple GUI 
    ==================================================================================  '''       
    def AppIsRunning(self):
        
        #Create the .run file
        RunPath = os.path.join(self.AlreadyRunning_Dir,self.RunFile)
        RunFile = open(RunPath, 'w')
        RunFile.write(str(os.getpid()))
        RunFile.close()

    ''' ==================================================================================
    FUNCTION Btn_Restore_Clicked: Restore the original default configuration
    ================================================================================== '''    
    def Btn_Restore_Clicked(self):
        
        return

    ''' ==================================================================================
    FUNCTION Btn_SaveDefault_Clicked: Saves the current configuration as default
    ================================================================================== '''    
    def Btn_SaveDefault_Clicked(self):
        
        return

    ''' ==================================================================================
    FUNCTION DisplayMessage: Display the message  
    ==================================================================================  '''    
    def DisplayMessage(self, msg, priority):
        
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
    FUNCTION Quit: Exit the application
    ==================================================================================  '''
    def Quit(self):
        
        self.Before_Quit()
        
        # Cannot quit while process is running
        if self.ProcessRunning is True and self.Run is not None:
            print("   The following process will be killed", self.Run)
            self.Run.kill()
            
        # Close any Wizard interface in Pymol if started
        if not self.ActiveWizard is None:
            if self.PyMOL:
                cmd.get_wizard().btn_Done()

        #Delete the .run file
        RunPath = os.path.join(self.AlreadyRunning_Dir,self.RunFile)
        if os.path.isfile(RunPath):
            try:
                os.remove(RunPath)
            except OSError:
                time.sleep(0.1)
                os.remove(RunPath)
                    
        self.Del_Trace()
        
        self.top.destroy()        

        print('   Closed ' + self.Name)

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes the trace of some variables
    ================================================================================== '''    
    def Del_Trace(self):
    
        for Tab in self.listTabs:
            Tab.Del_Trace()

