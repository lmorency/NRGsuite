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
import time

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
    FUNCTION Btn_SaveDefault_Clicked: Sets the default config
    ================================================================================== '''    
    def Btn_SaveDefault_Clicked(self):

        return

    ''' ==================================================================================
    FUNCTION Btn_Restore_Clicked: Sets back the original default config
    ================================================================================== '''    
    def Btn_Restore_Clicked(self):

        return

    ''' ==================================================================================
    FUNCTION SetActiveFrame: Switch up tabs in the uppper menu
    ================================================================================== '''    
    def SetActiveFrame(self, event, Frame):

        if not self.ActiveWizard is None:
            self.DisplayMessage("Cannot switch tab: A wizard is currently running...", 2)
            return

        if self.ProcessRunning:
            self.DisplayMessage("Cannot switch tab: A process is currently running...", 2)
            return

        if self.ActiveFrame != Frame:

            if not self.ActiveFrame is None:
                
                print "Validating all entries before switching frame"
                rv = self.ActiveFrame.Validate_Entries()
                print "rv=", rv
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


            print "Switching frame..."
            self.ActiveFrame = Frame
            
            print("Showing the frame", self.ActiveFrame.FrameName)
            self.ActiveFrame.Show()
            self.ActiveFrame.After_Show()
            self.ActiveFrame.Tab.config(bg=self.Color_Blue)

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
                self.ActiveWizard.btn_Done()
        
        #Delete the .run file
        RunPath = os.path.join(self.AlreadyRunning_Dir,self.RunFile)
        if os.path.isfile(RunPath):
            try:
                os.remove(RunPath)
            except OSError:
                print('   An error occured while clearing running file for ' + self.Name)                
                    
        self.Del_Trace()
        
        self.top.destroy()        

        print('   Closed ' + self.Name)

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes the trace of some variables
    ================================================================================== '''    
    def Del_Trace(self):
    
        for Tab in self.listTabs:
            Tab.Del_Trace()

