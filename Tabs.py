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

import General
import Queue

class Tab:

    # 100 ms
    TKINTER_UPDATE_INTERVAL = 100
    
    def __init__(self, top, PyMOL, FrameButton, FrameName, Vars):

        self.PyMOL = PyMOL

        self.top = top
        self.Tab = FrameButton
        self.FrameName = FrameName

        self.font_Text = self.top.font_Text
        self.font_Title = self.top.font_Title
        self.font_Title_H = self.top.font_Title_H

        self.Color_Black = self.top.Color_Black
        self.Color_Blue = self.top.Color_Blue
        self.Color_White = self.top.Color_White
        self.Color_Red = self.top.Color_Red
        
        self.DisplayMessage = self.top.DisplayMessage
        
        self.StateList = []
        self.Validator = []
        
        self.Vars = Vars
        self.Def_Vars()
        self.Init_Vars()
        
        self.fFrame = self.Frame()
        self.Trace()
                
    ''' ==================================================================================
    FUNCTION Def_Vars: Define extra variables for the tab
    ==================================================================================  '''  
    def Def_Vars(self):
    
        return

    ''' ==================================================================================
    FUNCTION Init_Vars: Sets the default values of the variables
    ==================================================================================  '''  
    def Init_Vars(self):
    
        return

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
    
        self.Load_Message()
        self.fFrame.pack(fill=BOTH, expand=True)
        
    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):
    
        return
        
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''    
    def Del_Trace(self):
    
        return

    ''' ==================================================================================
    FUNCTION Before_Kill_Frame: Actions related before killing a frame
    =================================================================================  '''    
    def Before_Kill_Frame(self):
        
        return True

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''
    def Kill_Frame(self):
        
        self.fFrame.pack_forget()

        return True
        
    ''' ==================================================================================
    FUNCTION Frame: Generate the frame in the the middle frame section of FlexAID root
    ==================================================================================  '''  
    def Frame(self):
    
        return
    
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):
    
        return

    ''' ==================================================================================
    FUNCTION Load_Session: Actions related to when a new session is loaded
    =================================================================================  '''    
    def Load_Session(self):

        return
    
    ''' ==================================================================================
    FUNCTION Check_Integrity: Compares the checksum of the physical files in the session
    =================================================================================  '''    
    def Check_Integrity(self):

        return

    ''' ==================================================================================
    FUNCTION Load_Message: Welcome message (frame-specific)
    =================================================================================  '''    
    def Load_Message(self):

        return
        
        ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    
    ''' ==================================================================================
    FUNCTION Validator_Fail: Runs actions when tab fails to validate
    =================================================================================  '''    
    def Validator_Fail(self):

        return

    ''' ==================================================================================
    FUNCTION Validate_Field: Validates an Entry Field in the FlexAID interface
    ==================================================================================  '''    
    def Validate_Field(self, *args, **kwargs):
        
        input = kwargs.pop('input')
        var = kwargs.pop('var')
        min = kwargs.pop('min')
        max = kwargs.pop('max')
        ndec = kwargs.pop('ndec')
        tag = kwargs.pop('tag')
        _type = kwargs.pop('_type')
        
        rv = -1
        
        if _type == float:

            if type(min) != float:
                try:
                    min = float(min.get())
                except:
                    min = 0

            if type(max) != float:
                try:
                    max = float(max.get())
                except:
                    max = 0
        
            rv = General.validate_Float(var.get(), min, max, ndec)
        
        elif _type == int:
        
            if type(min) != int:
                try:
                    min = int(min.get())
                except ValueError:
                    min = 0

            if type(max) != int:
                try:
                    max = int(max.get())
                except ValueError:
                    max = 100

            rv = General.validate_Integer(var.get(), min, max)
        
        # Return-value testing
        if rv == 0:
            input.config(bg=self.Color_White)
            return True
            
        elif rv == -1:
            return True
            
        else:
            if rv == 1:
                self.DisplayMessage("Value has erroneus format for field " + tag, 1)
            elif rv == 2:
                self.DisplayMessage("The number of decimals cannot exceed (" + str(ndec) + ") for field " + tag, 1)
            elif rv == 3:
                self.DisplayMessage("Value must be within the range[" + str(min) + "," + str(max) + "] for field " + tag, 1)

            input.config(bg=self.Color_Red)
            
            return False
    
        return True
        
    ''' ==================================================================================
    FUNCTION Validate_Fields: Validate all fields of a frame before switching
    ==================================================================================  '''    
    def Validate_Fields(self):
        
        for Validator in self.Validator:
            
            if not Validator[1] and \
                   Validator[2] != None and \
                   Validator[2]['state'] == 'normal' and \
                   Validator[2]['bg'] == self.Color_Red:
                   
                return Validator[0]
    
        return 0

    ''' ==================================================================================
    FUNCTION Disable_Frame: Disables all controls on main frame
    =================================================================================  '''    
    def Disable_Frame(self, *args):

        del self.StateList[:]
        General.saveState(self.fFrame, self.StateList)
        General.setState(self.fFrame)

        for arg in args:
            General.setState(arg,'normal')

    ''' ==================================================================================
    FUNCTION Enable_Frame: Enables all controls on main frame
    =================================================================================  '''    
    def Enable_Frame(self):

        General.backState(self.fFrame, self.StateList)

    ''' ==================================================================================
    FUNCTION Start_Update: Starts updating Tkinter with tasks in the queue
    ==================================================================================  '''               
    def Start_Update(self):
        
        self.queue = Queue.Queue()
        
        self.QueueParse = True
        self.Update_Tkinter()

    ''' ==================================================================================
    FUNCTION Condition_Update: Tests tab specific conditions to trigger stopping
    ==================================================================================  '''               
    def Condition_Update(self):

        if self.QueueParse:
            return True
        else:
            return False
        
    ''' ==================================================================================
    FUNCTION End_Update: Stops updating Tkinter interface using the queue
    ==================================================================================  '''               
    def End_Update(self):
        
        self.QueueParse = False
        
    ''' ==================================================================================
    FUNCTION After_Update: Executes tasks when done updating Tkinter
    ==================================================================================  '''               
    def After_Update(self):

        # override this method in base class
        return
        
    ''' ==================================================================================
    FUNCTION Update_Tkinter: update the tkinter interface (tasks queued from the worker)
    ==================================================================================  '''               
    def Update_Tkinter(self):
        
        # Check every 100 ms if there is something new in the queue.
        while self.queue.qsize():
            try:
                func = self.queue.get()
                func()
            except Queue.Empty:
                pass
        
        if self.Condition_Update():
            self.top.root.after(self.TKINTER_UPDATE_INTERVAL, self.Update_Tkinter)
        else:
            self.After_Update()
