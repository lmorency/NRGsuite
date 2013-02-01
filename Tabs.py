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

class Tab:

    def __init__(self, top, PyMOL, FrameButton, FrameName, Vars):

        self.PyMOL = PyMOL

        self.top = top
        self.Tab = FrameButton
        self.FrameName = FrameName

        self.font_Text = self.top.font_Text
        self.font_Title = self.top.font_Title

        self.Color_Black = self.top.Color_Black
        self.Color_Blue = self.top.Color_Blue
        self.Color_White = self.top.Color_White

        self.DisplayMessage = self.top.DisplayMessage

        self.Vars = Vars
        self.Update_Vars()
        self.Def_Vars()
        self.Init_Vars()

        self.fFrame = self.Frame()
        self.Trace()
        
        self.StateList = []
        self.Validator = []
        
    def Def_Vars(self):
    
        return

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
    FUNCTION Update_Vars: Update session variables when a session is loaded
    =================================================================================  '''    
    def Update_Vars(self):
    
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
    FUNCTION Load_Message: Welcome message (frame-specific)
    =================================================================================  '''    
    def Load_Message(self):

        return
        
        ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    def Validator_Fail(self):

        return

    ''' ==================================================================================
    FUNCTION Validate_Field: Validates an Entry Field in the FlexAID interface
    ==================================================================================  '''    
    def Validate_Field(self, args):

        Entry = args[0]
        StringVar = args[1]
        Min = args[2]
        Max = args[3]
        nDec = args[4]
        Validator = args[5]
        Label = args[6]
        Type = args[7]

        rv = -1

        if Type == 'float':
            # If the min-max value depends on another Widget Entry value
            if type(Min) != float:
                try:
                    Min = float(Min.get())
                except:
                    Min = 0.0

            if type(Max) != float:
                try:
                    Max = float(Max.get())
                except:
                    Max = 1.0

            rv = General.validate_Float(StringVar.get(), Min, Max, nDec)

        elif Type == 'int':
            # If the min-max value depends on another Widget Entry value
            if type(Min) != int:
                try:
                    Min = int(Min.get())
                except: 
                    Min = 1

            if type(Max) != int:
                try:
                    Max = int(Max.get())
                except:
                    Max = 100
        
            rv = General.validate_Integer(StringVar.get(), Min, Max)
                
        elif Type == 'str':
            rv = General.validate_String(StringVar.get())

        # Return-value testing
        if rv == 0:
            Entry.config(bg=self.Color_White)
            Validator[0] = True
 
        elif rv == 1:
            self.DisplayMessage("Value has erroneus format for field " + Label, 1)
            Validator[0] = False

        elif rv == 2:
            self.DisplayMessage("The number of decimals cannot exceed (" + str(nDec) + ") for field " + Label, 1)
            Validator[0] = False

        elif rv == 3:
            self.DisplayMessage("Value must be within the range[" + str(Min) + "," + str(Max) + "] for field " + Label, 1)
            Validator[0] = False

        elif rv == -1:
            print("Unknown data format")

        return
        
    ''' ==================================================================================
    FUNCTION Lost_Focus: Draws a red square in box if field was not validated
    ==================================================================================  '''    
    def Lost_Focus(self, event, args):

        self.Validate_Field(args)

        Entry = args[0]
        Validator = args[5]
        
        if not Validator[0]:
            Entry.config(bg=self.Color_Red)

        return "break"

    ''' ==================================================================================
    FUNCTION Validate_Entries: Validate all entries of a frame before switching
    ==================================================================================  '''    
    def Validate_Entries(self, list):
        
        for valid in list:
            if not valid[2] and valid[0] == False:

                if valid[3] != None:
                    valid[3].config(bg=self.Color_Red)

                return valid[1]

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

