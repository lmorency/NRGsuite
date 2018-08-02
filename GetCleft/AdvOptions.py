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

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
else:
    from tkinter import *

import functools
import General

'''
@title: FlexAID - IOFile tab - Interface

@summary: This is the IOFile tab interface of FlexAID application

@organization: Najmanovich Research Group
@creation date:  Nov. 24 2011
'''

class AdvOptions(object):

    def __init__(self,top):

        self.top = top
        self.FrameName = 'AdvOptions'
        self.Tab = self.top.Btn_Config

        self.DisplayMessage = self.top.DisplayMessage

        self.Def_Vars()
        self.Init_Vars()

        self.Frame()
        self.Trace()

    def Def_Vars(self):

        self.Entry_L = StringVar()
        self.Entry_U = StringVar()
        self.Entry_K = StringVar()
        self.Entry_R = StringVar()
        self.Entry_R_Value = IntVar()
        
        self.Entry_h = StringVar()
        self.Entry_h_Value = IntVar()

        self.Entry_H = StringVar()
        self.Entry_H_Value = IntVar()

        self.Entry_C = StringVar()

        self.ValidEntry_L = list()
        self.ValidEntry_U = list()
        self.ValidEntry_K = list()
        #self.ValidEntry_R = list()
        #self.ValidEntry_h = list()
        #self.ValidEntry_H = list()
        #self.ValidEntry_C = list()

        self.Validator = list()

    def Init_Vars(self):

        self.Entry_L.set('1.50')
        self.Entry_U.set('4.00')
        self.Entry_K.set('5.00')

        self.Entry_R.set('False')
        self.Entry_R_Value.set(0)

        self.Entry_h.set('False')
        self.Entry_h_Value.set(0)

        self.Entry_H.set('False')
        self.Entry_H_Value.set(0)

        self.Entry_C.set('ALL') 

        self.ValidEntry_L = [True, 1, 0, None]
        self.ValidEntry_U = [True, 1, 0, None]
        self.ValidEntry_K = [True, 1, 0, None]
        #self.ValidEntry_R = [True, 1, 0, None]
        #self.ValidEntry_h = [True, 1, 0, None]
        #self.ValidEntry_H = [True, 1, 0, None]
        #self.ValidEntry_C = [True, 1, 0, None]

        self.Validator = [ self.ValidEntry_L, self.ValidEntry_U, self.ValidEntry_K ]#, self.ValidEntry_C]

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
        
        self.fAdvOptions.pack_forget()
        #self.fAdvOptions.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    def Validator_Fail(self):

        return

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
        
        self.fAdvOptions.pack(fill=BOTH, expand=True)

        self.LoadMessage()

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        return

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        return

    ''' ==================================================================================
    FUNCTION Frame: Generate the Input / Output Files frame in the the middle 
                    frame section    
    ==================================================================================  '''  
    def Frame(self):

        self.fAdvOptions = Frame(self.top.fMiddle, relief=RIDGE)

        #==================================================================================
        #                           SPHERES RADII
        #==================================================================================                
        fSpheres = Frame(self.fAdvOptions)
        fSpheres.pack(side=TOP, fill=X, padx=5, pady=10)
        fSpheresLine1 = Frame(fSpheres)
        fSpheresLine1.pack(side=TOP, fill=X)
        fSpheresLine2 = Frame(fSpheres)
        fSpheresLine2.pack(side=TOP, fill=X)
        fSpheresLine3 = Frame(fSpheres)
        fSpheresLine3.pack(side=TOP, fill=X)

        Label(fSpheresLine1, text='Inserted spheres radii', font=self.top.font_Title).pack(side=LEFT)
        Label(fSpheresLine2, text='Minimum sphere radius:', width=20, font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Label(fSpheresLine3, text='Maximum sphere radius:', width=20, font=self.top.font_Text).pack(side=LEFT, anchor=W)

        EntryResult_L = Entry(fSpheresLine2, textvariable=self.Entry_L, background='white', width=6, justify=CENTER, font=self.top.font_Text)
        EntryResult_L.pack(side=RIGHT)
        EntryResult_U = Entry(fSpheresLine3, textvariable=self.Entry_U, background='white', width=6, justify=CENTER, font=self.top.font_Text)
        EntryResult_U.pack(side=RIGHT)

        args_list = [EntryResult_L, self.Entry_L, 0.1, EntryResult_U, 2, self.ValidEntry_L,'Min radius','float']
        EntryResult_L.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidEntry_L[3] = EntryResult_L

        args_list = [EntryResult_U, self.Entry_U, EntryResult_L, 5.0, 2, self.ValidEntry_U,'Max radius','float']
        EntryResult_U.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidEntry_U[3] = EntryResult_U

        #==================================================================================
        #                           OUTPUTTING OPTIONS
        #==================================================================================                
        # Option -r

        fOutput = Frame(self.fAdvOptions)
        fOutput.pack(side=TOP, fill=X, padx=5, pady=10)
        fOutputLine1 = Frame(fOutput)
        fOutputLine1.pack(side=TOP, fill=X)
        fOutputLine2 = Frame(fOutput)
        fOutputLine2.pack(side=TOP, fill=X)
        fOutputLine3 = Frame(fOutput)
        fOutputLine3.pack(side=TOP, fill=X)
        fOutputLine4 = Frame(fOutput)
        fOutputLine4.pack(side=TOP, fill=X)

        Label(fOutputLine1, text='Outputting options', font=self.top.font_Title).pack(side=LEFT, anchor=W)

        Label(fOutputLine2, text='Include all atoms of the residue', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Button(fOutputLine2, textvariable=self.Entry_R, width=6, bg='white', font=self.top.font_Text, command=self.Btn_Entry_R_Clicked).pack(side=RIGHT)
        
        #==================================================================================
        # Option -h

        Label(fOutputLine3, text = 'Output HET groups atoms', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Button(fOutputLine3, textvariable=self.Entry_h, width=6, bg='white', font=self.top.font_Text, command=self.Btn_Entry_h_Clicked).pack(side=RIGHT)

        #==================================================================================
        # Option -H

        Label(fOutputLine4, text = 'Output all atoms of HET groups', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        Button(fOutputLine4, textvariable=self.Entry_H, width=6, bg='white', font=self.top.font_Text, command=self.Btn_Entry_H_Clicked).pack(side=RIGHT)
                                        
        #==================================================================================
        #                           OUTPUTTING OPTIONS
        #==================================================================================                
        fOthers = Frame(self.fAdvOptions)
        fOthers.pack(side=TOP, fill=X, padx=5, pady=10)
        fOthersLine1 = Frame(fOthers)
        fOthersLine1.pack(side=TOP, fill=X)
        fOthersLine2 = Frame(fOthers)
        fOthersLine2.pack(side=TOP, fill=X)
        fOthersLine3 = Frame(fOthers)
        fOthersLine3.pack(side=TOP, fill=X)

        Label(fOthersLine1, text='Other options', font=self.top.font_Title).pack(side=LEFT, anchor=W)

        # Option -k

        Label(fOthersLine2, text='Threshold distance for contact definition', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        EntryResult_K = Entry(fOthersLine2, width=4, textvariable=self.Entry_K, background='white', justify=CENTER, font=self.top.font_Text)
        EntryResult_K.pack(side=RIGHT)
        args_list = [EntryResult_K, self.Entry_K, 0.1, 7.0, 2, self.ValidEntry_K,'Contact threshold','float']
        EntryResult_K.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidEntry_K[3] = EntryResult_K

        # Option -c
 
        Label(fOthersLine3, text = 'Chain ID to be considered', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        optionTuple = 'ALL',
        EntryResult_C = OptionMenu(*(fOthersLine3, self.Entry_C) + optionTuple)
        EntryResult_C.config(width=8, bg=self.top.Color_White, font=self.top.font_Text)
        EntryResult_C.pack(side=RIGHT)
        self.Entry_C.set('ALL')
        
        #==================================================================================
        '''                           --- BUTTONS AREA ---                              '''
        #================================================================================== 
        
        # Quit the advanced options menu
        Btn_Back = Button(self.fAdvOptions, text='Back', font=self.top.font_Text, relief=RIDGE, command=self.Btn_Back_Clicked)
        Btn_Back.pack(side=RIGHT)

    ''' ==================================================================================
    FUNCTION Btn_Back_Clicked: Goes back to the Default tab of GetCleft
    ==================================================================================  '''  
    def Btn_Back_Clicked(self):
        
        self.top.SetActiveFrame(self.top.Default)

    ''' ==================================================================================
    FUNCTION Btn_Entry_R_Clicked: (De)Activates Entry_R field
    ==================================================================================  '''  
    def Btn_Entry_R_Clicked(self):
        
        if self.Entry_R_Value.get():
            self.Entry_R.set('False')
            self.Entry_R_Value.set(0)
        else:
            self.Entry_R.set('True')
            self.Entry_R_Value.set(1)
            
    ''' ==================================================================================
    FUNCTION Btn_Entry_h_Clicked: (De)Activates Entry_h field
    ==================================================================================  '''  
    def Btn_Entry_h_Clicked(self):
        
        if self.Entry_h_Value.get():
            self.Entry_h.set('False')
            self.Entry_h_Value.set(0)
        else:
            self.Entry_h.set('True')
            self.Entry_h_Value.set(1)

    ''' ==================================================================================
    FUNCTION Btn_Entry_H_Clicked: (De)Activates Entry_H field
    ==================================================================================  '''  
    def Btn_Entry_H_Clicked(self):
        
        if self.Entry_H_Value.get():
            self.Entry_H.set('False')
            self.Entry_H_Value.set(0)
        else:
            self.Entry_H.set('True')
            self.Entry_H_Value.set(1)

    ''' ========================================================
                  Welcome message upon frame built
    ========================================================='''
    def Load_Message(self):

        self.DisplayMessage('', 0)
        self.DisplayMessage('  Opened the advanced options menu... ',0)        
