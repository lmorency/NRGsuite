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
@title: IsoCleft - Interface

@summary: This is the interface of IsoCleft application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Nov. 11, 2010
'''

from subprocess import Popen, PIPE
from Tkinter import *
from time import time

from pymol import cmd
from pymol import util

import tkFileDialog
import tkFont, os, shutil
import threading
import CfgIsoCleft

        
#=========================================================================================
'''                        ---   STARTING GETCLEFT  ---                                '''
#=========================================================================================     
class SimThread(threading.Thread):

    def __init__(self, top, value):
        threading.Thread.__init__(self)
        self.arg = value
        self.top = top        
        self.start()
        
    def run(self):        
        process = Popen(self.arg, shell=True, stderr=PIPE)
        process.wait()        
        
        if process.returncode != 0: # 0 = success, optional check
            
            the_outerr = process.stderr.read()

            duration = str(time() - self.top.timeStart)
            duration = self.top.FormatFloat(duration)
            
            self.top.DisplayMessage('  Process completed (' + duration + ' sec) ...', 0)
            self.top.DisplayMessage(the_outerr, 0)            
        else:

            self.top.Read_IsoCleftFile()            
            
            duration = str(time() - self.top.timeStart)
            duration = self.top.FormatFloat(duration)

            self.top.DisplayMessage('  Process completed (' + duration + ' sec) ...', 2)
            self.top.DisplayMessage('', 0)
        
        self.top.Btn_Start.config(state='normal')
            

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayIsoCleft:
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, top, UserPath, IsoCleftPath, GetCleftPath, RunDir_Path, OSid, ProtName):
        
        self.OSid = OSid                      
            
        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''
        self.RunDir_Path = RunDir_Path
        self.UserPath = UserPath
        self.UserIsoCleft = UserPath + 'IsoCleft/'
        self.UserGetCleft = UserPath + 'GetCleft/'
        self.SavePath = self.UserIsoCleft + 'Save/'
        self.TempPath = self.UserIsoCleft + 'Temp/'
        self.AdvOpt_Path = self.UserIsoCleft + 'AdvOption.cfg'
        
        self.ValidateFolders()        
             
        self.IsoCleft_Path = IsoCleftPath
        self.GetCleft_Path = GetCleftPath
        
        self.IsoCleftIsRunning()
        
        self.OutputPath = ''

        # CLEFT A and CLEFT B Area
        self.CleftPath_A = ''
        self.CleftName_A = ''
        self.CleftPath_B = ''
        self.CleftName_B = ''
        
        self.MidOrigin = 0
        self.ZoomVal = 0
        self.BiggestClf = ''
        
        # Superimposed Cleft A and B
        self.CleftName_Asi = ''
        self.CleftName_Bsi = ''
        self.CleftPath_Asi = ''     
        self.CleftPath_Bsi = ''
        
        self.CleftA_Coord = {}
        self.CleftB_Coord = {}
        
        self.Radius_A = 0.0
        self.Radius_B = 0.0
        
        self.IsoCleftRun = False
        self.CleftA_UsedAtoms = list()
        self.CleftB_UsedAtoms = list()
        
        self.Unused_AtomsCleftA = list()       
        self.Unused_AtomsCleftB = list()
        
        # GET THE INFORMATIONS NEEDED FOR THE ROTATION / TRANSLATION
        self.Center_CleftA_1 = list()
        self.Center_CleftB_1 = list()
        self.RotationMatrix_1 = list()
        
        self.Center_CleftA_2 = list()
        self.Center_CleftB_2 = list()
        self.RotationMatrix_2 = list()      
        
        # ADVEANCED OPTIONS Area
        self.DEFAULT_R = '5'
        self.DEFAULT_C = '3.50'
        self.DEFAULT_D = '4.00'
        self.DEFAULT_N = '4.00'
        self.DEFAULT_S = '4'
        self.DEFAULT_O = 'Test'
        
        self.OutputName = ''        
        
        self.Valid_ClfA = False
        self.Valid_ClfB = False       

        self.LAST_STEP = 1
        self.STEPS = IntVar()
        self.STEPS.set(self.LAST_STEP)         
        self.MsgLineCounter = 2
        
        self.ValidEntry_R = True
        self.ValidEntry_C = True
        self.ValidEntry_D = True
        self.ValidEntry_N = True
        self.ValidEntry_S = True
        self.ValidEntry_O = True
        
        self.Entry_R = StringVar()
        self.Entry_C = StringVar()        
        self.Entry_D = StringVar()
        self.Entry_N = StringVar()
        self.Entry_S = StringVar()
        self.Entry_O = StringVar()      # Output Filename
        
        if (os.path.isfile(self.AdvOpt_Path)):
            self.Read_AdvOptionFile()
        else:
            self.Entry_R.set(self.DEFAULT_R)
            self.Entry_C.set(self.DEFAULT_C)
            self.Entry_D.set(self.DEFAULT_D)
            self.Entry_N.set(self.DEFAULT_N)
            self.Entry_S.set(self.DEFAULT_S)
            self.Entry_O.set(self.DEFAULT_O)
       
        #================================================================================== 
        #                 SET the default fonts of the interface
        #==================================================================================
        FontType = CfgIsoCleft.GetFontType()
        FontSize = CfgIsoCleft.GetFontSize()     
        
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)        
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_U = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        
        self.Color_Gray = '#EDEDED'
        self.Color_White = '#FFFFFF'
        self.Color_Red = '#FF9999'
        self.Color_Green = '#CCFFCC'
        
        #================================================================================== 
        #                        Initialize the window
        #================================================================================== 
        self.top = top
        self.top.title('IsoCleft')
        self.top.geometry('400x420')   # Interface DIMENSIONS
        self.top.maxsize(400,420)
        self.top.minsize(400,420)
        self.top.protocol('WM_DELETE_WINDOW', self.Btn_Quit_Clicked)
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================
        self.frame = Frame(top)  # Main frame
        self.frame.pack(expand='true')
        
        #Add the Menu
        self.makeMenuBar()

        self.Frame_Main()
        

    ''' ==================================================================================
    FUNCTION makeMenuBar: Display a MenuBar in the GetCleft application    
    ==================================================================================  '''     
    def makeMenuBar(self): 
         self.menubar = Menu(self.top) 
         self.top.config(menu=self.menubar) 
         pulldown1 = Menu(self.menubar, tearoff=0) 
         pulldown1.add_command(label='Load Previous Results...', command=self.Btn_LoadResult) 
         pulldown1.add_command(label='Save Results...', command=self.Btn_SaveResult)
         pulldown1.add_separator() 
         pulldown1.add_command(label='Close', command=self.Btn_Quit_Clicked) 
         self.menubar.add_cascade(label='File', underline=0, menu=pulldown1)
         
         pulldown2 = Menu(self.menubar, tearoff=0) 
         pulldown2.add_command(label='Set Advanced Options ', command=self.Frame_AdvancedOptions)
         self.menubar.add_cascade(label='Edit', underline=0, menu=pulldown2)     


    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''        
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        self.fTop = Frame(self.frame, relief=RIDGE, border=0, width=400, height=330)        
        
        self.Frame_Default()
        
        self.fTop.pack(fill=X, expand=True, padx=10)
        self.fTop.pack_propagate(0) 
        
        #==================================================================================
        '''                  --- BOTTOM FRAME OF THE INTERFACE ---                          '''
        #==================================================================================
        
        fBottom = Frame(self.frame, relief=RIDGE, border=0, width=400, height=90)
        
        # Messages Box
        fMsg = Frame(fBottom, border=1, width=380, height=80, relief=SUNKEN)        

        scrollBar = Scrollbar(fMsg)
        self.TextMessage = Text(fMsg, width=380, height=5, background=self.Color_Gray)
        scrollBar.pack(side=RIGHT, fill=Y)
        self.TextMessage.pack(side=LEFT, anchor=S)
        scrollBar.config(command=self.TextMessage.yview)
        self.TextMessage['font'] = self.font_Text
        self.TextMessage.config(state='disabled', yscrollcommand=scrollBar.set)

        fMsg.pack(fill=X, expand=True, side=LEFT, anchor=S, padx=10, pady= 5)
        fMsg.pack_propagate(0)

        fBottom.pack(fill=X, expand=True, side=BOTTOM, anchor=S)
        fBottom.pack_propagate(0)
        

    ''' ==================================================================================
    FUNCTION Frame_Default: Default Displayed Interface 
    ==================================================================================  '''    
    def Frame_Default(self):
        
        # First Section (from Cleft to Start)
        self.fTopPDB = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=330)        
        fTopPDB1 = Frame(self.fTopPDB, relief=RIDGE, border=0, width=400, height=45)
        fTopPDB2 = Frame(self.fTopPDB, relief=RIDGE, border=0, width=400, height=40)
        fTopPDB3 = Frame(self.fTopPDB, relief=RIDGE, border=0, width=400, height=35)
        fTopPDB4 = Frame(self.fTopPDB, relief=RIDGE, border=0, width=400, height=45)       
        fTopOpt = Frame(self.fTopPDB, relief=SUNKEN, border=1, width=400, height=115)       
        
        lblClfA = Label(fTopPDB1, text='Cleft A:')
        lblClfA.pack(padx=5, pady=5, side=LEFT, anchor=SW)
        lblClfA['font'] = self.font_Title
        
        # Add the DROP DOWN LIST
        optionTuple = ('',)
        self.defaultOpt_A = StringVar()
        self.optionMenuWidget_A = apply(OptionMenu, (fTopPDB1, self.defaultOpt_A) + optionTuple)

        self.defaultOpt_A.set('')
        self.optionMenuWidget_A['width'] = 30
        self.optionMenuWidget_A['background'] = self.Color_Gray
        self.optionMenuWidget_A.pack(side=LEFT, anchor=SW, padx=5)       
        
        fTopPDB1.pack(fill=X, expand=True, side=TOP)
        fTopPDB1.pack_propagate(0)
        
        pad1 = Label(fTopPDB2, text='')
        pad1.pack(padx=10, side=LEFT, anchor=NW)
        
        lblLoadA = Label(fTopPDB2, text='Load from:')
        lblLoadA.pack(padx=5, pady=5, side=LEFT, anchor=NW)
        lblLoadA['font'] = self.font_Text_U

        # Refresh the list with the selections in Pymol
        Btn_LoadPymol_A = Button(fTopPDB2, text='Pymol Selections', width=15, command=self.Btn_LoadPymol_A_Clicked)
        Btn_LoadPymol_A.pack(side=LEFT, anchor=NW, padx=5)
        Btn_LoadPymol_A['font'] = self.font_Text
        
        # Refresh the list with a file selection
        Btn_LoadFile_A = Button(fTopPDB2, text='File', width=5, command=self.Btn_LoadFile_A_Clicked)
        Btn_LoadFile_A.pack(side=LEFT, anchor=NW)
        Btn_LoadFile_A['font'] = self.font_Text       
        
        fTopPDB2.pack(fill=X, expand=True, side=TOP)
        fTopPDB2.pack_propagate(0)

        lblClfB = Label(fTopPDB3, text='Cleft B:')
        lblClfB.pack(padx=5, pady=5, side=LEFT, anchor=SW)
        lblClfB['font'] = self.font_Title
        
        # Add the DROP DOWN LIST
        optionTuple = ('',)
        self.defaultOpt_B = StringVar()
        self.optionMenuWidget_B = apply(OptionMenu, (fTopPDB3, self.defaultOpt_B) + optionTuple)

        self.defaultOpt_B.set('')
        self.optionMenuWidget_B['width'] = 30
        self.optionMenuWidget_B['background'] = self.Color_Gray
        self.optionMenuWidget_B.pack(side=LEFT, anchor=SW, padx=5)       
        
        fTopPDB3.pack(fill=X, expand=True, side=TOP)
        fTopPDB3.pack_propagate(0)
        
        pad2 = Label(fTopPDB4, text='')
        pad2.pack(padx=10, side=LEFT, anchor=NW)

        lblLoadB = Label(fTopPDB4, text='Load from:')
        lblLoadB.pack(padx=5, pady=5, side=LEFT, anchor=NW)
        lblLoadB['font'] = self.font_Text_U
        
        # Refresh the list with the selections in Pymol
        Btn_LoadPymol_B = Button(fTopPDB4, text='Pymol Selections', width=15, command=self.Btn_LoadPymol_B_Clicked)
        Btn_LoadPymol_B.pack(side=LEFT, anchor=NW, padx=5)
        Btn_LoadPymol_B['font'] = self.font_Text
        
        # Refresh the list with a file selection
        Btn_LoadFile_B = Button(fTopPDB4, text='File', width=5, command=self.Btn_LoadFile_B_Clicked)
        Btn_LoadFile_B.pack(side=LEFT, anchor=NW)
        Btn_LoadFile_B['font'] = self.font_Text         
        
        fTopPDB4.pack(fill=X, expand=True, side=TOP)
        fTopPDB4.pack_propagate(0)       

        
        #==================================================================================
        '''          --- DISPLAY THE DISPLAY'S OPTIONS IN THE INTERFACE ---             '''
        #==================================================================================
        lblOptDisp = Label(fTopOpt, text='Display Steps')
        lblOptDisp.pack(padx=30, pady=5, side=TOP, anchor=W)
        lblOptDisp['font'] = self.font_Title
        
        self.RadioBtn_OP1 = Radiobutton(fTopOpt, text=' Step 1:   Starting Clefts', variable=self.STEPS, value=1, command=self.Step1_Clicked)
        self.RadioBtn_OP1.pack(padx=50, side=TOP, anchor=W)
        self.RadioBtn_OP1['font'] = self.font_Text            
        
        self.RadioBtn_OP2 = Radiobutton(fTopOpt, text=' Step 2:   Filtered Atoms (Unused in Grey)', variable=self.STEPS, value=2, command=self.Step2_Clicked)
        self.RadioBtn_OP2.pack(padx=50, pady=3, side=TOP, anchor=W)
        self.RadioBtn_OP2['font'] = self.font_Text      
        
        self.RadioBtn_OP3 = Radiobutton(fTopOpt, text=' Step 3:   Superimposed Atoms', variable=self.STEPS, value=3, command=self.Step3_Clicked)
        self.RadioBtn_OP3.pack(padx=50, side=TOP, anchor=W)
        self.RadioBtn_OP3['font'] = self.font_Text
        
        if not(self.IsoCleftRun): 
            self.RadioBtn_OP1.config(state='disable')
            self.RadioBtn_OP2.config(state='disable')      
            self.RadioBtn_OP3.config(state='disable')        
        
        fTopOpt.pack(fill=X, expand=True, side=TOP)
        fTopOpt.pack_propagate(0)
        
        #==================================================================================
        '''                           --- BUTTONS AREA ---                              '''
        #================================================================================== 
        
        fBtnOpt = Frame(self.fTopPDB, relief=RIDGE, border=0, width=400, height=50)
        
        # Refresh the list with a file selection
        self.Btn_Start = Button(fBtnOpt, text='Start', width=7, command=self.Btn_Start_Clicked)
        self.Btn_Start.pack(side=LEFT, anchor=W)
        self.Btn_Start['font'] = self.font_Text
        self.Btn_Start.config(state='disable')
        
        # Clear Pymol elements
        Btn_Clear = Button(fBtnOpt, text='Clear', width=7, command=self.Btn_Clear_Clicked)
        Btn_Clear.pack(side=LEFT, anchor=W, padx=5)
        Btn_Clear['font'] = self.font_Text
        
        # Quit the application
        Btn_Quit = Button(fBtnOpt, text='Close', width=7, command=self.Btn_Quit_Clicked)
        Btn_Quit.pack(side=RIGHT, anchor=E)
        Btn_Quit['font'] = self.font_Text
        
        fBtnOpt.pack(fill=X, expand=True, side=TOP)
        fBtnOpt.pack_propagate(0)
        
        self.fTopPDB.pack(fill=X, expand=True, side=TOP)
        self.fTopPDB.pack_propagate(0)
        
        
        # Update the informations if open after the Advanced Options Menu
        if self.CleftName_A != '':
            self.Refresh_DDL_A()
            
        if self.CleftName_B != '':
            self.Refresh_DDL_B()
            
        if self.Valid_ClfA and self.Valid_ClfB:
            self.Btn_Start.config(state='normal')
        
        
    ''' ==================================================================================
    FUNCTION Frame_AdvancedOptions: Display the Advanced Options in the Interface 
    ==================================================================================  '''    
    def Frame_AdvancedOptions(self):
        
        self.fTopPDB.pack_forget()
        self.Read_AdvOptionFile()
        
        #==================================================================================
        '''             --- DISPLAY OPTIONS FRAME OF THE INTERFACE ---                  '''
        #==================================================================================
        self.fGlbOpt = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=280)
        
        fGlbTitle = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=45)
        fGlb1 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=25)
        fGlb2 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=25)
        fGlb3 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=25)
        fGlb4 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=25)
        fGlb5 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=25)
        fGlb6 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=50)
        fGlb7 = Frame(self.fGlbOpt, relief=RIDGE, border=0, width=400, height=50)
        
        lblOutput = Label(fGlbTitle, text='Advanced Options [range]')
        lblOutput.pack(side=LEFT, anchor=SW, pady = 10)
        lblOutput['font'] = self.font_Title
        
        fGlbTitle.pack(fill=X, expand=True, side=TOP)
        fGlbTitle.pack_propagate(0)        
        
        # Option -r
        padOpt1L = Label(fGlb1, text='')
        padOpt1L.pack(padx=5, side=LEFT, anchor=W)
        padOpt1L['font'] = self.font_Text
        
        ChkOpt1 = Label(fGlb1, text = '  JTT rank Threshold [1 - 20]', height=5)
        ChkOpt1.pack(side=LEFT, anchor=W)
        ChkOpt1['font'] = self.font_Text
        
        padOpt1R = Label(fGlb1, text='')
        padOpt1R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt1R['font'] = self.font_Text
        
        self.EntryResult_R = Entry(fGlb1, width=8, textvariable=self.Entry_R, background='white', justify=CENTER)
        self.EntryResult_R.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_R['font'] = self.font_Text
        self.EntryResult_R.bind('<KeyRelease>', self.Validate_Entry_R)
        
        fGlb1.pack(fill=X, expand=True, side=TOP)
        fGlb1.pack_propagate(0)
        
        # Option -c
        padOpt2L = Label(fGlb2, text='')
        padOpt2L.pack(padx=5, side=LEFT, anchor=W)
        padOpt2L['font'] = self.font_Text
        
        ChkOpt2 = Label(fGlb2, text = '  C-alpha Node Distance Threshold [0.5 - 5.0]', height=5)
        ChkOpt2.pack(side=LEFT, anchor=W)
        ChkOpt2['font'] = self.font_Text
        
        padOpt2R = Label(fGlb2, text='')
        padOpt2R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt2R['font'] = self.font_Text
        
        self.EntryResult_C = Entry(fGlb2, width=8, textvariable=self.Entry_C, background='white', justify=CENTER)
        self.EntryResult_C.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_C['font'] = self.font_Text
        self.EntryResult_C.bind('<KeyRelease>', self.Validate_Entry_C)
        
        fGlb2.pack(fill=X, expand=True, side=TOP)
        fGlb2.pack_propagate(0)
        
        # Option -n
        padOpt3L = Label(fGlb3, text='')
        padOpt3L.pack(padx=5, side=LEFT, anchor=W)
        padOpt3L['font'] = self.font_Text
        
        ChkOpt3 = Label(fGlb3, text = '  Neighborhood sphere radius [0.5 - 6.0]', height=5)
        ChkOpt3.pack(side=LEFT, anchor=W)
        ChkOpt3['font'] = self.font_Text
        
        padOpt3R = Label(fGlb3, text='')
        padOpt3R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt3R['font'] = self.font_Text
        
        self.EntryResult_N = Entry(fGlb3, width=8, textvariable=self.Entry_N, background='white', justify=CENTER)
        self.EntryResult_N.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_N['font'] = self.font_Text
        self.EntryResult_N.bind('<KeyRelease>', self.Validate_Entry_N)
        
        fGlb3.pack(fill=X, expand=True, side=TOP)
        fGlb3.pack_propagate(0)
        
        # Option -d
        padOpt4L = Label(fGlb4, text='')
        padOpt4L.pack(padx=5, side=LEFT, anchor=W)
        padOpt4L['font'] = self.font_Text
        
        ChkOpt4 = Label(fGlb4, text='  Node Distance Threshold [0.5 - 5.0]', height=5)
        ChkOpt4.pack(side=LEFT, anchor=W)
        ChkOpt4['font'] = self.font_Text
        
        padOpt4R = Label(fGlb4, text='')
        padOpt4R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt4R['font'] = self.font_Text
        
        self.EntryResult_D = Entry(fGlb4, width=8, textvariable=self.Entry_D, background='white', justify=CENTER)
        self.EntryResult_D.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_D['font'] = self.font_Text
        self.EntryResult_D.bind('<KeyRelease>', self.Validate_Entry_D)
        
        fGlb4.pack(fill=X, expand=True, side=TOP)
        fGlb4.pack_propagate(0)
        
        # Option -s
        padOpt5L = Label(fGlb5, text='')
        padOpt5L.pack(padx=5, side=LEFT, anchor=W)
        padOpt5L['font'] = self.font_Text
        
        ChkOpt5 = Label(fGlb5, text = '  Level of BK search simplification [0 - 8]', height=5)
        ChkOpt5.pack(side=LEFT, anchor=W)
        ChkOpt5['font'] = self.font_Text
        
        padOpt5R = Label(fGlb5, text='')
        padOpt5R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt5R['font'] = self.font_Text

        self.EntryResult_S = Entry(fGlb5, width=8, textvariable=self.Entry_S, background='white', justify=CENTER)
        self.EntryResult_S.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_S['font'] = self.font_Text
        self.EntryResult_S.bind('<KeyRelease>', self.Validate_Entry_S)
                
        fGlb5.pack(fill=X, expand=True, side=TOP)
        fGlb5.pack_propagate(0)
        
        # Option -o
        padOpt6L = Label(fGlb6, text='')
        padOpt6L.pack(padx=5, side=LEFT, anchor=W)
        padOpt6L['font'] = self.font_Text
        
        ChkOpt6 = Label(fGlb6, text = '  Output Filename', height=5)
        ChkOpt6.pack(side=LEFT, anchor=W)
        ChkOpt6['font'] = self.font_Text
        
        padOpt6R = Label(fGlb6, text='')
        padOpt6R.pack(padx=5, side=RIGHT, anchor=E)
        padOpt6R['font'] = self.font_Text
        
        self.EntryResult_O = Entry(fGlb6, width=25, textvariable=self.Entry_O, background='white', justify=CENTER)
        self.EntryResult_O.pack(side=RIGHT, anchor=E, padx=5)
        self.EntryResult_O['font'] = self.font_Text
        self.EntryResult_O.bind('<KeyRelease>', self.Validate_Entry_O)
        
        fGlb6.pack(fill=X, expand=True, side=TOP)
        fGlb6.pack_propagate(0)
        
        fGlb7.pack(fill=X, expand=True, side=TOP)
        fGlb7.pack_propagate(0)        
        
        self.fGlbOpt.pack(fill=X, expand=True, side=TOP, padx=10)
        self.fGlbOpt.pack_propagate(0)
        
        #==================================================================================
        '''                           --- BUTTONS AREA ---                              '''
        #================================================================================== 
        
        self.fBtnOpt = Frame(self.fTop, relief=RIDGE, border=0, width=400, height=50)
        
        # Restore to the Defaults values
        Btn_Restore = Button(self.fBtnOpt, text='Restore Defaults', width=14, command=self.Btn_Restore_Clicked)
        Btn_Restore.pack(side=LEFT, anchor=W)
        Btn_Restore['font'] = self.font_Text
        
        # Apply the Changed Done
        Btn_Apply = Button(self.fBtnOpt, text='Apply', width=7, command=self.Btn_Apply_Clicked)
        Btn_Apply.pack(side=LEFT, anchor=W, padx=5)
        Btn_Apply['font'] = self.font_Text
        
        # Quit the advanced options menu
        Btn_Done = Button(self.fBtnOpt, text='Done', width=7, command=self.Btn_Done_Clicked)
        Btn_Done.pack(side=RIGHT, anchor=E)
        Btn_Done['font'] = self.font_Text
        
        self.fBtnOpt.pack(fill=X, expand=True, side=TOP)
        self.fBtnOpt.pack_propagate(0)


    ''' ==================================================================================
    FUNCTION Btn_Restore_Clicked: Restore the Fields to their default value. 
    ==================================================================================  '''         
    def Btn_Restore_Clicked(self):
        
        self.DisplayMessage('  Advanced Options restored to their Defaults values.', 0)
        
        self.Entry_R.set(self.DEFAULT_R)
        self.Entry_C.set(self.DEFAULT_C)
        self.Entry_D.set(self.DEFAULT_D)
        self.Entry_N.set(self.DEFAULT_N)
        self.Entry_S.set(self.DEFAULT_S)
        self.Entry_O.set(self.DEFAULT_O)
        
        if os.path.isfile(self.AdvOpt_Path):
            try:
                os.remove(self.AdvOpt_Path)
            except OSError:
                time.sleep(0.1)
                os.remove(self.AdvOpt_Path)     
        
        
    ''' ==================================================================================
    FUNCTION Btn_Apply_Clicked: Apply the Change done in the fields. 
    ==================================================================================  '''         
    def Btn_Apply_Clicked(self):
        
        if self.ValidEntry_R and self.ValidEntry_S and self.ValidEntry_D and \
           self.ValidEntry_N and self.ValidEntry_C and self.ValidEntry_O:
        
            self.DisplayMessage('  SAVED the changes done...', 0)        
            
            write_file = open(self.AdvOpt_Path, 'w')
            write_file.write('Entry_R: ' + self.Entry_R.get() + '\n')
            write_file.write('Entry_S: ' + self.Entry_S.get() + '\n')
            write_file.write('Entry_D: ' + self.Entry_D.get() + '\n')
            write_file.write('Entry_N: ' + self.Entry_N.get() + '\n')
            write_file.write('Entry_C: ' + self.Entry_C.get() + '\n')
            write_file.write('Entry_O: ' + self.Entry_O.get() + '\n')
            write_file.close()
            
        else:
            self.DisplayMessage('  Can''t apply the changed while their is Error(s) or Empty Field(s)!', 1)
        
        
    ''' ==================================================================================
    FUNCTION Read_AdvOptionFile: Read and Set the Fields to their Prefered Values. 
    ==================================================================================  '''   
    def Read_AdvOptionFile(self):
        
        if (os.path.isfile(self.AdvOpt_Path)):
        
            # Read the Output file
            file = open(self.AdvOpt_Path, 'r')
            txtFile = file.readlines()
            file.close()
            
            for line in txtFile:
                if (line[0:7] == 'Entry_R'):
                    self.Entry_R.set(line[9:].strip())
                elif (line[0:7] == 'Entry_S'):
                    self.Entry_S.set(line[9:].strip())
                elif (line[0:7] == 'Entry_D'):
                    self.Entry_D.set(line[9:].strip())
                elif (line[0:7] == 'Entry_N'):
                    self.Entry_N.set(line[9:].strip())
                elif (line[0:7] == 'Entry_C'):
                    self.Entry_C.set(line[9:].strip())
                elif (line[0:7] == 'Entry_O'):
                    self.Entry_O.set(line[9:].strip())
                             
        
    ''' ==================================================================================
    FUNCTION Btn_Done_Clicked: Close the Advanced Options Menu 
    ==================================================================================  '''         
    def Btn_Done_Clicked(self):
        
        if self.ValidEntry_R and self.ValidEntry_S and self.ValidEntry_D and \
           self.ValidEntry_N and self.ValidEntry_C and self.ValidEntry_O:
            pass
        else:
            self.DisplayMessage('  Values SET to DEFAULTS because there is Empty Field(s) or Error(s)... ',0)
            self.Btn_Restore_Clicked()
        
        self.fGlbOpt.pack_forget()
        self.fBtnOpt.pack_forget()
        
        self.Frame_Default()           

    ''' ==================================================================================
    FUNCTION Btn_LoadResult: Load a previous results from a folder 
    ==================================================================================  '''         
    def Btn_LoadResult(self):
        self.DisplayMessage('  Load a previous Results... ',0)


    ''' ==================================================================================
    FUNCTION Btn_SaveResult: Save the actual results in a folder 
    ==================================================================================  '''        
    def Btn_SaveResult(self):
        self.DisplayMessage('  Save the actual Results... ',0)
        
    
    ''' ==================================================================================
    FUNCTION Btn_LoadPymol_A_Clicked: Load the Cleft from Pymol. 
    ==================================================================================  '''    
    def Btn_LoadPymol_A_Clicked(self):
        
        self.DisplayMessage('',0)
        self.DisplayMessage('  Loading the Cleft from Pymol Selections (Cleft_A).',0)
        self.Refresh_DDL_A()
        if len(self.defaultOpt_A.get()) > 0:
            self.DisplayMessage('  Selecting the Cleft_A: ' + self.defaultOpt_A.get(),0)
            self.Copy_PDB_FileA(self.defaultOpt_A.get(), '')
            
            if self.Valid_ClfB:
                self.Btn_Start.config(state='normal')
                    
            if self.Valid_ClfA:
                self.Set_CleftCoordinates(self.CleftName_A, self.CleftPath_A, 'A')
                self.Set_CleftDisplay(self.CleftName_A, self.CleftPath_A, True)
                self.Refresh_DDL_A()
                self.DisplayMessage('  Load Completed (Cleft_A): ' + self.CleftName_A,0)
             
        

    ''' ==================================================================================
    FUNCTION Refresh_DDL_A: Refreshing the DropDownList items. 
    ==================================================================================  '''        
    def Refresh_DDL_A(self):
        
        #Clean the Drop Down List menu
        self.defaultOpt_A.set('')
        self.optionMenuWidget_A['menu'].delete(0, 'end')
                  
        # Get the Selection Names in Pymol        
        sess = cmd.get_names('all')
        
        for i in sess:
            if (i.count('_clf_') > 0) and (i.count('_S') == 0) and (i.count('clqa') == 0):
                
                # Display the model name              
                self.optionMenuWidget_A['menu'].add_command(label=i, command=lambda temp = i: self.NewDDL_SelectionA(temp))
                
                if self.CleftName_A != '':
                    self.defaultOpt_A.set(self.CleftName_A)
                else:                        
                    self.defaultOpt_A.set(i) 
            
        if len(self.defaultOpt_A.get()) == 0:
            self.DisplayMessage('  Attention: There is NO Cleft in the Pymol Selections.',1)
            self.Valid_ClfA = False
            self.Btn_Start.config(state='disabled')


    ''' ==================================================================================
    FUNCTION NewDDL_SelectionA: New selection in the Drop Down List. 
    ==================================================================================  '''           
    def NewDDL_SelectionA(self, temp):
        
        self.DisplayMessage('',0)
        
        if len(temp) > 0:
            self.DisplayMessage('  Selecting the Cleft_A: ' + temp,0)
            self.Copy_PDB_FileA(temp, '')
            
            if self.Valid_ClfB:
                self.Btn_Start.config(state='normal')
           
            if self.Valid_ClfA:
                self.Set_CleftCoordinates(self.CleftName_A, self.CleftPath_A, 'A')
                self.Set_CleftDisplay(self.CleftName_A, self.CleftPath_A, True)
                self.Refresh_DDL_A()
                self.DisplayMessage('  Load Completed (Cleft_A): ' + self.CleftName_A,0)
                            

    ''' ==================================================================================
    FUNCTION NewDDL_SelectionB: New selection in the Drop Down List. 
    ==================================================================================  '''           
    def NewDDL_SelectionB(self, temp):
        
        self.DisplayMessage('',0)
        
        if len(temp) > 0:
            self.DisplayMessage('  Selecting the Cleft_B: ' + temp,0)
            self.Copy_PDB_FileB(temp, '')
            
            if self.Valid_ClfA:
                self.Btn_Start.config(state='normal')
           
            if self.Valid_ClfB:
                self.Set_CleftCoordinates(self.CleftName_B, self.CleftPath_B, 'B')
                self.Set_CleftDisplay(self.CleftName_B, self.CleftPath_B, True)
                self.Refresh_DDL_B()
                self.DisplayMessage('  Load Completed (Cleft_B): ' + self.CleftName_B,0)
            
        
    ''' ==================================================================================
    FUNCTION Btn_LoadPymol_B_Clicked: Load the Cleft from Pymol. 
    ==================================================================================  '''    
    def Btn_LoadPymol_B_Clicked(self):

        self.DisplayMessage('',0)
        self.DisplayMessage('  Loading the Cleft from Pymol Selections (Cleft_B).',0)
        self.Refresh_DDL_B()
        
        if len(self.defaultOpt_B.get()) > 0:
            self.DisplayMessage('  Selecting the Cleft_B: ' + self.defaultOpt_B.get(),0)
            self.Copy_PDB_FileB(self.defaultOpt_B.get(), '')
            
            if self.Valid_ClfA:
                self.Btn_Start.config(state='normal')
           
            if self.Valid_ClfB:
                self.Set_CleftCoordinates(self.CleftName_B, self.CleftPath_B, 'B')
                self.Set_CleftDisplay(self.CleftName_B, self.CleftPath_B, True)
                self.Refresh_DDL_B()
                self.DisplayMessage('  Load Completed (Cleft_B): ' + self.CleftName_B,0)             
        

    ''' ==================================================================================
    FUNCTION Refresh_DDL_B: Refreshing the DropDownList items. 
    ==================================================================================  '''        
    def Refresh_DDL_B(self):
        
        #Clean the Drop Down List menu
        self.defaultOpt_B.set('')
        self.optionMenuWidget_B['menu'].delete(0, 'end')
                  
        # Get the Selection Names in Pymol        
        sess = cmd.get_names('all')
        NotFound = True
        
        for i in sess:
            if (i.count('_clf_') > 0) and (i.count('_S') == 0) and (i.count('clqa') == 0):
                # Display the model name              
                self.optionMenuWidget_B['menu'].add_command(label=i, command=lambda temp = i: self.NewDDL_SelectionB(temp))
                
                if self.CleftName_B != '':
                    self.defaultOpt_B.set(self.CleftName_B)
                else:                        
                    self.defaultOpt_B.set(i)
                               
        if len(self.defaultOpt_B.get()) == 0:
            self.DisplayMessage('  Attention: There is NO Cleft in the Pymol Selections.',1)
            self.Valid_ClfB = False
            self.Btn_Start.config(state='disabled')
        
    
    ''' ==================================================================================
    FUNCTION Btn_LoadFile_A_Clicked: Load the Cleft from a File. 
    ==================================================================================  '''    
    def Btn_LoadFile_A_Clicked(self):
        
        self.DisplayMessage('',0)
        self.DisplayMessage('  Loading the Cleft from a file (Cleft_A).',0)
        
        FilePath = tkFileDialog.askopenfilename(filetypes=[('_clf_ File','*.pdb')], initialdir=self.UserGetCleft, title='Select a CLEFT File to Import')
        
        if len(FilePath)>0:           
            self.CleftPath_A = FilePath
            self.CleftName_A = os.path.splitext(self.CleftPath_A)[0]
            self.CleftName_A = os.path.basename(self.CleftName_A)            
            
            if (self.CleftName_A.count('_clf_') > 0) and (self.CleftName_A.count('_S') == 0):         
            
                self.Copy_PDB_FileA(self.CleftName_A, self.CleftPath_A)
                self.Set_CleftCoordinates(self.CleftName_A, self.CleftPath_A, 'A')
                self.Set_CleftDisplay(self.CleftName_A, self.CleftPath_A, True) 
                self.Refresh_DDL_A()    # Refresh the Drop Down List of the Cleft A              
            
                if self.Valid_ClfB:
                    self.Btn_Start.config(state='normal')
                    
                if self.Valid_ClfA:
                    self.DisplayMessage('  Load Completed (Cleft_A): ' + self.CleftName_A,0)                
                
            else:
                self.Valid_ClfA = False
                self.Btn_Start.config(state='disabled')
                self.CleftPath_A = ''
                self.CleftName_A = ''
                
                self.DisplayMessage('  Cancelled the loading Cleft...',0)
                self.DisplayMessage('  Attention: This is NOT a valid Cleft File!',1)
            
        else:
            self.Valid_ClfA = False
            self.Btn_Start.config(state='disabled')
            self.DisplayMessage('  Cancelled the loading Cleft...',0)
        
        
    ''' ==================================================================================
    FUNCTION Btn_LoadFile_B_Clicked: Load the Cleft from a File. 
    ==================================================================================  '''    
    def Btn_LoadFile_B_Clicked(self):
        
        self.DisplayMessage('',0)
        self.DisplayMessage('  Loading the Cleft from a file (Cleft_B).',0)
        
        FilePath = tkFileDialog.askopenfilename(filetypes=[('_clf_ File','*.pdb')], initialdir=self.UserGetCleft, title='Select a CLEFT File to Import')
        
        if len(FilePath)>0:           
            self.CleftPath_B = FilePath
            self.CleftName_B = os.path.splitext(self.CleftPath_B)[0]
            self.CleftName_B = os.path.basename(self.CleftName_B)
            
            if (self.CleftName_B.count('_clf_')) > 0 and (self.CleftName_B.count('_S') == 0):      
                
                self.Copy_PDB_FileB(self.CleftName_B, self.CleftPath_B)
                self.Set_CleftCoordinates(self.CleftName_B, self.CleftPath_B, 'B')
                self.Set_CleftDisplay(self.CleftName_B, self.CleftPath_B, True)

                self.Refresh_DDL_B()     # Refresh the Drop Down List of the Cleft A
                
                if self.Valid_ClfA:
                    self.Btn_Start.config(state='normal')
                    
                if self.Valid_ClfB:
                    self.DisplayMessage('  Load Completed (Cleft_B): ' + self.CleftName_B,0)
                
            else:
                self.Valid_ClfB = False
                self.Btn_Start.config(state='disabled')
                self.CleftPath_B = ''
                self.CleftName_B = ''
                
                self.DisplayMessage('  Cancelled the loading Cleft...',0)
                self.DisplayMessage('  Attention: This is NOT a valid Cleft File!',1)
            
        else:
            self.Valid_ClfB = False
            self.Btn_Start.config(state='disabled')
            self.DisplayMessage('  Cancelled the loading Cleft...',0)
            
            
    ''' ==================================================================================
    FUNCTION Copy_PDB_FileA: Make a copy of the PDB File 
    ==================================================================================  '''   
    def Copy_PDB_FileA(self, PDB_Name, PDB_Path):
        
        if (len(PDB_Name)>0):

            if (PDB_Name.count('_r') > 0):
                self.CleftName_A = PDB_Name[:-1] + 'a'
                self.CleftName_Asi = PDB_Name[:-2] + 'clqa'      
            else:
                self.CleftName_A = PDB_Name + '_ra'
                self.CleftName_Asi = PDB_Name + '_clqa'
                    
            if self.CleftName_B == self.CleftName_A:       
                    self.CleftName_A = self.CleftName_A[:-3] + 'A_ra'
                    self.CleftName_Asi = self.CleftName_A[:-3] + 'A_clqa'
                
            self.CleftPath_A = self.TempPath + self.CleftName_A + '.pdb'
            self.CleftPath_Asi = self.TempPath + self.CleftName_Asi + '.pdb'
            
            self.CleanTempDirectory(self.TempPath)

            # Copy the Cleft File in the Temp Folder
            if PDB_Path != '':                
                shutil.copy(PDB_Path, self.CleftPath_A)          
            else:
                cmd.save(self.CleftPath_A, PDB_Name, 1, 'pdb')
            
            self.Valid_ClfA = True
            self.defaultOpt_A.set(self.CleftName_A)
        else:
            self.Valid_ClfA = False
            self.Btn_Start.config(state='disabled')
            
            
    ''' ==================================================================================
    FUNCTION Copy_PDB_FileB: Make a copy of the PDB File 
    ==================================================================================  '''   
    def Copy_PDB_FileB(self, PDB_Name, PDB_Path):
        
        if (len(PDB_Name)>0):
            
            if (PDB_Name.count('_r') > 0):
                self.CleftName_B = PDB_Name[:-1] + 'a'
                self.CleftName_Bsi = PDB_Name[:-2] + 'clqa'              
            else:
                self.CleftName_B = PDB_Name + '_ra'
                self.CleftName_Bsi = PDB_Name + '_clqa'
                    
            if self.CleftName_A == self.CleftName_B:       
                self.CleftName_B = self.CleftName_B[:-3] + 'B_ra'
                self.CleftName_Bsi = self.CleftName_B[:-3] + 'B_clqa'
                
            self.CleftPath_B = self.TempPath + self.CleftName_B + '.pdb'
            self.CleftPath_Bsi = self.TempPath + self.CleftName_Bsi + '.pdb'
            
            self.CleanTempDirectory(self.TempPath)
            
            # Copy the Cleft File in the Temp Folder
            if PDB_Path != '':                
                shutil.copy(PDB_Path, self.CleftPath_B)
            else:
                cmd.save(self.CleftPath_B, PDB_Name, 1, 'pdb')                
            
            self.Valid_ClfB = True
            self.defaultOpt_B.set(self.CleftName_B)
                            
        else:
            self.Valid_ClfB = False
            self.Btn_Start.config(state='disabled')            
            
            
    ''' ==================================================================================
    FUNCTION Read_IsoCleftFile: Read the Output file generated by IsoCleft
    ==================================================================================  '''            
    def Read_IsoCleftFile(self):        
        
        filePath = self.OutputPath + '.isd'
        
        if (os.path.isfile(filePath)):
        
            # Read the Output file
            file = open(filePath, 'r')
            txtFile = file.readlines()
            file.close()
            
            # Clear the Used Atoms Lists
            del self.CleftA_UsedAtoms[:]
            del self.CleftB_UsedAtoms[:]
            
            # Contain the Center Coordinates of the Clefts
            del self.Center_CleftA_1[:]
            del self.Center_CleftB_1[:]
            
            del self.Center_CleftA_2[:]
            del self.Center_CleftB_2[:]
            
            # Contain the Rotation Matrix Values
            del self.RotationMatrix_1[:]
            del self.RotationMatrix_2[:]
            
            nbLines = len(txtFile)
            
            for no in range(0, nbLines):
                line = txtFile[no]
                endChar = len(line) - 1
                self.DisplayMessage('  ' + line[0:endChar], 0)
                
                ValueLn = line[8:14].strip()
                
                if ValueLn == 'NODE':
                    self.CleftA_UsedAtoms.append(line[23:27].strip())
                    self.CleftB_UsedAtoms.append(line[29:33].strip())
                
                elif (ValueLn == 'Centre'):

                    coordX = float(line[18:26].strip())
                    coordY = float(line[27:34].strip())
                    coordZ = float(line[35:42].strip())                   
                    
                    if line[15:16].strip() == 'A':
                        if len(self.Center_CleftA_1) > 0:
                            self.Center_CleftA_2 = [coordX, coordY, coordZ]
                        else:
                            self.Center_CleftA_1 = [coordX, coordY, coordZ]   
                
                    elif line[15:16].strip() == 'B':                   
                        if len(self.Center_CleftB_1) > 0:
                            self.Center_CleftB_2 = [coordX, coordY, coordZ]
                        else:                   
                            self.Center_CleftB_1 = [coordX, coordY, coordZ]
                    
                elif (line[8:25].strip() == 'Rotation Matrix:'):
                    line = txtFile[no+1]
                    row1 = [float(line[8:16].strip()), float(line[17:24].strip()), float(line[25:32].strip())]
                    line = txtFile[no+2]                        
                    row2 = [float(line[8:16].strip()), float(line[17:24].strip()), float(line[25:32].strip())]
                    line = txtFile[no+3]
                    row3 = [float(line[8:16].strip()), float(line[17:24].strip()), float(line[25:32].strip())]
                    
                    if len(self.RotationMatrix_1) > 0:
                        self.RotationMatrix_2 = [row1, row2, row3]
                        self.CleftMove_A(self.Center_CleftA_2, self.Center_CleftB_2, self.RotationMatrix_2)
                        self.CleftMove_B(self.Center_CleftB_2, self.CleftPath_B)
                        
                        # Clear the Used Atoms Lists
                        del self.CleftA_UsedAtoms[:]
                        del self.CleftB_UsedAtoms[:]
                    else:
                        self.RotationMatrix_1 = [row1, row2, row3]
                        self.CleftMove_A(self.Center_CleftA_1, self.Center_CleftB_1, self.RotationMatrix_1)
                        self.CleftMove_B(self.Center_CleftB_1, self.CleftPath_B)
 
            self.Fill_UnusedList()
            
            # PDB Cleft A superimposed
            self.CreatePDB(self.CleftPath_A, self.CleftPath_Asi, self.CleftA_Coord)
            # PDB Cleft B superimposed
            self.CreatePDB(self.CleftPath_B, self.CleftPath_Bsi, self.CleftB_Coord)
            
            self.Set_CleftDisplay(self.CleftName_Asi, self.CleftPath_Asi, True)
            self.Set_CleftDisplay(self.CleftName_Bsi, self.CleftPath_Bsi, True)
            
            nbAtUsed = len(self.CleftA_UsedAtoms)

            for at in range(0,nbAtUsed):
                At_ClfA = self.CleftName_Asi + '  & ID ' + self.CleftA_UsedAtoms[at]
                At_ClfB = self.CleftName_Bsi + ' & ID ' + self.CleftB_UsedAtoms[at]                
                cmd.distance('Superimposed', At_ClfA , At_ClfB)
                
            cmd.hide('label', 'Superimposed')
            cmd.color('white', 'Superimposed')
            cmd.disable('Superimposed')
            
            self.RadioBtn_OP1.config(state='normal')   
            self.RadioBtn_OP2.config(state='normal')
            self.RadioBtn_OP3.config(state='normal')
                
        else:
            self.DisplayMessage('  The file <' + fileName + '> did''nt exist...',1) 
        

    ''' ==================================================================================
    FUNCTION CleftMove_A: Move the Cleft A
    ==================================================================================  '''     
    def CleftMove_A(self, CenterA, CenterB, RotMatrix):
        
        for key, value in self.CleftA_Coord.items():            
            
            CoordList = [0.0, 0.0, 0.0]
            for i in range(0,3):
                for j in range(0,3):
                    CoordList[i] += ((self.CleftA_Coord[key][j] - CenterA[j]) * RotMatrix[j][i])
                
                CoordList[i] = (CoordList[i] + CenterB[i])
            
            self.CleftA_Coord[key][0] = CoordList[0]
            self.CleftA_Coord[key][1] = CoordList[1]
            self.CleftA_Coord[key][2] = CoordList[2]            
            
    
    ''' ==================================================================================
    FUNCTION CleftMove_B: Move the Cleft B
    ==================================================================================  '''     
    def CleftMove_B(self, CenterB, cleft_Path):
        
        # Get the Mass Center
        MC = list()
        MC = self.Get_CenterOfMass(cleft_Path)
        
        # Translate
        for key, value in self.CleftB_Coord.items():            
            
            for i in range(0,3):
                self.CleftB_Coord[key][i] = (self.CleftB_Coord[key][i] - MC[i] + CenterB[i])
        
                
    ''' ==================================================================================
    FUNCTION CreatePDB: Create a PDB file based on the NEW coordinates
    ==================================================================================  '''                 
    def CreatePDB(self, inputPath, outputPath, dictClf):
        
        # Read the Input file
        readfile = open(inputPath, 'r')
        txtFile = readfile.readlines()
        readfile.close()
        
        # Write in the output File
        file = open(outputPath, 'w')
        
        for pdbLine in txtFile:
            type = pdbLine[0:6].strip()
            if (type == 'HETATM') or (type == 'ATOM'):

                NoAtom = int(pdbLine[7:11])
                atomX = self.FormatFloat(dictClf[NoAtom][0])
                atomY = self.FormatFloat(dictClf[NoAtom][1])
                atomZ = self.FormatFloat(dictClf[NoAtom][2])
                
                tmpLine = pdbLine[0:29]
                tmpLine += str(atomX).rjust(9, ' ')     # The atom X coordinate
                tmpLine += str(atomY).rjust(8, ' ')     # The atom Y coordinate
                tmpLine += str(atomZ).rjust(8, ' ')     # The atom Z coordinate
                tmpLine += pdbLine[54:]
                file.write(tmpLine)
            else:
                file.write(pdbLine)
            
        file.close()        
        
   
    ''' ==================================================================================
    FUNCTION IsDisplayed: Verify if the pdb is displayed 
    ==================================================================================  '''   
    def IsDisplayed(self, pdbName):
        
        # Get the Selection Names in Pymol        
        sess = cmd.get_names('all')
        NotFound = True
        
        for i in sess:
            if (i == pdbName):
                cmd.enable(i)
                NotFound = False
                break
        
        return NotFound
    

    ''' ==================================================================================
    FUNCTION Set_CleftCoordinates: Change the coordinates of the Cleft to Display the
                                   Cleft A and B next together. 
    ==================================================================================  '''    
    def Set_CleftCoordinates(self, Cleft_Name, Cleft_Path, CleftAorB):
        
        # MassCenterRadius
        MCR = list()
        MCR = self.Get_CenterOfMass(Cleft_Path)
        
        # Read the PDB file    
        file = open(Cleft_Path, 'r')
        txtFile = file.readlines()
        file.close()
         
        write_file = open(Cleft_Path, 'w') 
              
        for line in txtFile:
            type = line[0:6].strip()
            if (type == 'ATOM') or (type == 'HETATM'):

                coordX = float(line[30:38].strip()) - MCR[0]     # The atom X coordinate
                coordY = float(line[39:46].strip()) - MCR[1]     # The atom Y coordinate
                coordZ = float(line[47:54].strip()) - MCR[2]     # The atom Z coordinate
                
                if CleftAorB == 'A':
                    self.Radius_A = MCR[3]
                    coordX = coordX - MCR[3]
                else:
                    self.Radius_B = MCR[3]
                    coordX = coordX + MCR[3]
                
                tmpLine = line[0:29]
                tmpLine += str(self.FormatFloat(coordX)).rjust(9, ' ')     # The atom X coordinate
                tmpLine += str(self.FormatFloat(coordY)).rjust(8, ' ')     # The atom Y coordinate
                tmpLine += str(self.FormatFloat(coordZ)).rjust(8, ' ')     # The atom Z coordinate
                tmpLine += line[54:]
                
                write_file.write(tmpLine)
            else:
                write_file.write(line)
    

    #=======================================================================
    ''' Function centerOfMass: Get the coordinates of the center of mass 
                               of the protein, and the Radius            '''
    #=======================================================================       
    def Get_CenterOfMass(self, CleftPath):

        xMin = 0.0
        yMin = 0.0
        zMin = 0.0
        
        xMax = 0.0       
        yMax = 0.0       
        zMax = 0.0
        
        xMid = 0.0
        yMid = 0.0
        zMid = 0.0
        
        Range = 0.0
        count = 0
        
        # Read the PDB file    
        file = open(CleftPath, 'r')
        txtFile = file.readlines()
        file.close()         
                  
        firstIN = True
            
        for line in txtFile:
            type = line[0:6].strip()
            if (type == 'ATOM') or (type == 'HETATM'):
                count += 1                
                coordX = float(line[30:38].strip())     # The atom X coordinate
                coordY = float(line[39:46].strip())     # The atom Y coordinate
                coordZ = float(line[47:54].strip())     # The atom Z coordinate
                
                xMid += coordX
                yMid += coordY
                zMid += coordZ
            
                if firstIN:
                    xMin = xMax = coordX
                    yMin = yMax = coordY
                    zMin = zMax = coordZ
                    firstIN = False
                else:
                    if coordX < xMin:
                        xMin = coordX
                    elif coordX > xMax:
                        xMax = coordX
                    
                    if coordY < yMin:
                        yMin = coordY
                    elif coordY > yMax:
                        yMax = coordY
                        
                    if coordZ < zMin:
                        zMin = coordZ
                    elif coordZ > zMax:
                        zMax = coordZ

        xDist = abs(xMin - xMax)
        yDist = abs(yMin - yMax)
        zDist = abs(zMin - zMax)    
           
        if (xDist > yDist):
            if (xDist > zDist):
                Range = xDist
            else:
                Range = zDist
        elif (yDist > zDist):
            Range = yDist
        else:
            Range = zDist    
       
        return (xMid / float(count), yMid / float(count), zMid / float(count), Range)
 
    
    ''' ==================================================================================
    FUNCTION DisplayCleft: Display Cleft A 
    ==================================================================================  '''    
    def Set_CleftDisplay(self, CleftName, CleftPath, boolLines):
        
        if self.IsDisplayed(CleftName):
            #Add to the DropDownList            
            cmd.mstop()           # Stop the movie
            cmd.frame('1')        # Display the Frame no.1
            
            if boolLines:
                cmd.load(CleftPath, CleftName)    # Load the clf file in Pymol
            else:
                cmd.set('connect_mode', 1)
                cmd.load(CleftPath, CleftName)    # Load the clf file in Pymol
                cmd.set('connect_mode', 0)
            
        cmd.hide('surface', CleftName)
        cmd.hide('nonbonded', CleftName)
        cmd.show('sphere', CleftName)                    
        cmd.alter(CleftName,'vdw=0.40')
        
        if self.Radius_A != 0.0 and self.Radius_B != 0.0:
            self.MidOrigin = int(-1 * (self.Radius_A - self.Radius_B))
            cmd.origin(position=[self.MidOrigin,0,0])        
            cmd.center('origin')
            
            if  self.Radius_A > self.Radius_B:
                self.BiggestClf = self.CleftName_Asi
            else:
                self.BiggestClf = self.CleftName_Bsi   
                      
            self.ZoomVal = int((self.Radius_A + self.Radius_B)/1.2 )           
            cmd.zoom('origin', self.ZoomVal)

        util.cbag('all')
        
        self.DisplayONLY_ClfA_B()
        

    ''' ==================================================================================
    FUNCTION DisplayONLY_ClfA_B: Display ONLY the Cleft A and Cleft B. 
    ==================================================================================  '''   
    def DisplayONLY_ClfA_B(self):
        
        sess = cmd.get_names('all')
        for i in sess:
            if (i != self.CleftName_A) and (i != self.CleftName_B):
                cmd.disable(i)
    
    
    ''' ==================================================================================
    FUNCTION Step1_Clicked: Display the Starting Clefts. 
    ==================================================================================  '''    
    def Step1_Clicked(self):
        
        if self.STEPS.get() != self.LAST_STEP:
            if self.LAST_STEP == 3:
                cmd.enable(self.CleftName_A)
                cmd.enable(self.CleftName_B)
                cmd.disable(self.CleftName_Asi)
                cmd.disable(self.CleftName_Bsi)
                cmd.disable('Superimposed')          
 
                cmd.origin(position=[self.MidOrigin,0,0])
                cmd.center('origin')
                cmd.zoom('origin', self.ZoomVal)
                
            self.LAST_STEP = 1
            self.DisplayMessage('  Display STEP 1: The Starting Clefts.',0)
            util.cbag('all')
        
    
    ''' ==================================================================================
    FUNCTION Step2_Clicked: Display the Filtered Atoms. 
    ==================================================================================  '''    
    def Step2_Clicked(self):
        
        if self.STEPS.get() != self.LAST_STEP:
            if self.LAST_STEP == 3:
                cmd.enable(self.CleftName_A)
                cmd.enable(self.CleftName_B)
                cmd.disable(self.CleftName_Asi)
                cmd.disable(self.CleftName_Bsi)
                cmd.disable('Superimposed')
                
                cmd.origin(position=[self.MidOrigin,0,0])
                cmd.center('origin')
                cmd.zoom('origin', self.ZoomVal)           
           
            if  self.LAST_STEP == 1: 
                # Put in Grey the Atoms that are'nt used for the SuperImposed 
                if (len(self.Unused_AtomsCleftA) > 0):
                    self.Set_AtomsInGrey(self.CleftName_A, self.Unused_AtomsCleftA)
                    
                if (len(self.Unused_AtomsCleftB) > 0):
                    self.Set_AtomsInGrey(self.CleftName_B, self.Unused_AtomsCleftB)            

            self.LAST_STEP = 2            
            self.DisplayMessage('  Display STEP 2: The Filtered Atoms (Unused).',0)
        

    ''' ==================================================================================
    FUNCTION Set_AtomsInGrey: Display in Grey the Atoms that are'nt used for the Superimposed. 
    ==================================================================================  '''       
    def Set_AtomsInGrey(self, CleftName, Unused_AtomsCleft):
        
        for elem in Unused_AtomsCleft:
            cmd.color('grey20','(id %s)'%elem)

                    
    ''' ==================================================================================
    FUNCTION Fill_UnusedList: Fill the Unused Atoms List of each Cleft. 
    ==================================================================================  '''                    
    def Fill_UnusedList(self):
        
        # CLEFT A
        file = open(self.CleftPath_A, 'r')
        clfA_File = file.readlines()
        file.close()
        
        del self.Unused_AtomsCleftA[:]      # Clear the list
        
        for lineA in clfA_File:
            type = lineA[0:6].strip()
            if (type == 'ATOM'):             
                index = lineA[7:11].strip()   # The atom number
                NotFound = True
                for elem in self.CleftA_UsedAtoms:
                    if index == elem:
                        NotFound = False
                        break
                if NotFound:
                    self.Unused_AtomsCleftA.append(index)                
                
        # CLEFT B 
        file = open(self.CleftPath_B, 'r')
        clfB_File = file.readlines()
        file.close()
        
        del self.Unused_AtomsCleftB[:]      # Clear the list
        
        for lineB in clfB_File:
            type = lineB[0:6].strip()
            if (type == 'ATOM'):             
                index = lineB[7:11].strip()   # The atom number
                NotFound = True
                for elem in self.CleftB_UsedAtoms:
                    if index == elem:
                        NotFound = False
                        break
                if NotFound:
                    self.Unused_AtomsCleftB.append(index)         
        
        
    ''' ==================================================================================
    FUNCTION Step3_Clicked: Display the Superimposed Atoms. 
    ==================================================================================  '''    
    def Step3_Clicked(self):
        
        if self.STEPS.get() != self.LAST_STEP:
            
            cmd.disable(self.CleftName_A)
            cmd.disable(self.CleftName_B)
            cmd.enable(self.CleftName_Asi)
            cmd.enable(self.CleftName_Bsi)
            cmd.enable('Superimposed')            
        
            if  self.LAST_STEP == 1:  
                # Put in Grey the Atoms that are'nt used for the SuperImposed 
                if (len(self.Unused_AtomsCleftA) > 0):
                    self.Set_AtomsInGrey(self.CleftName_A, self.Unused_AtomsCleftA)
                    
                if (len(self.Unused_AtomsCleftB) > 0):
                    self.Set_AtomsInGrey(self.CleftName_B, self.Unused_AtomsCleftB)
                    
            cmd.center(self.CleftName_Asi)
            cmd.zoom(self.BiggestClf)
                   
            self.LAST_STEP = 3
            self.DisplayMessage('  Display STEP 3: The Superimposed Atoms.',0)
        

    ''' ==================================================================================
    FUNCTION GetCleftIsRunning: Update or Create the Running File to BLOCK multiple GUI 
    ==================================================================================  '''       
    def IsoCleftIsRunning(self):
        
        #Create the .run file
        RunPath = self.RunDir_Path + '.irun'
        RunFile = open(RunPath, 'w')
        RunFile.write(str(os.getpid()))
        RunFile.close()
  

    ''' ==================================================================================
    FUNCTION FormatFloat: Format a string float number to 3 numbers after the . 
    ==================================================================================  '''     
    def FormatFloat(self, strFloat):
        
        noFloat = str(strFloat)

        if noFloat.count('.') > 0:                        
            dotIndex = noFloat.find('.') + 4
            nbChar = len(noFloat)
            if nbChar <  dotIndex:
                diffZero = dotIndex - nbChar
                noFloat = noFloat[0:nbChar]
                for i in range(0, diffZero):
                    noFloat += '0'
            else:
                # round the value if necessary
                if ((nbChar - dotIndex) > 0):
                    if (int(noFloat[dotIndex]) > 4):                    
                        value = float(strFloat) + 0.001
                        noFloat = str(value)
                
                noFloat = noFloat[0:dotIndex]
        else:
            noFloat = noFloat + '.000'

        return noFloat
        
    
    ''' ==================================================================================
    FUNCTION Btn_Clear_Clicked: Clear the selections in Pymol 
    ==================================================================================  '''
    def Btn_Clear_Clicked(self):
     
        self.DisplayMessage('',0)
        self.DisplayMessage('  Clear ALL the selections...',0) 
        
        #self.Radius_A = 0.0
        #self.Radius_B = 0.0
        
        
    ''' ==================================================================================
    FUNCTION Btn_Quit_Clicked: Exit the application 
    ==================================================================================  '''
    def Btn_Quit_Clicked(self):
        
        self.DisplayMessage('  *** Closing IsoCleft... ***',0)
        
        #Delete the .run file
        RunPath = self.RunDir_Path + '.irun'
        if (os.path.isfile(RunPath)):
            try:
                os.remove(RunPath)
            except OSError:
                time.sleep(0.1)
                os.remove(RunPath)
                    
        self.top.destroy()
        


    ''' ==================================================================================
    FUNCTION Btn_Start_Clicked: Run IsoCleft and display the result in Pymol 
    ==================================================================================  '''
    def Btn_Start_Clicked(self):
        
        self.timeStart = time()
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  Starting IsoCleft...', 0)
        self.IsoCleftRun = True
        
        self.RadioBtn_OP1.config(state='disable')   
        self.RadioBtn_OP2.config(state='disable')
        self.RadioBtn_OP3.config(state='disable')
        
        # Display ONLY the CleftA and CleftB
        self.DisplayONLY_ClfA_B()          
        
        Argument = self.IsoCleft_Path + 'WRK/IsoCleft -i ' + self.CleftPath_A + ' ' + self.CleftPath_B
        
        # Extract all the Advanced Options
        
        if self.Entry_R.get() != self.DEFAULT_R:
            Argument += ' -r ' + self.Entry_R.get()
            
        if self.Entry_C.get() != self.DEFAULT_C:
            Argument += ' -c ' + self.FormatFloat(self.Entry_C.get())
            
        if self.Entry_D.get() != self.DEFAULT_D:
            Argument += ' -d ' + self.FormatFloat(self.Entry_D.get())

        if self.Entry_N.get() != self.DEFAULT_N:
            Argument += ' -n ' + self.FormatFloat(self.Entry_N.get())
        
        if self.Entry_S.get() != self.DEFAULT_S:
            Argument += ' -s ' + self.Entry_S.get()                
            
        self.OutputPath = self.TempPath + self.Entry_O.get()        
        
        Argument += ' -o ' + self.OutputPath
        
        self.CleftA_Coord.clear()
        self.CleftB_Coord.clear()
        
        self.CleftA_Coord = self.Get_CleftCoordinates(self.CleftPath_A)
        self.CleftB_Coord = self.Get_CleftCoordinates(self.CleftPath_B)  
                    
        StartSim = SimThread(self, Argument)
        self.Btn_Start.config(state='disable')
        

    ''' ==================================================================================
    FUNCTION Get_CleftCoordinates: Get the coordinates of each atoms in each cleft pdb
    ==================================================================================  '''    
    def Get_CleftCoordinates(self, CleftPath):
        
        CoordAtom = {}

        file = open(CleftPath)
        pdbFile = file.readlines()
        file.close()
        
        # Creation of a dictionary containing the 3 neighbors of each atoms
        for line in pdbFile:
            if (line.startswith('ATOM'))  or (line.startswith('HETATM')):
                noLine = int(line[7:11])
                CoordAtom[noLine] = [float(line[30:38]),float(line[39:46]),float(line[47:54])]
               
        return CoordAtom
                 
    
    ''' ==================================================================================
    FUNCTION CleanTempDirectory: Clean the TEMP directory
    ==================================================================================  '''    
    def CleanTempDirectory(self, path):        
        
        files=os.listdir(path)

        for x in files:
            fullpath=os.path.join(path, x)
            if (os.path.isfile(fullpath)):
                if (fullpath != self.CleftPath_A) and (fullpath != self.CleftPath_B):
                    try:
                        os.remove(fullpath)
                    except OSError:
                        time.sleep(0.1)
                        os.remove(fullpath)
                        
                        
    ''' ==================================================================================
                            VALIDATION FUNCTION REGION 
    ================================================================================== '''

    ''' ==================================================================================
    FUNCTION Validate_Entry_R: Validate the Entry R
    ================================================================================== '''   
    def Validate_Entry_R(self, event):        
        
        term = self.Entry_R.get().upper()
        
        self.ValidEntry_R = False
        
        if len(term) > 0:        
            if term.isdigit():
                value = int(term)
                if (value >= 1) and (value <= 20):
                    self.EntryResult_R.config(bg=self.Color_White)                    
                    self.ValidEntry_R = True
                else:
                    self.DisplayMessage('  ERROR: The JTT rank Threshold must be between [1 - 20]!', 1)
                    self.EntryResult_R.config(bg=self.Color_Red)
            else:
                self.DisplayMessage('  ERROR: The JTT rank Threshold must be an INTEGER!', 1)
                self.EntryResult_R.config(bg=self.Color_Red)
        else:
            self.EntryResult_R.config(bg=self.Color_White) 
            

    ''' ==================================================================================
    FUNCTION Validate_Entry_C: Validate the Entry C
    ================================================================================== '''   
    def Validate_Entry_C(self, event):        
        
        term = self.Entry_C.get().upper()
        
        self.ValidEntry_C = False
        
        if len(term) > 0:
            try:
                valFloat = float(term)

                if (valFloat >= 0.5) and (valFloat <= 5.0): 
                    self.EntryResult_C.config(bg=self.Color_White)
                    self.ValidEntry_C = True
                else:
                    self.DisplayMessage('  ERROR: The C-alpha Node Distance must be between 0.5 and 5.0!!', 1)
                    self.EntryResult_C.config(bg=self.Color_Red)
            except:
                self.DisplayMessage('  ERROR: The C-alpha Node Distance must be an INTEGER or a DECIMAL!', 1)
                self.EntryResult_C.config(bg=self.Color_Red)
        else:
            self.EntryResult_C.config(bg=self.Color_White)   
           

    ''' ==================================================================================
    FUNCTION Validate_Entry_D: Validate the Entry D
    ================================================================================== '''   
    def Validate_Entry_D(self, event):        
        
        term = self.Entry_D.get().upper()
        
        self.ValidEntry_D = False
        
        if len(term) > 0:
            try:
                valFloat = float(term)

                if (valFloat >= 0.5) and (valFloat <= 6.0):     
                    self.EntryResult_D.config(bg=self.Color_White)
                    self.ValidEntry_N = True
                else:
                    self.DisplayMessage('  ERROR: The Neighborhood Sphere Radius must be between 0.5 and 6.0!', 1)
                    self.EntryResult_D.config(bg=self.Color_Red)
            except:
                self.DisplayMessage('  ERROR: The Neighborhood Sphere Radius must be an INTEGER or a DECIMAL!', 1)
                self.EntryResult_D.config(bg=self.Color_Red)
        else:
            self.EntryResult_D.config(bg=self.Color_White)
            

    ''' ==================================================================================
    FUNCTION Validate_Entry_N: Validate the Entry N
    ================================================================================== '''   
    def Validate_Entry_N(self, event):        
        
        term = self.Entry_N.get().upper()
        
        self.ValidEntry_N = False
        
        if len(term) > 0:                       
            try:
                valFloat = float(term)

                if (valFloat >= 0.5) and (valFloat <= 5.0):     
                    self.EntryResult_N.config(bg=self.Color_White)
                    self.ValidEntry_N = True
                else:
                    self.DisplayMessage('  ERROR: The Node Distance Threshold must be between 0.5 and 5.0!!', 1)
                    self.EntryResult_N.config(bg=self.Color_Red)
            except:
                self.DisplayMessage('  ERROR: The Node Distance Threshold must be an INTEGER or a DECIMAL!', 1)
                self.EntryResult_N.config(bg=self.Color_Red)
        else:
            self.EntryResult_N.config(bg=self.Color_White)           


    ''' ==================================================================================
    FUNCTION Validate_Entry_S: Validate the Entry S
    ================================================================================== '''   
    def Validate_Entry_S(self, event):        
        
        term = self.Entry_S.get().upper()
        
        self.ValidEntry_S = False
        
        if len(term) > 0:          
            if term.isdigit():
                value = int(term)
                if (value >= 0) and (value <= 8):
                    self.EntryResult_S.config(bg=self.Color_White)
                    self.ValidEntry_S = True
                else:
                    self.DisplayMessage('  ERROR: The Level of BK must be between [0 - 8]!', 1)
                    self.EntryResult_S.config(bg=self.Color_Red)
            else:
                self.DisplayMessage('  ERROR: The Level of BK must be an INTEGER!', 1)
                self.EntryResult_S.config(bg=self.Color_Red)
        else:
            self.EntryResult_S.config(bg=self.Color_White)            


    ''' ==================================================================================
    FUNCTION Validate_Entry_O: Validate the Entry O - Output Filename
    ================================================================================== '''   
    def Validate_Entry_O(self, event):        
        
        term = self.Entry_O.get().upper()
        
        self.ValidEntry_O = False
        
        if len(term) > 0:
            self.ValidEntry_O = True
            
            
    ''' ==================================================================================
    FUNCTION ValidateFolders: Be sure the folders Exists 
    ==================================================================================  '''    
    def ValidateFolders(self):       
       
        if not(os.path.isdir(self.UserIsoCleft)):
            os.mkdir(self.UserIsoCleft)
            
        if not(os.path.isdir(self.UserGetCleft)):
            os.mkdir(self.UserGetCleft)
        
        if not(os.path.isdir(self.SavePath)):
            os.mkdir(self.SavePath)
            
        if not(os.path.isdir(self.TempPath)):
            os.mkdir(self.TempPath)                

               
    ''' ==================================================================================
    FUNCTION DisplayMessage: Display the message  
    ==================================================================================  '''    
    def DisplayMessage(self, msg, priority):
        
        NBCHARMAX = 45
        self.TextMessage.config(state='normal')
 
        lineNo = str(self.MsgLineCounter)
        insNo = lineNo + '.0'      

        if msg != '':
            # Verify if the message need to be split
            NbChar = len(msg)
            if NbChar > NBCHARMAX:                
              
                CharStart = 0
                CharEnd = 0
                words = list()
                words = msg.split()
                NbWords = len(words)
                line = ''
                
                for i in range(0, NbWords):
                   
                    # Building a line...
                    line += words[i] + ' '
                    if (len(line) > NBCHARMAX) or (i == NbWords-1):
                        Size = len(line)
       
                        self.TextMessage.insert(INSERT, '\n  ' + line)
                        CharStart += NBCHARMAX
                        
                        if priority == 1:
                            self.TextMessage.tag_add('warn', lineNo + '.0', lineNo + '.' + str(Size + 2))
                            self.TextMessage.tag_config('warn', foreground='red')
                        elif priority == 2:
                            self.TextMessage.tag_add('notice', lineNo + '.0', lineNo + '.' + str(Size + 2))
                            self.TextMessage.tag_config('notice', foreground='blue')   
                            
                        self.TextMessage.yview(INSERT)
                        self.MsgLineCounter += 1
                        lineNo = str(self.MsgLineCounter)
                        
                        line = ''
            else:
                #self.TextMessage.config(font='red')
                self.TextMessage.insert(INSERT, '\n' + msg)
                
                if priority == 1:
                    self.TextMessage.tag_add('warn', lineNo + '.0', lineNo + '.' + str(NbChar))
                    self.TextMessage.tag_config('warn', foreground='red')
                elif priority == 2:
                    self.TextMessage.tag_add('notice', lineNo + '.0', lineNo + '.' + str(NbChar))
                    self.TextMessage.tag_config('notice', foreground='blue')   
                    
                self.TextMessage.yview(INSERT)
                self.MsgLineCounter += 1
        else:
            # Skip a line...
            self.TextMessage.insert(INSERT, '\n')
            self.TextMessage.yview(INSERT)
            self.MsgLineCounter += 1        
        
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


