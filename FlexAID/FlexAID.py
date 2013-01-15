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
@title: FlexAID - Interface

@summary: This is the interface of FlexAID application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Aug. 25, 2010
'''

from Tkinter import *

import os
import time
import tkFont
import tkFileDialog

import Prefs
import General
import MessageBox
import IOFile
import Config1
import Config2
import Config3
import GAParam
import Simulate        

if __debug__:
	from pymol import cmd
	from pymol.cgo import *
	from pymol.vfont import plain

	import General_cmd

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayFlexAID:
    
    ''' ==================================================================================
    FUNCTION __init__ : Initialization of the variables of the interface
    ================================================================================== '''
    def __init__(self, top, ActiveWizard, Project_Dir, Install_Dir, AlreadyRunning_Dir, OSid, PyMOL):
        
        #print("New instance of FlexAID")
        self.PyMOL = PyMOL

        self.WINDOWWIDTH = 700
        self.WINDOWHEIGHT = 600

        #================================================================================== 
        ''' USER SAVE PATH DIRECTORY  '''
        #================================================================================== 
        self.AlreadyRunning_Dir = AlreadyRunning_Dir

        self.Install_Dir = Install_Dir
        self.FlexAIDInstall_Dir = os.path.join(self.Install_Dir,'FlexAID')

        self.OSid = OSid
        if self.OSid == 'WIN':
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID.exe')
        else:
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID')

        self.Project_Dir = Project_Dir
        self.FlexAIDProject_Dir = os.path.join(self.Project_Dir,'FlexAID')
        self.GetCleftProject_Dir = os.path.join(self.Project_Dir,'GetCleft')
        self.GetCleftSaveProject_Dir = os.path.join(self.GetCleftProject_Dir,'Save')

        self.CleftProject_Dir = os.path.join(self.Project_Dir,'Cleft')
        self.ProteinProject_Dir = os.path.join(self.Project_Dir,'Target')
        self.BindingSiteProject_Dir = os.path.join(self.Project_Dir,'Binding_Site')

        self.FlexAIDLigandProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Ligand')
        self.FlexAIDSimulationProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Simulation')
        self.FlexAIDTargetFlexProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Target_Flexibility')

        self.LOGFILE = os.path.join(self.FlexAIDProject_Dir,'logfile.txt')
        #self.LOGCOMMANDS = os.path.join(self.UserFlexAID,'logcmd.pml')
        
        self.ValidateFolders()      

        #================================================================================== 
        ''' ROOT PATH TO THE FLEXAID DIRECTORY  '''
        #================================================================================== 
        self.OSid = OSid

        self.BindingSiteDisplay = 'BINDING_SITE_AREA'
        self.SphereDisplay = 'SPHERE_AREA'
                
        self.FlexAIDIsRunning()
        
        #self.MsgLineCounter = 2
        
        self.Color_Green = '#CCFFCC'
        self.Color_Grey = '#EDEDED'
        self.Color_Blue = '#6699FF'
        self.Color_Red = '#FF9999'
        self.Color_White = '#FFFFFF'
        self.Color_Black = 'black'                

        #print("FlexAID: initializing window")
        # Initialize the window
        self.top = top        
        self.top.title('FlexAID')

        #print("FlexAID: center window")
        General.CenterWindow(self.top,self.WINDOWWIDTH,self.WINDOWHEIGHT)


        #self.top.geometry()   # Interface DIMENSIONS
        #self.top.maxsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        #self.top.minsize(self.WINDOWWIDTH,self.WINDOWHEIGHT)
        self.top.protocol('WM_DELETE_WINDOW', self.Btn_Quit_Clicked)       

        #================================================================================== 
        #                 SET the default fonts of the interface
        #==================================================================================
        FontType = Prefs.GetFontType()
        FontSize = Prefs.GetFontSize()
        
        self.font_Title = tkFont.Font(family=FontType,size=FontSize, weight=tkFont.BOLD)        
        self.font_Text = tkFont.Font(family=FontType,size=FontSize)
        self.font_Text_I = tkFont.Font(family=FontType,size=FontSize, slant=tkFont.ITALIC)
        self.font_Text_U = tkFont.Font(family=FontType,size=FontSize, underline=True)       
        self.Img_Msg = list()        
        
        #================================================================================== 
        #                       FRAMES Settings and startup
        #==================================================================================
                
        self.frame = Frame(top)
        self.frame.pack(expand=True)
        self.curCursor = self.frame['cursor']

        self.Frame_Main()
	#self.DisplayMessage = self.MessageBox.add_message

        #print "Main frame generated"

        self.ActiveFrame = None
        self.ActiveWizard = ActiveWizard
        self.WizardError = False
        self.WizardResult = 0

        self.ProcessRunning = False
        self.Run = None

        # Build class objects of each tab
        self.IOFile = IOFile.IOFile(self, self.PyMOL)        
        self.Config1 = Config1.Config1(self, self.PyMOL)
        self.Config2 = Config2.Config2(self, self.PyMOL)
        self.Config3 = Config3.Config3(self, self.PyMOL)
        self.GAParam = GAParam.GAParam(self, self.PyMOL)
        self.Simulate = Simulate.Simulate(self, self.PyMOL)

        self.listBtnTabs = [self.Btn_IOFiles,self.Btn_Config1,self.Btn_Config2,self.Btn_Config3,self.Btn_GAParam,self.Btn_Simulate]

        #print "Created instances of all tab classes"

        # default Tab (IOFile)
        self.Btn_IOFiles_Clicked()

    #=====================================================================================
    '''                     --- FRAME DISPLAY SETTINGS ---                             '''
    #=====================================================================================
    
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''  
    def Frame_Main(self):
   
        #==================================================================================
        '''                  --- TOP MENU OF THE INTERFACE ---                          '''
        #==================================================================================

        fTop = Frame(self.frame, relief=RIDGE, border=4, height=50)
        fTop.pack(fill=BOTH, expand=True)#, padx=10, pady=10, ipady=10, ipadx=10, side=TOP)
        fTop.pack_propagate(0)
        
        self.Btn_IOFiles = Button(fTop, text='Input Files', bg=self.Color_White, command=self.Btn_IOFiles_Clicked, font=self.font_Text)
        self.Btn_IOFiles.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_IOFiles.config(state='normal')

        self.Btn_Config1 = Button(fTop, text='Target Cfg', bg=self.Color_Grey, command=self.Btn_Config1_Clicked, font=self.font_Text)
        self.Btn_Config1.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config1.config(state='disabled')

        self.Btn_Config2 = Button(fTop, text='Ligand Cfg', bg=self.Color_Grey, command=self.Btn_Config2_Clicked, font=self.font_Text)
        self.Btn_Config2.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config2.config(state='disabled')

        self.Btn_Config3 = Button(fTop, text='Scoring Cfg', bg=self.Color_Grey, command=self.Btn_Config3_Clicked, font=self.font_Text)
        self.Btn_Config3.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config3.config(state='disabled')
        
        self.Btn_GAParam = Button(fTop, text='GA Param', bg=self.Color_Grey, command=self.Btn_GAParam_Clicked, font=self.font_Text)
        self.Btn_GAParam.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_GAParam.config(state='disabled')    

        self.Btn_Simulate = Button(fTop, text='Simulate', bg=self.Color_Grey, command=self.Btn_Simulate_Clicked, font=self.font_Text)
        self.Btn_Simulate.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Simulate.config(state='disabled') 
        
        #==================================================================================
        '''                MIDDLE DISPLAY SECTION OF THE INTERFACE                      '''
        #==================================================================================
           
        self.fMiddle = Frame(self.frame)#, height=355)
                
        self.fMiddle.pack(fill=X, expand=True, padx=10, side=TOP)
        #self.fMiddle.config(takefocus=1)
        self.fMiddle.bind('<FocusIn>', lambda e: None)
        #self.fMiddle.pack_propagate(0)

        #==================================================================================
        '''                 BOTTOM DISPLAY SECTION OF THE INTERFACE                     '''
        #==================================================================================
 
        fBottom = Frame(self.frame, border=1, relief=SUNKEN)
        fBottom.pack(fill=BOTH, expand=True, padx=10, pady=10, ipadx=10, ipady=10, side=TOP)

        fBottomRight = Frame(fBottom)
        fBottomRight.pack(fill=Y, side=RIGHT)

        Btn_Default = Button(fBottomRight, text='Default', command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Default.pack(side=TOP, fill=X)

        Btn_SaveDefault = Button(fBottomRight, text='Save as default', command=self.Btn_SaveDefault_Clicked, font=self.font_Text)
        Btn_SaveDefault.pack(side=TOP, fill=X)

        Btn_Restore = Button(fBottomRight, text='Restore', command=self.Btn_Restore_Clicked, font=self.font_Text)
        Btn_Restore.pack(side=TOP, fill=X)

        Btn_Quit = Button(fBottomRight, text='Close', command=self.Btn_Quit_Clicked, font=self.font_Text)

        Btn_Quit.pack(side=BOTTOM, fill=X)

        fBottomLeft = Frame(fBottom)
        fBottomLeft.pack(side=LEFT, fill=Y, ipadx=20, ipady=20)

        #self.MessageBox = MessageBox.scrollableContainer(fBottomLeft, bd=2, bg="black")
        #self.MessageBox.pack(fill=BOTH, expand=True)

        scrollBar = Scrollbar(fBottomLeft)
        scrollBar.pack(side=RIGHT, fill=Y)

        self.TextMessage = Text(fBottomLeft, border=1, background=self.Color_Grey, font=self.font_Text)
        self.TextMessage.pack(side=RIGHT, fill=BOTH, expand=True)

        scrollBar.config(command=self.TextMessage.yview)
        self.TextMessage.config(state='disabled', yscrollcommand=scrollBar.set)                                       
    
    ''' ==================================================================================
    FUNCTION FlexAIDIsRunning: Update or Create the Running File to BLOCK multiple GUI 
    ==================================================================================  '''       
    def FlexAIDIsRunning(self):
        
        #Create the .run fileBtn_DelResidu_Clicked
        RunPath = os.path.join(self.AlreadyRunning_Dir,'.frun')
        RunFile = open(RunPath, 'w')
        RunFile.write(str(os.getpid()))
        RunFile.close()
        
                
    ''' ==================================================================================
    FUNCTION Btn_*_Clicked: Display the Tab options menu
    ================================================================================== '''    
    def Btn_IOFiles_Clicked(self):
        
        self.SetActiveFrame(self.IOFile)

    def Btn_Config1_Clicked(self):
        
        self.SetActiveFrame(self.Config1)

    def Btn_Config2_Clicked(self):
        
        self.SetActiveFrame(self.Config2)

    def Btn_Config3_Clicked(self):
        
        self.SetActiveFrame(self.Config3)

    def Btn_GAParam_Clicked(self):
        
        self.SetActiveFrame(self.GAParam)

    def Btn_Simulate_Clicked(self):
        
        self.SetActiveFrame(self.Simulate)

    ''' ============================================================================= '''

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
    def SetActiveFrame(self,Frame):

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

                rv = self.Validate_Entries(self.ActiveFrame.Validator)
                if rv > 0:
                    if rv == 1:
                        self.DisplayMessage("Cannot switch tab: Not all fields are validated", 2)
                    elif rv == 2:
                        self.ActiveFrame.Validator_Fail()
                    return

                if not self.ActiveFrame.Kill_Frame():
                    self.DisplayMessage("Cannot switch tab: Not all fields are validated", 2)
                    return

                self.fMiddle.update_idletasks()

                #self.ActiveFrame.Del_Trace()
                self.ActiveFrame.Tab.config(bg=self.Color_White)
                #print "Killed Frame " + self.ActiveFrame.FrameName


            self.ActiveFrame = Frame
            #print "New active frame " + self.ActiveFrame.FrameName
            self.ActiveFrame.Show()
            self.ActiveFrame.Tab.config(bg=self.Color_Blue)

	    self.fMiddle.update_idletasks()


    ''' ==================================================================================
    FUNCTION Go_Step1: Enables/Disables buttons for step 1
    ================================================================================== '''    
    def Go_Step1(self):

        #print "Setting Tab buttons to Step 1"        
        for Btn in self.listBtnTabs:
            Btn.config(state='disabled',bg=self.Color_Grey)

        self.Btn_IOFiles.config(state='normal')

    ''' ==================================================================================
    FUNCTION Go_Step2: Enables/Disables buttons for step 2
    ================================================================================== '''    
    def Go_Step2(self):
        
        #print "Setting Tab buttons to Step 2"        
        for Btn in self.listBtnTabs:
            Btn.config(state='normal',bg=self.Color_White)
        
    ''' ==================================================================================
    FUNCTION Reset_Step2: Reset ALL the parameters
    ==================================================================================  '''    
    def Reset_Step2(self):
        
        self.Config1.Init_Vars()
        self.Config2.Init_Vars()
        self.Config3.Init_Vars()
        self.GAParam.Init_Vars()
        self.Simulate.Init_Vars()        


    #=====================================================================================
    '''                            --- BUTTONS EVENT ---                               '''
    #=====================================================================================     

    ''' ==================================================================================
    FUNCTION Btn_LoadCfg_Clicked: Load the Configuration from a Text File 
    ==================================================================================  '''
    def Btn_LoadCfg_Clicked(self):
        # Validate the selection of a PDB name and a Ligand
        filename = tkFileDialog.askopenfilename(initialdir=os.path.join(self.path,'Config'), title='Load the Configuration File', filetypes=[('Configuration File','*.cfg')])
            
        if len(filename ) > 0:
            text_file = open(os.path.join(self.path,'Config/write_it.txt'), 'r')
            #print text_file.read()
            text_file.close()
            

    ''' ==================================================================================
    FUNCTION Btn_Quit_Clicked: Exit the application 
    ==================================================================================  '''
    def Btn_Quit_Clicked(self):

        # Cannot quit while process is running
        if self.ProcessRunning is True:
            self.DisplayMessage('  GUI cannot be closed while a process is currently running', 0)
            return
        
        # Close any Wizard interface in Pymol if started
        if not self.ActiveWizard is None:
            if self.PyMOL:
                cmd.get_wizard().btn_Done()

        #Delete the .run file
        RunPath = os.path.join(self.AlreadyRunning_Dir,'.frun')
        
        if os.path.isfile(RunPath):
            try:
                os.remove(RunPath)
            except OSError:
                time.sleep(0.1)
                os.remove(RunPath)       

        if self.PyMOL:
            cmd.set_wizard()
            cmd.set_wizard()
        
        # Kill main application window
        self.top.destroy()
   
        
    ''' ==================================================================================
    FUNCTION ValidateResiduValue: Validate the residue entered 
    ================================================================================== '''
    def ValidateResiduValue(self, event):        
        
        term = self.ResiduValue.get().upper()
        self.ResiduValue.set(term)
        
        NbChar = len(term)
        self.ValidResn = False
        
        if NbChar > 0:
            # Be sure the number of characters do not exceed 10
            if NbChar > 10:
                term = term[0:10]
                self.ResiduValue.set(term)
        
            Notfound = True            
            
            self.NoResn = 0 
            for elem in self.listResidu:                
                if elem.startswith(term):                    
                    if (NbChar == len(elem)) or (elem.startswith(term + '-')):
                        self.ValidResn = True
                    Notfound = False
                    break
                self.NoResn += 1
                
            if Notfound:
                self.EntryResidu.config(bg=self.Color_Red)
                self.Btn_AddResidu.config(state='disabled')
                
            else:
                if self.ValidResn:
                    self.EntryResidu.config(bg=self.Color_Green)
                    self.Btn_AddResidu.config(state='normal')
                else:
                    self.EntryResidu.config(bg=self.Color_White)
                    self.Btn_AddResidu.config(state='disabled')
        else:
            self.EntryResidu.config(bg=self.Color_White)
            
    ''' ==================================================================================
    FUNCTION ValidateFolders: Be sure the folders Exists 
    ==================================================================================  '''    
    def ValidateFolders(self):
       
        if not(os.path.isdir(self.FlexAIDProject_Dir)):
            os.makedirs(self.FlexAIDProject_Dir)
            
        if not(os.path.isdir(self.FlexAIDLigandProject_Dir)):
            os.makedirs(self.FlexAIDLigandProject_Dir)
            
        if not(os.path.isdir(self.ProteinProject_Dir)):
            os.makedirs(self.ProteinProject_Dir)
            
        if not(os.path.isdir(self.FlexAIDSimulationProject_Dir)):
            os.makedirs(self.FlexAIDSimulationProject_Dir)
            
        if not(os.path.isdir(self.BindingSiteProject_Dir)):
            os.makedirs(self.BindingSiteProject_Dir)
            
        if not(os.path.isdir(self.FlexAIDTargetFlexProject_Dir)):
            os.makedirs(self.FlexAIDTargetFlexProject_Dir)
                                      
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
            print "Unknown data format"

    ''' ==================================================================================
    FUNCTION Validate_Entries: Validate all entries of a frame before switching
    ==================================================================================  '''    
    def Validate_Entries(self, list):
        
        #print list
        for valid in list:
            if not valid[2] and valid[0] == False:

                if valid[3] != None:
                    valid[3].config(bg=self.Color_Red)

                return valid[1]

        return 0

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

