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
@title: FlexAID - Interface

@summary: This is the interface of FlexAID application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Aug. 25, 2010
'''

import sys
if sys.version_info[0] < 3:
    from Tkinter import *
    import tkFileDialog
else:
    from tkinter import *
    import tkinter.filedialog as tkFileDialog

import os
import sys
import pickle
import shutil
import random
import string

import Base
import General
import IOFile
import Config1
import Config2
import Config3
import GAParam
import Simulate

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayFlexAID(Base.Base):
    
    def Init_Vars(self):
        
        self.SessionPath = ''
        self.SaveSessionFile = ''
        
    def After_Init(self):
                
        # default Tab (IOFile)
        self.Btn_IOFiles_Clicked()
        
        # By default hide advanced tabs
        if self.Prefs.AlwaysShowAdvancedView == 1:
            self.bAdvancedView = True
        else:
            self.bAdvancedView = False
        self.Btn_Toggle_AdvView()

        # set files to copy when saving a session
        self.SessionVars = [ self.IOFile.Vars.LigandPath, self.IOFile.Vars.ProcessedLigandPath,
                             self.IOFile.Vars.ProcessedLigandINPPath, self.IOFile.Vars.ProcessedLigandICPath,
                             self.IOFile.Vars.TargetPath, self.IOFile.Vars.ProcessedTargetPath ]
        
    def Build_Tabs(self):
        
        # Build class objects of each tab
        self.IOFile = IOFile.IOFile(self, self.PyMOL, self.Btn_IOFiles, 'IOFile', IOFile.IOFileVars(), self.Prefs)
        self.Config1 = Config1.Config1(self, self.PyMOL, self.Btn_Config1, 'Config1', Config1.Config1Vars(), self.Prefs)
        self.Config2 = Config2.Config2(self, self.PyMOL, self.Btn_Config2, 'Config2', Config2.Config2Vars(), self.Prefs)
        self.Config3 = Config3.Config3(self, self.PyMOL, self.Btn_Config3, 'Config3', Config3.Config3Vars(), self.Prefs)
        self.GAParam = GAParam.GAParam(self, self.PyMOL, self.Btn_GAParam, 'GAParam', GAParam.GAParamVars(), self.Prefs)
        self.Simulate = Simulate.Simulate(self, self.PyMOL, self.Btn_Simulate, 'Simulate', Simulate.SimulateVars(), self.Prefs)

        self.listTabs = [self.IOFile, self.Config1, self.Config2, self.Config3, self.GAParam, self.Simulate]
        
        self.listBtnTabs = [self.Btn_IOFiles, self.Btn_Config1, self.Btn_Config2, 
                            self.Btn_Config3, self.Btn_GAParam, self.Btn_Simulate]
        
        return

    def Clean(self):
        
        files = os.listdir(self.FlexAIDTempProject_Dir)

        for file in files:
            try:
                os.remove(os.path.join(self.FlexAIDTempProject_Dir,file))
            except OSError:
                pass
    
    def Set_Folders(self):

        self.FlexAIDInstall_Dir = os.path.join(self.Install_Dir,'FlexAID')

        if self.OSid == 'WIN':
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID.exe')
            self.Process_LigandExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','Process_Ligand.exe')
        else:
            self.FlexAIDExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','FlexAID')
            self.Process_LigandExecutable = os.path.join(self.FlexAIDInstall_Dir,'WRK','Process_Ligand')

        self.FlexAIDProject_Dir = os.path.join(self.Project_Dir,'FlexAID')
        self.GetCleftProject_Dir = os.path.join(self.Project_Dir,'GetCleft')
        self.GetCleftSaveProject_Dir = os.path.join(self.GetCleftProject_Dir,'.Save')
        
        self.CleftProject_Dir = os.path.join(self.Project_Dir,'Cleft')
        self.TargetProject_Dir = os.path.join(self.Project_Dir,'Target')

        self.FlexAIDLigandProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Ligand')
        self.FlexAIDSimulationProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Simulation')
        self.FlexAIDSessionProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Session')
        self.FlexAIDTempProject_Dir = os.path.join(self.FlexAIDProject_Dir,'.Temp')
        self.FlexAIDSaveProject_Dir = os.path.join(self.FlexAIDProject_Dir,'.Save')
        self.FlexAIDResultsProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Results')
        self.FlexAIDBindingSiteProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Binding_Site')
        self.FlexAIDTargetFlexProject_Dir = os.path.join(self.FlexAIDProject_Dir,'Target_Flexibility')

        self.Folders.extend( [  self.FlexAIDProject_Dir, self.GetCleftProject_Dir, self.GetCleftSaveProject_Dir,
                                self.CleftProject_Dir, self.TargetProject_Dir, self.FlexAIDLigandProject_Dir,
                                self.FlexAIDSimulationProject_Dir, self.FlexAIDSessionProject_Dir, self.FlexAIDResultsProject_Dir,
                                self.FlexAIDBindingSiteProject_Dir, self.FlexAIDTargetFlexProject_Dir,
                                self.FlexAIDSaveProject_Dir, self.FlexAIDTempProject_Dir ] )
        
    ''' ==================================================================================
    FUNCTION MakeMenuBar: Builds the menu on the upper left corner    
    ==================================================================================  '''        
    def MakeMenuBar(self):
        
        self.menubar = Menu(self.root)
        
        loadmenu = Menu(self.menubar, tearoff=0)
        loadmenu.add_command(label="Load session", command=self.Btn_Load_Session)
        self.menubar.add_cascade(label="Load", menu=loadmenu)

        savemenu = Menu(self.menubar, tearoff=0)        
        savemenu.add_command(label="Save session", command=self.Btn_Save_Session)
        savemenu.add_command(label="Save session as...", command=self.Btn_SaveAs_Session)
        self.menubar.add_cascade(label="Save", menu=savemenu)

        viewmenu = Menu(self.menubar, tearoff=0)
        viewmenu.add_command(label="Advanced view", command=self.Btn_Toggle_AdvView)
        self.menubar.add_cascade(label="View", menu=viewmenu)
        
        self.root.config(menu=self.menubar)
    
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes the trace of some variables
    ================================================================================== '''    
    def Del_Trace(self):
    
        for Tab in self.listTabs:
            Tab.Del_Trace()
    
    ''' ==================================================================================
    FUNCTION Btn_Toggle_AdvView: Hides/shows the advanced tabs (scoring + ga)
    ==================================================================================  '''  
    def Btn_Toggle_AdvView(self):
    
        if not self.bAdvancedView:
            self.Btn_Config3.pack_forget()
            self.Btn_GAParam.pack_forget()
            self.bAdvancedView = True
            
            if self.ActiveFrame == self.Config3 or self.ActiveFrame == self.GAParam:
                self.Btn_IOFiles_Clicked()

            self.menubar.entryconfig(3,label="Show")
            
        else:
            self.Btn_Simulate.pack_forget()
            self.Btn_Config3.pack(side=LEFT, fill=BOTH, expand=True)
            self.Btn_GAParam.pack(side=LEFT, fill=BOTH, expand=True)
            self.Btn_Simulate.pack(side=LEFT, fill=BOTH, expand=True)
            self.bAdvancedView = False
            
            self.menubar.entryconfig(3,label="Hide")
            
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''  
    def Frame_Main(self):
   
        #==================================================================================
        '''                  --- TOP MENU OF THE INTERFACE ---                          '''
        #==================================================================================

        self.fTop = Frame(self.fMain, relief=RIDGE, border=4, height=50)
        self.fTop.pack(fill=BOTH, expand=True)#, padx=10, pady=10, ipady=10, ipadx=10, side=TOP)
        self.fTop.pack_propagate(0)

        self.Btn_IOFiles = Button(self.fTop, text='Input Files', bg=self.Color_White, command=self.Btn_IOFiles_Clicked, font=self.font_Text)
        self.Btn_IOFiles.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_IOFiles.config(state='normal')

        self.Btn_Config1 = Button(self.fTop, text='Target Cfg', bg=self.Color_Grey, command=self.Btn_Config1_Clicked, font=self.font_Text)
        self.Btn_Config1.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config1.config(state='disabled')

        self.Btn_Config2 = Button(self.fTop, text='Ligand Cfg', bg=self.Color_Grey, command=self.Btn_Config2_Clicked, font=self.font_Text)
        self.Btn_Config2.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config2.config(state='disabled')

        self.Btn_Config3 = Button(self.fTop, text='Scoring Cfg', bg=self.Color_Grey, command=self.Btn_Config3_Clicked, font=self.font_Text)
        self.Btn_Config3.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Config3.config(state='disabled')
        
        self.Btn_GAParam = Button(self.fTop, text='GA Param', bg=self.Color_Grey, command=self.Btn_GAParam_Clicked, font=self.font_Text)
        self.Btn_GAParam.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_GAParam.config(state='disabled')    

        self.Btn_Simulate = Button(self.fTop, text='Simulate', bg=self.Color_Grey, command=self.Btn_Simulate_Clicked, font=self.font_Text)
        self.Btn_Simulate.pack(side=LEFT, fill=BOTH, expand=True)
        self.Btn_Simulate.config(state='disabled') 
                
        #==================================================================================
        '''                MIDDLE DISPLAY SECTION OF THE INTERFACE                      '''
        #==================================================================================
           
        self.fMiddle = Frame(self.fMain)#, height=355)
                
        self.fMiddle.pack(fill=X, expand=True, padx=10, side=TOP)
        #self.fMiddle.config(takefocus=1)
        #self.fMiddle.bind('<FocusIn>', lambda e: None)
        #self.fMiddle.pack_propagate(0)
        
        #==================================================================================
        '''                 BOTTOM DISPLAY SECTION OF THE INTERFACE                     '''
        #==================================================================================
 
        fBottom = Frame(self.fMain, border=1, relief=SUNKEN)
        fBottom.pack(fill=BOTH, expand=True, padx=10, pady=10, ipadx=10, ipady=10, side=TOP)

        fBottomRight = Frame(fBottom)
        fBottomRight.pack(fill=Y, side=RIGHT)

        Btn_Default = Button(fBottomRight, text='Default', width=15, command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Default.pack(side=TOP, fill=X)

        #Btn_SaveDefault = Button(fBottomRight, text='Save as default', command=self.Btn_SaveDefault_Clicked, font=self.font_Text)
        #Btn_SaveDefault.pack(side=TOP, fill=X)

        #Btn_Restore = Button(fBottomRight, text='Restore', command=self.Btn_Restore_Clicked, font=self.font_Text)
        #Btn_Restore.pack(side=TOP, fill=X)

        Btn_Quit = Button(fBottomRight, text='Close', width=15, command=self.Quit, font=self.font_Text)
        Btn_Quit.pack(side=BOTTOM, fill=X)

        fBottomLeft = Frame(fBottom)
        fBottomLeft.pack(side=LEFT, fill=Y, ipadx=20, ipady=20)

        scrollBar = Scrollbar(fBottomLeft)
        scrollBar.pack(side=RIGHT, fill=Y)

        self.TextMessage = Text(fBottomLeft, border=1, background=self.Color_Grey, font=self.font_Text)
        self.TextMessage.pack(side=RIGHT, fill=BOTH, expand=True)

        scrollBar.config(command=self.TextMessage.yview)
        self.TextMessage.config(state='disabled', yscrollcommand=scrollBar.set)                                       
    
    ''' ==================================================================================
    FUNCTION Btn_Load_Session: Loads a previously saved session
    ==================================================================================  '''        
    def Btn_Load_Session(self):

        if self.ValidateProcessRunning() or self.ValidateWizardRunning() or \
            self.ValidateWindowRunning():
            return

        LoadFile = tkFileDialog.askopenfilename(initialdir=self.FlexAIDSessionProject_Dir,
                                                filetypes=[('NRG FlexAID Session','*.nrgfs')],
                                                title='Select the Session to load')
        # LoadFile = self.root.master.splitlist(LoadFile)
        if len(LoadFile) > 0:
            
            self.Btn_IOFiles_Clicked()
            
            LoadFile = os.path.normpath(LoadFile)
            
            try:
                in_ = open(LoadFile, 'rb')
                self.SessionPath = pickle.load(in_)
                
                for Tab in self.listTabs:
                    Tab.Vars = pickle.load(in_)
                    Tab.Vars.refresh()
                    
                    if Tab.Check_Integrity():
                        self.Reset_All()
                        self.DisplayMessage("  ERROR: The loading of the session failed the integrity check.", 2)
                        return

                    Tab.Load_Session()
                    Tab.Def_Vars()
                        
                in_.close()
                
                self.SaveSessionFile = LoadFile
                self.DisplayMessage("  The session '" + os.path.split(LoadFile)[1] + "' was loaded successfully.", 2)
                
                return

            except pickle.UnpicklingError:
                self.DisplayMessage("  ERROR: Could not properly load the session: " + str(sys.exc_info()), 2)
            except:
                self.DisplayMessage("  Unexpected error: " + str(sys.exc_info()), 2)

            self.SessionPath = ''
            self.SaveSessionFile = ''
                
    ''' ==================================================================================
    FUNCTION ValidateSaveProject: Validates if a file was saved in the right folder
    ==================================================================================  '''    
    def ValidateSaveProject(self, file, objtype):
    
        basepath, filename = os.path.split(file)
        
        if objtype == 'Ligand' and not self.FlexAIDLigandProject_Dir in basepath:
            return 1
        elif objtype == 'Target' and not self.TargetProject_Dir in basepath:
            return 1
        elif objtype == 'Session' and not self.FlexAIDSessionProject_Dir in basepath:
            return 1
        elif objtype == 'Results' and not self.FlexAIDResultsProject_Dir in basepath:
            return 1
        elif objtype == 'TargetFlex' and not self.FlexAIDTargetFlexProject_Dir in basepath:
            return 1
        elif objtype == 'BindingSite' and not self.FlexAIDBindingSiteProject_Dir in basepath:
            return 1
        
        return 0

    ''' ==================================================================================
    FUNCTION Save_SessionFile: Saves the current session
    ==================================================================================  '''        
    def Save_SessionFile(self, SaveFile):

        if len(SaveFile) > 0:
            
            SaveFile = os.path.normpath(SaveFile)
            
            if General.validate_String(SaveFile, '.nrgfs', True, True, False):
                self.DisplayMessage("  ERROR: Could not save the file because you entered an invalid name.", 2)
                return
            
            if self.ValidateSaveProject(SaveFile, 'Session'):
                self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                return
            
            try:
                out = open(SaveFile, 'wb')
                pickle.dump(self.SessionPath, out)
                
                for Tab in self.listTabs:
                    pickle.dump(Tab.Vars, out)
                out.close()
                
                self.SaveSessionFile = SaveFile
                self.DisplayMessage("  The session '" + os.path.split(SaveFile)[1] + "' was saved successfully.", 2)
                
            except pickle.PicklingError:
                self.DisplayMessage("  ERROR: Could not properly save the session: " + str(sys.exc_info()), 2)
            except:
                self.DisplayMessage("  Unexpected error: " + str(sys.exc_info()), 2)
                
    ''' ==================================================================================
    FUNCTION Btn_Save_Session: Saves the current session
    ==================================================================================  '''        
    def Btn_Save_Session(self):

        if self.ValidateProcessRunning() or self.ValidateWizardRunning() or \
            self.ValidateWindowRunning():
            return
            
        if self.SaveSessionFile == '':
            self.Btn_SaveAs_Session()
        else:
            self.Copy_SessionFiles()
            self.Save_SessionFile(self.SaveSessionFile)

    ''' ==================================================================================
    FUNCTION Btn_SaveAs_Session: Saves the current session under a new name
    ==================================================================================  '''        
    def Btn_SaveAs_Session(self):

        if self.ValidateProcessRunning() or self.ValidateWizardRunning() or \
            self.ValidateWindowRunning():
            return
        
        SaveFile = tkFileDialog.asksaveasfilename(initialdir=self.FlexAIDSessionProject_Dir,
                                  title='Save the Session file', initialfile='default_session',
                                  filetypes=[('NRG FlexAID Session','*.nrgfs')], defaultextension='.nrgfs')
        
        if len(SaveFile):
            
            while True:
                Session = ''.join(random.choice(string.ascii_lowercase + string.digits) for x in range(30))
                self.SessionPath = os.path.join(self.FlexAIDSaveProject_Dir,Session)
                
                if not os.path.isdir(self.SessionPath):
                    os.makedirs(self.SessionPath)
                    break

            self.Copy_SessionFiles()
            self.Save_SessionFile(SaveFile)
        
    ''' ==================================================================================
    FUNCTION Copy_SessionFiles: Moves the files that are project-specific
    ==================================================================================  '''        
    def Copy_SessionFiles(self):
        
        for Var in self.SessionVars:
            if Var.get():
                try:
                    shutil.copy(Var.get(), self.SessionPath)
                    Var.set(os.path.join(self.SessionPath,os.path.split(Var.get())[1]))
                except shutil.Error:
                    self.DisplayMessage("The following session file was not copied: " + Var.get(), 1)
            
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

    ''' ==================================================================================
    FUNCTION Go_Step1: Enables/Disables buttons for step 1
    ================================================================================== '''    
    def Go_Step1(self):

        #print "Setting Tab buttons to Step 1"        
        for Btn in self.listBtnTabs:
            Btn.config(state='disabled',bg=self.Color_Grey)

        self.Btn_IOFiles.config(state='normal',bg=self.Color_Blue)
    
    ''' ==================================================================================
    FUNCTION Go_Step2: Enables/Disables buttons for step 2
    ================================================================================== '''    
    def Go_Step2(self):
        
        #print "Setting Tab buttons to Step 2"        
        for Btn in self.listBtnTabs:
            Btn.config(state='normal',bg=self.Color_White)

        self.Btn_IOFiles.config(bg=self.Color_Blue)
        
    ''' ==================================================================================
    FUNCTION Reset_Step2: Reset all except IOFiles
    ==================================================================================  '''    
    def Reset_Step2(self):
        
        for Tab in self.listTabs:
            if Tab != self.IOFile:
                Tab.Init_Vars()
    
    ''' ==================================================================================
    FUNCTION Reset_All: Reset ALL the parameters
    ==================================================================================  '''    
    def Reset_All(self):
        
        for Tab in self.listTabs:
            Tab.Init_Vars()
    
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
            
