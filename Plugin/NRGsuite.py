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
@title: Plugin_NRG suite (for Pymol)

@summary: Plugin that add NRG suite applications access via the pymol application

@organization: Najmanovich Research Group
@creation date:  Aug. 25, 2010
'''

import tkSimpleDialog
import sys, os, stat
import platform
import tkMessageBox
import threading
import time

# IMPORTANT!!! NO UNDERSCORE or MINUS in the directories names
HomeDir = os.path.expanduser('~')
   
# DETECT the operating system
OS = platform.system().upper()
OSid = 'UNKNOWN'

# SET the configuration location
if OS == 'LINUX':
    OSid = 'LINUX'
elif OS == 'DARWIN':       
    OSid = 'MAC'    
elif OS == 'WINDOWS' or OS == 'MICROSOFT' or OS == 'WIN32':
    OSid = 'WIN'
else:
    OSid = 'UNKNOWN'


if OSid == 'WIN':
    ProgramFiles = os.getenv('PROGRAMFILES')
    if ProgramFiles:
        RootDir = os.path.join(ProgramFiles,'NRGsuite')
    else:
        RootDir = os.path.join('C:\\','Program Files','NRGsuite')
else:
    RootDir = os.path.join('/usr','local','NRGsuite')


if OSid != 'UNKNOWN':
    
    if os.path.isdir(RootDir):        
        
        RunDir_Path = os.path.join(HomeDir,'Documents','NRGsuite')
        Project_Path = os.path.join(RootDir,'Project')
        FlexAID_Path = os.path.join(RootDir,'FlexAID')
        GetCleft_Path = os.path.join(RootDir,'GetCleft')
        IsoCleft_Path = os.path.join(RootDir,'IsoCleft')
        
        if os.path.isdir(Project_Path) and os.path.isdir(FlexAID_Path) and os.path.isdir(GetCleft_Path) and os.path.isdir(IsoCleft_Path):
            
            if os.path.isfile(os.path.join(Project_Path,'NewProject.py'))  and \
                    os.path.isfile(os.path.join(FlexAID_Path,'FlexAID.py')) and \
                    os.path.isfile(os.path.join(GetCleft_Path,'GetCleft.py')):# and \
                    #os.path.isfile(IsoCleft_Path + 'IsoCleft.py'):
        
                sys.path.append(RootDir)

                sys.path.append(Project_Path)
                import NewProject
                import LoadProject
        
                sys.path.append(FlexAID_Path)
                import FlexAID

                sys.path.append(GetCleft_Path)
                import GetCleft
                
                sys.path.append(IsoCleft_Path)
                import IsoCleft
                
                from Tkinter import *
                from pymol import cmd
                
                #------------------------------------------------------------------#
                #                        SETTING THE MENUS                         #
                #------------------------------------------------------------------# 
                
                def __init__(self):
                    
                    # To Handle the Creation and load of Project                     
                    self.StopThread = False
                    self.EnableMenu = False
                    
                    self.ProjectName = ''
                    self.UserPath = ''
                    self.UserFlexAID = ''
                    self.UserGetCleft = ''
                    self.UserIsoCleft = ''
                    
                    self.RootFlexAID = None
                    self.RootGetCleft = None
                    self.RootIsoCleft = None
                    self.ActiveWizard = None

                    # ADD a new main menu named Project to the menu bar
                    #self.menuBar.addmenu('NRGsuite', 'Najmanovich Research Group Suite')
                    self.menuBar.addcascademenu('Plugin',
                                'NRGsuite',
                                'NRGsuite',
                                label = 'NRG Suite')

                    
                    #------------------------------------------------------------------
                    # ADD a menu named Open New Project to the NRGsuite main menu
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   New Project...',
                                      command = lambda s=self : projThread(s, 'NewProject'))
                    
                    #------------------------------------------------------------------
                    # ADD a menu named Open New Project to the NRGsuite main menu
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Load Project...',
                                      command = lambda s=self : projThread(s, 'LoadProject'))
                    
                    #------------------------------------------------------------------
                    # ADD a menu named Open New Project to the NRGsuite main menu
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Close Project',
                                      command = lambda s=self : StartCloseProject(s))
                    
                    # ADD a seperator to the main menu 
                    self.menuBar.addmenuitem('NRGsuite','separator')
                
                    #------------------------------------------------------------------
                    # ADD a menu named Preferences
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Preferences',
                                      command = lambda s=self : LoadPreferences(s))
                    
                    # ADD a seperator to the main menu 
                    self.menuBar.addmenuitem('NRGsuite','separator')

                    #------------------------------------------------------------------
                    # ADD a menu named Open FlexAID to the NRGsuite main menu
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Open FlexAID...',
                                      command = lambda s=self : StartFlexAID(s))                    
                    
                 
                    #------------------------------------------------------------------
                    # ADD a menu named Open GetCleft to the NRGsuite main menu
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Open GetCleft...',
                                      command = lambda s=self : StartGetCleft(s))
                    
                    #------------------------------------------------------------------
                    # ADD a menu named Open IsoCleft to the NRGsuite main menuimport threading
                    self.menuBar.addmenuitem('NRGsuite', 'command',
                                      'CfgFile',
                                      label='   Open IsoCleft...',
                                      command = lambda s=self : StartIsoCleft(s))
                    
                    EnableDisableMenu(self, False, True, '', '', '')
                 
                 
                #================================================================================
                # STARTING GetCleft Conditional... 
                #================================================================================  
                def StartGetCleft(self):
                    
                    Stopped = True
                    
                    #Detect is GetCleft is already running...
                    RunPath = os.path.join(RunDir_Path,'.grun')
                    if os.path.isfile(RunPath):
                        file = open(RunPath)
                        Runfile = file.readline()
                        file.close()                       
                        
                        ProcID = str(os.getpid())
                        
                        if Runfile.startswith(ProcID):
                            Stopped = False
                            print('\n   *** GetCleft is ALREADY started! ***')
                    
                    if Stopped:
                        print('   Starting GetCleft...')
                        self.RootGetCleft = Toplevel(self.root)
                        GetCleft.displayGetCleft(self.RootGetCleft, self.ActiveWizard, self.UserPath, RootDir, RunDir_Path, OSid, True)
                        
                        
                #================================================================================
                # STARTING IsoCleft Conditional... 
                #================================================================================  
                def StartIsoCleft(self):
                    
                    Stopped = True
                    
                    #Detect is IsoCleft is already running...
                    RunPath = os.path.join(RunDir_Path,'.irun')
                    if os.path.isfile(RunPath):
                        file = open(RunPath)
                        Runfile = file.readline()
                        file.close()                       
                        
                        ProcID = str(os.getpid())
                        
                        if Runfile.startswith(ProcID):
                            Stopped = False
                            print('\n   *** IsoCleft is ALREADY started! ***')
                    
                    if Stopped:
                        print('   Starting IsoCleft...')
                        self.RootIsoCleft = Toplevel(self.root)
                        IsoCleft.displayIsoCleft(self.RootIsoCleft, self.UserPath, IsoCleft_Path, RootDir, RunDir_Path, OSid)
                        
                        
                #================================================================================
                # STARTING FlexAID Conditional... 
                #================================================================================  
                def StartFlexAID(self):
                    
                    Stopped = True

                    #Detect is FlexAID is already running...
                    RunPath = os.path.join(RunDir_Path,'.frun')
                    if os.path.isfile(RunPath):
                        file = open(RunPath)
                        Runfile = file.readline()
                        file.close()                       
                        
                        ProcID = str(os.getpid())
                        
                        if Runfile.startswith(ProcID):
                            Stopped = False
                            print('\n   *** FlexAID is ALREADY started! ***')
                    
                    if Stopped:
                        print('   Starting FlexAID...')
                        self.RootFlexAID = Toplevel(self.root)
                        FlexAID.displayFlexAID(self.RootFlexAID, self.ActiveWizard, self.UserPath, RootDir, RunDir_Path, OSid, True)
                        
                
                #================================================================================        
                # Set the NRGsuite Menu Options to Enable when a Project is created or loaded.
                # Set to disable at start up or when a project is closed.
                #================================================================================
                def EnableDisableMenu(self, boolEnable, boolProj, SelProject, UserPath, ProjectName):
                    
                    if SelProject == 'NewProject':
                        if boolProj:
                            self.menuBar.component('NRGsuite-menu').entryconfig(1, state='normal')
                        else:
                            self.menuBar.component('NRGsuite-menu').entryconfig(1, state='disabled')
                            
                    elif SelProject == 'LoadProject':
                        if boolProj:
                            self.menuBar.component('NRGsuite-menu').entryconfig(0, state='normal')
                        else:
                            self.menuBar.component('NRGsuite-menu').entryconfig(0, state='disabled')            
                    
                    if boolEnable:
                        self.ProjectName = ProjectName
                        UpdatePath(self, UserPath)
                        self.menuBar.component('NRGsuite-menu').entryconfig(0, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(1, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(2, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(4, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(6, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(7, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(8, state='disabled')                       
                    else:
                        self.menuBar.component('NRGsuite-menu').entryconfig(0, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(1, state='normal')
                        self.menuBar.component('NRGsuite-menu').entryconfig(2, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(4, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(6, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(7, state='disabled')
                        self.menuBar.component('NRGsuite-menu').entryconfig(8, state='disabled')
                        
                
                #================================================================================
                # UpdatePath the General Path
                #================================================================================
                def UpdatePath(self, UserPath):
                    
                    self.UserPath = UserPath
                    self.UserFlexAID = os.path.join(UserPath,'FlexAID')
                    self.UserGetCleft = os.path.join(UserPath,'GetCleft')
                    self.UserIsoCleft = os.path.join(UserPath,'IsoCleft')

                #================================================================================
                # Loads the preferences menu to set the different options
                #================================================================================
                def LoadPreferences(self):
                    
                    return
                
                #================================================================================
                # Close the actual Project
                #================================================================================
                def StartCloseProject(self):

                    if not self.ActiveWizard is None:
                        print "The Wizard must be closed before you can Close this project"
                        return
                    
                    EnableDisableMenu(self, False, True, '', '', '')
                    
                    ProcID = str(os.getpid())
                    
                    # Close all the Opened Interfaces
                    # Close FlexAID
                    FlexRunPath = os.path.join(RunDir_Path,'.frun')
                    if os.path.isfile(FlexRunPath):
                       file = open(FlexRunPath)
                       RunfileFlex = file.readline()
                       file.close()                       
                       
                       if RunfileFlex.startswith(ProcID):
                           try:
                                os.remove(FlexRunPath)
                           except OSError:
                                time.sleep(0.1)
                                os.remove(FlexRunPath)
                                           
                           self.RootFlexAID.destroy()
                           print('\n   *** CLOSING FlexAID! ***')
                           
                    # Close GetCleft
                    GetCleftRunPath = os.path.join(RunDir_Path,'.grun')
                    if os.path.isfile(GetCleftRunPath):
                       file = open(GetCleftRunPath)
                       RunfileFlex = file.readline()
                       file.close()
                       
                       if RunfileFlex.startswith(ProcID):
                           try:
                                os.remove(GetCleftRunPath)
                           except OSError:
                                time.sleep(0.1)
                                os.remove(GetCleftRunPath)
                                                
                           self.RootGetCleft.destroy()
                           print('\n   *** CLOSING GetCleft! ***')
                           
                    # Close IsoCleft
                    IsoCleftRunPath = os.path.join(RunDir_Path,'.irun')
                    if os.path.isfile(IsoCleftRunPath):
                       file = open(IsoCleftRunPath)
                       RunfileFlex = file.readline()
                       file.close()
                       
                       if RunfileFlex.startswith(ProcID):
                           try:
                                os.remove(IsoCleftRunPath)
                           except OSError:
                                time.sleep(0.1)
                                os.remove(IsoCleftRunPath)
                
                           self.RootIsoCleft.destroy()
                           print('\n   *** CLOSING IsoCleft! ***')
                           
                    print('\n   CLOSING Project: ' + self.ProjectName)
            
            
            #================================================================================
            # ERRORS!!
            #================================================================================
            else:

                if not os.path.isfile(os.path.join(GetCleft_Path,'GetCleft.py')):
                    tkMessageBox.showerror('GetCleft.py not found',
                                          'The GetCleft.py file need to be present at this location: ' + 
                                           os.path.join(GetCleft_Path,'GetCleft.py') + '\n\nInstallation CANCELLED...')

                elif not os.path.isfile(os.path.join(FlexAID_Path,'FlexAID.py')):
                    tkMessageBox.showerror('FlexAID.py not found',
                                          'The FlexAID.py file need to be present at this location: ' + 
                                           os.path.join(FlexAID_Path,'FlexAID.py') + '\n\nInstallation CANCELLED...')

                #elif not os.path.isfile(os.path.join(IsoCleft_Path,'IsoCleft.py')):
                #    tkMessageBox.showerror('IsoCleft.py not found',
                #                          'The IsoCleft.py file need to be present at this location: ' + 
                #                           os.path.join(IsoCleft_Path,'IsoCleft.py') + '\n\nInstallation CANCELLED...')
           
        else:
            tkMessageBox.showerror('Programs MISSING',
                                   'The FlexAID, GetCleft, and IsoCleft applications ' + 
                                   'need to be installed at this location: ' +
                                   RootDir + '\n\nInstallation CANCELLED...')

    else:
        tkMessageBox.showerror('NRGsuite not found',
                               'The NRGsuite applications need to be installed at this location: ' + 
                               RootDir + '\n\nInstallation CANCELLED...')
else:
    tkMessageBox.showerror('System Not Recognized',
                           'The plugin installation was interrupted because the operating system used is unknown...') 
    
    
#=========================================================================================
'''                        ---   STARTING Project  ---                                '''
#=========================================================================================     
class projThread(threading.Thread):   

    def __init__(self, top, Selection):
       
        threading.Thread.__init__(self)
        
        self.top = top 
        
        self.Selection = Selection
        self.EnableMenu = False
        self.StopThread = False
        self.UserPath = ''
        self.ProjectName = ''
        
        EnableDisableMenu(self.top, False, False, self.Selection, '', '')
        
        if Selection == 'NewProject':
            self.top.menuBar.component('NRGsuite-menu').entryconfig(1, state='disabled')
            self.StartNewProject()

        elif Selection == 'LoadProject':
            self.top.menuBar.component('NRGsuite-menu').entryconfig(0, state='disabled')   
            self.StartLoadProject()
            
        
    #================================================================================
    # Start the Thread
    #================================================================================       
    def run(self):
           
        while(1):
            time.sleep(0.1)
            
            if self.StopThread:
                break
        
        if self.EnableMenu:    
            EnableDisableMenu(self.top, True, True, self.Selection, self.UserPath, self.ProjectName)
        else:
            EnableDisableMenu(self.top, False, True, self.Selection, self.UserPath, self.ProjectName)
            

    #================================================================================
    # Create a New Project
    #================================================================================
    def StartNewProject(self):
        
        self.EnableMenu = False
        self.StopThread = False
        
        Stopped = True
        
        #Detect if Project is already running...
        RunPath = os.path.join(RunDir_Path,'.nrun')
        if os.path.isfile(RunPath):
            file = open(RunPath)
            Runfile = file.readline()
            file.close()                       
            
            ProcID = str(os.getpid())
            
            if Runfile.startswith(ProcID):
                Stopped = False
                print('\n   *** New Project is ALREADY started! ***')
        
        if Stopped:
            print('\n   Create a new NRG suite Project...')
            NewProject.displayNewProject(self.top.root, self, Project_Path, HomeDir, RunDir_Path, OSid)
            self.start()
            
            
    #================================================================================
    # Load an existing Project
    #================================================================================
    def StartLoadProject(self):
        
        self.EnableMenu = False
        self.StopThread = False        

        Stopped = True
        
        #Detect if Project is already running...
        RunPath = os.path.join(RunDir_Path,'.lrun')
        if os.path.isfile(RunPath):
            file = open(RunPath)
            Runfile = file.readline()
            file.close()                       
            
            ProcID = str(os.getpid())
            
            if Runfile.startswith(ProcID):
                Stopped = False
                print('\n   *** Load Project is ALREADY started! ***')
        
        if Stopped:
            print('\n   Load an existing NRG suite Project...')
            LoadProject.displayLoadProject(self.top.root, self, Project_Path, HomeDir, RunDir_Path, OSid)
            self.start()                   
