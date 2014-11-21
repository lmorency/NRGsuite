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
@title: Plugin_NRG suite (for Pymol)

@summary: Plugin that add NRG suite applications access via the pymol application

@organization: Najmanovich Research Group
@creation date:  Aug. 25, 2010
'''

import sys, os, stat
import platform
import tkMessageBox
import time

'''=================================================================================================
FUNCTION get_default_path_for_OSid: returns the default (and preferred) installation path for any OS
================================================================================================='''
def get_default_path_for_OSid():
    # DETECT the operating system
    OS = platform.system().upper()
    OSid = 'UNKNOWN'

    # SET the configuration location
    if OS == 'LINUX' or OS == 'BSD':
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
            Install_Dir = os.path.join(ProgramFiles,'NRGsuite')
        else:
            Install_Dir = os.path.join('C:\\','Program Files','NRGsuite')
            
    elif OSid == 'MAC':
        Install_Dir = os.path.join('/Applications','NRGsuite')

    else:
        Install_Dir = os.path.join('/usr','local','NRGsuite')

    return Install_Dir

# get the installation path of the NRGsuite with the ENV variable NRGSUITE or
# sets Install_Dir to the default (and preferred) installation directory
Install_Dir = os.environ.get('NRGSUITE_INSTALLATION', get_default_path_for_OSid())
 
if Install_Dir is '' or not os.path.isdir(Install_Dir):
    Install_Dir = get_default_path_for_OSid()

if os.path.isdir(Install_Dir):
    
    NRGsuite_Path = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite')
    Plugin_Path = os.path.join(Install_Dir,'Plugin')
    Project_Path = os.path.join(Install_Dir,'Project')
    FlexAID_Path = os.path.join(Install_Dir,'FlexAID')
    GetCleft_Path = os.path.join(Install_Dir,'GetCleft')
    About_Path = os.path.join(Install_Dir,'About')
    
    if os.path.isdir(Project_Path) and os.path.isdir(FlexAID_Path) and \
       os.path.isdir(GetCleft_Path) and os.path.isdir(About_Path):
    
        sys.path.append(Install_Dir)
        import Prefs
        
        sys.path.append(Plugin_Path)
        
        sys.path.append(Project_Path)
        import NewProject
        import LoadProject
        
        sys.path.append(FlexAID_Path)
        import FlexAID
        
        sys.path.append(GetCleft_Path)
        import GetCleft
        
        sys.path.append(About_Path)
        import About
        

        from Tkinter import *
        
        #------------------------------------------------------------------#
        #                        SETTING THE MENUS                         #
        #------------------------------------------------------------------# 
        
        def __init__(self):
            
            self.ProjectName = ''
            self.Project_Dir = ''
            
            self.RootFlexAID = None
            self.RootGetCleft = None

            self.RootPrefs = Prefs.Prefs()

            # ADD a new main menu named Project to the menu bar
            #self.menuBar.addmenu('NRGsuite', 'Najmanovich Research Group Suite')
            self.menuBar.addcascademenu('Plugin',
                        'NRGsuite',
                        'NRGsuite',
                        label = 'NRGsuite')

            #------------------------------------------------------------------
            # ADD a menu named Open New Project to the NRGsuite main menu
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   New Project...',
                                     command = lambda s=self, menuindex=0 : StartNewProject(s, menuindex))
            
            #------------------------------------------------------------------
            # ADD a menu named Open New Project to the NRGsuite main menu
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   Load Project...',
                                     command = lambda s=self, menuindex=1 : StartLoadProject(s, menuindex))
            
            #------------------------------------------------------------------
            # ADD a menu named Open New Project to the NRGsuite main menu
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   Close Project',
                                     command = lambda s=self, menuindex=2 : StartCloseProject(s, menuindex))
            
            # ADD a seperator to the main menu 
            self.menuBar.addmenuitem('NRGsuite','separator')
        
            #------------------------------------------------------------------
            # ADD a menu named Preferences
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   Preferences',
                                     command = lambda s=self, menuindex=4 : StartPreferences(s, menuindex))
            
            # ADD a seperator to the main menu 
            self.menuBar.addmenuitem('NRGsuite','separator')

            #------------------------------------------------------------------
            # ADD a menu named Open FlexAID to the NRGsuite main menu
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   Open FlexAID...',
                                     command = lambda s=self, menuindex=6 : StartFlexAID(s, menuindex))
            
            #------------------------------------------------------------------
            # ADD a menu named Open GetCleft to the NRGsuite main menu
            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   Open GetCleft...',
                                     command = lambda s=self, menuindex=7 : StartGetCleft(s, menuindex))
            
            # ADD a seperator to the main menu 
            self.menuBar.addmenuitem('NRGsuite','separator')

            self.menuBar.addmenuitem('NRGsuite', 'command',
                                     'CfgFile',
                                     label='   About',
                                     command = lambda s=self, menuindex=9 : StartAbout(s, menuindex))
            
            # Make the Load  - Create projects button clickable only
            EnableDisableMenu(self, ['normal','normal','disabled','normal','disabled','disabled','normal'] )
        
        #================================================================================
        # STARTING GetCleft Conditional... 
        #================================================================================  
        def StartGetCleft(self, menuindex):
            
            self.RootGetCleft = Toplevel(self.root)
            GetCleft.displayGetCleft(self.RootGetCleft, self, menuindex, self.Project_Dir, Install_Dir,
                                     NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - GetCleft', 500, 550, self.RootPrefs)
                
        #================================================================================
        # STARTING FlexAID Conditional... 
        #================================================================================  
        def StartFlexAID(self, menuindex):
            
            self.RootFlexAID = Toplevel(self.root)
            FlexAID.displayFlexAID(self.RootFlexAID, self, menuindex, self.Project_Dir, Install_Dir, 
                                   NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - FlexAID', 800, 600, self.RootPrefs)
                
        #================================================================================
        # Load an existing Project
        #================================================================================
        def StartLoadProject(self, menuindex):
            
            LoadProject.displayLoadProject(Toplevel(self.root), self, menuindex, self.Project_Dir, Install_Dir,
                                           NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - Load Project', 460, 350, self.RootPrefs)
            
        #================================================================================
        # Create a New Project
        #================================================================================
        def StartNewProject(self, menuindex):
                    
            NewProject.displayNewProject(Toplevel(self.root), self, menuindex, self.Project_Dir, Install_Dir,
                                         NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - New Project', 420, 300, self.RootPrefs)
            
        #================================================================================
        # Loads the preferences menu to set the different options
        #================================================================================
        def StartPreferences(self, menuindex):
            Prefs.displayPrefs(Toplevel(self.root), self, menuindex, self.Project_Dir, Install_Dir,
                               NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - Preferences', 450, 250, self.RootPrefs)
            
        #================================================================================
        # Loads the about menu to see versionning
        #================================================================================
        def StartAbout(self, menuindex):
            
            About.displayAbout(Toplevel(self.root), self, menuindex, self.Project_Dir, Install_Dir,
                               NRGsuite_Path, self.RootPrefs.OSid, True, 'NRGsuite - About', 400, 350, self.RootPrefs)
        
        #================================================================================
        # Close the actual Project
        #================================================================================
        def StartCloseProject(self, menuitem):
                            
            if self.RootFlexAID is not None:
                self.RootFlexAID.destroy()
                print('  FlexAID interface successfully closed.')
            
            if self.RootGetCleft is not None:
                self.RootGetCleft.destroy()
                print('  GetCleft interface successfully closed.')
                   
            EnableDisableMenu(self, ['normal','normal','disabled','disabled','disabled','disabled','normal'] )
            
            print('\n   The project \'' + self.ProjectName + '\' was closed.')

        #================================================================================        
        # Set the NRGsuite Menu Options to Enable when a Project is created or loaded.
        # Set to disable at start up or when a project is closed.
        #================================================================================
        def EnableDisableMenu(self, states):
            
            self.menuBar.component('NRGsuite-menu').entryconfig(0, state=states.pop(0))    # new
            self.menuBar.component('NRGsuite-menu').entryconfig(1, state=states.pop(0))    # load
            self.menuBar.component('NRGsuite-menu').entryconfig(2, state=states.pop(0))    # close
            self.menuBar.component('NRGsuite-menu').entryconfig(4, state=states.pop(0))    # preferences
            self.menuBar.component('NRGsuite-menu').entryconfig(6, state=states.pop(0))    # flexaid
            self.menuBar.component('NRGsuite-menu').entryconfig(7, state=states.pop(0))    # getcleft
            self.menuBar.component('NRGsuite-menu').entryconfig(9, state=states.pop(0))    # about
            
    else:
        tkMessageBox.showerror('Programs MISSING',
                               'The FlexAID and GetCleft applications ' + 
                               'need to be installed at this location: ' +
                               Install_Dir + '\n\nInstallation CANCELLED...')

else:
    tkMessageBox.showerror('NRGsuite not found',
                           'The NRGsuite applications need to be installed at this location: ' + 
                           Install_Dir + '\n\nInstallation CANCELLED...')
