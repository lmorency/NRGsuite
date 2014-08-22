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
@title: Prefs.py

@summary: Module that define the preferences of the NRG suite applications.

@organization: Najmanovich Research Group
@creation date:  feb. 18, 2011
@modification date : aug. 19 2104
'''

from Tkinter import *

import os
import cPickle as pickle

import Base
# import NRGsuite


class Prefs(object):

    # FONT SETTINGS
    DefaultFontType = 'helvetica'
    DefaultFontSize = 10

    def __init__(self, FontType = None, FontSize = 0, ToggleAllFlexibleBonds = 0,PreferenceFilePath = None):
        self.FontType = self.DefaultFontType
        self.FontSize = self.DefaultFontSize
        self.ToggleAllFlexibleBonds = 0
        self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')
        self.Load_User_Prefs()

    ''' ==================================================================================
    FUNCTION Load_User_Prefs: Read the user preferences from the Preference file
    =================================================================================  '''
    def Load_User_Prefs(self):
        if os.path.isfile(self.PreferenceFilePath):
            try:
                # loading pickle object directly from user's preference file
                Preferences = pickle.loads(open(self.PreferenceFilePath,'r'))
                self.FontType = Preferences.FontType
                self.FontSize = Preferences.FontSize
                self.ToggleAllFlexibleBonds = Preferences.ToggleAllFlexibleBonds
                self.PreferenceFilePath = Preferences.PreferenceFilePath
            except:
                # Putting back default's values
                print 'exception entered in Prefs.Load_User_Prefs'
                # self.FontType = self.DefaultFontType
                # self.FontSize = self.DefaultFontSize
                # self.ToggleAllFlexibleBonds = 0;
                # self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')
    
    ''' ==================================================================================
    FUNCTION GetFontType: Returns the font type preferred by the user
    =================================================================================  '''
    def GetFontType(self):
        return self.FontType

    ''' ==================================================================================
    FUNCTION GetFontSize: Returns the font size preferred by the user
    =================================================================================  '''
    def GetFontSize(self):
        return self.FontSize

    ''' ==================================================================================
    FUNCTION Write_User_Prefs: Save & Write the user's preferences into Preference file
    =================================================================================  '''
    def Write_User_Prefs(self):
        try:
            prefs_string = pickle.dumps(self)
            file = open(self.PreferenceFilePath,'w')
            file.write(prefs_string)
            file.close()
        except:
            # error code || error message required "unable to write preferences"
            print pickle.dumps(self)

    ''' ==================================================================================
    FUNCTION Write_User_Prefs: Save & Write the user's preferences into Preference file
    =================================================================================  '''
    def Restore_Default_Prefs(self):
        # self.
        pass

    def Update_ToggleAllFlexibleBonds(self):
        if self.ToggleAllFlexibleBonds == 1:
            self.ToggleAllFlexibleBonds = 0
        elif self.ToggleAllFlexibleBonds == 0:
            self.ToggleAllFlexibleBonds = 1


class displayPrefs(Base.Base):

    # ''' ==================================================================================
    # FUNCTION SaveDefault: Save & Write the user's preferences into Preference file
    # =================================================================================  '''
    # def SaveDefault(self):
    #     self.Write_User_Prefs()
    
    ''' ==================================================================================
    FUNCTION Frame_Main: 
    =================================================================================  '''    
    def Frame_Main(self):
        fTop = Frame(self.fMain, height=100, border=1)
        fTop.pack(fill=X, side=TOP, pady=10)

        fText = Frame(self.fMain, border=1)
        fText.pack(fill=X, side=TOP, pady=10)

        Title = Label(fText, text='NRGsuite Preferences Panel', height=3, font=self.font_Title)
        Title.pack(side=TOP, anchor=N)

        FontType = Label(fText, text='Preferred Font Type')
        FontType.pack(side=TOP)
        FontSize = Label(fText, text='Preferred Font Size')
        FontSize.pack(side=TOP)

        fSep = Frame(fText, height=20, border=1) #, bg='green')
        fSep.pack(side=TOP, fill=X)

        ToggleAllFlexibleBonds_Var = IntVar()
        if self.Prefs.ToggleAllFlexibleBonds:
            ToggleAllFlexibleBonds_Var.set(1)

        ToggleAllFlexibleBonds = Checkbutton(fText, text='Automatically consider all rotable bonds of the ligand as flexible during the simulation',variable=self.Prefs.ToggleAllFlexibleBonds)#,command=self.Prefs.Update_ToggleAllFlexibleBonds())
        ToggleAllFlexibleBonds.pack(side=TOP)

        fButtons = Frame(fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=50)

        Btn_Save = Button(fButtons, text='Save Preferences', width=12, command=self.Btn_Save_Clicked, font=self.font_Text)
        Btn_Save.pack(side=RIGHT,anchor=SE,pady=3)

        Btn_Restore_Default = Button(fButtons, text='Restore Defaults', width=12, command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Restore_Default.pack(side=RIGHT, anchor=SE,pady=3)

        Btn_Cancel = Button(fButtons, text='Cancel', width=12, command=self.Btn_Cancel_Clicked, font=self.font_Text)
        Btn_Cancel.pack(side=RIGHT, anchor=SE, pady=3)

        fButtons.pack(side=RIGHT, fill=X, expand=True)
        fButtons.pack_propagate(0)

    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        # NRGsuite.EnableDisableMenu(self.top, [ 'disabled', 'disabled', 'normal', 'normal', 'normal', 'normal' ])
        return
    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Save_Clicked(self):
        self.Prefs.Write_User_Prefs()
        return

    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Default_Clicked(self):
        self.Prefs.Restore_Default_Prefs()

    ''' ==================================================================================
    FUNCTION Btn_Cancel_Clicked: Cancel the project creation then quit the application.
    ==================================================================================  '''        
    def Btn_Cancel_Clicked(self):
        self.Quit()


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
        
        # Prepend text at end of control text
        self.TextMessage.mark_set(INSERT, END)

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
        