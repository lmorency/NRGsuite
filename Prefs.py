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


""" *********************************************************************************
    CLASS Prefs: Main Class used to define NRGsuite preferences 
    *********************************************************************************  """ 
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
                with open(self.PreferenceFilePath, "rb") as f:
                    Preferences = pickle.load(f)

                self.FontType = Preferences.FontType
                self.FontSize = Preferences.FontSize
                if Preferences.ToggleAllFlexibleBonds == 1:
                    self.ToggleAllFlexibleBonds = 1
                if os.path.isfile(Preferences.PreferenceFilePath):
                    self.PreferenceFilePath = Preferences.PreferenceFilePath
                else:
                    self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')

            except Exception, e:
                # Catch exceptions 
                print 'exception entered in Prefs.Load_User_Prefs'

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
            with open(self.PreferenceFilePath,"wb") as f:
                pickle.dump(self,f)
        except Exception, e:
            # error code || error message required "unable to write preferences"
            print 'exception in Prefs.Write_User_Prefs'
            print self.PreferenceFilePath

    ''' ==================================================================================
    FUNCTION Write_User_Prefs: Save & Write the user's preferences into Preference file
    =================================================================================  '''
    def Restore_Default_Prefs(self):
        self.FontType = self.DefaultFontType;
        self.FontSize = self.DefaultFontSize;
        self.ToggleAllFlexibleBonds = 0
        self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')
        self.Write_User_Prefs()



class displayPrefs(Base.Base):

    def Def_Vars(self):
        # IntVar() used for the ToggleAllFlexibleBonds Checkbutoton()
        self.ToggleAllFlexibleBonds_Var = IntVar(0)
        # IntVar() used for the FontSize OptionMenu()
        self.FontSize_IntVar = IntVar()
        self.FontSize_IntVar.set(self.Prefs.DefaultFontSize)
        # StringVar() used for the FontType OptionMenu()
        self.FontType_StringVar = StringVar()
        self.FontType_StringVar.set(self.Prefs.DefaultFontType)

    def Init_Vars(self):
        # ToggleAllFlexibleBonds preferred value set
        if self.Prefs.ToggleAllFlexibleBonds == 1:
            self.ToggleAllFlexibleBonds_Var.set(1)
        # FontType preferred value set
        if self.Prefs.FontType != self.FontType_StringVar.get():
            self.FontType_StringVar.set(self.Prefs.FontType)
        # FontSize preferred value set
        if self.Prefs.FontSize != self.FontSize_IntVar.get():
            self.FontSize_IntVar.set(self.Prefs.FontSize)

    ''' ====================================================================================================
    FUNCTION Update_ToggleAllFlexibleBonds: Update the Prefs class with current ToggleAllFlexibleBonds value
    ========================================================================================================  '''    
    def Update_ToggleAllFlexibleBonds(self):
        # if self.ToggleAllFlexibleBonds_Var.get() == 1 and self.Prefs.ToggleAllFlexibleBonds == 1:
        if self.Prefs.ToggleAllFlexibleBonds == 1:
            self.ToggleAllFlexibleBonds_Var.set(0)
            self.Prefs.ToggleAllFlexibleBonds = 0
            return 0
        # elif self.ToggleAllFlexibleBonds_Var.get() == 0 and self.Prefs.ToggleAllFlexibleBonds == 0:
        elif self.Prefs.ToggleAllFlexibleBonds == 0:
            self.ToggleAllFlexibleBonds_Var.set(1)
            self.Prefs.ToggleAllFlexibleBonds = 1
            return 1

    ''' ====================================================================================================
    FUNCTION Update_FontSize: Update the Prefs class with current FontSize value
    ========================================================================================================  '''    
    def Update_FontSize(self, val):
        self.Prefs.FontSize = val
        self.FontSize_IntVar.set(val)

    ''' ====================================================================================================
    FUNCTION Update_FontTYpe: Update the Prefs class with current FontType value
    ========================================================================================================  '''    
    def Update_FontType(self, val):
        self.Prefs.FontType = val
        self.FontType_StringVar.set(val)

    ''' ==================================================================================
    FUNCTION SaveDefault: Save & Write the user's preferences into Preference file
    =================================================================================  '''
    def SaveDefault(self):
        self.Prefs.Write_User_Prefs()
    
    ''' ==================================================================================
    FUNCTION Frame_Main: 
    =================================================================================  '''    
    def Frame_Main(self):
        self.Prefs.Load_User_Prefs()
        fTop = Frame(self.fMain, height=100, border=1)#, bg='blue')

        fText = Frame(self.fMain, border=1)#, bg='red')
        fText.pack(fill=X, side=TOP, pady=10)

        Title = Label(fText, text='NRGsuite Preferences Panel', height=3, font=(self.Prefs.FontType,18))
        Title.pack(side=TOP, anchor=N)
        Title.pack_propagate(0)
        # fSep = Frame(fText, height=20, border=1) #, bg='green')
        # fSep.pack(side=TOP, fill=X)
        # fSep.pack_propagate(0)

        # FonType OptionMenu widget
        fFont_options = Frame(fText)
        fFont_options.pack(side=TOP,fill=X,padx=5,pady=5)
        
        fFontType_Label = Label(fFont_options, text='Preferred Font Type : ', font=(self.Prefs.FontType,16))
        fFontType_Label.pack(side=LEFT)
        # fFontType_Label.pack_propagate(0)

        fontypes = ["arial", "tahoma", "helvetica", "times"]
        fFontType_OptionMenu = OptionMenu(fFont_options, self.FontType_StringVar,*fontypes,command=self.Update_FontType)
        fFontType_OptionMenu.pack(side=LEFT)
        # fFontType_OptionMenu.pack_propagate(0)

        fFontSize_Label = Label(fFont_options, text=' : Preferred Font Size', font=(self.Prefs.FontType,16))
        fFontSize_Label.pack(side=RIGHT)

        fontsizes = [10,11,12,13,14]
        fFontSize_OptionMenu = OptionMenu(fFont_options, self.FontSize_IntVar,*fontsizes,command=self.Update_FontSize)
        fFontSize_OptionMenu.pack(side=RIGHT)
        

        ToggleAllFlexibleBonds = Checkbutton(fTop,variable=self.ToggleAllFlexibleBonds_Var,command=self.Update_ToggleAllFlexibleBonds,text='Automatically consider all rotable bonds of the ligand as flexible during the simulation')
        ToggleAllFlexibleBonds.pack(side=TOP)

        fButtons = Frame(fTop, relief=RIDGE, border=0, width=self.WINDOWWIDTH, height=50)#, bg='green')

        Btn_Save = Button(fButtons, text='Save Preferences', width=12, command=self.Btn_Save_Clicked, font=self.font_Text)
        Btn_Save.pack(side=RIGHT,anchor=SE,pady=3)

        Btn_Restore_Default = Button(fButtons, text='Restore Defaults', width=12, command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Restore_Default.pack(side=RIGHT, anchor=SE,pady=3)

        Btn_Cancel = Button(fButtons, text='Cancel', width=12, command=self.Btn_Cancel_Clicked, font=self.font_Text)
        Btn_Cancel.pack(side=RIGHT, anchor=SE, pady=3)

        fButtons.pack(side=RIGHT, fill=X, expand=True)
        fButtons.pack_propagate(0)
        fTop.pack(fill=X, side=TOP, pady=10)

    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        return
    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Save_Clicked(self):
        self.SaveDefault()
        self.Quit()

    ''' ==================================================================================
    FUNCTION Btn_Default_Clicked: Sets back the default config
    ================================================================================== '''    
    def Btn_Default_Clicked(self):
        self.Prefs.Restore_Default_Prefs()
        self.Quit()

    ''' ==================================================================================
    FUNCTION Btn_Cancel_Clicked: Cancel the project creation then quit the application.
    ==================================================================================  '''        
    def Btn_Cancel_Clicked(self):
        self.Quit()

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
        