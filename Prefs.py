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
    ********************************************************************************* """ 
class Prefs(object):

    # FONT SETTINGS
    DefaultFontType = 'Helvetica'
    DefaultFontSize = 12

    def __init__(self, FontType = None, FontSize = 0, ToggleAllFlexibleBonds = 0, PreferenceFilePath = None, AlwaysShowAdvancedView = 0):

        self.FontType = self.DefaultFontType
        self.FontSize = self.DefaultFontSize
        self.ToggleAllFlexibleBonds = ToggleAllFlexibleBonds
        self.AlwaysShowAdvancedView = AlwaysShowAdvancedView
        self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')
        self.Load_User_Prefs()

    ''' ==================================================================================
    FUNCTION Load_User_Prefs: Read the user preferences from the Preference file
    =================================================================================  '''
    def Load_User_Prefs(self):

        if os.path.isfile(self.PreferenceFilePath):
            try:
                # loading pickle object directly from user's preference file
                f = open(self.PreferenceFilePath, "rb")
                Preferences = pickle.load(f)
                f.close()
                self.FontType = Preferences.FontType
                self.FontSize = Preferences.FontSize

                if Preferences.ToggleAllFlexibleBonds == 1:
                    self.ToggleAllFlexibleBonds = 1
                
                if Preferences.AlwaysShowAdvancedView == 1:
                    self.AlwaysShowAdvancedView = 1

                if os.path.isfile(Preferences.PreferenceFilePath):
                    self.PreferenceFilePath = Preferences.PreferenceFilePath
                else:
                    self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')


            except Exception, e:
                # Catch exceptions
                shutil.rmtree(self.PreferenceFilePath)
                self.Load_User_Prefs()
                print 'Exception caught in Prefs.Load_User_Prefs(). Old Preferences file removed.'
                print 'Reloading User\'s preferences.'

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
            f = open(self.PreferenceFilePath,"wb")
            pickle.dump(self,f)
            f.close()
        except Exception, e:
            # error code || error message required "unable to write preferences"
            print 'exception caught in Prefs.Write_User_Prefs()'
            pass

    ''' ==================================================================================
    FUNCTION Write_User_Prefs: Save & Write the user's preferences into Preference file
    =================================================================================  '''
    def Restore_Default_Prefs(self):
        self.FontType = self.DefaultFontType
        self.FontSize = self.DefaultFontSize
        self.ToggleAllFlexibleBonds = 0
        self.AlwaysShowAdvancedView = 0
        self.PreferenceFilePath = os.path.join(os.path.expanduser('~'),'Documents','NRGsuite','.NRGprefs')
        self.Write_User_Prefs()


''' 
==================================================================================
CLASS displayPrefs :     
==================================================================================
'''
class displayPrefs(Base.Base):

    def Def_Vars(self):
        # IntVar() used for the ToggleAllFlexibleBonds and AlwaysShowAdvacnedView Checkbuttons()
        self.ToggleAllFlexibleBonds_Var = IntVar(0)
        self.AlwaysShowAdvancedView_Var = IntVar(0)
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
        if self.Prefs.AlwaysShowAdvancedView == 1:
            self.AlwaysShowAdvancedView_Var.set(1)
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

        if self.Prefs.ToggleAllFlexibleBonds == 1:
            self.ToggleAllFlexibleBonds_Var.set(0)
            self.Prefs.ToggleAllFlexibleBonds = 0
            return 0

        elif self.Prefs.ToggleAllFlexibleBonds == 0:
            self.ToggleAllFlexibleBonds_Var.set(1)
            self.Prefs.ToggleAllFlexibleBonds = 1
            return 1

    ''' ====================================================================================================
    FUNCTION Update_AlwaysShowAdvancedView: Update the Prefs class with current ToggleAllFlexibleBonds value
    ========================================================================================================  '''    
    def Update_AlwaysShowAdvancedView(self):

        if self.Prefs.AlwaysShowAdvancedView == 1:
            self.AlwaysShowAdvancedView_Var.set(0)
            self.Prefs.AlwaysShowAdvancedView = 0
            return 0

        elif self.Prefs.AlwaysShowAdvancedView == 0:
            self.AlwaysShowAdvancedView_Var.set(1)
            self.Prefs.AlwaysShowAdvancedView = 1
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
        try:
            self.Prefs.Write_User_Prefs()
        except:
            print "Exception caught in Prefs.SaveDefault()"
            pass
    
    ''' ==================================================================================
    FUNCTION Frame_Main: 
    =================================================================================  '''    
    def Frame_Main(self):

        # Load User Preferences
        self.Prefs.Load_User_Prefs()

        fTop = Frame(self.fMain)#,bg='orange')

        ######################## Font Options

        fText = Frame(self.fMain)#, bg='red')
        fText.pack(side=TOP, fill=BOTH, expand=True, padx=5, pady=2)

        Title_Font = Label(fText, text='Font Options', font=self.font_Title_H)
        Title_Font.pack(side=TOP, anchor=W, padx=5, pady=2)
        Title_Font.pack_propagate(0)

        fFont_options = Frame(fText)#,bg='green')
        fFont_options.pack(side=TOP,fill=BOTH,padx=5,pady=0)
        fFont_options2 = Frame(fText)#,bg='purple')
        fFont_options2.pack(side=TOP,fill=BOTH,padx=5,pady=0)

        fFontType_Label = Label(fFont_options, text='Preferred Font Type : ', font=self.font_Text)
        fFontType_Label.pack(side=LEFT,anchor=W)
        fFontType_Label.pack_propagate(0)

        fontypes = ["Arial", "Lucida", "Helvetica", "Tahoma", "Times", "System", "Fixed", "Courrier"]
        fontypes.sort()
        fFontType_OptionMenu = OptionMenu(fFont_options, self.FontType_StringVar,command=self.Update_FontType, *fontypes)
        fFontType_OptionMenu.configure(font=self.font_Text)#, bg='white')#,width=12)
        fFontType_OptionMenu['menu'].config(font=self.font_Text)#, bg='white')
        fFontType_OptionMenu.pack(side=RIGHT,anchor=E)#,fill=BOTH,expand=True)
        fFontType_OptionMenu.pack_propagate(0)
        
        fFontSize_Label = Label(fFont_options2, text='Preferred Font Size : ', font=self.font_Text)
        fFontSize_Label.pack(side=LEFT,anchor=W)

        fontsizes = [8,9,10,11,12,13,14]
        fFontSize_OptionMenu = OptionMenu(fFont_options2, self.FontSize_IntVar,command=self.Update_FontSize, *fontsizes)
        fFontSize_OptionMenu.configure(font=self.font_Text)#, bg='white')#, width=8)
        fFontSize_OptionMenu['menu'].config(font=self.font_Text)#, bg='white')
        fFontSize_OptionMenu.pack(side=RIGHT,anchor=E)#,fill=BOTH,expand=True)
        fFontSize_OptionMenu.pack_propagate(0)

        ######################## FlexAID Options

        fOptions = Frame(self.fMain)#, bg='yellow')

        Title_FlexAID_Options = Label(fOptions, text='FlexAID Options', font=self.font_Title_H)
        Title_FlexAID_Options.pack(side=TOP, anchor=W, padx=5, pady=2)
        Title_FlexAID_Options.pack_propagate(0)

        ToggleAllFlexibleBonds = Checkbutton(fOptions, variable=self.ToggleAllFlexibleBonds_Var, command=self.Update_ToggleAllFlexibleBonds, text='Automatically consider all rotable bonds of the ligand as flexible', font=self.font_Text)
        ToggleAllFlexibleBonds.pack(side=TOP,anchor=W,padx=5, pady=2)
        ToggleAllFlexibleBonds.pack_propagate(0)

        AlwaysShowAdvancedView = Checkbutton(fOptions, variable=self.AlwaysShowAdvancedView_Var, command=self.Update_AlwaysShowAdvancedView, text='Always show Advanced View in FlexAID', font=self.font_Text)
        AlwaysShowAdvancedView.pack(side=TOP,anchor=W,padx=5, pady=2)
        AlwaysShowAdvancedView.pack_propagate(0)
        fOptions.pack(side=TOP, fill=BOTH, padx=5, pady=2)

        fButtons = Frame(fOptions, relief=RIDGE, border=0, width=self.WINDOWWIDTH)#, bg='blue')

        Btn_Save = Button(fButtons, text='Save Preferences', width=12, command=self.Btn_Save_Clicked, font=self.font_Text)
        Btn_Save.pack(side=RIGHT,anchor=S, expand=True, fill=BOTH)

        Btn_Restore_Default = Button(fButtons, text='Restore Defaults', width=12, command=self.Btn_Default_Clicked, font=self.font_Text)
        Btn_Restore_Default.pack(side=RIGHT, anchor=S, expand=True, fill=BOTH)

        Btn_Cancel = Button(fButtons, text='Cancel', width=12, command=self.Btn_Cancel_Clicked, font=self.font_Text)
        Btn_Cancel.pack(side=RIGHT, anchor=S,  expand=True, fill=BOTH)

        fButtons.pack(side=RIGHT, fill=X, expand=True,pady=5, padx=2)

        fTop.pack(side=TOP, fill=BOTH, expand=True)
        fTop.pack_propagate(0)

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
        
