# -*- coding: utf-8 -*-
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
@title: NewProject - Interface

@summary: This is the interface of the New Project application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec 8, 2010
'''

from Tkinter import *

import os
import webbrowser
import Base
import NRGsuite

#=========================================================================================
'''                           ---   PARENT WINDOW  ---                                 '''
#========================================================================================= 

class displayAbout(Base.Base):
    
    ''' ==================================================================================
    FUNCTION Set_Folders: Builds the list of folders
    ==================================================================================  '''  
    def Set_Folders(self):
        
        self.AboutInstall_Dir = os.path.join(self.Install_Dir,'About')
        
    ''' ==================================================================================
    FUNCTION Frame_Main: Generate the Main interface containing ALL the Frames    
    ==================================================================================  '''
    def Frame_Main(self):
        
        #==================================================================================
        '''                  --- TOP FRAME OF THE INTERFACE ---                         '''
        #==================================================================================
        
        fTop = Frame(self.fMain, height=100, border=1) #, bg='red')
        fTop.pack(fill=X, side=TOP,pady=10)
        #fTop.pack_propagate(0)
        
        image_path = os.path.join(self.AboutInstall_Dir,'images','banner-crop-small.gif')
        self.Banner = PhotoImage(file=image_path)        
        btnBanner = Button(fTop, image=self.Banner, command=self.Banner_Clicked)
        btnBanner.pack()
        
        fText = Frame(self.fMain, border=1) #, bg='yellow')
        fText.pack(fill=X,side=TOP, pady=10)
        
        Title = Label(fText, text='NRGsuite', height=3, font=(self.Prefs.FontType, 16))
        Title.pack(side=TOP, anchor=N)
        
        Development = Label(fText, text='Developped by Francis Gaudreault', font=(self.Prefs.FontType,self.Prefs.FontSize))
        Development.pack(side=TOP)
        Supervision = Label(fText, text='Supervised by Dr. Rafael Najmanovich', font=(self.Prefs.FontType,self.Prefs.FontSize))
        Supervision.pack(side=TOP)
        Thanks = Label(fText, text='Special thanks to Eugène Morin and Dominic Duchêne', font=(self.Prefs.FontType,self.Prefs.FontSize))
        Thanks.pack(side=TOP)

        fSep = Frame(fText, height=20, border=1) #, bg='green')
        fSep.pack(side=TOP, fill=X)

        Version = Label(fText, text='Version ' + self.VERSION, height=3, font=(self.Prefs.FontType,self.Prefs.FontSize))
        Version.pack(side=TOP)
        
        Copyright = Label(fText, text='Copyright © 2010-2014 Najmanovich Research Group', font=(self.Prefs.FontType,self.Prefs.FontSize))
        Copyright.pack(side=TOP)

        Rights = Label(fText, text='All rights reserved', font=(self.Prefs.FontType,self.Prefs.FontSize))
        Rights.pack(side=TOP)
        
    ''' ==================================================================================
    FUNCTION After_Quit: Do some tasks after killing a frame
    ==================================================================================  '''
    def After_Quit(self):
        
        NRGsuite.EnableDisableMenu(self.top, [ 'disabled', 'disabled', 'normal', 'normal', 'normal', 'normal' ] )

    ''' ==================================================================================
    FUNCTION Banner_Clicked: Opens up the NRG website
    ==================================================================================  '''
    def Banner_Clicked(self):
        
        webbrowser.open_new_tab(self.WEBSITE)
        
        