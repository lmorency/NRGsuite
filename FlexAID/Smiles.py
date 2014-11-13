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
@title: FlexAID - Smiles.py

@summary: Class that handle the input of a smiles string.
          This class inherits from IOFile

@organization: Najmanovich Research Group
@creation date:  Feb. 9, 2013
'''

from Tkinter import *

import General

class Smiles(object):
    
    def __init__(self, top, SmilesString, SmilesName):
    
        self.top = top
        
        self.inputWindow = Toplevel()
        self.inputWindow.title = 'SMILES'
        self.inputWindow.protocol('WM_DELETE_WINDOW', self.Quit)
        
        WINDOWWIDTH = 400
        WINDOWHEIGHT = 130

        self.inputWindow.maxsize(WINDOWWIDTH,WINDOWHEIGHT)
        self.inputWindow.minsize(WINDOWWIDTH,WINDOWHEIGHT)
        
        General.CenterWindow(self.inputWindow,WINDOWWIDTH,WINDOWHEIGHT)

        self.SmilesString = SmilesString
        self.SmilesName = SmilesName
        
        self.Frame()
    
    def Frame(self):
    
        Label(self.inputWindow, text='Input your SMILES string:', font=self.top.font_Title).pack(side=TOP, pady=5, anchor=W, padx=10)
        Entry(self.inputWindow, font=self.top.font_Text, textvariable=self.SmilesString, width=100).pack(side=TOP, fill=X, padx=10)

        Label(self.inputWindow, text='Name the object:', font=self.top.font_Text).pack(side=TOP, pady=5, anchor=W, padx=10)
        Entry(self.inputWindow, font=self.top.font_Text, textvariable=self.SmilesName, width=30).pack(side=TOP, padx=10, anchor=W)

        Button(self.inputWindow, font=self.top.font_Text, text='Cancel', command=self.Quit).pack(side=RIGHT, padx=10)
        Button(self.inputWindow, font=self.top.font_Text, text='Enter', command=self.Enter).pack(side=RIGHT)
        
    def Quit(self):
        
        self.inputWindow.destroy()
        self.top.SmilesRunning(False, False)
        
    def Enter(self):

        if General.validate_String(self.SmilesName.get(), '', False, False, True):
            self.top.DisplayMessage("  ERROR: Could not process the smiles string because you entered an invalid name.", 2)
            return
            
        #Test : CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO
        self.inputWindow.destroy()
        self.top.SmilesRunning(False, True)
    