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

from Tkinter import *

import tkTable
import functools
import Grid

class EstimateVolume:

    def __init__(self, top):

        self.top = top
        self.FrameName = 'Estimate Volume'
        self.Tab = self.top.Btn_Volume

        self.DisplayMessage = self.top.DisplayMessage

        self.TempBindingSite = self.top.Default.TempBindingSite

        self.Def_Vars()
        self.Init_Vars()

        self.Frame()
        self.Trace()

    def Def_Vars(self):

        self.Precision = StringVar()
        self.CleftVolume = DoubleVar()
        
        self.ValidPrecision = list()
        self.Validator = list()

    def Init_Vars(self):

        self.ProcessError = False
        self.Precision.set('1.0')

        self.ValidPrecision = [True, 1, 0, None]
        self.Validator = [  ]

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
                
        self.fVolume.pack_forget()
        #self.fVolume.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        return

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        return

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
        
        self.fVolume.pack(fill=BOTH, expand=True)

        #self.LoadMessage()
        
        self.Init_Table()

        self.top.Default.Update_TempBindingSite()

    ''' ==================================================================================
    FUNCTION Frame_CalcVolume: Display the Calculate Volume Option in the Interface 
    ==================================================================================  '''    
    def Frame(self):
                        
        self.fVolume = Frame(self.top.fMiddle)

        #==================================================================================
        '''                           --- BUTTONS AREA ---                              '''
        #==================================================================================                
        
        fButtons = Frame(self.fVolume)
        fButtons.pack(side=TOP, fill=X, padx=5, pady=5)
        fButtonsLine1 = Frame(fButtons)
        fButtonsLine1.pack(side=TOP, fill=X)
        fButtonsLine2 = Frame(fButtons)
        fButtonsLine2.pack(side=TOP, fill=X)

        Label(fButtonsLine1, text='Volume estimator', font=self.top.font_Title).pack(side=LEFT)

        self.Btn_Selected = Button(fButtonsLine2, text='Selected', font=self.top.font_Text,command=self.Btn_Selected_Clicked)
        self.Btn_Selected.pack(side=LEFT)

        self.Btn_Remaining = Button(fButtonsLine2, text='Remaining', font=self.top.font_Text, command=self.Btn_Remaining_Clicked)
        self.Btn_Remaining.pack(side=LEFT)

        self.Btn_ALL = Button(fButtonsLine2, text='ALL', font=self.top.font_Text, command=self.Btn_ALL_Clicked)
        self.Btn_ALL.pack(side=LEFT)

        #fTitle = Frame(self.fVolume, relief=RIDGE, border=0, width=400, height=30)
        #fTitle.pack(fill=X, expand=True, side=TOP)

        #fProbe = Frame(self.fVolume, relief=RIDGE, border=0, width=400, height=30)
        #fProbe.pack(fill=X, expand=True, side=TOP)

        fList = Frame(self.fVolume, relief=SUNKEN, border=1, width=400, height=300)
        fList.pack(fill=X, expand=True, side=TOP, pady=20)
                
        #EntryPrecision = Entry(fProbe, width=6, background='white', justify=CENTER, font=self.top.font_Text, textvariable=self.Precision)
        #EntryPrecision.pack(side=RIGHT, anchor=SE, padx=5)
        #args_list = [EntryPrecision, self.Precision, 0.25, 2.00, 2, self.ValidPrecision,'Precision','float']
        #EntryPrecision.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        #self.ValidPrecision[3] = EntryPrecision


        #lblProbeSpace = Label(fProbe, text='Precision:', font=self.top.font_Text)
        #lblProbeSpace.pack(side=RIGHT, anchor=SE)
                

        self.Table = tkTable.Table(fList, 3,
                                   [ 'Cleft Name', 'Estimated', 'Volume' ],
                                   [ 180, 70, 130 ],
                                   [ 1, 4, 10 ],
                                   [ True, True, True ],
                                   self.top.font_Text,
                                   self.top.Color_Blue)

        self.Table.Draw()
        
        self.Selected = self.Table.Columns['Cleft Name']['StringVar']
        
        # Quit the advanced options menu
        Btn_Back = Button(self.fVolume, text='Back', font=self.top.font_Text, command=self.Btn_Back_Clicked)
        Btn_Back.pack(side=BOTTOM, anchor=E)

    ''' ==================================================================================
    FUNCTION Validate_Precision: Validates the value in Entry before starting Grid
    ==================================================================================  '''                 
    def Validate_Precision(self):
        
        self.fVolume.focus_set()
        self.fVolume.update_idletasks()

        if self.ValidPrecision[0]:
            return 0
        
        return 1

    ''' ==================================================================================
    FUNCTION Btn_Back_Clicked: Goes back to default frame
    ==================================================================================  '''                 
    def Btn_Back_Clicked(self):
        
        #self.fVolume.focus_set()
        #self.fVolume.update_idletasks()
        
        self.top.SetActiveFrame(self.top.Default)


    ''' ==================================================================================
    FUNCTION Calculates the volume of the selected cleft only 
    ==================================================================================  '''    
    def Btn_Selected_Clicked(self):

        if not self.Validate_Precision():

            Selected = self.Selected.get()

            if Selected != '':
                self.Cleft = Selected

                try:
                    CleftFile = self.dictTempClefts[self.Cleft][2]

                    self.top.ProcessRunning = True
                    self.GenGridRunning(True)

                    Process = Grid.Grid(self, CleftFile, '', self.Precision.get(), True)

                except:
                    self.DisplayMessage("The cleft '" + self.Cleft + "' no longer exists", 2)
                    return
            
    ''' ==================================================================================
    FUNCTION Calculates the volume of all clefts in the list
    ==================================================================================  '''    
    def Btn_Remaining_Clicked(self):

        if not self.Validate_Precision():

            for item in self.Table.Columns['Cleft Name']['List'].get(0, END):

                self.Cleft = item.lstrip()

                # skip already Estimated
                if self.dictTempClefts[self.Cleft][6]:
                    continue

                #if not already Estimated
                try:
                    CleftFile = self.dictTempClefts[self.Cleft][2]

                    self.top.ProcessRunning = True
                    self.GenGridRunning(True)

                    self.DisplayMessage("Calculating volume of " + self.Cleft + "...", 0)
                    Process = Grid.Grid(self, CleftFile, '', self.Precision.get(), True)
                    Process.join()

                except:
                    self.DisplayMessage("The cleft '" + self.Cleft + "' no longer exists", 2)
                    continue

    ''' ==================================================================================
    FUNCTION Calculates the volume of all clefts in the list
    ==================================================================================  '''    
    def Btn_ALL_Clicked(self):

        if not self.Validate_Precision():

            for item in self.Table.Columns['Cleft Name']['List'].get(0, END):

                self.Cleft = item.lstrip()

                #if not already Estimated
                try:
                    CleftFile = self.dictTempClefts[self.Cleft][2]

                    self.top.ProcessRunning = True
                    self.GenGridRunning(True)

                    self.DisplayMessage("Calculating volume of " + self.Cleft + "...", 0)
                    Process = Grid.Grid(self, CleftFile, '', self.Precision.get(), True)
                    Process.join()

                except:
                    self.DisplayMessage("The cleft '" + self.Cleft + "' no longer exists", 2)
                    continue

    ''' ==================================================================================
    FUNCTION Init_Table: updates the list of cleft in the table
    ==================================================================================  '''    
    def Init_Table(self):

        for Cleft in iter(self.TempBindingSite.listClefts):
            
            self.Table.Add( [ Cleft.CleftName, 'False', Cleft.Volume ],
                            [ None, None, None ] )

    ''' ==================================================================================
    FUNCTION GenGridRunning: Actives/Deactives controls when a process is running
    ==================================================================================  '''    
    def GenGridRunning(self, boolRun):

        if boolRun:
            self.Btn_Selected.config(state='disabled')
            self.Btn_Remaining.config(state='disabled')
            self.Btn_ALL.config(state='disabled')

            self.Btn_Selected.update_idletasks()
            self.Btn_Remaining.update_idletasks()
            self.Btn_ALL.update_idletasks()

        else:
            self.Btn_Selected.config(state='normal')
            self.Btn_Remaining.config(state='normal')
            self.Btn_ALL.config(state='normal')

            # Catch errors
            if not self.ProcessError:

                return
                
            #self.dictTempClefts[self.Cleft][6] = True
            #    self.dictTempClefts[self.Cleft][7] = '%.3f' % self.CleftVolume.get()

            #    self.Table.Set(self.Cleft, 'Cleft Name', True, 'Estimated')
            #    self.Table.Set(self.Cleft, 'Cleft Name', 
            #                   self.dictTempClefts[self.Cleft][7], 'Volume')
                            
            self.top.ProcessRunning = False

