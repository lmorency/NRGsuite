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
from subprocess import Popen, PIPE

import Tabs
import tkTable
import functools
import threading

if __debug__:
    import General_cmd

class RunVolume(threading.Thread):

    def __init__(self, top, Cleft, Iterations):
        
        threading.Thread.__init__(self)

        self.top = top
        self.GetCleft = self.top.top
        self.Cleft = Cleft
        
        self.cmdline  = '"' + self.GetCleft.VolumeExecutable + '"'
        self.cmdline += ' -i "' + self.Cleft.CleftFile + '"'
        self.cmdline += ' -t ' + Iterations
        #print self.cmdline

        self.start()
        
        
    def run(self):
    
        self.GetCleft.ProcessRunning = True
        
        if self.GetCleft.OSid == 'WIN':
            self.GetCleft.Run = Popen(self.cmdline, shell=False, stdout=PIPE)
        else:
            self.GetCleft.Run = Popen(self.cmdline, shell=True, stdout=PIPE)
            
        self.GetCleft.Run.wait()        
        
        
        if self.GetCleft.Run == None:
            self.top.DisplayMessage("  Thread timed-out for " + self.Cleft.CleftName, 0)
                        
        elif self.GetCleft.Run.returncode != 0:
            self.top.ProcessError = True
            self.top.DisplayMessage("  ERROR: An error occured while executing volume_calc.", 1)
            
        else:
            self.top.DisplayMessage("Thread terminated for " + self.Cleft.CleftName, 0)
            
            Line = self.GetCleft.Run.stdout.readline()
            if Line.startswith('Volume'):
                self.Cleft.Volume = float(Line[8:].strip())

        self.GetCleft.Run = None

        self.GetCleft.ProcessRunning = False
        
        

class EstimateVolume(Tabs.Tab):

    def Def_Vars(self):

        self.Iterations = StringVar()
        self.CleftVolume = DoubleVar()
        
        self.ValidIterations = list()
        self.Validator = list()

    def Init_Vars(self):

        self.ProcessError = False
        self.Iterations.set('3')

        self.ValidIterations = [True, 1, 0, None]
        self.Validator = [ self.ValidIterations ]

    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''  
    def After_Show(self):

        self.top.Default.Update_TempBindingSite()
        
        self.Init_Table()

    ''' ==================================================================================
    FUNCTION Frame_CalcVolume: Display the Calculate Volume Option in the Interface 
    ==================================================================================  '''    
    def Frame(self):
                        
        self.fVolume = Frame(self.top.fMiddle)

        '''===========================================================================
                                    Number of Iterations
        ==========================================================================='''
        fIterations = Frame(self.fVolume, relief=RIDGE, border=0, width=400, height=30)
        fIterations.pack(fill=X, expand=True, side=TOP)
                
        EntryIterations = Entry(fIterations, width=6, background='white', justify=CENTER, font=self.top.font_Text, textvariable=self.Iterations)
        EntryIterations.pack(side=RIGHT, anchor=SE, padx=5)
        args_list = [EntryIterations, self.Iterations, 1, 5, -1, self.ValidIterations,'Iterations','int']
        EntryIterations.bind('<FocusOut>', functools.partial(self.Lost_Focus,args=args_list))
        self.ValidIterations[3] = EntryIterations

        lblIterations = Label(fIterations, text='Iterations:', font=self.top.font_Text)
        lblIterations.pack(side=RIGHT, anchor=SE)
                

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
        self.Btn_Selected.pack(side=LEFT, padx=2)

        self.Btn_Remaining = Button(fButtonsLine2, text='Remaining', font=self.top.font_Text, command=self.Btn_Remaining_Clicked)
        self.Btn_Remaining.pack(side=LEFT, padx=2)

        self.Btn_ALL = Button(fButtonsLine2, text='ALL', font=self.top.font_Text, command=self.Btn_ALL_Clicked)
        self.Btn_ALL.pack(side=LEFT, padx=2)


        fList = Frame(self.fVolume, relief=SUNKEN, border=1, width=375, height=200)
        fList.pack(fill=X, expand=True, side=TOP, pady=20)
        fList.pack_propagate(0)
        
        self.Table = tkTable.Table(fList, 3,
                                   [ 'Cleft Name', 'Estimated', 'Volume' ],
                                   [ 180, 70, 130 ],
                                   [ 1, 6, 13 ],
                                   [ True, True, True ],
                                   self.top.font_Text,
                                   self.top.Color_Blue)

        self.Table.Draw()
        
        self.SelectedCleft = self.Table.Columns['Cleft Name']['StringVar']
        
        
        fButtons2 = Frame(self.fVolume)
        fButtons2.pack(side=TOP, fill=X, padx=5, pady=5)
        fButtons2Line1 = Frame(fButtons2)
        fButtons2Line1.pack(side=TOP, fill=X)

        Button(fButtons2Line1, text='Save volumes', font=self.top.font_Text).pack(side=RIGHT, padx=2)
        Button(fButtons2Line1, text='Refresh clefts', font=self.top.font_Text, command=self.Init_Table).pack(side=RIGHT, padx=2)
        
        return self.fVolume
        
    ''' ==================================================================================
    FUNCTION Validate_Iterations: Validates the value in Entry before starting Grid
    ==================================================================================  '''                 
    def Validate_Iterations(self):
        
        self.fVolume.focus_set()
        self.fVolume.update_idletasks()

        if self.ValidIterations[0]:
            return 0
        
        return 1

    ''' ==================================================================================
    FUNCTION Calculates the volume of the selected cleft only 
    ==================================================================================  '''    
    def Btn_Selected_Clicked(self):

        if not self.Validate_Iterations():

            Selected = self.SelectedCleft.get()

            if Selected != '':
                self.Cleft = self.top.Default.TempBindingSite.Get_CleftName(Selected)
                
                self.VolumeRunning(True)

                if self.Cleft != None:            
                    try:
                        Process = RunVolume(self,self.Cleft, self.Iterations.get())
                        Process.join(30.0)
                        
                        self.Init_Table()
                        
                    except:
                        self.DisplayMessage("The cleft file for '" + self.Cleft + "' no longer exists", 2)
                else:
                    self.DisplayMessage("The cleft object '" + Selected + "' no longer exists", 2)                    
                
            
                self.VolumeRunning(False)
            
    ''' ==================================================================================
    FUNCTION Calculates the volume of remaining clefts (volume=0) in the list
    ==================================================================================  '''    
    def Btn_Remaining_Clicked(self):

        if not self.Validate_Iterations():

            self.VolumeRunning(True)

            for item in self.Table.Columns['Cleft Name']['List'].get(0, END):

                self.Cleft = self.top.Default.TempBindingSite.Get_CleftName(item.lstrip())
                
                if self.Cleft != None:            
                    try:
                        Process = RunVolume(self,self.Cleft, self.Iterations.get())
                        Process.join(30.0)
                        
                        self.Init_Table()

                    except:
                        self.DisplayMessage("The cleft file for '" + self.Cleft.CleftName + "' no longer exists", 2)
                else:
                    self.DisplayMessage("The cleft object '" + item.lstrip() + "' no longer exists", 2)                    
                
            self.VolumeRunning(False)

    ''' ==================================================================================
    FUNCTION Calculates the volume of all clefts in the list
    ==================================================================================  '''    
    def Btn_ALL_Clicked(self):

        if not self.Validate_Iterations():

            self.VolumeRunning(True)
            
            for item in self.Table.Columns['Cleft Name']['List'].get(0, END):

                self.Cleft = self.top.Default.TempBindingSite.Get_CleftName(item.lstrip())
                
                if self.Cleft != None:            
                    try:
                        Process = RunVolume(self,self.Cleft, self.Iterations.get())
                        Process.join(30.0)
                        
                        self.Init_Table()
                        
                    except:
                        self.DisplayMessage("The cleft file for '" + self.Cleft.CleftName + "' no longer exists", 2)
                else:
                    self.DisplayMessage("The cleft object '" + item.lstrip() + "' no longer exists", 2)
                                        
            self.VolumeRunning(False)

    ''' ==================================================================================
    FUNCTION Init_Table: updates the list of cleft in the table
    ==================================================================================  '''    
    def Init_Table(self):

        self.Table.Clear()
        
        for CleftName in self.top.Default.TempBindingSite.Get_SortedCleftNames():
            
            self.Table.Add( [ CleftName, 'False', self.top.Default.TempBindingSite.Get_CleftName(CleftName).Volume ],
                            [ None, None, None ] )

    ''' ==================================================================================
    FUNCTION VolumeRunning: Actives/Deactives controls when a process is running
    ==================================================================================  '''    
    def VolumeRunning(self, boolRun):
                
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
                
