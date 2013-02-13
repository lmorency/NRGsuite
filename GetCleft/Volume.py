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
import tkMessageBox
import MultiList
import threading

if __debug__:
    import General_cmd

class RunVolume(threading.Thread):

    def __init__(self, top, Clefts, Iterations):
        
        threading.Thread.__init__(self)

        self.top = top
        self.GetCleft = self.top.top
        self.Clefts = Clefts
        self.Iterations = Iterations
        
        self.start()
    
    def run(self):        
    
        self.GetCleft.ProcessRunning = True

        for Cleft in self.Clefts:
        
            cmdline  = '"' + self.GetCleft.VolumeExecutable + '"'
            cmdline += ' -i "' + Cleft.CleftFile + '"'
            cmdline += ' -t ' + self.Iterations
            print(cmdline)
        
            if self.GetCleft.OSid == 'WIN':
                self.GetCleft.Run = Popen(cmdline, shell=False, stdout=PIPE)
            else:
                self.GetCleft.Run = Popen(cmdline, shell=True, stdout=PIPE)
                
            # Wait for task to end in thread
            (out,err) = self.GetCleft.Run.communicate()
            
            if self.GetCleft.Run.returncode != 0:
                self.top.ProcessError = True
                self.top.DisplayMessage("  ERROR: An error occured while executing volume_calc.", 1)
                
            else:
                #self.top.DisplayMessage("Thread terminated for " + self.Cleft.CleftName, 0)
                
                Lines = out.splitlines()
                for Line in Lines:
                    if Line.startswith('Volume'):
                        Cleft.Volume = float(Line[8:].strip())
                        self.top.Init_Table()
                        break
        
        self.GetCleft.Run = None
        self.GetCleft.ProcessRunning = False
        
        self.top.VolumeRunning(False)
    

class EstimateVolume(Tabs.Tab):

    def Def_Vars(self):

        self.Iterations = StringVar()
        self.CleftVolume = DoubleVar()
        
    def Init_Vars(self):

        self.ProcessError = False
        self.Iterations.set('3')
        
    def Trace(self):

        try:
            self.IterationsTrace = self.Iterations.trace('w', lambda *args, **kwargs:
                                                              self.Validate_Field(input=self.EntryIterations, var=self.Iterations, min=1,
                                                              max=5, ndec=-1, tag='Number of iterations', _type=int))
        except:
            pass
            
    def Del_Trace(self):
    
        try:
            self.Iterations.trace_vdelete('w', self.IterationsTrace)
        except:
            pass
    
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
                
        self.EntryIterations = Entry(fIterations, width=6, background='white', justify=CENTER, font=self.top.font_Text, textvariable=self.Iterations)
        self.EntryIterations.pack(side=RIGHT, anchor=SE, padx=5)
        self.ValidIterations = [1, False, self.EntryIterations]

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
        
        self.Table = MultiList.Table(fList, 3,
                                   [ 'Color', 'Cleft object', 'Volume' ],
                                   [ 40, 200, 130 ],
                                   [ 0, 6, 6 ],
                                   [ False, True, True ],
                                   self.top.font_Text,
                                   self.top.Color_Blue)
        
        self.Table.Draw()
        
        self.SelectedCleft = self.Table.Columns['Cleft object']['StringVar']
        
        fButtons2 = Frame(self.fVolume)
        fButtons2.pack(side=TOP, fill=X, padx=5, pady=5)
        fButtons2Line1 = Frame(fButtons2)
        fButtons2Line1.pack(side=TOP, fill=X)

        Button(fButtons2Line1, text='Save volumes', font=self.top.font_Text, command=self.Save_Volume).pack(side=RIGHT, padx=2)
        Button(fButtons2Line1, text='Refresh clefts', font=self.top.font_Text, command=self.Init_Table).pack(side=RIGHT, padx=2)

        self.Validator = [ self.ValidIterations ]
        
        return self.fVolume
        
    ''' ==================================================================================
    FUNCTION Calculates the volume of the selected cleft only 
    ==================================================================================  '''    
    def Btn_Selected_Clicked(self):

        if not self.Validate_Fields():
        
            Selected = self.SelectedCleft.get().lstrip()
            if Selected != '':
            
                Cleft = self.top.Default.TempBindingSite.Get_CleftName(Selected)
                if Cleft != None:
                    try:
                        self.VolumeRunning(True)
                        
                        Process = RunVolume(self, [ Cleft ], self.Iterations.get())
                    except:
                        self.DisplayMessage("The cleft file for '" + self.Cleft.CleftFile + "' no longer exists", 2)
                else:
                    self.DisplayMessage("The cleft object '" + Selected + "' no longer exists", 2)                    
    
    ''' ==================================================================================
    FUNCTION Calculates the volume of remaining clefts (volume=0) in the list
    ==================================================================================  '''    
    def Btn_Remaining_Clicked(self):

        if not self.Validate_Fields():
        
            Clefts = []
            for item in self.Table.Columns['Cleft object']['List'].get(0, END):

                Cleft = self.top.Default.TempBindingSite.Get_CleftName(item.lstrip())
                if Cleft != None and Cleft.Volume == 0.0:
                    Clefts.append(Cleft)
            
            try:
                self.VolumeRunning(True)
                Process = RunVolume(self, Clefts, self.Iterations.get())
            except:
                pass
    
    ''' ==================================================================================
    FUNCTION Calculates the volume of all clefts in the list
    ==================================================================================  '''    
    def Btn_ALL_Clicked(self):

        if not self.Validate_Fields():
                    
            Clefts = []
            for item in self.Table.Columns['Cleft object']['List'].get(0, END):

                Cleft = self.top.Default.TempBindingSite.Get_CleftName(item.lstrip())
                if Cleft != None:
                    Clefts.append(Cleft)                

            try:
                self.VolumeRunning(True)
                Process = RunVolume(self, Clefts, self.Iterations.get())
            except:
                pass
    
    ''' ==================================================================================
    FUNCTION Init_Table: updates the list of cleft in the table
    ==================================================================================  '''    
    def Init_Table(self):
        
        self.Table.Clear()
        
        for CleftName in self.top.Default.TempBindingSite.Get_SortedCleftNames():
            self.Table.Add( [ '', CleftName, str(self.top.Default.TempBindingSite.Get_CleftName(CleftName).Volume) ],
                            [ self.top.Default.TempBindingSite.Get_CleftName(CleftName).Color, None, None ] )

    ''' ==================================================================================
    FUNCTION VolumeRunning: Actives/Deactives controls when a process is running
    ==================================================================================  '''    
    def VolumeRunning(self, boolRun):
                
        if boolRun:
            self.Disable_Frame()
        else:
            self.Enable_Frame()
            
            self.Init_Table()

    ''' ==================================================================================
    FUNCTION SaveVolumes: Saves the calculated volume of the clefts
    ==================================================================================  '''    
    def Save_Volume(self):
        
        if not self.top.CopySession:
            if tkMessageBox.askquestion("Question", message="The cleft(s) first need to be saved before you can save the volumes. Do you still want to proceed?",
                                        icon='warning') == 'yes':
                self.top.Default.Btn_Save_Clefts()
        else:
            self.top.Default.Btn_Save_Clefts()
            
