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

import General
import functools

class GAParam:

    def __init__(self,top,PyMOL):
        
        #print "New instance of GAParam Class"
        self.PyMOL = PyMOL
        
        self.top = top
        self.Tab = self.top.Btn_GAParam
        self.FrameName = 'GAParam'

        self.Def_Vars()
        self.Init_Vars()

        self.Frame()
        self.Trace()

    def Def_Vars(self):
        
        self.NbTopChrom = StringVar()
        self.NbGenFreq = StringVar()
        self.NbGen = StringVar()
        self.NbChrom = StringVar()
        self.CrossRate = StringVar()
        self.MutaRate = StringVar()
        self.FitModel = StringVar()
        self.LastFitModel = StringVar()
        self.FitAlpha = StringVar()
        self.FitPeak = StringVar()
        self.FitScale = StringVar()
        self.RepModel = StringVar()
        self.LastRepModel = StringVar()
        self.RepSS = StringVar()
        self.RepB = StringVar()
        self.RepDup = IntVar()
        self.UseAGA = IntVar()
        self.LastUseAGA = IntVar()
        self.AGAk1 = StringVar()
        self.AGAk2 = StringVar()
        self.AGAk3 = StringVar()
        self.AGAk4 = StringVar()

        self.ValidNbGen = list()
        self.ValidNbGenFreq = list()
        self.ValidNbChrom = list()
        self.ValidNbTopChrom = list()
        self.ValidCrossRate = list()
        self.ValidMutaRate = list()
        self.ValidFitAlpha = list()
        self.ValidFitPeak = list()
        self.ValidFitScale = list()
        self.ValidRepSS = list()
        self.ValidRepB = list()
        self.ValidAGAk1 = list()
        self.ValidAGAk2 = list()
        self.ValidAGAk3 = list()
        self.ValidAGAk4 = list()

        self.Validator = list()


    def Init_Vars(self):
                
        self.NbTopChrom.set('5')
        self.NbGenFreq.set('10')
        self.NbGen.set('300')
        self.NbChrom.set('300')
        self.CrossRate.set('0.900')
        self.MutaRate.set('0.025')
        self.FitModel.set('PSHARE')
        self.FitAlpha.set('4.0')
        self.FitPeak.set('5.0')
        self.FitScale.set('10.0')
        self.RepModel.set('BOOM')
        self.RepDup.set(0)
        self.RepSS.set('500')
        self.RepB.set('1.0')
        self.UseAGA.set(0)
        self.AGAk1.set('0.95')
        self.AGAk2.set('0.10')
        self.AGAk3.set('0.95')
        self.AGAk4.set('0.10')

        # [ Validated, Callback_Validator_Fail(2), Skip Validation, Widget reference ]
        self.ValidNbGen = [True, 1, 0, None]
        self.ValidNbGenFreq = [True, 1, 0, None]
        self.ValidNbChrom = [True, 1, 0, None]
        self.ValidNbTopChrom = [True, 1, 0, None]
        self.ValidCrossRate = [True, 1, 0, None]
        self.ValidMutaRate = [True, 1, 0, None]
        self.ValidFitAlpha = [True, 1, 0, None]
        self.ValidFitPeak = [True, 1, 0, None]
        self.ValidFitScale = [True, 1, 0, None]
        self.ValidRepSS = [True, 1, 0, None]
        self.ValidRepB = [True, 1, 0, None]
        self.ValidAGAk1 = [True, 1, 0, None]
        self.ValidAGAk2 = [True, 1, 0, None]
        self.ValidAGAk3 = [True, 1, 0, None]
        self.ValidAGAk4 = [True, 1, 0, None]
        
        self.Validator = [self.ValidNbGen, self.ValidNbGenFreq, self.ValidNbChrom, self.ValidNbTopChrom,
                          self.ValidCrossRate, self.ValidMutaRate, self.ValidFitAlpha, self.ValidFitPeak,
                          self.ValidFitScale, self.ValidRepSS, self.ValidRepB, 
                          self.ValidAGAk1, self.ValidAGAk2, self.ValidAGAk3, self.ValidAGAk4]

    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
        
        if int(self.NbTopChrom.get()) > int(self.NbChrom.get()):
            self.NbTopChrom.set(self.NbChrom.get())

        if int(self.NbGenFreq.get()) > int(self.NbGen.get()):
            self.NbGenFreq.set(self.NbGen.get())

        self.fParam.pack_forget()
        #self.fParam.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Validator_Fail: Triggers visual events upon validation failure
    =================================================================================  '''    
    def Validator_Fail(self):

        pass

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):
        
        self.LoadMessage()

        self.fParam.pack(fill=BOTH, expand=True)

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        self.AGATrace = self.UseAGA.trace('w',self.AGA_Toggle)
        self.RepModelTrace = self.RepModel.trace('w',self.RepModel_Toggle)
        self.FitModelTrace = self.FitModel.trace('w',self.FitModel_Toggle)

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''    
    def Del_Trace(self):

        self.UseAGA.trace_vdelete('w',self.AGATrace)
        self.RepModel.trace_vdelete('w',self.RepModelTrace)
        self.FitModel.trace_vdelete('w',self.FitModelTrace)

    ''' ==================================================================================
    FUNCTION Frame: Generate the Parameter Options frame in the the middle 
                    frame section
    =================================================================================  '''    
    def Frame(self):
        
        self.fParam = Frame(self.top.fMiddle)
        fPLeft = Frame(self.fParam)
        fPLeft.pack(side=LEFT, fill=BOTH, expand=True)
        fPRight = Frame(self.fParam)
        fPRight.pack(side=RIGHT, fill=BOTH, expand=True)
        
        ''' ========================================================================== '''

        fGenetic = Frame(fPLeft)
        fGenetic.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        fGeneticLine1 = Frame(fGenetic)
        fGeneticLine1.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)
        fGeneticLine2 = Frame(fGenetic)
        fGeneticLine2.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)
        fGeneticLine3 = Frame(fGenetic)
        fGeneticLine3.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)

        Label(fGeneticLine1, text='Genetic parameters', font=self.top.font_Title).pack(side=LEFT)
        Label(fGeneticLine2, text='Number of chromosomes:', font=self.top.font_Text).pack(side=LEFT)
        inputNbChr = Entry(fGeneticLine2, width=5, background='white', justify=CENTER, textvariable=self.NbChrom, font=self.top.font_Text)
        inputNbChr.pack(side=RIGHT, anchor=NW)
        args_list = [inputNbChr, self.NbChrom, 1, 100000, -1, self.ValidNbChrom,'Chromosomes','int']
        inputNbChr.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidNbChrom[3] = inputNbChr

        Label(fGeneticLine3, text='Number of generations:', font=self.top.font_Text).pack(side=LEFT)
        inputNbGen = Entry(fGeneticLine3, width=5, background='white', justify=CENTER, textvariable=self.NbGen, font=self.top.font_Text)
        inputNbGen.pack(side=RIGHT, anchor=NW)
        args_list = [inputNbGen, self.NbGen, 1, 100000, -1, self.ValidNbGen,'Generations','int']
        inputNbGen.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidNbGen[3] = inputNbGen

        ''' ========================================================================== '''

        fOperators = Frame(fPLeft)
        fOperators.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        fOperatorsLine1 = Frame(fOperators)
        fOperatorsLine1.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)
        fOperatorsLine2 = Frame(fOperators)
        fOperatorsLine2.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)
        fOperatorsLine3 = Frame(fOperators)
        fOperatorsLine3.pack(side=TOP, anchor=W, fill=X, padx=5, pady=2)

        Label(fOperatorsLine1, text='Genetic operators', font=self.top.font_Title).pack(side=LEFT)
        Label(fOperatorsLine2, text='Crossover rate:', font=self.top.font_Text).pack(side=LEFT, anchor=W)
        self.inputCR = Entry(fOperatorsLine2, width=5, background='white', justify=CENTER, textvariable=self.CrossRate, font=self.top.font_Text)
        self.inputCR.pack(side=RIGHT, anchor=NW)
        args_list = [self.inputCR, self.CrossRate, 0.000, 1.000, 3, self.ValidCrossRate,'Crossover Rate','float']
        self.inputCR.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidCrossRate[3] = self.inputCR
        
        Label(fOperatorsLine3, text='Mutation rate:', font=self.top.font_Text).pack(side=LEFT, anchor=W) 
        self.inputMR = Entry(fOperatorsLine3, width=5, background='white', justify=CENTER, textvariable=self.MutaRate, font=self.top.font_Text)
        self.inputMR.pack(side=RIGHT, anchor=W)
        args_list = [self.inputMR, self.MutaRate, 0.001, 1.000, 3, self.ValidMutaRate,'Mutation Rate','float']
        self.inputMR.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidMutaRate[3] = self.inputMR

        ''' ========================================================================== '''

        fAGA = Frame(fPLeft)
        fAGA.pack(side=TOP, fill=X, padx=5, pady=5)
        fAGALine1 = Frame(fAGA)
        fAGALine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fAGALine2 = Frame(fAGA)
        fAGALine2.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fAGALine1, text='Adaptive operators:', font=self.top.font_Text).pack(side=LEFT)
        Checkbutton(fAGALine1, text='Use', variable=self.UseAGA, font=self.top.font_Text).pack(side=RIGHT)

        self.entAGAk4 = Entry(fAGALine2, width=5, background='white', justify=CENTER, textvariable=self.AGAk4, font=self.top.font_Text)
        self.entAGAk4.pack(side=RIGHT, anchor=W)
        args_list = [self.entAGAk4, self.AGAk4, 0.001, 1.000, 3, self.ValidAGAk4,'Adative GA Constant k4','float']
        self.entAGAk4.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidAGAk4[3] = self.entAGAk4
        Label(fAGALine2, text='k4:',font=self.top.font_Text).pack(side=RIGHT)

        self.entAGAk3 = Entry(fAGALine2, width=5, background='white', justify=CENTER, textvariable=self.AGAk3, font=self.top.font_Text)
        self.entAGAk3.pack(side=RIGHT, anchor=W)
        args_list = [self.entAGAk3, self.AGAk3, 0.001, 1.000, 3, self.ValidAGAk3,'Adative GA Constant k3','float']
        self.entAGAk3.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidAGAk3[3] = self.entAGAk3
        Label(fAGALine2, text='k3:',font=self.top.font_Text).pack(side=RIGHT)
        
        self.entAGAk2 = Entry(fAGALine2, width=5, background='white', justify=CENTER, textvariable=self.AGAk2, font=self.top.font_Text)
        self.entAGAk2.pack(side=RIGHT, anchor=W)
        args_list = [self.entAGAk2, self.AGAk2, 0.001, 1.000, 3, self.ValidAGAk2,'Adative GA Constant k2','float']
        self.entAGAk2.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidAGAk2[3] = self.entAGAk2
        Label(fAGALine2, text='k2:',font=self.top.font_Text).pack(side=RIGHT)

        self.entAGAk1 = Entry(fAGALine2, width=5, background='white', justify=CENTER, textvariable=self.AGAk1, font=self.top.font_Text)
        self.entAGAk1.pack(side=RIGHT, anchor=W)
        args_list = [self.entAGAk1, self.AGAk1, 0.001, 1.000, 3, self.ValidAGAk1,'Adative GA Constant k1','float']
        self.entAGAk1.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidAGAk1[3] = self.entAGAk1
        Label(fAGALine2, text='k1:',font=self.top.font_Text).pack(side=RIGHT)

        ''' ========================================================================== '''

        fVisual = Frame(fPLeft)
        fVisual.pack(side=TOP, anchor=W, fill=X, padx=5, pady=5)
        fVisualLine1 = Frame(fVisual)
        fVisualLine1.pack(side=TOP, anchor=W, fill=X)
        fVisualLine2 = Frame(fVisual)
        fVisualLine2.pack(side=TOP, anchor=W, fill=X)
        fVisualLine3 = Frame(fVisual)
        fVisualLine3.pack(side=TOP, anchor=W, fill=X)

        Label(fVisualLine1, text='Visual display', font=self.top.font_Title).pack(side=LEFT)
        Label(fVisualLine2, text='Number of TOP complexes:', font=self.top.font_Text, justify=LEFT).pack(side=LEFT)
        inputTopChr = Entry(fVisualLine2, width=5, background='white', justify=CENTER, textvariable=self.NbTopChrom, font=self.top.font_Text)
        inputTopChr.pack(side=RIGHT)
        args_list = [inputTopChr, self.NbTopChrom, 0, inputNbChr, -1, self.ValidNbTopChrom,'TOP Complexes','int']
        inputTopChr.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidNbTopChrom[3] = inputTopChr

        Label(fVisualLine3, text='Refresh interval:', font=self.top.font_Text).pack(side=LEFT)
        inputGenFq = Entry(fVisualLine3, width=5, background='white', justify=CENTER, textvariable=self.NbGenFreq, font=self.top.font_Text)
        inputGenFq.pack(side=RIGHT)
        args_list = [inputGenFq, self.NbGenFreq, 1, inputNbGen, -1, self.ValidNbGenFreq,'Generations Interval','int']
        inputGenFq.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidNbGenFreq[3] = inputGenFq

        ''' ========================================================================== '''        

        fFitness = Frame(fPRight)#, borderwidth=2, relief=SUNKEN)
        fFitness.pack(side=TOP, fill=X, padx=5, pady=5)
        fFitnessLine1 = Frame(fFitness)
        fFitnessLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fFitnessLine2 = Frame(fFitness)
        fFitnessLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fFitnessLine3 = Frame(fFitness)
        fFitnessLine3.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fFitnessLine1, text='Fitness model', font=self.top.font_Title).pack(side=LEFT)
        Radiobutton(fFitnessLine2, text='Linear', value='LINEAR', variable=self.FitModel, font=self.top.font_Text).pack(side=LEFT)
        Radiobutton(fFitnessLine3, text='Share', value='PSHARE', variable=self.FitModel, font=self.top.font_Text).pack(side=LEFT)

        self.entScale = Entry(fFitnessLine3, width=4, background='white', justify=CENTER, textvariable=self.FitScale, font=self.top.font_Text)
        self.entScale.pack(side=RIGHT)
        args_list = [self.entScale, self.FitScale, 0.0, 100.0, 2, self.ValidFitScale,'Fitness scale','float']
        self.entScale.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidFitScale[3] = self.entScale
        Label(fFitnessLine3, text='Scale:', font=self.top.font_Text).pack(side=RIGHT)

        self.entPeak = Entry(fFitnessLine3, width=4, background='white', justify=CENTER, textvariable=self.FitPeak, font=self.top.font_Text)
        self.entPeak.pack(side=RIGHT)
        args_list = [self.entPeak, self.FitPeak, 0.0, 100.0, 2, self.ValidFitPeak,'Fitness peak','float']
        self.entPeak.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidFitPeak[3] = self.entPeak
        Label(fFitnessLine3, text='Peak:', font=self.top.font_Text).pack(side=RIGHT)

        self.entAlpha = Entry(fFitnessLine3, width=4, background='white', justify=CENTER, textvariable=self.FitAlpha, font=self.top.font_Text)
        self.entAlpha.pack(side=RIGHT)
        args_list = [self.entAlpha, self.FitAlpha, 0.0, 100.0, 2, self.ValidFitAlpha,'Fitness alpha','float']
        self.entAlpha.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidFitAlpha[3] = self.entAlpha
        Label(fFitnessLine3, text='alpha:', font=self.top.font_Text).pack(side=RIGHT)
        
        ''' ========================================================================== '''                

        fRepModel = Frame(fPRight)
        fRepModel.pack(side=TOP, fill=X, padx=5, pady=5)
        fRepModelLine1 = Frame(fRepModel)
        fRepModelLine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fRepModelLine2 = Frame(fRepModel)
        fRepModelLine2.pack(side=TOP, fill=X, padx=5, pady=2)
        fRepModelLine3 = Frame(fRepModel)
        fRepModelLine3.pack(side=TOP, fill=X, padx=5, pady=2)
        fRepModelLine4 = Frame(fRepModel)
        fRepModelLine4.pack(side=TOP, fill=X, padx=5, pady=2)
        
        Label(fRepModelLine1, text='Reproduction model', font=self.top.font_Title).pack(side=LEFT)
        Radiobutton(fRepModelLine2, text='Elitism (Steady-state):', value='STEADY', variable=self.RepModel, font=self.top.font_Text).pack(side=LEFT)
        self.entSS = Entry(fRepModelLine2, width=5, background='white', justify=CENTER, textvariable=self.RepSS, font=self.top.font_Text)
        self.entSS.pack(side=RIGHT)
        args_list = [self.entSS, self.RepSS, 1, 1000, -1, self.ValidRepSS,'Reproduction Steady-State','int']
        self.entSS.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidRepSS[3] = self.entSS
        Label(fRepModelLine2, text='Chr.:', font=self.top.font_Text).pack(side=RIGHT)

        Radiobutton(fRepModelLine3, text='Population boom:', value='BOOM', variable=self.RepModel, font=self.top.font_Text).pack(side=LEFT)
        self.entB = Entry(fRepModelLine3, width=5, background='white', justify=CENTER, textvariable=self.RepB, font=self.top.font_Text)
        self.entB.pack(side=RIGHT)
        args_list = [self.entB, self.RepB, 0.01, 5.00, 2, self.ValidRepB,'Reproduction PopBoom','float']
        self.entB.bind('<FocusOut>', functools.partial(self.top.Lost_Focus,args=args_list))
        self.ValidRepB[3] = self.entB
        Label(fRepModelLine3, text='Ratio.:', font=self.top.font_Text).pack(side=RIGHT)
        
        Checkbutton(fRepModelLine4, text='Allow duplicates', variable=self.RepDup, font=self.top.font_Text).pack(side=LEFT)

    ''' ==================================================================================
    FUNCTION RepModel_Toggle: Disables/enables the edit box according to checkstate
    =================================================================================  '''    
    def RepModel_Toggle(self, *args):
        
        if self.RepModel.get() == 'STEADY':
            self.entSS.config(state='normal')
            self.entB.config(state='disabled')
            
        elif self.RepModel.get() == 'BOOM':
            self.entSS.config(state='disabled')
            self.entB.config(state='normal')
            
    ''' ==================================================================================
    FUNCTION FitModel_Toggle: Disables/enables the edit box according to checkstate
    =================================================================================  '''    
    def FitModel_Toggle(self, *args):
        
        if self.FitModel.get() == 'PSHARE':
            self.entAlpha.config(state='normal')
            self.entPeak.config(state='normal')
            self.entScale.config(state='normal')
        else:
            self.entAlpha.config(state='disabled')
            self.entPeak.config(state='disabled')
            self.entScale.config(state='disabled')

    ''' ==================================================================================
    FUNCTION AGA_Toggle: Disables/enables the edit box according to checkstate
    =================================================================================  '''    
    def AGA_Toggle(self, *args):
        
        if self.UseAGA.get() == 1:
            self.inputCR.config(state='disabled')
            self.inputMR.config(state='disabled')

            self.entAGAk1.config(state='normal')
            self.entAGAk2.config(state='normal')
            self.entAGAk3.config(state='normal')
            self.entAGAk4.config(state='normal')

        else:
            self.inputCR.config(state='normal')
            self.inputMR.config(state='normal')

            self.entAGAk1.config(state='disabled')
            self.entAGAk2.config(state='disabled')
            self.entAGAk3.config(state='disabled')
            self.entAGAk4.config(state='disabled')

    ''' ==================================================================================
    FUNCTION MenuLoadMessage: Display the message based on the menu selected
    ================================================================================== '''        
    def LoadMessage(self):
        
        self.top.DisplayMessage('', 0)
        self.top.DisplayMessage('  FlexAID < GA Parameters > Menu.', 2)
        self.top.DisplayMessage('  INFO:   Set the different options of the Genetic Algorithm.', 2)

    
