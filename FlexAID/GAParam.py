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

from Tkinter import *

import General
import Vars
import Tabs 

class GAParamVars(Vars.Vars):
    
    NbTopChrom = StringVar()
    NbGenFreq = StringVar()
    NbGen = StringVar()
    NbChrom = StringVar()
    CrossRate = StringVar()
    MutaRate = StringVar()
    FitModel = StringVar()
    FitAlpha = StringVar()
    FitPeak = StringVar()
    FitScale = StringVar()
    RepModel = StringVar()
    RepSS = StringVar()
    #RepB = StringVar()
    RepDup = IntVar()
    UseAGA = IntVar()
    AGAk1 = StringVar()
    AGAk2 = StringVar()
    AGAk3 = StringVar()
    AGAk4 = StringVar()

class GAParam(Tabs.Tab):

    def Def_Vars(self):
        
        self.NbTopChrom = self.Vars.NbTopChrom
        self.NbGenFreq = self.Vars.NbGenFreq
        self.NbGen = self.Vars.NbGen
        self.NbChrom = self.Vars.NbChrom
        self.CrossRate = self.Vars.CrossRate
        self.MutaRate = self.Vars.MutaRate
        self.FitModel = self.Vars.FitModel
        self.FitAlpha = self.Vars.FitAlpha
        self.FitPeak = self.Vars.FitPeak
        self.FitScale = self.Vars.FitScale
        self.RepModel = self.Vars.RepModel
        self.RepSS = self.Vars.RepSS
        #self.RepB = self.Vars.RepB
        self.RepDup = self.Vars.RepDup
        self.UseAGA = self.Vars.UseAGA
        self.AGAk1 = self.Vars.AGAk1
        self.AGAk2 = self.Vars.AGAk2
        self.AGAk3 = self.Vars.AGAk3
        self.AGAk4 = self.Vars.AGAk4

    def Init_Vars(self):
                
        self.NbTopChrom.set('5')
        self.NbGenFreq.set('25')
        self.NbGen.set('500')
        self.NbChrom.set('500')
        self.CrossRate.set('0.900')
        self.MutaRate.set('0.025')
        self.FitModel.set('PSHARE')
        self.FitAlpha.set('4.0')
        self.FitPeak.set('5.0')
        self.FitScale.set('10.0')
        self.RepModel.set('BOOM')
        self.RepDup.set(0)
        self.RepSS.set('0.10')
        #self.RepB.set('500')
        self.UseAGA.set(1)
        self.AGAk1.set('0.95')
        self.AGAk2.set('0.10')
        self.AGAk3.set('0.95')
        self.AGAk4.set('0.10')

    ''' ==================================================================================
    FUNCTION Before_Kill_Frame: Actions related before killing a frame
    =================================================================================  '''    
    def Before_Kill_Frame(self):
        
        if int(self.NbTopChrom.get()) > int(self.NbChrom.get()):
            self.NbTopChrom.set(self.NbChrom.get())

        if int(self.NbGenFreq.get()) > int(self.NbGen.get()):
            self.NbGenFreq.set(self.NbGen.get())

        # Set k3 = k1 and k4 = k2 because now using max cr (k1) and max mr (k2)
        self.AGAk3.set(self.AGAk1.get())
        self.AGAk4.set(self.AGAk2.get())
        
        return True

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''    
    def Trace(self):

        try:
            self.AGATrace = self.UseAGA.trace('w',self.AGA_Toggle)
            self.RepModelTrace = self.RepModel.trace('w',self.RepModel_Toggle)
            self.FitModelTrace = self.FitModel.trace('w',self.FitModel_Toggle)

            self.NbChromTrace = self.NbChrom.trace('w', lambda *args, **kwargs:
                                                   self.Validate_Field(input=self.inputNbChr, var=self.NbChrom, min=1, max=100000,
                                                                       ndec=-1, tag='Chromosomes', _type=int))

            self.NbGenTrace = self.NbGen.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.inputNbGen, var=self.NbGen, min=1, 
                                                                   max=100000, ndec=-1, tag='Generations', _type=int))
            
            self.CrossRateTrace = self.CrossRate.trace('w', lambda *args, **kwargs:
                                                       self.Validate_Field(input=self.inputCR, var=self.CrossRate, min=0.000, 
                                                                           max=1.000, ndec=3, tag='Crossover rate', _type=float))
            
            self.MutaRateTrace = self.MutaRate.trace('w', lambda *args, **kwargs:
                                                     self.Validate_Field(input=self.inputMR, var=self.MutaRate, min=0.000,
                                                                         max=1.000, ndec=3, tag='Mutation rate', _type=float))

            self.AGAk1Trace = self.AGAk1.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.entAGAk1, var=self.AGAk1, min=0.000,
                                                                   max=1.000, ndec=3, tag='Adaptive constant k1', _type=float))
            
            self.AGAk2Trace = self.AGAk2.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.entAGAk2, var=self.AGAk2, min=0.000,
                                                                   max=1.000, ndec=3, tag='Adaptive constant k2', _type=float))

            self.AGAk3Trace = self.AGAk3.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.entAGAk3, var=self.AGAk3, min=0.000,
                                                                   max=1.000, ndec=3, tag='Adaptive constant k3', _type=float))

            self.AGAk4Trace = self.AGAk4.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.entAGAk4, var=self.AGAk4, min=0.000,
                                                                   max=1.000, ndec=3, tag='Adaptive constant k4', _type=float))

            self.NbTopChromTrace = self.NbTopChrom.trace('w', lambda *args, **kwargs:
                                                         self.Validate_Field(input=self.inputTopChr, var=self.NbTopChrom, min=0,
                                                                             max=self.inputNbChr, ndec=-1, tag='TOP complexes', _type=int))

            self.NbGenFreqTrace = self.NbGenFreq.trace('w', lambda *args, **kwargs:
                                                       self.Validate_Field(input=self.inputGenFq, var=self.NbGenFreq, min=1,
                                                                           max=self.inputNbGen, ndec=-1, tag='Generations interval', _type=int))

            self.FitScaleTrace = self.FitScale.trace('w', lambda *args, **kwargs:
                                                     self.Validate_Field(input=self.entScale, var=self.FitScale, min=0.00,
                                                                         max=100.00, ndec=2, tag='Fitness scale', _type=float))

            self.FitPeakTrace = self.FitPeak.trace('w', lambda *args, **kwargs:
                                                   self.Validate_Field(input=self.entPeak, var=self.FitPeak, min=0.00,
                                                                       max=100.00, ndec=2, tag='Fitness peak', _type=float))

            self.FitAlphaTrace = self.FitAlpha.trace('w', lambda *args, **kwargs:
                                                     self.Validate_Field(input=self.entAlpha, var=self.FitAlpha, min=0.00,
                                                                         max=100.00, ndec=2, tag='Fitness alpha', _type=float))
            
            self.RepSSTrace = self.RepSS.trace('w', lambda *args, **kwargs:
                                               self.Validate_Field(input=self.entSS, var=self.RepSS, min=0.00,
                                                                   max=1.00, ndec=2, tag='Reproduction steady-state', _type=float))

            #self.RepBTrace = self.RepB.trace('w', lambda *args, **kwargs:
            #                                      self.Validate_Field(input=self.entB, var=self.RepB, min=0.10,
            #                                      max=5.00, ndec=2, tag='Reproduction population boom', _type=float))


        except:
            pass
    
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''    
    def Del_Trace(self):

        try:
            self.UseAGA.trace_vdelete('w',self.AGATrace)
            self.RepModel.trace_vdelete('w',self.RepModelTrace)
            self.FitModel.trace_vdelete('w',self.FitModelTrace)
            self.NbChrom.trace_vdelete('w',self.NbChromTrace)
            self.NbGen.trace_vdelete('w',self.NbGenTrace)
            self.CrossRate.trace_vdelete('w',self.CrossRateTrace)
            self.MutaRate.trace_vdelete('w',self.MutaRateTrace)
            self.AGAk1.trace_vdelete('w',self.AGAk1Trace)
            self.AGAk2.trace_vdelete('w',self.AGAk2Trace)
            self.AGAk3.trace_vdelete('w',self.AGAk3Trace)
            self.AGAk4.trace_vdelete('w',self.AGAk4Trace)
            self.NbTopChrom.trace_vdelete('w',self.NbTopChromTrace)
            self.NbGenFreq.trace_vdelete('w',self.NbGenFreqTrace)
            self.FitScale.trace_vdelete('w',self.FitScaleTrace)
            self.FitPeak.trace_vdelete('w',self.FitPeakTrace)
            self.FitAlpha.trace_vdelete('w',self.FitAlphaTrace)
            self.RepSS.trace_vdelete('w',self.RepSSTrace)
            #self.RepB.trace_vdelete('w',self.RepBTrace)
        except:
            pass

    ''' ==================================================================================
    FUNCTION Frame: Generate the Parameter Options frame in the the middle 
                    frame section
    =================================================================================  '''    
    def Frame(self):
        
        self.fGAParam = Frame(self.top.fMiddle)
        
        fPLeft = Frame(self.fGAParam)
        fPLeft.pack(side=LEFT, fill=BOTH, expand=True)
        fPRight = Frame(self.fGAParam)
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
        self.inputNbChr = Entry(fGeneticLine2, width=5, background='white', justify=CENTER, textvariable=self.NbChrom, font=self.top.font_Text)
        self.inputNbChr.pack(side=RIGHT, anchor=NW)
        self.ValidNbChrom = [1, False, self.inputNbChr]

        Label(fGeneticLine3, text='Number of generations:', font=self.top.font_Text).pack(side=LEFT)
        self.inputNbGen = Entry(fGeneticLine3, width=5, background='white', justify=CENTER, textvariable=self.NbGen, font=self.top.font_Text)
        self.inputNbGen.pack(side=RIGHT, anchor=NW)
        self.ValidNbGen = [1, False, self.inputNbGen]

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
        self.inputCR = Entry(fOperatorsLine2, width=5, background='white', justify=CENTER,
                             textvariable=self.CrossRate, font=self.top.font_Text, state='disabled')
        self.inputCR.pack(side=RIGHT, anchor=NW)
        self.ValidCrossRate = [1, False, self.inputCR]
        
        Label(fOperatorsLine3, text='Mutation rate:', font=self.top.font_Text).pack(side=LEFT, anchor=W) 
        self.inputMR = Entry(fOperatorsLine3, width=5, background='white', justify=CENTER,
                             textvariable=self.MutaRate, font=self.top.font_Text, state='disabled')
        self.inputMR.pack(side=RIGHT, anchor=W)
        self.ValidMutaRate = [1, False, self.inputMR]

        ''' ========================================================================== '''

        fAGA = Frame(fPLeft)
        fAGA.pack(side=TOP, fill=X, padx=5, pady=5)
        fAGALine1 = Frame(fAGA)
        fAGALine1.pack(side=TOP, fill=X, padx=5, pady=2)
        fAGALine2 = Frame(fAGA)
        fAGALine2.pack(side=TOP, fill=X, padx=5, pady=2)

        Label(fAGALine1, text='Adaptive operators:', font=self.top.font_Text).pack(side=LEFT)
        Checkbutton(fAGALine1, text=' Use', variable=self.UseAGA, font=self.top.font_Text).pack(side=RIGHT)

        self.entAGAk2 = Entry(fAGALine2, width=5, background='white', justify=CENTER, 
                              textvariable=self.AGAk2, font=self.top.font_Text)
        self.entAGAk2.pack(side=RIGHT, anchor=W)
        Label(fAGALine2, text='Max M.R.:',font=self.top.font_Text).pack(side=RIGHT)
        self.ValidAGAk2 = [1, False, self.entAGAk2]

        self.entAGAk1 = Entry(fAGALine2, width=5, background='white', justify=CENTER,
                              textvariable=self.AGAk1, font=self.top.font_Text)
        self.entAGAk1.pack(side=RIGHT, anchor=W)
        Label(fAGALine2, text='Max C.R.:',font=self.top.font_Text).pack(side=RIGHT)
        self.ValidAGAk1 = [1, False, self.entAGAk1]
        
        self.entAGAk4 = Entry(fAGALine2, width=5, background='white', justify=CENTER,
                              textvariable=self.AGAk4, font=self.top.font_Text)
        #self.entAGAk4.pack(side=RIGHT, anchor=W)
        #Label(fAGALine2, text='k4:',font=self.top.font_Text).pack(side=RIGHT)
        self.ValidAGAk4 = [1, False, self.entAGAk4]

        self.entAGAk3 = Entry(fAGALine2, width=5, background='white', justify=CENTER,
                              textvariable=self.AGAk3, font=self.top.font_Text)
        #self.entAGAk3.pack(side=RIGHT, anchor=W)
        #Label(fAGALine2, text='k3:',font=self.top.font_Text).pack(side=RIGHT)
        self.ValidAGAk3 = [1, False, self.entAGAk3]
        
        #self.entAGAk2 = Entry(fAGALine2, width=5, background='white', justify=CENTER,textvariable=self.AGAk2, font=self.top.font_Text)
        #self.entAGAk2.pack(side=RIGHT, anchor=W)
        #Label(fAGALine2, text='k2:',font=self.top.font_Text).pack(side=RIGHT)
        #self.ValidAGAk2 = [1, False, self.entAGAk2]

        #self.entAGAk1 = Entry(fAGALine2, width=5, background='white', justify=CENTER, textvariable=self.AGAk1, font=self.top.font_Text)
        #self.entAGAk1.pack(side=RIGHT, anchor=W)
        #Label(fAGALine2, text='k1:',font=self.top.font_Text).pack(side=RIGHT)
        #self.ValidAGAk1 = [1, False, self.entAGAk1]

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
        self.inputTopChr = Entry(fVisualLine2, width=5, background='white', justify=CENTER, textvariable=self.NbTopChrom, font=self.top.font_Text)
        self.inputTopChr.pack(side=RIGHT)
        self.ValidNbTopChrom = [1, False, self.inputTopChr]
        
        Label(fVisualLine3, text='Refresh interval:', font=self.top.font_Text).pack(side=LEFT)
        self.inputGenFq = Entry(fVisualLine3, width=5, background='white', justify=CENTER, textvariable=self.NbGenFreq, font=self.top.font_Text)
        self.inputGenFq.pack(side=RIGHT)
        self.ValidNbGenFreq = [1, False, self.inputGenFq]

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
        #self.entScale.pack(side=RIGHT)
        #Label(fFitnessLine3, text='Scale:', font=self.top.font_Text).pack(side=RIGHT)
        self.ValidFitScale = [1, False, self.entScale]
        
        self.entPeak = Entry(fFitnessLine3, width=4, background='white', justify=CENTER, textvariable=self.FitPeak, font=self.top.font_Text)
        #self.entPeak.pack(side=RIGHT)
        #Label(fFitnessLine3, text='Peak:', font=self.top.font_Text).pack(side=RIGHT)
        self.ValidFitPeak = [1, False, self.entPeak]

        self.entAlpha = Entry(fFitnessLine3, width=4, background='white', justify=CENTER, textvariable=self.FitAlpha, font=self.top.font_Text)
        #self.entAlpha.pack(side=RIGHT)
        #Label(fFitnessLine3, text='alpha:', font=self.top.font_Text).pack(side=RIGHT)
        self.ValidFitAlpha = [1, False, self.entAlpha]

        
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
        Radiobutton(fRepModelLine2, text='Elitism (Steady-state)', value='STEADY', variable=self.RepModel, font=self.top.font_Text).pack(side=LEFT)
        self.entSS = Entry(fRepModelLine2, width=5, background='white', justify=CENTER,
                           textvariable=self.RepSS, font=self.top.font_Text, state='disabled')
        self.entSS.pack(side=RIGHT)
        Label(fRepModelLine2, text='Fraction:', font=self.top.font_Text).pack(side=RIGHT)
        self.ValidRepSS = [1, False, self.entSS]

        Radiobutton(fRepModelLine3, text='Population boom', value='BOOM', variable=self.RepModel, font=self.top.font_Text).pack(side=LEFT)
        #self.entB = Entry(fRepModelLine3, width=5, background='white', justify=CENTER, textvariable=self.RepB, font=self.top.font_Text)
        #self.entB.pack(side=RIGHT)
        #Label(fRepModelLine3, text='Chr.:', font=self.top.font_Text).pack(side=RIGHT)
        #self.ValidRepB = [1, False, self.entB]
        
        Checkbutton(fRepModelLine4, text=' Allow duplicates', variable=self.RepDup, font=self.top.font_Text).pack(side=LEFT)
        
        self.Validator = [self.ValidNbGen, self.ValidNbGenFreq, self.ValidNbChrom, self.ValidNbTopChrom,
                          self.ValidCrossRate, self.ValidMutaRate, self.ValidFitAlpha, self.ValidFitPeak,
                          self.ValidFitScale, self.ValidRepSS, #self.ValidRepB, 
                          self.ValidAGAk1, self.ValidAGAk2, self.ValidAGAk3, self.ValidAGAk4]
        
        return self.fGAParam

    ''' ==================================================================================
    FUNCTION RepModel_Toggle: Disables/enables the edit box according to checkstate
    =================================================================================  '''    
    def RepModel_Toggle(self, *args):
        
        if self.RepModel.get() == 'STEADY':
            self.entSS.config(state='normal')
            #self.entB.config(state='disabled')
        elif self.RepModel.get() == 'BOOM':
            self.entSS.config(state='disabled')
            #self.entB.config(state='normal')
        
        return
 
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
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''        
    def Load_Message(self):
        
        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < GA Parameters > Menu.', 2)
        self.DisplayMessage('  INFO:   Set the different options of the Genetic Algorithm.', 2)
