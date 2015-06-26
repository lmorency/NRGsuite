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

import os
import tkFileDialog
import Vars
import Tabs
import Queue
import Constants
import General
import SphereObj
import CleftObj
import BindingSite
import TargetFlex
import pickle

if __debug__:
    from pymol import cmd
    from pymol.wizard import Wizard

    import Sphere
    import FlexSideChain
    import General_cmd


class Config1Vars(Vars.Vars):

    RngOpt = StringVar()
    ResidueValue = StringVar()
    defOptSphere = StringVar()
    defOptCleft = StringVar()
    defOptResidue = StringVar()
    SphereSize = DoubleVar()
    FlexSCStatus = StringVar()
    TargetFlexName = StringVar()
    BindingSiteName = StringVar()

    def __init__(self):

        self.BindingSite = BindingSite.BindingSite()
        self.TargetFlex = TargetFlex.TargetFlex()

class Config1(Tabs.Tab):

    HIGHLIGHT_RELIEF = RAISED
    HIGHLIGHT_RELIEF = SUNKEN
    HIGHLIGHT_WIDTH = 2

    BindingSiteDisplay = 'BINDINGSITE_AREA__'
    SphereDisplay = 'SPHERE__'
    ScaleResolution = 0.50
    
    def Def_Vars(self):

        # class instance objects
        self.listResidues = list()

        # vars class objects
        self.RngOpt = self.Vars.RngOpt
        self.ResidueValue = self.Vars.ResidueValue
        self.defOptSphere = self.Vars.defOptSphere
        self.defOptCleft = self.Vars.defOptCleft
        self.defOptResidue = self.Vars.defOptResidue
        self.SphereSize = self.Vars.SphereSize
        self.FlexSCStatus = self.Vars.FlexSCStatus
        self.TargetFlexName = self.Vars.TargetFlexName
        self.BindingSiteName = self.Vars.BindingSiteName

    def Init_Vars(self):

        self.RngOpt.set('')

        self.defOptSphere.set('')
        self.defOptCleft.set('')
        self.defOptResidue.set('')
        self.SphereSize.set(0.5)

        self.ResidueValue.set('')

        self.TargetFlexName.set('')
        self.BindingSiteName.set('')

        self.FlexSCStatus.set('No flexible side-chain(s) set')

        self.Vars.TargetFlex.Clear_SideChain()
        self.Vars.BindingSite.Clear()

        self.CleftTmpPath = os.path.join(self.top.FlexAIDBindingSiteProject_Dir,'tmp.pdb')
        self.TargetName = self.top.IOFile.TargetName.get()

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback to StringVars
    =================================================================================  '''
    def Trace(self):

        try:
            self.RngOptTrace = self.RngOpt.trace('w',self.RngOpt_Toggle)
            self.defOptCleftTrace = self.defOptCleft.trace('w',self.defOptCleft_Toggle)
        except:
            pass

    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    =================================================================================  '''
    def Del_Trace(self):

        try:
            self.RngOpt.trace_vdelete('w',self.RngOptTrace)
            self.defOptCleft.trace_vdelete('w',self.defOptCleftTrace)
        except:
            pass

    ''' ==================================================================================
    FUNCTION Before_Kill_Frame: Kills the main frame window
    =================================================================================  '''
    def Before_Kill_Frame(self):

        self.Highlight_SelectedCleft('')

        return True

    ''' ==================================================================================
    FUNCTION After_Show: Actions related after showing the frame
    ==================================================================================  '''
    def After_Show(self):

        self.Btn_OptSphRefresh_Clicked()
        #self.Update_FlexSideChain_DDL()
        self.Update_Clefts_DDL()

    ''' ==================================================================================
    FUNCTION Frame: Generate the Configuration Options frame in the the middle
                    frame section
    =================================================================================  '''
    def Frame(self):

        self.fConfig = Frame(self.top.fMiddle)

        #==================================================================================
        '''                        ---  FRAME -  TOP   SIDE  ---                        '''
        #==================================================================================

        self.fBindingSite = Frame(self.fConfig, relief=RAISED, border=1)
        self.fBindingSite.pack(fill=BOTH, side=TOP, padx=5, pady=10)

        Label(self.fBindingSite, text='Binding-site definition', font=self.top.font_Title_H).pack(side=TOP, fill=X, pady=3)

        fConfBS = Frame(self.fBindingSite)
        fConfBS.pack(side=TOP, fill=X, expand=True, padx=5, pady=5)

        Label(fConfBS, text='Pre-configured binding-site:', width=30, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfBS, text='Load', command=self.Btn_LoadConfBS_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfBS, text='Save', command=self.Btn_SaveConfBS_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Entry(fConfBS, textvariable=self.BindingSiteName, font=self.top.font_Text, state='disabled', disabledbackground=self.top.Color_White,
                        disabledforeground=self.top.Color_Black, justify=CENTER).pack(side=LEFT, fill=X, expand=True)

        #=================================================================================
        # Range Optimization fields
        #==================================================================================
        fRangeOpt = Frame(self.fBindingSite)
        fRangeOpt.pack(fill=BOTH, expand=True, padx=5, pady=5)

        Label(fRangeOpt, text='Choose binding-site type', font=self.top.font_Title).pack(side=TOP, fill=X)

        self.fRangeOptLeft = Frame(fRangeOpt)
        self.fRangeOptLeft.pack(side=LEFT, padx=3, expand=True, fill=BOTH)
        fRangeOptInLeft = Frame(self.fRangeOptLeft)
        fRangeOptInLeft.pack(side=LEFT, padx=5, expand=True, fill=BOTH, pady=5)

        fRangeOptLeftLine1 = Frame(fRangeOptInLeft)
        fRangeOptLeftLine1.pack(side=TOP, fill=X)
        fRangeOptLeftLine2 = Frame(fRangeOptInLeft)
        fRangeOptLeftLine2.pack(side=TOP, fill=X)
        fRangeOptLeftLine3 = Frame(fRangeOptInLeft)
        fRangeOptLeftLine3.pack(side=TOP, fill=X)

        self.fRangeOptRight = Frame(fRangeOpt)
        self.fRangeOptRight.pack(side=LEFT, padx=3, expand=True, fill=BOTH)
        fRangeOptInRight = Frame(self.fRangeOptRight)
        fRangeOptInRight.pack(side=LEFT, padx=5, expand=True, fill=BOTH, pady=5)

        fRangeOptRightLine1 = Frame(fRangeOptInRight)
        fRangeOptRightLine1.pack(side=TOP, fill=X)
        fRangeOptRightLine2 = Frame(fRangeOptInRight)
        fRangeOptRightLine2.pack(side=TOP, fill=X)
        fRangeOptRightLine3 = Frame(fRangeOptInRight)
        fRangeOptRightLine3.pack(side=TOP, fill=X)

        #=====================================================================================
        #SPHERE Section
        self.RadioBtn_Sphere = Radiobutton(fRangeOptLeftLine1, text='SPHERE', variable=self.RngOpt, value='LOCCEN', font=self.top.font_Text, disabledforeground=self.top.Color_White)
        self.RadioBtn_Sphere.pack(side=TOP, anchor=N, pady=3)

        self.lblCenter = Label(fRangeOptLeftLine2, text='Center with:', width=10, font=self.top.font_Text)
        self.lblCenter.pack(side=LEFT)
        optionTuple = '',
        self.OptMenuSphere = apply(OptionMenu, (fRangeOptLeftLine2, self.defOptSphere) + optionTuple)
        self.OptMenuSphere.config(state='disabled', width=15)
        self.OptMenuSphere.pack(side=LEFT, fill=X, expand=True)

        self.Btn_OptSphRefresh = Button(fRangeOptLeftLine2, text='Refresh', command=self.Btn_OptSphRefresh_Clicked, font=self.top.font_Text)
        self.Btn_OptSphRefresh.pack(side=LEFT)
        self.Btn_OptSphRefresh.config(state='disabled')

        self.lblRadius = Label(fRangeOptLeftLine3, text='Radius:', font=self.top.font_Text, justify=RIGHT)
        self.lblRadius.pack(side=LEFT, anchor=SE)

        self.sclResizeSphere = Scale(fRangeOptLeftLine3, showvalue=0, from_ = 0.5, to = 0.5, orient=HORIZONTAL, length=150,
                                     resolution=self.ScaleResolution, command=self.ResizeSphere, variable=self.SphereSize)
        self.sclResizeSphere.pack(side=LEFT, fill=X, expand=True)
        self.sclResizeSphere.config(state='disabled')

        self.Btn_EditSphere = Button(fRangeOptLeftLine3, text='Edit', command=self.Btn_EditSphere_Clicked, font=self.top.font_Text)
        self.Btn_EditSphere.pack(side=LEFT, anchor=S)
        self.Btn_EditSphere.config(state='disabled')

        #=====================================================================================
        # CLEFT Section
        #=====================================================================================
        self.RadioBtn_Cleft = Radiobutton(fRangeOptRightLine1, text='CLEFT', variable=self.RngOpt, value='LOCCLF', font=self.top.font_Text, disabledforeground=self.top.Color_White)
        self.RadioBtn_Cleft.pack(side=TOP, anchor=N, pady=3)

        #Label(fRangeOptRightLine2, text='Clefts:', width=10, font=self.top.font_Text).pack(side=LEFT)

        self.Btn_ImportCleft = Button(fRangeOptRightLine2, text='Import clefts', command=self.Btn_Import_Clefts, font=self.top.font_Text, width=47)
        self.Btn_ImportCleft.pack(side=TOP, fill=X, expand=True)
        self.Btn_ImportCleft.config(state='disabled')

        self.Btn_ClearCleft = Button(fRangeOptRightLine3, text='Clear', command=self.Btn_ClearCleft_Clicked, font=self.top.font_Text)
        self.Btn_ClearCleft.pack(side=RIGHT)
        self.Btn_ClearCleft.config(state='disabled')

        self.Btn_DeleteCleft = Button(fRangeOptRightLine3, text='Delete', command=self.Btn_DeleteCleft_Clicked, font=self.top.font_Text)
        self.Btn_DeleteCleft.pack(side=RIGHT)
        self.Btn_DeleteCleft.config(state='disabled')

        self.Btn_DeleteOthersCleft = Button(fRangeOptRightLine3, text='Delete others', command=self.Btn_DeleteOthersCleft_Clicked, font=self.top.font_Text)
        self.Btn_DeleteOthersCleft.pack(side=RIGHT)
        self.Btn_DeleteOthersCleft.config(state='disabled')

        optionTuple = '',
        self.OptMenuCleft = apply(OptionMenu, (fRangeOptRightLine3, self.defOptCleft) + optionTuple)
        self.OptMenuCleft.config(state='disabled')
        self.OptMenuCleft.pack(side=RIGHT, fill=X, expand=True)

        #==================================================================================
        '''                        ---  FRAME - BOTTOM SIDE  ---                        '''
        #==================================================================================

        fTargetFlex = Frame(self.fConfig, border=1, relief=RAISED)
        fTargetFlex.pack(fill=BOTH, side=TOP, padx=5, pady=10)

        Label(fTargetFlex, text='Target flexibility', font=self.top.font_Title_H).pack(side=TOP, fill=X, pady=3)

        #==================================================================================
        # Pre-configured protein flexibility
        #==================================================================================
        fConfFlex = Frame(fTargetFlex)
        fConfFlex.pack(side=TOP, fill=X, expand=True, padx=5, pady=5)

        Label(fConfFlex, text='Pre-configured target flexibility:', width=30, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfFlex, text='Load', command=self.Btn_LoadConfFlex_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Button(fConfFlex, text='Save', command=self.Btn_SaveConfFlex_Clicked, font=self.top.font_Text).pack(side=LEFT)
        Entry(fConfFlex, textvariable=self.TargetFlexName, font=self.top.font_Text, state='disabled', disabledbackground=self.top.Color_White,
                        disabledforeground=self.top.Color_Black, justify=CENTER).pack(side=LEFT, fill=X, expand=True)

        #==================================================================================
        # Flexible Side-Chains
        #==================================================================================
        self.fFlexSC = Frame(fTargetFlex)
        self.fFlexSC.pack(side=TOP, fill=X, padx=5, pady=5)

        fFlexSCLine1 = Frame(self.fFlexSC)
        fFlexSCLine1.pack(fill=X, side=TOP)

        self.fFlexSCLeft = Frame(self.fFlexSC)
        self.fFlexSCLeft.pack(side=LEFT, padx=5, expand=True, fill=BOTH)
        self.fFlexSCRight = Frame(self.fFlexSC)
        self.fFlexSCRight.pack(side=RIGHT, padx=5, expand=True, fill=BOTH)

        fFlexSCLeftLine1 = Frame(self.fFlexSCLeft)
        fFlexSCLeftLine1.pack(fill=X, side=TOP)
        fFlexSCLeftLine2 = Frame(self.fFlexSCLeft)
        fFlexSCLeftLine2.pack(fill=X, side=TOP)

        fFlexSCRightLine1 = Frame(self.fFlexSCRight)
        fFlexSCRightLine1.pack(fill=X, side=TOP)
        fFlexSCRightLine2 = Frame(self.fFlexSCRight)
        fFlexSCRightLine2.pack(fill=X, side=TOP)

        Label(fFlexSCLine1, text='Side-chain flexibility', font=self.top.font_Title).pack(side=TOP,fill=X)

        Button(fFlexSCLeftLine1, text='Add/Delete flexible side-chains', font=self.top.font_Text,
               command=self.SelectFlexibleSideChains, width=40).pack(side=TOP, anchor=NW, expand=True, fill=X)

        Entry(fFlexSCLeftLine2, text='', state='disabled', textvariable=self.FlexSCStatus, font=self.top.font_Text,
                disabledforeground=self.top.Color_Black, width=40, disabledbackground=self.top.Color_White,justify=CENTER).pack(side=LEFT,anchor=NW,expand=True, fill=X)

        #Button(fFlexSCLine3, text='Clear', command=self.Btn_DelAllResidu_Clicked, font=self.top.font_Text).pack(side=RIGHT)
        #Button(fFlexSCLine3, text='Delete', command=self.Btn_DelResidu_Clicked, font=self.top.font_Text).pack(side=RIGHT)

        Button(fFlexSCRightLine1, text='Enter', command=self.Btn_EnterResidue_Clicked, font=self.top.font_Text).pack(side=RIGHT)

        self.EntryResidu = Entry(fFlexSCRightLine1, textvariable=self.ResidueValue, background='white',
                                 width=10, justify=CENTER, font=self.top.font_Text)
        self.EntryResidu.pack(side=RIGHT, anchor=E)

        Label(fFlexSCRightLine1, text='Residue code (e.g. ALA31A):', font=self.top.font_Text, justify=LEFT).pack(side=RIGHT, anchor=W, fill=X)

        '''
        #==================================================================================
        # Normal-modes
        #==================================================================================
        fNormalModes = Frame(fTargetFlex)
        #fNormalModes.pack(fill=X, side=TOP, padx=5, pady=5)

        fNMrow1 = Frame(fNormalModes)
        fNMrow1.pack(fill=X, side=TOP)

        fNMrow2 = Frame(fNormalModes)
        fNMrow2.pack(fill=X, side=TOP)

        Label(fNMrow1, text='Normal Modes', font=self.top.font_Title).pack(side=LEFT)
        Label(fNMrow2, text='Upcoming feature...', width=30, font=self.top.font_Text_I).pack(side=LEFT, anchor=W)
        '''

        return self.fConfig

    ''' ==================================================================================
    FUNCTION Get_TargetFlexPath: Returns the default path of target flexibility
    =================================================================================  '''
    def Get_TargetFlexPath(self):
        TARGETNAME = self.TargetName.upper()
        TargetFlexPath = os.path.join(self.top.FlexAIDTargetFlexProject_Dir,TARGETNAME)

        return TargetFlexPath

    ''' ==================================================================================
    FUNCTION Highlight_SelectedCleft: Highlights the selected cleft
    =================================================================================  '''
    def Highlight_SelectedCleft(self, CleftName):

        if self.top.ActiveFrame == self and self.Vars.BindingSite.Type == 2:

            try:
                cmd.color('purpleblue', self.BindingSiteDisplay)
                cmd.refresh()

                if CleftName != '':
                    Cleft = self.Vars.BindingSite.Get_CleftName(CleftName)
                    cmd.color('oxygen', self.BindingSiteDisplay + ' & resi ' + str(Cleft.Index))
                    cmd.refresh()

                    self.defOptCleft.set(CleftName)
            except:
                pass

    ''' ==================================================================================
    FUNCTION defOptCleft_Toggle: Highlights the selected cleft
    =================================================================================  '''
    def defOptCleft_Toggle(self,*args):

        self.Highlight_SelectedCleft(self.defOptCleft.get())

    ''' ==================================================================================
    FUNCTION RngOpt_Toggle: Changes the optimization range
    =================================================================================  '''
    def RngOpt_Toggle(self,*args):

        if self.RngOpt.get() == 'LOCCLF':
            # Cleft controls
            self.Btn_ImportCleft.config(state='normal')
            self.Btn_DeleteCleft.config(state='normal')
            self.Btn_DeleteOthersCleft.config(state='normal')
            self.Btn_ClearCleft.config(state='normal')
            self.OptMenuCleft.config(state='normal')

            # Sphere controls
            self.Btn_EditSphere.config(state='disabled')
            self.Btn_OptSphRefresh.config(state='disabled')
            self.OptMenuSphere.config(state='disabled')

            self.Vars.BindingSite.Set_Cleft()
            self.Display_BindingSite()

        elif self.RngOpt.get() == 'LOCCEN':
            # Cleft controls
            self.Btn_ImportCleft.config(state='disabled')
            self.Btn_DeleteCleft.config(state='disabled')
            self.Btn_DeleteOthersCleft.config(state='disabled')
            self.Btn_ClearCleft.config(state='disabled')
            self.OptMenuCleft.config(state='disabled')

            # Sphere controls
            self.Btn_EditSphere.config(state='normal')
            self.Btn_OptSphRefresh.config(state='disabled')
            self.OptMenuSphere.config(state='disabled')

            if self.Vars.BindingSite.Sphere == None:
                self.Create_NewSphere()

            self.Vars.BindingSite.Set_Sphere()
            self.Display_BindingSite()

        else:
            # Cleft controls
            self.Btn_ImportCleft.config(state='disabled')
            self.Btn_DeleteCleft.config(state='disabled')
            self.Btn_DeleteOthersCleft.config(state='disabled')
            self.Btn_ClearCleft.config(state='disabled')
            self.OptMenuCleft.config(state='disabled')

            # Sphere controls
            self.Btn_EditSphere.config(state='disabled')
            self.Btn_OptSphRefresh.config(state='disabled')
            self.OptMenuSphere.config(state='disabled')

            self.Vars.BindingSite.Unset()
            self.Delete_BindingSite()

        self.Highlight_RngOpt()
        self.top.fMiddle.update_idletasks()

    ''' ==================================================================================
    FUNCTION Create_NewSphere: Creates a new sphere using the target as center of geometry
    ==================================================================================  '''
    def Create_NewSphere(self):

        Center = General_cmd.Get_CenterOfMass2(self.TargetName, cmd.get_state())
        Width = General_cmd.Get_MaxWidth(self.TargetName, cmd.get_state())

        if len(Center) > 0 and Width != -1:
            self.Vars.BindingSite.Sphere = SphereObj.SphereObj(Width/4.0,Width/2.0,Center)
            self.sclResizeSphere.config(from_=0.5,to=self.Vars.BindingSite.Sphere.MaxRadius)
            self.SphereSize.set(self.Vars.BindingSite.Sphere.Radius)
        else:
            self.DisplayMessage("  ERROR: Could not display the default sphere", 1)
            self.RngOpt.set('')

    ''' ==================================================================================
    FUNCTION SelectFlexibleSideChains: Starts the flexible side-chain wizard
    ==================================================================================  '''
    def SelectFlexibleSideChains(self):

        self.FlexSCRunning(True)
        self.top.ActiveWizard = FlexSideChain.flexSC(self,self.queue)
        cmd.set_wizard(self.top.ActiveWizard)
        cmd.refresh()
        self.top.ActiveWizard.Start()

    ''' ==================================================================================
    FUNCTION FlexSCRunning: Disables/Enables controls related to Flexible side-chains
    ==================================================================================  '''
    def FlexSCRunning(self, boolRun):

        if boolRun:
            self.Start_Update()
            self.Disable_Frame()
        else:
            self.End_Update()
            self.Enable_Frame()

            #self.Update_FlexSideChain_DDL()

            if self.top.WizardResult:
                Status = '(' + str(self.top.WizardResult) + ') flexible side-chain(s) set'
            else:
                Status = 'No flexible side-chain(s) set'

            self.FlexSCStatus.set(Status)

    ''' ==================================================================================
    FUNCTION Update_FlexSideChain_DDL: Updates the drop-down-list meny of flex sc.
    ==================================================================================  '''
    def Update_FlexSideChain_DDL(self):

        self.optMenuWidgetRES["menu"].delete(0, END)

        Residue = ''
        for res in sorted(self.Vars.TargetFlex.listSideChain, key=str.lower):
            self.optMenuWidgetRES["menu"].add_command(label=res, command=lambda temp = res: self.optMenuWidgetRES.setvar(self.optMenuWidgetRES.cget("textvariable"), value = temp))
            Residue = res

        self.defOptResidue.set(Residue)

    ''' ==================================================================================
    FUNCTION Update_Clefts_DDL: Updates the drop-down-list of clefts
    ==================================================================================  '''
    def Update_Clefts_DDL(self):

        self.OptMenuCleft["menu"].delete(0, END)

        CleftName = ''
        for clf in self.Vars.BindingSite.Get_SortedCleftNames():
            self.OptMenuCleft["menu"].add_command(label=clf, command=lambda temp = clf: self.Highlight_SelectedCleft(temp))
            CleftName = clf

        self.defOptCleft.set(CleftName)

    ''' ==================================================================================
    FUNCTION Btn_LoadConfFlex_Clicked: Loads a preconfigured protein flexibility
    ==================================================================================  '''
    def Btn_LoadConfFlex_Clicked(self):

        TargetFlexPath = self.Get_TargetFlexPath()

        if not os.path.isdir(TargetFlexPath):
            self.DisplayMessage("  Could not find a TargetFlex folder for your target:", 2)
            self.DisplayMessage("  The default TargetFlex folder is used.", 2)

            TargetFlexPath = self.top.FlexAIDTargetFlexProject_Dir

        LoadPath = tkFileDialog.askopenfilename(filetypes=[('NRG TargetFlex file','*.nrgtf')],
                                                initialdir=TargetFlexPath, title='Load a TargetFlex file')

        if len(LoadPath) > 0:

            LoadPath = os.path.normpath(LoadPath)

            try:
                in_ = open(LoadPath, 'rb')
                TargetFlex = pickle.load(in_)
                in_.close()

                TargetFlex.listSideChain[:] = [ res for res in TargetFlex.listSideChain if self.listResidues.count(res) ]

                if TargetFlex.Count_SideChain() > 0: # or normal modes
                    self.Vars.TargetFlex = TargetFlex
                    #self.Update_FlexSideChain_DDL()

                    self.TargetFlexName.set(os.path.basename(os.path.splitext(LoadPath)[0]))
                else:
                    self.DisplayMessage("  ERROR: The TargetFlex '" + LoadPath + "' does not match the selected target", 2)

            except:
                self.DisplayMessage("  ERROR: Could not read the TargetFlex", 2)

    ''' ==================================================================================
    FUNCTION Btn_SaveConfFlex_Clicked: Saves the current protein flexibility
    ==================================================================================  '''
    def Btn_SaveConfFlex_Clicked(self):

        TargetFlexPath = self.Get_TargetFlexPath()
        #DefaultName = self.Test

        if self.Vars.TargetFlex.Count_SideChain() > 0:

            if not os.path.isdir(TargetFlexPath):
                os.makedirs(TargetFlexPath)

            SaveFile = tkFileDialog.asksaveasfilename(initialdir=TargetFlexPath,
                                                      filetypes=[('NRG TargetFlex file','*.nrgtf')],
                                                      title='Save the TargetFlex file', initialfile='default_targetflex',
                                                      defaultextension='.nrgtf')

            if len(SaveFile) > 0:

                SaveFile = os.path.normpath(SaveFile)

                if General.validate_String(SaveFile, '.nrgtf', True, True, False):
                    self.DisplayMessage("  ERROR: Could not save the file because you entered an ` name.", 2)
                    return

                if self.top.ValidateSaveProject(SaveFile, 'TargetFlex'):
                    self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                    return

                try:
                    out = open(SaveFile, 'wb')
                    pickle.dump(self.Vars.TargetFlex, out)
                    out.close()

                    self.TargetFlexName.set(os.path.basename(os.path.splitext(SaveFile)[0]))

                    self.DisplayMessage("  Successfully saved '" + SaveFile + "'", 0)
                except:
                    self.DisplayMessage("  ERROR: Could not save target flexibility configuration", 0)

    ''' ==================================================================================
    FUNCTION Btn_LoadConfBS_Clicked: Asks for user to load binding-site
    ==================================================================================  '''
    def Btn_LoadConfBS_Clicked(self):

        if self.top.ValidateWizardRunning():
            return

        BindingSitePath = self.Get_BindingSitePath()

        if not os.path.isdir(BindingSitePath):
            self.DisplayMessage("  Could not find a BindingSite folder for your target:", 2)
            self.DisplayMessage("  The default BindingSite folder is used.", 2)

            BindingSitePath = self.top.FlexAIDBindingSiteProject_Dir

        LoadPath = tkFileDialog.askopenfilename(filetypes=[('NRG BindingSite','*.nrgbs')],
                                                initialdir=BindingSitePath, title='Load a BindingSite file')

        if len(LoadPath) > 0:

            LoadPath = os.path.normpath(LoadPath)

            try:
                in_ = open(LoadPath, 'rb')
                TmpBindingSite = pickle.load(in_)
                in_.close()

                if TmpBindingSite.Count_Cleft() > 0 or TmpBindingSite.Sphere != None:
                    self.Vars.BindingSite = TmpBindingSite

                    if self.Vars.BindingSite.Type == 1:
                        self.sclResizeSphere.config(from_=0.5,to=self.Vars.BindingSite.Sphere.MaxRadius)
                        self.SphereSize.set(self.Vars.BindingSite.Sphere.Radius)
                        self.RngOpt.set('LOCCEN')
                    elif self.Vars.BindingSite.Type == 2:
                        self.RngOpt.set('LOCCLF')

                    self.Update_Clefts_DDL()
                    self.BindingSiteName.set(os.path.basename(os.path.splitext(LoadPath)[0]))

                    self.Update_Vars()

                else:
                    self.DisplayMessage("  ERROR: The BindingSite file has unknown format", 2)

            except:
                self.DisplayMessage("  ERROR: Could not read the BindingSite file", 2)

    ''' ==================================================================================
    FUNCTION Btn_SaveConfBS_Clicked: Asks for user to save binding-site
    ==================================================================================  '''
    def Btn_SaveConfBS_Clicked(self):

        if self.top.ValidateWizardRunning():
            return

        BindingSitePath = self.Get_BindingSitePath()

        if self.Vars.BindingSite.Type == 1 or (self.Vars.BindingSite.Type == 2 and self.Vars.BindingSite.Count_Cleft() > 0):

            if not os.path.isdir(BindingSitePath):
                os.makedirs(BindingSitePath)

            SaveFile = tkFileDialog.asksaveasfilename(initialdir=BindingSitePath,
                                                      filetypes=[('NRG BindingSite','*.nrgbs')],
                                                      title='Save the BindingSite file', initialfile='default_bindingsite',
                                                      defaultextension='.nrgbs')
            if len(SaveFile) > 0:
                SaveFile = os.path.normpath(SaveFile)

                if General.validate_String(SaveFile, '.nrgbs', True, True, False):
                    self.DisplayMessage("  ERROR: Could not save the file because you entered an invalid name.", 2)
                    return

                if self.top.ValidateSaveProject(SaveFile, 'BindingSite'):
                    self.DisplayMessage("  ERROR: The file can only be saved at its default location", 2)
                    return

                try:
                    out = open(SaveFile, 'wb')
                    pickle.dump(self.Vars.BindingSite, out)
                    out.close()

                    self.BindingSiteName.set(os.path.basename(os.path.splitext(SaveFile)[0]))
                    self.DisplayMessage("  Successfully saved '" + SaveFile + "'", 0)

                except:
                    self.DisplayMessage("  ERROR: Could not save binding-site configuration", 1)

        else:
            self.DisplayMessage("  No clefts or sphere to save as 'binding-site'", 2)

    ''' ==================================================================================
    FUNCTION Validate_EnterResidue: Validates the residue code
    ==================================================================================  '''
    def Validate_EnterResidue(self, residue):

        res = residue[0:3]
        num = residue[3:(len(residue)-1)]
        chn = residue[(len(residue)-1):len(residue)]

        selString = ' resn ' + res
        selString += ' & resi ' + num

        if chn != '-':
            selString += ' & chain ' + chn
        else:
            selString += ' & chain \'\''

        selString += ' & ! name H*'
        selString += ' & ! name OXT'

        selString += ' & ' + self.top.IOFile.Target

        try:
            n = cmd.count_atoms(selString, state=1)
        except:
            return "An unexpected error occured while validating the residue."

        if res == 'GLY' or res == 'ALA' or res == 'PRO':
            return "Could not validate the residue: The residue does not have flexible bonds."
        elif res not in Constants.nAtoms:
            return "Could not validate the residue: The residue is unknown."
        elif n == 0:
            return "Could not validate the residue: No such residue."
        elif n < Constants.nAtoms[res]:
            return "Could not validate the residue: The residue is missing atoms."
        elif n > Constants.nAtoms[res]:
            return "Could not validate the residue: The residue has too many atoms."

        return ""

    ''' ==================================================================================
    FUNCTION Btn_EnterResidue_Clicked: Add a residue name in the Flexible side chain ddl.
    ==================================================================================  '''
    def Btn_EnterResidue_Clicked(self):

        Residue = self.ResidueValue.get()

        Error = self.Validate_EnterResidue(Residue)

        if Error:
            self.DisplayMessage(Error, 2)
            self.EntryResidu.config(bg=self.top.Color_Red)
        else:
            self.Vars.TargetFlex.Add_SideChain(Residue)

            self.ResidueValue.set('')
            self.EntryResidu.config(bg=self.top.Color_White)

            Status = '(' + str(self.Vars.TargetFlex.Count_SideChain()) + ') flexible side-chain(s) set'
            self.FlexSCStatus.set(Status)

    ''' ==================================================================================
    FUNCTION Highlight_RngOpt: Highlight with a grey bar the RngOpt selected
    ==================================================================================  '''
    def Highlight_RngOpt(self):

        if self.RngOpt.get() == 'LOCCEN':
            self.fRangeOptRight.config(bd=0)
            self.fRangeOptLeft.config(bd=self.HIGHLIGHT_WIDTH, relief=self.HIGHLIGHT_RELIEF)

        elif self.RngOpt.get() == 'LOCCLF':
            self.fRangeOptLeft.config(bd=0)
            self.fRangeOptRight.config(bd=self.HIGHLIGHT_WIDTH, relief=self.HIGHLIGHT_RELIEF)

        else:
            self.fRangeOptRight.config(bd=0)
            self.fRangeOptLeft.config(bd=0)

    ''' ==================================================================================
    FUNCTION Btn_DeleteCleft_Clicked: Deletes the selected item in the DDL
    ==================================================================================  '''
    def Btn_DeleteCleft_Clicked(self):

        if self.defOptCleft.get() != '':
            self.Vars.BindingSite.Remove_CleftName(self.defOptCleft.get())

            self.Display_BindingSite()
            self.Update_Clefts_DDL()

    ''' ==================================================================================
    FUNCTION Btn_DeleteOthersCleft_Clicked: Deletes all clefts except the selected one
    ==================================================================================  '''
    def Btn_DeleteOthersCleft_Clicked(self):

        if self.defOptCleft.get() != '':

            self.Vars.BindingSite.listClefts = [ Cleft for Cleft in self.Vars.BindingSite.listClefts if Cleft.CleftName == self.defOptCleft.get() ]

            self.Display_BindingSite()
            self.Update_Clefts_DDL()

    ''' ==================================================================================
    FUNCTION Btn_ClearCleft_Clicked: Clears the DDL of the clefts
    ==================================================================================  '''
    def Btn_ClearCleft_Clicked(self):

        self.Vars.BindingSite.Clear_Cleft()
        self.Display_BindingSite()
        self.Update_Clefts_DDL()

    ''' ==================================================================================
                                  Edit Sphere radius-center
    ================================================================================== '''
    def Btn_EditSphere_Clicked(self):

        self.Sphere_Clicked()

    ''' ==================================================================================
                                  Sphere binding-site: Starts the wizard
    ================================================================================== '''
    def Sphere_Clicked(self):

        self.SphereRunning(True)
        self.top.ActiveWizard = Sphere.Sphere(self, self.queue, self.Vars.BindingSite.Sphere, self.SphereDisplay, self.SphereSize, '')
        cmd.set_wizard(self.top.ActiveWizard)
        cmd.refresh()

        self.top.ActiveWizard.Start()

    # Deletes the binding-site object
    def Delete_BindingSite(self):
        try:
            cmd.delete(self.BindingSiteDisplay)
            cmd.refresh()
        except:
            pass

    # Disable/Enable widgets
    def SphereRunning(self, boolRun):

        if boolRun:
            self.Start_Update()
            self.Disable_Frame(self.Btn_OptSphRefresh,self.OptMenuSphere,self.sclResizeSphere,self.lblCenter,self.lblRadius)
            self.Delete_BindingSite()

        else:
            self.End_Update()
            self.Enable_Frame()

            # Reset RngOpt to None if there was an error
            if self.top.WizardError:
                self.RngOpt.set('')

            self.SphereSize.set(self.Vars.BindingSite.Sphere.Radius)
            self.Display_BindingSite()

    # ResizeSphere = Get the scale value, then resize the sphere
    def ResizeSphere(self, *args):

        if not self.top.WizardRunning():
            return

        self.top.ActiveWizard.ResizeSphere()

    # Center the sphere based on the selection
    def CenterSphere(self, sel):

        if not self.top.WizardRunning():
            return

        Center = General_cmd.Get_CenterOfMass2(sel, cmd.get_state())

        if len(Center) > 0:
            self.top.ActiveWizard.SphereView.Set_Center(Center)
            self.top.ActiveWizard.DisplaySphere()
            self.defOptSphere.set(sel)

    # Listbox Menu Refresh values
    def Btn_OptSphRefresh_Clicked(self):

        exc = [ self.BindingSiteDisplay, self.SphereDisplay ]
        General_cmd.Refresh_DDL(self.OptMenuSphere,self.defOptSphere,exc,self.CenterSphere)

    ''' ==================================================================================
    FUNCTION Get_BindingSitePath: Retrieves the default path of the bindingsite
    ==================================================================================  '''
    def Get_BindingSitePath(self):

        TARGETNAME = self.TargetName.upper()
        BindingSitePath = os.path.join(self.top.FlexAIDBindingSiteProject_Dir,TARGETNAME)

        return BindingSitePath

    ''' ==================================================================================
    FUNCTION Get_BindingSitePath: Retrieves the default path of clefts
    ==================================================================================  '''
    def Get_CleftPath(self):

        TARGETNAME = self.TargetName.upper()
        CleftPath = os.path.join(self.top.CleftProject_Dir,TARGETNAME)

        return CleftPath

    ''' ==================================================================================
    FUNCTION Btn_Import_Clefts: Asks for user to load clefts
    ==================================================================================  '''
    def Btn_Import_Clefts(self):

        CleftPath = self.Get_CleftPath()

        if not os.path.isdir(CleftPath):
            self.DisplayMessage("  Could not find a Cleft folder for your target:", 2)
            self.DisplayMessage("  The default Cleft folder is used.", 2)

            CleftPath = os.path.normpath(self.top.CleftProject_Dir)

        LoadFiles = tkFileDialog.askopenfilename(filetypes=[('NRG Cleft files','*.nrgclf')],
                                                 initialdir=CleftPath, title='Select cleft file(s) to load',
                                                 multiple=1)

        if len(LoadFiles) > 0:

            for LoadFile in iter(LoadFiles):

                LoadFile = os.path.normpath(LoadFile)

                if os.path.exists(LoadFile):
                    try:
                        in_ = open(LoadFile, 'rb')
                    except (Exception, IOError) as IOerr:
                        print "Error while loading cleft file :"
                        print "   filename: ", LoadFile
                        print "   err no: ", IOerr.errno
                        print "   err code: ", IOerr.errorcode[IOerr.errno]
                        print "   err message: ", os.strerror(IOerr.errno)
                        print "   err: ", IOerr
                        # look if pass is necessary here because of the prints
                        pass

                    try:
                        Cleft = pickle.load(in_)
                    except (Exception, pickle.PickleError) as PicklingErr:
                        print "error while unpickling the cleft file :"
                        print "   err no: ", PicklingErr.errno
                        print "   err code: ", PicklingErr.errorcode[PicklingErr.errno]
                        print "   err message: ", os.strerror(PicklingErr.errno)
                        print "   err: ", PicklingErr
                        # look if pass is necessary here because of the prints
                        pass

                    in_.close()

                    self.Vars.BindingSite.Add_Cleft(Cleft)

            self.Display_BindingSite()
            self.Update_Clefts_DDL()

    ''' ==================================================================================
    FUNCTION Display_BindingSite: Displays the binding-site of the target
    ================================================================================== '''
    def Display_BindingSite(self):

        auto_zoom = cmd.get("auto_zoom")
        cmd.set("auto_zoom", 0)

        self.Delete_BindingSite()

        try:
            if self.Vars.BindingSite.Type == 1:

                cmd.pseudoatom(self.BindingSiteDisplay,
                               pos=self.Vars.BindingSite.Sphere.Center,
                               vdw=self.Vars.BindingSite.Sphere.Radius,
                               state=1)
                cmd.refresh()

                cmd.hide('everything', self.BindingSiteDisplay)
                cmd.refresh()

                cmd.show('spheres', self.BindingSiteDisplay)
                cmd.refresh()

                cmd.set('sphere_transparency', 0.7, self.BindingSiteDisplay)
                cmd.color('purpleblue', self.BindingSiteDisplay)
                cmd.refresh()

            elif self.Vars.BindingSite.Type == 2 and self.Vars.BindingSite.listClefts:
                self.Generate_CleftBindingSite()

                cmd.load(self.CleftTmpPath, self.BindingSiteDisplay, state=1)
                cmd.refresh()

                cmd.hide('everything', self.BindingSiteDisplay)
                cmd.refresh()

                cmd.show('mesh', self.BindingSiteDisplay)
                cmd.refresh()

                cmd.alter(self.BindingSiteDisplay,'vdw=2.00')
                cmd.rebuild(self.BindingSiteDisplay)
                cmd.color('purpleblue', self.BindingSiteDisplay)
                cmd.refresh()
            # cmd.color('purpleblue', self.BindingSiteDisplay)
            # cmd.refresh()

        except:
            self.DisplayMessage("  ERROR: while displaying the binding-site", 2)

        cmd.set("auto_zoom", auto_zoom)

    ''' ==================================================================================
    FUNCTION Generate_CleftBindingSite: Generate a file combining all clefts together
    ================================================================================== '''
    def Generate_CleftBindingSite(self):

        #self.listSpheres = []

        out = file(self.CleftTmpPath, 'w')
        # print self.Vars.BindingSite.listClefts
        for Cleft in self.Vars.BindingSite.listClefts:
            in_ = open(Cleft.CleftFile, 'r')
            lines = in_.readlines()
            in_.close()

            #0         1         2         3         4         5         6         7
            #0123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM     16  C   SPH Z   1      11.271   0.268  -8.282  1.00  2.17
            for line in lines:
                if line.startswith('ATOM  '):
                    outline = line[0:22]
                    outline += '%4d' % Cleft.Index
                    outline += line[26:]

                    '''
                    sphere = SphereObj.SphereObj()
                    sphere.Set_Radius(float(line[60:66].strip()))

                    sphere.Set_Center( [ float(line[30:38].strip()),
                                         float(line[38:46].strip()),
                                         float(line[46:54].strip()) ] )

                    self.listSpheres.append(sphere)
                    '''

                else:
                    outline = line

                out.write(outline)

        out.close()

    ''' ==================================================================================
    FUNCTION Load_Message: Display the message based on the menu selected
    ================================================================================== '''
    def Load_Message(self):

        self.DisplayMessage('', 0)
        self.DisplayMessage('  FlexAID < Target Cfg > Menu.', 2)
        self.DisplayMessage('  INFO: Configure the target molecule', 2)
        self.DisplayMessage('          1) Define the binding-site using a sphere or one or more cleft(s)', 2)
        self.DisplayMessage('          2) Include flexibility in the target, if necessary', 2)
