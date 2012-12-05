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
@title: CropCleft - Interface

@summary: This is an optional interface of GetCleft application, which is accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Dec. 6, 2010
'''

from Tkinter import *

import os
import General
import Geometry
import time
import tkFileDialog

if __debug__:
    from pymol import cmd
    import Sphere

#=========================================================================================
'''                        ---  GETCLEFT's CROPCLEFT FRAME  ---                        '''
#=========================================================================================     
class CropCleft:

    def __init__(self, top, PyMOL):
        
        self.PyMOL = PyMOL

        self.top = top
        self.FrameName = 'CropCleft'
        self.Tab = self.top.Btn_CropCleft

        self.DisplayMessage = self.top.DisplayMessage

        self.Def_Vars()
        self.Init_Vars()


    def Def_Vars(self):
        
        self.Img_Button = list()
        self.Filename = list()
        
        self.Step1Selection = StringVar()
        self.LastStep1Selection = StringVar()
        self.Step2Selection = StringVar()
        self.Step3Output = StringVar()
        self.Step3Check = IntVar()
        
        self.SphereList = dict()
        self.SphereSize = DoubleVar()
        self.SphereCoord = list()

        self.Validator = list()

    def Init_Vars(self):
        
        self.Step = 1
        self.Filename = []

        del self.Img_Button[:]
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/stop.gif')))
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/down.gif')))
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/ok.gif')))

        self.Step1Selection.set('')
        self.LastStep1Selection.set('')
        self.Step2Selection.set('')
        self.Step3Output.set('')
        self.Step3Check.set(1)

        self.Color_Default = '#d9d9d9'
        self.Color_Green = self.top.Color_Green
        self.TempPartition = os.path.join(self.top.GetCleftProject_Dir,'tmppt.pdb')
 
        self.ScaleResolution = 0.25
        self.PartitionDisplay = 'PARTITION_AREA__'
        self.SphereDisplay = 'SPHERE_PT_AREA__'
        self.SphereCoord = []
        self.SphereSize.set(0.0)

        self.Validator = []
    
    ''' ==================================================================================
    FUNCTION Kill_Frame: Kills the main frame window
    =================================================================================  '''    
    def Kill_Frame(self):
                
        self.fCrop.pack_forget()
        self.fCrop.destroy()

        return True

    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        self.Step1SelectionTrace = self.Step1Selection.trace('w', self.Toggle_Step1)
        self.Step2SelectionTrace = self.Step2Selection.trace('w', self.Toggle_Step2)

        self.Step3CheckTrace = self.Step3Check.trace('w', self.Toggle_Step3)
        
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        self.Step1Selection.trace_vdelete('w',self.Step1SelectionTrace)
        self.Step2Selection.trace_vdelete('w',self.Step2SelectionTrace)

        self.Step3Check.trace_vdelete('w', self.Step3CheckTrace)

    ''' ==================================================================================
    FUNCTION Show: Displays the frame onto the middle main frame
    ==================================================================================  '''  
    def Show(self):

        self.Frame()
        
        self.Trace()
        
        # Reset to defaults
        #self.Init_Vars()
        
        # Highlight
        self.activate_Step(self.Step)
        self.highlight_Step(self.Step)

        
        # update list of cleft spheres
        self.update_Step1_DDL()
        # update list of spheres added
        self.update_Spheres()

    ''' ==================================================================================
    FUNCTION Frame: Displays the Default Frame
    ==================================================================================  '''          
    def Frame(self):
              
        self.fCrop = Frame(self.top.fMiddle, relief=RIDGE, border=0, width=400, height=460)
        self.fCrop.pack(fill=X, expand=True)
        self.fCrop.pack_propagate(0)
  
        fCropTitle = Frame(self.fCrop, relief=RIDGE, border=0, width=400, height=40)
        fCropTitle.pack(fill=X, expand=True, side=TOP)
        fCropTitle.pack_propagate(0)

        self.fCropStep1 = Frame(self.fCrop, relief=SUNKEN, border=1, width=400, height=100)
        self.fCropStep1.pack(fill=X, expand=True, side=TOP)
        self.fCropStep1.pack_propagate(0)

        self.fCropStep2 = Frame(self.fCrop, relief=SUNKEN, border=1, width=400, height=180)
        self.fCropStep2.pack(fill=X, expand=True, side=TOP)
        self.fCropStep2.pack_propagate(0)

        self.fCropStep3 = Frame(self.fCrop, relief=SUNKEN, border=1, width=400, height=80)
        self.fCropStep3.pack(fill=X, expand=True, side=TOP)
        self.fCropStep3.pack_propagate(0)

        self.fCropBtn = Frame(self.fCrop, width=400, height=60)
        self.fCropBtn.pack(fill=X, expand=True, side=TOP)
        self.fCropBtn.pack_propagate(0)

        Title = Label(fCropTitle, text = 'Cleft Partitioner', height=5, font=self.top.font_Title)
        Title.pack(side=LEFT, anchor=W, padx=5)
                
        #==================================================================================
        '''            --- Selection of the Cleft to be partioned ---                   '''
        #==================================================================================
        fCropStep1a = Frame(self.fCropStep1, relief=RIDGE, border=0, width=400, height=30)
        fCropStep1a.pack(fill=X, expand=True, side=TOP)
        fCropStep1a.pack_propagate(0)
        
        Step1aTitle = Label(fCropStep1a, text = 'STEP 1:', font=self.top.font_Title)
        Step1aTitle.pack(side=LEFT, anchor=W, padx=5)
        
        Step1aTxt = Label(fCropStep1a, text = 'Choose the Cleft to Partition', font=self.top.font_Text)
        Step1aTxt.pack(side=LEFT, anchor=W)
        
        fCropStep1a.pack(fill=X, expand=True, side=TOP)
        fCropStep1a.pack_propagate(0)

        #--------------------------------------------------------------------------------
        
        fCropStep1b = Frame(self.fCropStep1, relief=RIDGE, border=0, width=400, height=40)
        fCropStep1b.pack(fill=X, expand=True, side=TOP)
        fCropStep1b.pack_propagate(0)
        
        Step1bPad = Label(fCropStep1b, text = '')
        Step1bPad.pack(side=LEFT, anchor=W, padx=1)
        
        Step1bTxt = Label(fCropStep1b, text = ' Cleft Spheres objects:', height=5, font=self.top.font_Text)
        Step1bTxt.pack(side=LEFT, anchor=W, padx=5)
        
        Step1b2Pad = Label(fCropStep1b, text = '')
        Step1b2Pad.pack(side=RIGHT, anchor=E, padx=3)        
        
        optionTuple = '',
        self.OptCleftStep1 = apply(OptionMenu, (fCropStep1b, self.Step1Selection) + optionTuple)
        self.OptCleftStep1.config(width=20, bg=self.top.Color_White)
        self.OptCleftStep1.pack(side=LEFT, anchor=W)
                        
        self.Step1bBtn = Button(fCropStep1b, width=20, height=20, image=self.Img_Button[1], command=self.Step1_Next) 
        self.Step1bBtn.pack(side=RIGHT, anchor=E, padx=7)
        self.Step1bBtn.config(state='disabled')
                
        
        #--------------------------------------------------------------------------------

        #==================================================================================
        '''            --- Add spheres to partition the cleft ---                       '''
        #==================================================================================
        fCropStep2a = Frame(self.fCropStep2, relief=RIDGE, border=0, width=400, height=30)
        fCropStep2a.pack(fill=X, expand=True, side=TOP)
        fCropStep2a.pack_propagate(0)

        Step2aTitle = Label(fCropStep2a, text = 'STEP 2:', font=self.top.font_Title)
        Step2aTitle.pack(side=LEFT, anchor=S, padx=5)
        
        Step2aTxt = Label(fCropStep2a, text = 'Add Sphere(s) to Partition the Cleft', font=self.top.font_Text)
        Step2aTxt.pack(side=LEFT, anchor=S)
                
        #--------------------------------------------------------------------------------

        fCropStep2b = Frame(self.fCropStep2, relief=RIDGE, border=0, width=400, height=15)
        fCropStep2b.pack(fill=X, expand=True, side=TOP)
        fCropStep2b.pack_propagate(0)

        Step2bPad = Label(fCropStep2b, text = '')
        Step2bPad.pack(side=LEFT, anchor=N, padx=28)
        
        Step2bTxt = Label(fCropStep2b, text = '', font=self.top.font_Text)
        Step2bTxt.pack(side=LEFT, anchor=N)
        
        #--------------------------------------------------------------------------------
        
        fCropStep2c = Frame(self.fCropStep2, relief=RIDGE, border=0, width=400, height=40)
        fCropStep2c.pack(fill=X, expand=True, side=TOP)
        fCropStep2c.pack_propagate(0)
        
        Step2cPad = Label(fCropStep2c, text = '')
        Step2cPad.pack(side=LEFT, anchor=W, padx=1)
                        
        Step2c3Pad = Label(fCropStep2c, text = '')
        Step2c3Pad.pack(side=LEFT, anchor=E, padx=10)
        
        Step2cTxt = Label(fCropStep2c, text = 'Sphere(s) list:', height=5, font=self.top.font_Text)
        Step2cTxt.pack(side=LEFT, anchor=E)
        
        optionTuple = '',
        self.OptCleftStep2 = apply(OptionMenu, (fCropStep2c, self.Step2Selection) + optionTuple)
        self.OptCleftStep2.config(width=20, bg=self.top.Color_White)
        self.OptCleftStep2.pack(side=LEFT, anchor=E)

        #--------------------------------------------------------------------------------

        fCropStep2d = Frame(self.fCropStep2, relief=RIDGE, border=0, width=400, height=25)
        fCropStep2d.pack(fill=X, expand=True, side=TOP)
        fCropStep2d.pack_propagate(0)
                
        Step2dTxt = Label(fCropStep2d, text = '', height=5, font=self.top.font_Text)
        Step2dTxt.pack(side=LEFT, anchor=W)
        
        Step2d2Pad = Label(fCropStep2d, text = '')
        Step2d2Pad.pack(side=LEFT, anchor=E, padx=50)
        
        self.Step2d1Btn = Button(fCropStep2d, text='Add', width=5, command=self.Step2_Add)
        self.Step2d1Btn.pack(side=LEFT, anchor=E)         

        self.Step2d2Btn = Button(fCropStep2d, text='Delete', width=5, command=self.Step2_Delete)
        self.Step2d2Btn.pack(side=LEFT, anchor=E, padx=4)

        self.Step2d3Btn = Button(fCropStep2d, text='Edit', width=5, command=self.Step2_Edit)
        self.Step2d3Btn.pack(side=LEFT, anchor=E)
                        
        Step2dPad = Label(fCropStep2d, text = '')
        Step2dPad.pack(side=LEFT, anchor=E, padx=35)

        #--------------------------------------------------------------------------------
        
        fCropStep2e = Frame(self.fCropStep2, relief=RIDGE, border=0, width=400, height=65)
        fCropStep2e.pack(fill=X, expand=True, side=TOP)
        fCropStep2e.pack_propagate(0)       

        Step2ePad = Label(fCropStep2e, text = '')
        Step2ePad.pack(side=LEFT, anchor=S, padx=40)
        
        Step2eTxt = Label(fCropStep2e, text = 'Radius:', font=self.top.font_Text_I)
        Step2eTxt.pack(side=LEFT, anchor=S, pady=20)
        
        Step2e2Pad = Label(fCropStep2e, text = '')
        Step2e2Pad.pack(side=LEFT, anchor=W, padx=5)
        
        self.ResizeSphere = Scale(fCropStep2e, from_=0.0, to=0.0, resolution=self.ScaleResolution, variable=self.SphereSize, orient=HORIZONTAL, length=120, command=self.MoveResizeSphere, highlightthickness=0)
        self.ResizeSphere.pack(side=LEFT, anchor=NW, pady=5)
        self.ResizeSphere.config(state='disabled')
        
        Step2e2Pad = Label(fCropStep2e, text = '')
        Step2e2Pad.pack(side=RIGHT, anchor=E, padx=5)

        Step2eBtn = Button(fCropStep2e, width=20, height=20, image=self.Img_Button[1], command=self.Step2_Next)
        Step2eBtn.pack(side=RIGHT, anchor=E, padx=2)             

        Step2eBtn2 = Button(fCropStep2e, width=20, height=20, image=self.Img_Button[0], command=self.Step2_Prev)
        Step2eBtn2.pack(side=RIGHT, anchor=E, padx=2)             

        #--------------------------------------------------------------------------------

        #==================================================================================
        '''                --- Save the partition ---                       '''
        #==================================================================================

        fCropStep3a = Frame(self.fCropStep3, relief=RIDGE, border=0, width=400, height=30)
        fCropStep3a.pack(fill=X, expand=True, side=TOP)
        fCropStep3a.pack_propagate(0)

        Step3aTitle = Label(fCropStep3a, text = 'STEP 3:', height=5, font=self.top.font_Title)
        Step3aTitle.pack(side=LEFT, anchor=W, padx=5)
        
        Step3aTxt = Label(fCropStep3a, text = 'Name the newly created Cleft Partition', height=5, font=self.top.font_Text)
        Step3aTxt.pack(side=LEFT, anchor=W)
                
        #--------------------------------------------------------------------------------

        fCropStep3b = Frame(self.fCropStep3, relief=RIDGE, border=0, width=400, height=40)
        fCropStep3b.pack(fill=X, expand=True, side=TOP)
        fCropStep3b.pack_propagate(0)  

        Step3bPad = Label(fCropStep3b, text = '')
        Step3bPad.pack(side=LEFT, anchor=W, padx=8)
        
        self.Step3bChk = Checkbutton(fCropStep3b, text=' Default', variable=self.Step3Check, highlightthickness=0, font=self.top.font_Text)#, command=self.Step3_Checked )
        self.Step3bChk.pack(side=LEFT, anchor=W)
        self.Step3bChk.config(state='disabled')
        
        Step3b2Pad = Label(fCropStep3b, text = '')
        Step3b2Pad.pack(side=LEFT, anchor=E, padx=3)
                
        self.Step3bEntry = Entry(fCropStep3b, width=22, textvariable=self.Step3Output, disabledforeground = 'black', justify=CENTER, font=self.top.font_Text)
        self.Step3bEntry.pack(side=LEFT, anchor=E, padx=4)
        self.Step3bEntry.config(state='disabled')
        
        Step3b3Pad = Label(fCropStep3b, text = '')
        Step3b3Pad.pack(side=RIGHT, anchor=E, padx=5)

        Step3bBtn = Button(fCropStep3b, width=20, height=20, image=self.Img_Button[1], command=self.Step3_Next)
        Step3bBtn.pack(side=RIGHT, anchor=E, padx=2)             

        Step3bBtn2 = Button(fCropStep3b, width=20, height=20, image=self.Img_Button[0], command=self.Step3_Prev)
        Step3bBtn2.pack(side=RIGHT, anchor=E, padx=2)             


        #--------------------------------------------------------------------------------
        # Quit the partition options menu
        Btn_Back = Button(self.fCropBtn, text='Back', width=7, command=self.Btn_Back_Clicked, font=self.top.font_Text)
        Btn_Back.pack(side=RIGHT, anchor=E)


    ''' ==================================================================================
    FUNCTION Btn_Back_Clicked: Goes back to default frame
    ==================================================================================  '''                 
    def Btn_Back_Clicked(self):
        
        self.top.SetActiveFrame(self.top.Default)

    ''' ==================================================================================
    FUNCTION SphereRunning: Disables/Enables controls related to Sphere Wizard
    ==================================================================================  '''                 
    def SphereRunning(self, boolRun):

        if boolRun:
            self.ResizeSphere.config(state='normal')

        else:
            self.ResizeSphere.config(state='disabled')

            if not self.top.WizardError:

                self.SphereList[self.Sphere] = [ self.SphereSize.get(), self.SphereCoord ]
                
                self.update_Spheres()                
                self.write_Partition(self.Step1Selection.get())
                self.displayPartition(self.Step1Selection.get())

    ''' ==========================================================
    highlight_Step: highlights active step in Cropping
    ==========================================================='''           
    def highlight_Step(self, Step):

        if Step == 1:
            self.color_Step(self.fCropStep1, self.Color_Green)
            self.color_Step(self.fCropStep2, self.Color_Default)
            self.color_Step(self.fCropStep3, self.Color_Default)

        elif Step == 2:
            self.color_Step(self.fCropStep1, self.Color_Default)
            self.color_Step(self.fCropStep2, self.Color_Green)
            self.color_Step(self.fCropStep3, self.Color_Default)

        elif Step == 3:
            self.color_Step(self.fCropStep1, self.Color_Default)
            self.color_Step(self.fCropStep2, self.Color_Default)
            self.color_Step(self.fCropStep3, self.Color_Green)

    ''' ==========================================================
    activate_Step: activates the controls of the active step
    ==========================================================='''           
    def activate_Step(self, Step):

        if Step == 1:
            General.setState(self.fCropStep1, 'normal')
            General.setState(self.fCropStep2)
            General.setState(self.fCropStep3)
            self.Toggle_Step1()
            
        elif Step == 2:
            General.setState(self.fCropStep1)
            General.setState(self.fCropStep2, 'normal')
            General.setState(self.fCropStep3)
            self.Toggle_Step2()

        elif Step == 3:
            General.setState(self.fCropStep1)
            General.setState(self.fCropStep2)
            General.setState(self.fCropStep3, 'normal')
            self.Toggle_Step3()

    ''' ==========================================================
    color_Step: Colors recursively widgets
    ==========================================================='''           
    def color_Step(self, widget, bg):

        try:
            Class = str(widget.winfo_class())
 #Class != 'Scale' and 
            if Class != 'Button' and Class != 'Menubutton' and \
               Class != 'Entry' and \
               Class != 'Menu' and Class != 'OptionMenu':

                widget.config(bg=bg)
        except:
            pass
        for child in widget.winfo_children():
            self.color_Step(child, bg)

    ''' ==========================================================
    Flash: Flashes the selected cleft in the PyMOL Viewer
    ==========================================================='''           
    def Flash(self, Sel):

        for i in range(0, 8):
            General.Oscillate(Sel, 0.03)

    ''' ==========================================================
    Toggle_Step1: Disable/Enable the 'Next Step' Button
    ==========================================================='''           
    def Toggle_Step1(self, *args):
        
        Sel = self.Step1Selection.get()

        if Sel != '':
            if General.object_Exists(Sel):
                self.Flash(Sel)
                self.Step1bBtn.config(state='normal')
            else:
                self.DisplayMessage("The cleft '" + Sel + "'no longer exists.", 2)
                self.Step1bBtn.config(state='disabled')
                self.Step1Selection.set('')

        else:
            self.Step1bBtn.config(state='disabled')

        if self.Step1Selection.get() != self.LastStep1Selection.get():
            self.SphereList.clear()
            self.update_Spheres()

        self.LastStep1Selection.set(Sel)
        
    ''' ==========================================================
    Step1_Next: Enables the controls for step 2 and highlight it
    ==========================================================='''
    def Step1_Next(self):
        
        self.activate_Step(2)
        self.highlight_Step(2)
        self.Step = 2

        self.displaySphere(self.Step2Selection.get())

    ''' ==========================================================
    Step2_Next: Enables the controls for step 3 and highlight it
    ==========================================================='''
    def Step2_Next(self):
        
        if self.top.ActiveWizard is None:
            
            self.activate_Step(3)
            self.highlight_Step(3)
            self.Step = 3

            cmd.delete(self.SphereDisplay)

            self.Filename = self.top.Manage.update_Filename(self.top.dictTempClefts[self.Step1Selection.get()])
            self.Step3Output.set(self.Filename[1])

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step2_Prev: Enables the controls for step 1 and highlight it
    ==========================================================='''           
    def Step2_Prev(self):

        if self.top.ActiveWizard is None:
        
            self.activate_Step(1)
            self.highlight_Step(1)
            self.Step = 1

            cmd.delete(self.PartitionDisplay)
            cmd.delete(self.SphereDisplay)

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step3_Next: Enables the controls for step 3 and highlight it
    ==========================================================='''           
    def Step3_Next(self):
        
        Filename = os.path.join(self.Filename[0],self.Filename[1])
        Prefix = self.Filename[1].replace('.pdb','')

        self.top.Manage.copy_TempPartition(self.TempPartition, Filename)
        
        # No hard-link
        self.top.dictTempClefts[Prefix] = [ 1, 1, Filename, 0, 'NULL', self.Step1Selection.get(), False, 0.000 ]
        self.top.load_Clefts(self.top.Default.ColorList)

        # Put the partition cleft over its parent
        General.Oscillate(self.Step1Selection.get(), 0.0)

        self.top.Default.show_Spheres(True)
        self.top.Default.show_Atoms(True)
        
        self.SphereList.clear()
        self.Step1Selection.set('')
        self.Step2Selection.set('')
        self.Step3Output.set('')

        cmd.delete(self.PartitionDisplay)
        cmd.delete(self.SphereDisplay)

        # Go back to Step 1
        self.Step = 1
        self.activate_Step(1)
        self.highlight_Step(1)
        
        self.update_Step1_DDL()

    ''' ==========================================================
    Step3_Prev: Enables the controls for step 2 and highlight it
    ==========================================================='''           
    def Step3_Prev(self):

        self.activate_Step(2)
        self.highlight_Step(2)
        self.Step = 2

        self.displaySphere(self.Step2Selection.get())

    ''' ==========================================================
    update_Step1_DDL: Updates the Drop-Down-List of the Spheres objects
    ==========================================================='''           
    def update_Step1_DDL(self):

        self.OptCleftStep1['menu'].delete(0, 'end')
        
        for key in sorted(iter(self.top.dictTempClefts)):
            if self.top.dictTempClefts[key][0]:
                self.OptCleftStep1['menu'].add_command(label=key, command=lambda temp = key: self.OptCleftStep1.setvar(self.OptCleftStep1.cget("textvariable"), value = temp))

    ''' ==========================================================
    update_Spheres: Updates the Drop-Down-List of Spheres
    ==========================================================='''           
    def update_Spheres(self):

        self.OptCleftStep2['menu'].delete(0, 'end')

        Sphere = ''
        for key in sorted(iter(self.SphereList)):
            self.OptCleftStep2['menu'].add_command(label=key, command=lambda temp = key: self.OptCleftStep2.setvar(self.OptCleftStep2.cget("textvariable"), value = temp))
            Sphere = key
        
        self.Step2Selection.set(Sphere)

    ''' ==========================================================
    Step2_Add: Adds a new sphere for partitionning
    ==========================================================='''           
    def Step2_Add(self):

        if self.top.ActiveWizard is None:

            self.Sphere = 'SPHERE_' + str(len(self.SphereList) + 1)
            self.SphereRunning(True)
                
            self.top.ActiveWizard = Sphere.Sphere(self, self.Step1Selection.get(), 1.0, False,0,[])
            cmd.set_wizard(self.top.ActiveWizard)

            self.top.ActiveWizard.Start()

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step2_Delete: Deletes a sphere for partitionning
    ==========================================================='''           
    def Step2_Delete(self):

        if self.top.ActiveWizard is None:
            del self.SphereList[self.Step2Selection.get()]

            self.update_Spheres()
            self.write_Partition(self.Step1Selection.get())
            self.displayPartition(self.Step1Selection.get())

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step2_Edit: Edits a sphere for partitionning
    ==========================================================='''           
    def Step2_Edit(self):

        if self.top.ActiveWizard is None:
            
            self.Sphere = self.Step2Selection.get()
            self.SphereRunning(True)

            self.top.ActiveWizard = Sphere.Sphere(self, self.Step1Selection.get(), 1.0, False,
                                                  self.SphereList[self.Sphere][0],
                                                  self.SphereList[self.Sphere][1])

            cmd.set_wizard(self.top.ActiveWizard)

            self.top.ActiveWizard.Start()

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Toggle_Step2: Disable/Enable the 'Edit' Button
    ==========================================================='''           
    def Toggle_Step2(self, *args):

        if self.Step2Selection.get() != '':
            self.Step2d2Btn.config(state='normal')
            self.Step2d3Btn.config(state='normal')
        else:
            self.Step2d2Btn.config(state='disabled')
            self.Step2d3Btn.config(state='disabled')

        self.displaySphere(self.Step2Selection.get())

    ''' ==========================================================
    Toggle_Step3: Disable/Enable the 'Ouput entry' Button
    ==========================================================='''           
    def Toggle_Step3(self, *args):

        if self.Step3Check.get():
            self.Filename = self.top.Manage.update_Filename(self.top.dictTempClefts[self.Step1Selection.get()][2])
            self.Step3Output.set(self.Filename[1])

        else:
            FilePath = tkFileDialog.asksaveasfilename(filetypes=[('GetCleft Spheres file','*.pdb')], initialdir=self.Filename[0], title='Choose the Cleft Spheres filename to Save', initialfile = self.Filename[1])

            if len(FilePath) > 0:
                self.Filename = os.path.split(FilePath)
                self.Step3Output.set(self.Filename[1])
            else:
                self.Step3Check.set(1)

        self.Step3bEntry.config(state='disabled')

    # Shows the sphere in transparency
    def displaySphere(self, Sel):

        cmd.delete(self.SphereDisplay)

        if Sel != '':
            # Only display sphere in current state
            cmd.pseudoatom(self.SphereDisplay, pos=self.SphereList[Sel][1], vdw=self.SphereList[Sel][0], color='purpleblue')#, state=cmd.get_state())
            cmd.hide('nonbonded',self.SphereDisplay)
            cmd.show('spheres', self.SphereDisplay)
            cmd.set('sphere_transparency', 0.7, self.SphereDisplay)
            cmd.rebuild()

    # Shows the sphere in transparency
    def displayPartition(self, Sel):

        View = cmd.get_view()

        cmd.delete(self.PartitionDisplay)

        cmd.load(self.TempPartition, self.PartitionDisplay)
        cmd.hide('everything', self.PartitionDisplay)
        cmd.color('grey60', self.PartitionDisplay)
        cmd.show('surface', self.PartitionDisplay)

        cmd.set_view(View)

        General.Oscillate(Sel, 0.0)

        
    # MoveResizeSphere = Get the scale value, then resize the sphere          
    def MoveResizeSphere(self, arg):

        if not self.top.ActiveWizard is None:
            self.top.ActiveWizard.resize()


    ''' ==================================================================================
    FUNCTION write_Partition: Writes the new Partition cleft files
    ==================================================================================  '''        
    def write_Partition(self, key):
        
        listNoAtom = list()
        FromFile = self.top.dictTempClefts[key][2]

        file = open(FromFile, 'r')
        SPHLines = file.readlines()
        file.close()
        
        # Write in the PDB file
        TMPFile = open(self.TempPartition, 'w')
        TMPFile.write('#REMARK  PARENTFILE  ' + FromFile + '\n')  

        for Line in SPHLines:
            if Line.startswith('ATOM  '):

                index = Line[7:11].strip()              # The atom number
                coordX = float(Line[30:38].strip())     # The atom X coordinate
                coordY = float(Line[39:46].strip())     # The atom Y coordinate
                coordZ = float(Line[47:54].strip())     # The atom Z coordinate
                radius = float(Line[61:].strip())       # The radius of the sphere                    
                    
                if listNoAtom.count(index) == 0:
                    
                    for sph in sorted(iter(self.SphereList)):
                        
                        sqrrad  = self.SphereList[sph][0]**2
                        sqrdist = Geometry.sqrdistance( [ coordX, coordY, coordZ ], self.SphereList[sph][1] )

                        #if sqrdist <= (self.SphereList[sph][0]+radius)*(self.SphereList[sph][0]+radius):

                        # ** NEW
                        # The center of the sphere needs to be inside the 'inserted Spheres'
                        if sqrdist <= sqrrad:

                            listNoAtom.append(index)
                            TMPFile.write(Line)
                            break
                                    
        TMPFile.close()
