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
import Tabs
import time
import tkMessageBox

if __debug__:
    from pymol import cmd
    
    import Sphere
    import SphereObj
    import CleftObj
    
    import General_cmd

#=========================================================================================
'''                        ---  GETCLEFT's CROPCLEFT FRAME  ---                        '''
#=========================================================================================     
class CropCleft(Tabs.Tab):

    def Def_Vars(self):
        
        self.Img_Button = list()
        self.Filename = list()
        
        self.Step1Selection = StringVar()
        self.LastStep1Selection = StringVar()
        self.Step2Selection = StringVar()
        self.Step3Output = StringVar()
        self.Step3Check = IntVar()
        self.dictSpheres = dict()
        self.SphereSize = DoubleVar()
        self.SphereCoord = list()

    def Init_Vars(self):
        
        self.Step = 1
        self.Filename = []

        del self.Img_Button[:]
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/stop.gif')))
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/down.gif')))
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/up.gif')))
        self.Img_Button.append(PhotoImage(file=os.path.join(self.top.GetCleftInstall_Dir,'images/ok.gif')))

        self.Step1Selection.set('')
        self.LastStep1Selection.set('')
        self.Step2Selection.set('')
        self.Step3Output.set('')
        self.Step3Check.set(1)

        self.Color_Default = '#d9d9d9'
        self.Color_Green = self.top.Color_Green
        self.TempPartition = os.path.join(self.top.GetCleftProject_Dir,'tmppart.pdb')
 
        self.ScaleResolution = 0.25
        self.PartitionDisplay = 'PARTITION_AREA__'
        
        self.SphereDisplay = 'SPHERE_PT_AREA__'
        self.SphereCoord = []
        self.SphereSize.set(0.0)

    
    ''' ==================================================================================
    FUNCTION Trace: Adds a callback function to StringVars
    ==================================================================================  '''  
    def Trace(self):

        try:
            self.Step1SelectionTrace = self.Step1Selection.trace('w', self.Toggle_Step1)
            self.Step2SelectionTrace = self.Step2Selection.trace('w', self.Toggle_Step2)
        except:
            pass
        
    ''' ==================================================================================
    FUNCTION Del_Trace: Deletes observer callbacks
    ==================================================================================  '''  
    def Del_Trace(self):

        try:
            self.Step1Selection.trace_vdelete('w',self.Step1SelectionTrace)
            self.Step2Selection.trace_vdelete('w',self.Step2SelectionTrace)
        except:
            pass

    ''' ==================================================================================
    FUNCTION After_Show: Displays the frame
    ==================================================================================  '''  
    def After_Show(self):
        
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
              
        self.fCrop = Frame(self.top.fMiddle)
                        
        #==================================================================================
        '''            --- Selection of the Cleft to be partioned ---                   '''
        #==================================================================================
        
        self.fCropStep1 = Frame(self.fCrop, relief=SUNKEN, border=1)
        self.fCropStep1.pack(fill=X, expand=True, side=TOP)        
        fCropStep1Line1 = Frame(self.fCropStep1, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep1Line1.pack(fill=X, expand=True, side=TOP)
        fCropStep1Line2 = Frame(self.fCropStep1, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep1Line2.pack(fill=X, expand=True, side=TOP)
        fCropStep1Line3 = Frame(self.fCropStep1, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep1Line3.pack(fill=X, expand=True, side=TOP)
        
        Label(fCropStep1Line1, text = 'STEP 1: Select a parent cleft to partition', font=self.top.font_Title).pack(side=LEFT)
        Label(fCropStep1Line2, text = 'Cleft objects:', width=15, font=self.top.font_Text).pack(side=LEFT)

        optionTuple = '',
        self.OptCleftStep1 = apply(OptionMenu, (fCropStep1Line2, self.Step1Selection) + optionTuple)
        self.OptCleftStep1.config(width=20, bg=self.top.Color_White)
        self.OptCleftStep1.pack(side=LEFT, anchor=W)

        Button(fCropStep1Line2, text='Refresh', command=self.update_Step1_DDL).pack(side=RIGHT, anchor=E)
                        
        self.Step1bBtn = Button(fCropStep1Line3, width=20, height=20, image=self.Img_Button[1], command=self.Step1_Next)
        self.Step1bBtn.pack(side=RIGHT, anchor=E, padx=2)
        self.Step1bBtn.config(state='disabled')
                
        
        #--------------------------------------------------------------------------------

        #==================================================================================
        '''            --- Add spheres to partition the cleft ---                       '''
        #==================================================================================
        
        self.fCropStep2 = Frame(self.fCrop, relief=SUNKEN, border=1)
        self.fCropStep2.pack(fill=X, expand=True, side=TOP)        
        fCropStep2Line1 = Frame(self.fCropStep2, relief=SUNKEN, border=0, pady=8, padx=10)
        fCropStep2Line1.pack(fill=X, expand=True, side=TOP)
        fCropStep2Line2 = Frame(self.fCropStep2, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep2Line2.pack(fill=X, expand=True, side=TOP)
        fCropStep2Line3 = Frame(self.fCropStep2, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep2Line3.pack(fill=X, expand=True, side=TOP)
        
        Label(fCropStep2Line1, text = 'STEP 2: Add one or more spheres to partition', font=self.top.font_Title).pack(side=LEFT)
        Label(fCropStep2Line2, text = 'Sphere objects:', width=15,font=self.top.font_Text).pack(side=LEFT)

        Step2eBtn2 = Button(fCropStep2Line1, width=20, height=20, image=self.Img_Button[2], command=self.Step2_Prev)
        Step2eBtn2.pack(side=RIGHT, anchor=E, padx=2)
        
        optionTuple = '',
        self.OptCleftStep2 = apply(OptionMenu, (fCropStep2Line2, self.Step2Selection) + optionTuple)
        self.OptCleftStep2.config(width=20, bg=self.top.Color_White)
        self.OptCleftStep2.pack(side=LEFT, anchor=W)

        self.Step2d3Btn = Button(fCropStep2Line2, text='Edit', command=self.Step2_Edit)
        self.Step2d3Btn.pack(side=RIGHT, anchor=E)
        self.Step2d2Btn = Button(fCropStep2Line2, text='Delete', command=self.Step2_Delete)
        self.Step2d2Btn.pack(side=RIGHT, anchor=E, padx=3)
        self.Step2d1Btn = Button(fCropStep2Line2, text='Add', command=self.Step2_Add)
        self.Step2d1Btn.pack(side=RIGHT, anchor=E)         

        #--------------------------------------------------------------------------------

        Label(fCropStep2Line3, text = 'Radius:', width=15, font=self.top.font_Text).pack(side=LEFT)
        self.ResizeSphere = Scale(fCropStep2Line3, from_=0.0, to=0.0, resolution=self.ScaleResolution, variable=self.SphereSize, orient=HORIZONTAL, length=120, command=self.MoveResizeSphere, highlightthickness=0)
        self.ResizeSphere.pack(side=LEFT, anchor=NW)
        self.ResizeSphere.config(state='disabled')

        Step2eBtn = Button(fCropStep2Line3, width=20, height=20, image=self.Img_Button[1], command=self.Step2_Next)
        Step2eBtn.pack(side=RIGHT, anchor=E, padx=2)

        #--------------------------------------------------------------------------------
        
        #==================================================================================
        '''                --- Create the partition ---                       '''
        #==================================================================================

        self.fCropStep3 = Frame(self.fCrop, relief=SUNKEN, border=1)
        self.fCropStep3.pack(fill=X, expand=True, side=TOP)        
        fCropStep3Line1 = Frame(self.fCropStep3, relief=SUNKEN, border=0, pady=8, padx=10)
        fCropStep3Line1.pack(fill=X, expand=True, side=TOP)
        fCropStep3Line2 = Frame(self.fCropStep3, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep3Line2.pack(fill=X, expand=True, side=TOP)
        fCropStep3Line3 = Frame(self.fCropStep3, relief=SUNKEN, border=0, pady=3, padx=10)
        fCropStep3Line3.pack(fill=X, expand=True, side=TOP)

        Label(fCropStep3Line1, text = 'STEP 3: Name the partitionned child cleft', font=self.top.font_Title).pack(side=LEFT)

        Step3bBtn2 = Button(fCropStep3Line1, width=20, height=20, image=self.Img_Button[2], command=self.Step3_Prev)
        Step3bBtn2.pack(side=RIGHT, anchor=E, padx=2)
        
        Label(fCropStep3Line2, text='Cleft name:', font=self.top.font_Text, width=15).pack(side=LEFT)
                
        self.Step3bEntry = Entry(fCropStep3Line2, width=25,textvariable=self.Step3Output, disabledforeground = 'black', justify=CENTER, font=self.top.font_Text)
        self.Step3bEntry.pack(side=LEFT, anchor=E)
        self.Step3bEntry.config(state='disabled')

        self.Step3bBtnCreate = Button(fCropStep3Line2, disabledforeground = 'black', text='Create', justify=CENTER, font=self.top.font_Text, command=self.Btn_CreatePartition)
        self.Step3bBtnCreate.pack(side=RIGHT, anchor=E)

        return self.fCrop

    ''' ==================================================================================
    FUNCTION SphereRunning: Disables/Enables controls related to Sphere Wizard
    ==================================================================================  '''                 
    def SphereRunning(self, boolRun):

        if boolRun:

            try:
                cmd.delete(self.SphereDisplay)
            except:
                pass

            self.ResizeSphere.config(state='normal')
            
            self.ResizeSphere.config(from_=0.0,to=self.Sphere.MaxRadius)
            self.SphereSize.set(self.Sphere.Radius)

        else:
            self.ResizeSphere.config(state='disabled')
            
            if not self.top.WizardError and self.top.WizardResult == 2:
                
                self.dictSpheres[self.SphereID] = self.Sphere.Copy()
                
                self.update_Spheres()
                self.Vertex = self.write_Partition()
                self.displayPartition()

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

    ''' ==========================================================
    color_Step: Colors recursively widgets
    ==========================================================='''           
    def color_Step(self, widget, bg):

        try:
            Class = str(widget.winfo_class())
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
            General_cmd.Oscillate(Sel, 0.03)

    ''' ==========================================================
    Toggle_Step1: Disable/Enable the 'Next Step' Button
    ==========================================================='''           
    def Toggle_Step1(self, *args):
        
        Sel = self.Step1Selection.get()

        if Sel != '':
            if General_cmd.object_Exists(Sel):
                self.Flash(Sel)
                self.Step1bBtn.config(state='normal')
                
                self.Center = General_cmd.Get_CenterOfMass2(Sel, state=1)
                self.Width = General_cmd.Get_MaxWidth(Sel, state=1)
                self.Cleft = self.top.Default.TempBindingSite.Get_CleftName(Sel)

            else:
                self.DisplayMessage("The cleft '" + Sel + "'no longer exists.", 2)
                self.Step1bBtn.config(state='disabled')
                
                self.update_Step1_DDL()                
                self.Step1Selection.set('')

        else:
            self.Step1bBtn.config(state='disabled')

        if self.Step1Selection.get() != self.LastStep1Selection.get():
            self.dictSpheres.clear()
            self.update_Spheres()

        self.LastStep1Selection.set(Sel)
        
    ''' ==========================================================
    Step1_Next: Enables the controls for step 2 and highlight it
    ==========================================================='''
    def Step1_Next(self):
        
        self.activate_Step(2)
        self.highlight_Step(2)
        self.Step = 2

        self.displaySphere()
        
        self.update_Spheres()
        self.Vertex = self.write_Partition()
        self.displayPartition()

    ''' ==========================================================
    Step2_Next: Enables the controls for step 3 and highlight it
    ==========================================================='''
    def Step2_Next(self):
        
        if self.top.ActiveWizard is None:
            
            self.activate_Step(3)
            self.highlight_Step(3)
            self.Step = 3

            try:
                cmd.delete(self.SphereDisplay)
                #cmd.delete(self.PartitionDisplay)
            except:
                pass

            CleftName = os.path.basename(os.path.splitext(self.Cleft.CleftName)[0]) + '_pt'
            self.Step3Output.set(CleftName)

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
    def Btn_CreatePartition(self):

        DestFile = ''
        Output = self.Step3Output.get()
        
        if Output != '':
            
            if not self.Vertex:
                self.top.DisplayMessage("  ERROR: Cannot save an empty cleft. Add spheres that include volume of the parent.", 1)
                return
            
            if General_cmd.object_Exists(Output):
                if tkMessageBox.askquestion("Question", message="An object with that name already exists. Would you like to overwrite it?",icon='warning') == 'no':
                    return
                    
            try:
                View = cmd.get_view()
                
                DestFile = os.path.join(self.top.GetCleftTempProject_Dir, Output + '.pdb')
                self.top.Manage.copy_TempPartition(self.TempPartition, DestFile)
                
                cmd.load(DestFile, Output, format='pdb', state=1)
                cmd.hide('everything', Output)
                cmd.show('surface', Output)
                
                self.top.Default.SetPartitionColor(self.Cleft.CleftName)
                cmd.color('partition', Output)
            
                cmd.delete(self.PartitionDisplay)
                cmd.set_view(View)            
                
            except:
                self.top.DisplayMessage("  ERROR: Could not create partition object: File not found.", 1)
                return

        else:
            self.top.DisplayMessage("  ERROR: Could not create partition object: Output is null.", 1)
            return
            
        
        Cleft = CleftObj.CleftObj()
        Cleft.CleftFile = DestFile
        Cleft.CleftName = self.Step3Output.get()
        Cleft.PartitionParent = self.Cleft
        Cleft.Partition = True
        Cleft.UTarget = self.Cleft.UTarget
        Cleft.Set_CleftMD5()
                
        self.top.Default.TempBindingSite.Add_Cleft(Cleft)
        self.top.CopySession = False

        # Put the partition cleft over its parent
        General_cmd.Oscillate(self.Cleft.CleftName, 0.0)
        
        # reset everything
        self.dictSpheres.clear()
        self.Step1Selection.set('')
        self.Step2Selection.set('')
        self.Step3Output.set('')

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

        self.displaySphere()
        
        self.update_Spheres()
        self.Vertex = self.write_Partition()
        self.displayPartition()

    ''' ==========================================================
    update_Step1_DDL: Updates the Drop-Down-List of the Spheres objects
    ==========================================================='''           
    def update_Step1_DDL(self):

        self.top.Default.Update_TempBindingSite()
        
        self.OptCleftStep1['menu'].delete(0, 'end')
        
        for CleftName in self.top.Default.TempBindingSite.Get_SortedCleftNames():
            self.OptCleftStep1['menu'].add_command(label=CleftName, command=lambda temp = CleftName: self.OptCleftStep1.setvar(self.OptCleftStep1.cget("textvariable"), value = temp))

    ''' ==========================================================
    update_Spheres: Updates the Drop-Down-List of Spheres
    ==========================================================='''           
    def update_Spheres(self):

        self.OptCleftStep2['menu'].delete(0, 'end')

        self.Step2Selection.set('')
        for key in sorted(self.dictSpheres.keys(), key=str.lower):
            self.OptCleftStep2['menu'].add_command(label=key, command=lambda temp = key: self.OptCleftStep2.setvar(self.OptCleftStep2.cget("textvariable"), value = temp))
            self.Step2Selection.set(key)
        
    ''' ==========================================================
    Get_SphereID: Returns a valid SphereID in dictionary
    ==========================================================='''           
    def Get_SphereID(self):

        ID = 1
        while self.dictSpheres.get('SPHERE_' + str(ID), ''):
            ID = ID + 1

        return 'SPHERE_' + str(ID)
        
    ''' ==========================================================
    Step2_Add: Adds a new sphere for partitionning
    ==========================================================='''           
    def Step2_Add(self):

        if not self.PyMOL:
            return

        if self.top.ActiveWizard is None:
            
            self.SphereID = self.Get_SphereID()
            self.Sphere = SphereObj.SphereObj(self.Width/1.5, self.Width, self.Center)
            
            self.SphereRunning(True)

            self.top.ActiveWizard = Sphere.Sphere(self, self.Sphere, self.SphereID, self.SphereSize)
            cmd.set_wizard(self.top.ActiveWizard)
            self.top.ActiveWizard.Start()

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step2_Delete: Deletes a sphere for partitionning
    ==========================================================='''           
    def Step2_Delete(self):

        if self.top.ActiveWizard is None:
            del self.dictSpheres[self.Step2Selection.get()]

            try:
                cmd.delete(self.SphereDisplay)
            except:
                pass
                
            self.update_Spheres()
            self.Vertex = self.write_Partition()
            self.displayPartition()

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Step2_Edit: Edits a sphere for partitionning
    ==========================================================='''           
    def Step2_Edit(self):

        if self.top.ActiveWizard is None and self.Sphere is not None:
            
            self.SphereRunning(True)
            
            self.top.ActiveWizard = Sphere.Sphere(self, self.Sphere, self.Step2Selection.get(), self.SphereSize)
            cmd.set_wizard(self.top.ActiveWizard)
            self.top.ActiveWizard.Start()

        else:
            self.top.DisplayMessage("Could not execute task: A Wizard is active")

    ''' ==========================================================
    Toggle_Step2: Disable/Enable the 'Edit' Button
    ==========================================================='''           
    def Toggle_Step2(self, *args):

        Sel = self.Step2Selection.get()
        
        if Sel != '':
            self.Step2d2Btn.config(state='normal')
            self.Step2d3Btn.config(state='normal')
            
            self.Sphere = self.dictSpheres[Sel]
        else:
            self.Step2d2Btn.config(state='disabled')
            self.Step2d3Btn.config(state='disabled')
            
            self.Sphere = None

        self.Step3Output.set('')
        self.displaySphere()
        

    # Shows the sphere in transparency
    def displaySphere(self):
        
        try:
            cmd.delete(self.SphereDisplay)
        except:
            pass
        
        if self.Sphere != None:
            # Only display sphere in current state
            cmd.pseudoatom(self.SphereDisplay, pos=self.Sphere.Center, vdw=self.Sphere.Radius, color='purpleblue', state=cmd.get_state())
            cmd.hide('nonbonded',self.SphereDisplay)
            cmd.show('spheres', self.SphereDisplay)
            cmd.set('sphere_transparency', 0.7, self.SphereDisplay)
            cmd.rebuild()

    # Shows the sphere in transparency
    def displayPartition(self):

        View = cmd.get_view()

        try:
            cmd.delete(self.PartitionDisplay)
        except:
            pass
            
        try:
            cmd.load(self.TempPartition, self.PartitionDisplay, format='pdb')
            cmd.hide('everything', self.PartitionDisplay)
            cmd.color('grey60', self.PartitionDisplay)
            cmd.show('surface', self.PartitionDisplay)
        except:
            self.top.DisplayMessage("  ERROR: An error occured while displaying the partitionned cleft", 1)
            return
            
        cmd.set_view(View)

        General_cmd.Oscillate(self.Cleft.CleftName, 0.0)

        
    # MoveResizeSphere = Get the scale value, then resize the sphere          
    def MoveResizeSphere(self, arg):

        if self.top.ActiveWizard is not None:
            self.top.ActiveWizard.ResizeSphere()


    ''' ==================================================================================
    FUNCTION write_Partition: Writes the new Partition cleft files
    ==================================================================================  '''        
    def write_Partition(self):
        
        listNoAtom = list()
        FromFile = self.Cleft.CleftFile

        try:
            file = open(FromFile, 'r')
            SPHLines = file.readlines()
            file.close()
        except:
            self.top.DisplayMessage("  ERROR: Could not find the parent cleft file.", 1)
            return
            
        # Write in the PDB file
        TMPFile = open(self.TempPartition, 'w')
        TMPFile.write('REMARK  PARENTFILE  ' + FromFile + '\n')            

        Vertex = 0
        for Line in SPHLines:
            if Line.startswith('ATOM  '):

                index = Line[7:11].strip()              # The atom number
                coordX = float(Line[30:38].strip())     # The atom X coordinate
                coordY = float(Line[39:46].strip())     # The atom Y coordinate
                coordZ = float(Line[47:54].strip())     # The atom Z coordinate
                radius = float(Line[61:].strip())       # The radius of the sphere                    
                    
                if listNoAtom.count(index) == 0:
                    
                    for sph in sorted(self.dictSpheres, key=str.lower):
                        
                        sqrrad  = self.dictSpheres[sph].Radius ** 2
                        sqrdist = Geometry.sqrdistance( [ coordX, coordY, coordZ ], self.dictSpheres[sph].Center )

                        #if sqrdist <= (self.dictSpheres[sph].Radius+radius)*(self.dictSpheres[sph].Radius+radius):

                        # ** NEW
                        # The center of the sphere needs to be inside the 'inserted Spheres'
                        if sqrdist <= sqrrad:

                            listNoAtom.append(index)
                            TMPFile.write(Line)
                            Vertex = Vertex + 1
                            break
                                    
        TMPFile.close()

        return Vertex
        
