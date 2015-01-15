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
from pymol.wizard import Wizard
from pymol import cmd

import pymol
import SphereObj
import General_cmd

class Sphere(Wizard):

    def __init__(self, top, queue, Sphere, SphereDisplay, SphereSize, ExtraPrompt):

        #print "New instance of flexsphere Wizard"

        Wizard.__init__(self)

        self.top = top
        self.queue = queue
        self.App = self.top.top

        self.App.WizardError = False
        self.App.WizardResult = 0

        self.Sphere = Sphere
        self.SphereView = self.Sphere.Copy()

        self.SphereSize = SphereSize
        self.SphereDisplay = SphereDisplay

        self.ErrorCode = 0

        self.Prompt = [ "Press Shift + Mouse3 (Wheel Click) to move the sphere.",
                        "Use the scaler in the " + self.App.Name + " interface to edit its radius." ]
                        
        if ExtraPrompt:
            self.Prompt.append( ExtraPrompt )

    #=======================================================================
    ''' Execute the first steps of the wizard '''
    #=======================================================================    
    def Start(self):

        cmd.window('hide')
        self.queue.put(lambda: cmd.window('show'))
        cmd.refresh_wizard()

        self.ErrorCode = 1

        try:
            self.State = cmd.get_state()
            self.config_mouse = General_cmd.get_config_mouse()
            cmd.config_mouse('three_button_editing')
        
            self.exc = [ self.SphereDisplay ]
            General_cmd.mask_Objects(self.exc)

            # success
            self.ErrorCode = 0

        except:
            self.queue.put(lambda: self.App.DisplayMessage("  ERROR: Could not start the Sphere wizard", 1))
            self.queue.put(lambda: self.App.DisplayMessage("         The wizard will abort prematurely", 1))
            self.Quit_Wizard()
            return
            
        if self.DisplaySphere():
            self.queue.put(lambda: self.App.DisplayMessage("  ERROR: Could not display the Sphere", 1))
            self.queue.put(lambda: self.App.DisplayMessage("         The wizard will abort prematurely", 1))
            self.Quit_Wizard()
            return
            
    #=======================================================================
    ''' Quits the wizard '''
    #=======================================================================    
    def Quit_Wizard(self):

        try:
            # Delete Sphere object
            General_cmd.unmask_Objects(self.exc)
            cmd.config_mouse(self.config_mouse)

            cmd.delete(self.SphereDisplay)
            cmd.refresh()
        except:
            pass

        # Catch error in App
        if self.ErrorCode > 0:
            self.App.WizardError = True

        # Re-enable controls
        self.queue.put(lambda: self.top.SphereRunning(False))
        self.queue.put(lambda: self.top.top.root.lift())
        self.App.ActiveWizard = None

        cmd.set_wizard()
        cmd.refresh()
        
    #=======================================================================
    ''' Display the Sphere in Pymol '''
    #=======================================================================    
    def DisplaySphere(self):
   
        try:
            # Display the Sphere
            cmd.delete(self.SphereDisplay)
            cmd.refresh()
        except:
            pass

        try:
            cmd.pseudoatom(self.SphereDisplay,
                           pos=self.SphereView.Center,
                           vdw=self.SphereView.Radius,
                           state=self.State)
            cmd.refresh()

            cmd.color('oxygen', self.SphereDisplay)
            cmd.refresh()

            cmd.hide('everything', self.SphereDisplay)
            cmd.refresh()

            cmd.show('spheres', self.SphereDisplay)
            cmd.refresh()
            
        except:
            self.ErrorCode = 1

        return self.ErrorCode
    
    #=======================================================================
    ''' Resize the sphere '''
    #=======================================================================      
    def ResizeSphere(self):
        
        self.SphereView.Set_Radius(self.SphereSize.get())

        try:
            cmd.alter(self.SphereDisplay,'vdw=' + str(self.SphereView.Radius))
            cmd.rebuild(self.SphereDisplay)
            cmd.refresh()
        except:
            self.queue.put(lambda: self.App.DisplayMessage("  ERROR: Could not resize the Sphere", 1))
            self.queue.put(lambda: self.App.DisplayMessage("         The wizard will abort prematurely", 1))
            self.Quit_Wizard()
    
    def get_prompt(self):
        
        return self.Prompt

    def reset(self):

        self.SphereView = self.Sphere.Copy()
        self.SphereSize.set(self.SphereView.Radius)

        if self.DisplaySphere():
            self.queue.put(lambda: self.App.DisplayMessage("  ERROR: Could not display the Sphere", 1))
            self.queue.put(lambda: self.App.DisplayMessage("         The wizard will abort prematurely", 1))
            self.Quit_Wizard()
            return

    def cancel(self):

        self.App.WizardResult = 1
        self.Quit_Wizard()
            
    def btn_Done(self):

        Center = General_cmd.Get_CenterOfMass2(self.SphereDisplay, self.State)
        if len(Center) > 0:
            self.SphereView.Set_Center(Center)

        self.Sphere.Set_Radius(self.SphereView.Radius)
        self.Sphere.Set_Center(self.SphereView.Center)

        self.App.WizardResult = 2
        self.Quit_Wizard()

    def get_panel(self):
        return [
        [ 1, '* Sphere Controls *',''],
        [ 2, 'Reset','cmd.get_wizard().reset()'],
        [ 2, 'Cancel','cmd.get_wizard().cancel()'],
        [ 2, 'Done','cmd.get_wizard().btn_Done()'],         
        ]
         
