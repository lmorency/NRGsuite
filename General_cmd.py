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
from pymol import cmd

import pymol
import time

''' ==========================================================
  is_ATOM: Determines whether a residue is an ATOM or HETATM  
==========================================================='''
def is_ATOM(Residue,Prot):
    
    resn = Residue[0:3]
    resi = Residue[3:len(Residue)-1]
  
    try:
        cmd.select("rtmp1__", "resi " + resi + " & resn " + resn + " & " + Prot)
        cmd.select("rtmp2__", "resi " + resi + " & resn " + resn + " & " + Prot + " & not hetatm")

        if cmd.count_atoms("rtmp1__") == cmd.count_atoms("rtmp2__"):
            cmd.delete("rtmp*__")
            return 1

        cmd.delete("rtmp*__")
        
    except:
        return -1

    return 0
    

''' ==========================================================
  mask_Objects: Mask ALL objects in list except those in exclude list  
==========================================================='''
def mask_Objects(lstexc):
    
    lstobj = cmd.get_names('objects')

    for obj in lstobj:
        if lstexc.count(str(obj)) == 0:
            try:
                type = cmd.get_type(str(obj))
                
                if type == 'object:molecule' or type == 'selection':
                    cmd.mask(str(obj))
            except:
                continue

''' ==========================================================
  unmask_Objects: UnMask ALL objects in list except those in exclude list  
==========================================================='''
def unmask_Objects(lstexc):
    
    lstobj = cmd.get_names('objects')
    
    for obj in lstobj:
        if lstexc.count(str(obj)) == 0:
            try:
                type = cmd.get_type(str(obj))
                
                if type == 'object:molecule' or type == 'selection':
                    cmd.mask(str(obj))
            except:
                continue

''' ==========================================================
object_Exists: Check in the list of selection/object if the object exists  
==========================================================='''           
def object_Exists(object):
    
    found = False
    sess = cmd.get_names('objects')
    for obj in sess:
        if str(obj).upper() == str(object).upper():
            found = True
            break
        
    return found

''' ==========================================================
get_config_mouse: gets the configuration of the mouse
==========================================================='''           
def get_config_mouse():
    
    config_mouse = ''

    try:
        name = cmd.get("button_mode_name")
        
        if name[0] == '1':
            config_mouse += 'one'
        elif name[0] == '2':
            config_mouse += 'two'
        elif name[0] == '3':
            config_mouse += 'three'

        config_mouse += '_button'

        if name[0] != '1':
            if name[9:] == 'Viewing':
                config_mouse += '_viewing'
            elif name[9:] == 'Editing':
                config_mouse += '_editing'
            elif name[9:] == 'Motions':
                config_mouse += '_motions'

        return config_mouse

    except:
        return 'three_button_viewing'

''' ==========================================================
Refresh_DDL: Refresh drop-down-list with all objects/selections
==========================================================='''               
def Refresh_DDL(widget, var, exc, fun):

    try:
        cmd.unpick()
        list = cmd.get_names('all')
    
        # Delete all entries in DDL
        widget['menu'].delete(0, END)

        if len(list) > 0:
            for item in sorted(list, key=str.lower):
                type = cmd.get_type(str(item))
                
                if (type == 'object:molecule' or type == 'selection') and exc.count(str(item)) == 0:
                    if fun != None:
                        widget['menu'].add_command(label=item, command=lambda temp = item: fun(temp))
                    else:
                        widget['menu'].add_command(label=item, command=lambda temp = item: widget.setvar(widget.cget('textvariable'), value = temp))

                    var.set(item)
                        
        else:            
            # Dummy value
            widget['menu'].add_command(label='', command=lambda temp = '': widget.setvar(widget.cget('textvariable'), value = ''))
            var.set('')

    except:
        pass

''' ==========================================================
Oscillate: Oscillates an object between 2 states
==========================================================='''           
def Oscillate(sel, interval):

    # Determine if selection is enabled/disabled
    Visible = 0
    if cmd.count_atoms('v. & ' + sel) > 0:
        Visible = 1
    
    if Visible:
        cmd.disable(sel)
        time.sleep(interval)
        cmd.enable(sel)
        time.sleep(interval)

    else:
        cmd.enable(sel)
        time.sleep(interval)
        cmd.disable(sel)
        time.sleep(interval)

''' ==================================================================================
FUNCTION repeat: repeats a character N number of times
==================================================================================  '''            
def repeat(string, length):
    L = len(string)
    return string * (length // L) + string[:length % L]


''' ==================================================================================
FUNCTION get_ID: retrieves atom number from an index and selection
==================================================================================  '''            
def get_ID(index, sel):
        
    return int(cmd.id_atom('index ' + str(index) + ' & ' + sel + ' & present'))

''' ==================================================================================
FUNCTION Get_CenterOfMass2: Gets the center of mass of an object
==================================================================================  '''            
def Get_CenterOfMass2(selection, state):

    try:
        nAtoms = cmd.count_atoms(selection, state=state)

        if nAtoms > 0:
            MinMax = cmd.get_extent(selection, state=state)

            MaxWidth = 0.0
            for i in range(0,3):
                if (MinMax[1][i]-MinMax[0][i]) > MaxWidth:
                    MaxWidth = (MinMax[1][i]-MinMax[0][i])

            return [ (MinMax[1][0]+MinMax[0][0])/2.0,
                     (MinMax[1][1]+MinMax[0][1])/2.0,
                     (MinMax[1][2]+MinMax[0][2])/2.0 ]

    except:
        return []

''' ==================================================================================
FUNCTION Get_MaxWidth: Returns the max width from the center of an object
==================================================================================  '''            
def Get_MaxWidth(selection,state):

    try:
        nAtoms = cmd.count_atoms(selection, state=state)

        if nAtoms > 0:
            MinMax = cmd.get_extent(selection, state=state)

            MaxWidth = 0.0
            for i in range(0,3):
                if (MinMax[1][i]-MinMax[0][i]) > MaxWidth:
                    MaxWidth = (MinMax[1][i]-MinMax[0][i])

            return MaxWidth

    except:
        return -1
