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

''' ==========================================================
  validate_Float: validate the entry field (float value)  
==========================================================='''
def validate_Float(sVal, Min, Max, nDec):

    try:
        fVal = float(sVal)
    except:
        return 1

    if len(str(fVal).split('.')[1]) > nDec:
        return 2

    elif fVal < Min or fVal > Max:
        return 3

    return 0
    
    
''' ==========================================================
  validate_Integer: validate the entry field (int value)  
==========================================================='''
def validate_Integer(sVal, Min, Max):

    try:
        iVal = int(sVal)
    except:
        return 1

    if iVal < Min or iVal > Max:
        return 3

    return 0
    
    
''' ==========================================================
  validate_String: validate the entry field (str value)  
==========================================================='''
def validate_String(sVal):

    return 0

''' ==========================================================
setState: set states for all child widget  
==========================================================='''           
def setState(widget, state='disabled'):
    try:
        widget.configure(state=state)
    except:
        pass
    for child in widget.winfo_children():
        setState(child, state=state)

''' ==========================================================
backState: back states from list for all child widget  
==========================================================='''           
def backState(widget, li):
    try:
        s = widget.cget('state')
        state = li.pop(0)
        widget.configure(state=state)
    except:
        pass

    for child in widget.winfo_children():
        backState(child, li)

''' ==========================================================
saveState: save states list for all child widget  
==========================================================='''           
def saveState(widget, li):
    try:
        state = widget.cget('state')
        li.append(state)
    except:
        pass

    for child in widget.winfo_children():
        saveState(child, li)

    return li

''' ==========================================================
CenterWindow: Centers the window in the middle of the monitor
==========================================================='''               
def CenterWindow(top, w, h):

    # Get screen width and height
    ws = top.winfo_screenwidth()
    hs = top.winfo_screenheight()

    # Calculate position x,y
    x = (ws/2) - (w/2)
    y = (hs/2) - (h/2)

    top.geometry('%dx%d+%d+%d' % (w,h,x,y))

#=======================================================================
''' Get center of geometry of molecule '''
#=======================================================================   
def get_CenterGeometry(CG,PDBFile):

    del CG[:]
    CG.append(0.0)
    CG.append(0.0)
    CG.append(0.0)

    tot = 0

    try:
        # Read the PDB file
        file = open(PDBFile, 'r')
        PDBLines = file.readlines()
        file.close()

        # Read every ATOM/HETATM line
        for line in PDBLines:

            if line.startswith('ATOM  ') or \
               line.startswith('HETATM'):
            
                #ATOM      3  C   ALA A  13      24.276   5.552   9.942  1.00 43.61           C  
                CG[0] += float(line[30:38])
                CG[1] += float(line[38:46])
                CG[2] += float(line[46:54])
                tot += 1
                
    except:
        return -1

    CG[0] /= float(tot)
    CG[1] /= float(tot)
    CG[2] /= float(tot)
    
    return tot


#=======================================================================
''' Store Residue List from PDB file '''
#=======================================================================   
def store_Residues(listResidues, PDBFile, HETATM):

    del listResidues[:]

    max = 0

    try:
        # Read the PDB file
        file = open(PDBFile, 'r')
        PDBLines = file.readlines()
        file.close()

        # Get the possibles Residues
        for line in PDBLines:

            if line.startswith('ATOM  ') or \
               line.startswith('HETATM'):
            
                index = int(line[6:11].strip())
                resn = line[17:20].strip()
                chain = line[21:22].strip()
                resi = line[23:26].strip() 

                if index > max:
                    max = index

                if chain == ' ':
                    chain = '-'

                residue = resn + resi + chain

                if line.startswith('HETATM') and not HETATM:
                    continue

                if listResidues.count(residue) == 0:
                    listResidues.append(residue)

    except:
        return -1


    return max

''' ==================================================================================
FUNCTION repeat: repeats a character N number of times
==================================================================================  '''            
def repeat(string, length):
    L = len(string)
    return string * (length // L) + string[:length % L]

