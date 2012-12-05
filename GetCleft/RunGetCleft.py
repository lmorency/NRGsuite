#!/usr/bin/python

from Tkinter import *

import tkFileDialog, tkMessageBox
import tkFont, os, sys

sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite/PLUGIN')
sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite')

import GetCleft

Project_Dir = '/Users/francisgaudreault/Documents/NRGsuite/Default'
Install_Dir = '/usr/local/NRGsuite'
AlreadyRunning_Dir = '/Users/francisgaudreault/Documents/NRGsuite'

root = Tk()
GetCleft = GetCleft.displayGetCleft(root, None, Project_Dir, Install_Dir, AlreadyRunning_Dir, 'Mac', 0)
root.mainloop()
