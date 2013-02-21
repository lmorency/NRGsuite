#!/usr/bin/python

from Tkinter import *

import tkFileDialog, tkMessageBox
import tkFont, os, sys

sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite/PLUGIN')
sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite')

Project_Dir = '/Users/francisgaudreault/Documents/NRGsuite/Default'
Install_Dir = '/usr/local/NRGsuite'
NRGsuite_Dir = '/Users/francisgaudreault/Documents/NRGsuite'

root = Tk()

import GetCleft

GetCleft = GetCleft.displayGetCleft(root, None, -1, Project_Dir, Install_Dir, NRGsuite_Dir, 'MAC', 0, 'GetCleft', 500, 550)
root.mainloop()
