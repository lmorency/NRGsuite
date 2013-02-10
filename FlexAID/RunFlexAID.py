#!/usr/bin/python

from Tkinter import *

import tkFileDialog, tkMessageBox
import tkFont, os, sys

sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite/PLUGIN')
sys.path.append('/Users/francisgaudreault/Development/NRGsuite/NRGsuite')

Project_Dir = '/Users/francisgaudreault/Documents/NRGsuite/Default'
Install_Dir = '/usr/local/NRGsuite'
AlreadyRunning_Dir = '/Users/francisgaudreault/Documents/NRGsuite'

root = Tk()

import FlexAID

FlexAID = FlexAID.displayFlexAID(root, None, Project_Dir, Install_Dir, AlreadyRunning_Dir, 'MAC', 0,
                                 'FlexAID', '.frun', 700, 600)
root.mainloop()
