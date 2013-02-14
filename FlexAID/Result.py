'''
    NRGsuite: PyMOL molecular tools interface
    Copyright (C) 2013 Gaudreault, F., Morin, E. & Najmanovich, R.

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

class CF:

    def __init__(self):
    
        self.sum = 0.0
        self.app = 0.0
        self.com = 0.0
        self.wal = 0.0
        self.sas = 0.0
        self.con = 0.0
    

class Result:
    
    def __init__(self):
    
        self.CF = CF()
        self.ResultFile = ''
    
    # Reads the physical file to retrieve the CF information
    def get_CF_info(self):
    
        return


class ResultsContainer:
    
    # Result PDB files
    Results = []
    
    # File containing the parameters allowing to continue the simulation
    ResultParams = ''
    
    # If these results were generated from other results this will be set
    ParentResult = None
    

    