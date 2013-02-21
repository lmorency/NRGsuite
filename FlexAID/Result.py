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

import re


class CF:

    def __init__(self, ResNumC):
        
        self.rnc = ResNumC
        
        self.com = 0.0
        self.wal = 0.0
        self.sas = 0.0
        self.con = 0.0
    

class Result:
    
    def __init__(self):
    
        self.CF = 'N/A'
        self.CFapp = 'N/A'
        
        self.RMSD = 'N/A'
        
        self.Optimizable = []
        
        self.ResultID = 0
        self.ResultFile = ''
    
    def __init__(self, ResultFile, ResultID):
    
        self.CF = 'N/A'
        self.CFapp = 'N/A'
        
        self.RMSD = 'N/A'
        
        self.Optimizable = []
        
        self.ResultID = ResultID
        self.ResultFile = ResultFile
        
        self.get_CF_info()
        
    # Reads the physical file to retrieve the CF information
    def get_CF_info(self):
    
        if self.ResultFile:
    
            try:
                fh = open(self.ResultFile,'r')
            except IOError:
                return
            
            Lines = fh.readlines()
            fh.close()
            
            '''
            REMARK CF.app=-845.79439
            REMARK optimizable residue LIG   9999
            REMARK CF.com=370.63038
            REMARK CF.sas=-1196.13076
            REMARK CF.wal=20.29402
            REMARK CF.con= 0.00000
            '''
            
            for Line in Lines:
                
                if re.match('REMARK', Line):
                    
                    m = re.search('optimizable residue (.{3}) (.) (.{4})', Line)
                    if m:
                        Res = m.group(1)
                        C = m.group(2)
                        Num = m.group(3)
                        
                        if C == ' ':
                            C = '-'
                        
                        Res = Res.replace(' ','-')
                        Num = Num.replace(' ','')
                        ResNumC = Res + Num + C
                        
                        self.Opt = CF(ResNumC)
                        self.Optimizable.append(self.Opt)
                        continue

                    m = re.search('CF=\s*(\S+)', Line)
                    if m:
                        self.CF = float(m.group(1))
                        continue

                    m = re.search('CF\.app=\s*(\S+)', Line)
                    if m:
                        self.CFapp = float(m.group(1))
                        continue

                    m = re.search('CF\.con=\s*(\S+)', Line)
                    if m:
                        self.Opt.con = float(m.group(1))
                        continue
                        
                    m = re.search('CF\.wal=\s*(\S+)', Line)
                    if m:
                        self.Opt.wal = float(m.group(1))
                        continue
                    
                    m = re.search('CF\.sas=\s*(\S+)', Line)
                    if m:
                        self.Opt.sas = float(m.group(1))
                        continue
                        
                    m = re.search('CF\.com=\s*(\S+)', Line)
                    if m:
                        self.Opt.com = float(m.group(1))
                        continue
                    
                    #REMARK  7.32206 RMSD to ref. structure
                    m = re.search('(\S+) RMSD to ref\. structure', Line)
                    if m:
                        self.RMSD = float(m.group(1))
                        continue
                    
                elif re.match('ATOM  ', Line):
                    break
                    

class ResultsContainer:
    
    def __init__(self):

        # Result PDB files
        self.Results = list()

        # File containing the parameters allowing to continue the simulation
        self.ResultParams = ''
        
        # If these results were generated from other results this will be set
        self.ParentResult = None
        
