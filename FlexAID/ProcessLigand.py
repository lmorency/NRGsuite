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
@title: FlexAID - ProcessLigand.py

@summary: Class that create the Protein folder in WRK folder
          and add the required files for the simulation.

@contain: 

@organization: Najmanovich Research Group
@creation date:  oct. 1, 2010
'''

from shutil import copy, Error

import os
import time
import shutil
import threading
import General

from subprocess import Popen, PIPE


#class ProcLig(threading.Thread):
class ProcLig():

    def __init__(self, top, MoleculeFile, StartAtomIndex, AtomTypes, AnchorAtom, 
                 ConvertOnly, ProcessOnly, Gen3D, Target, Callback):
        
        #threading.Thread.__init__(self)

        self.top = top
        self.FlexAID = self.top.top
        
        self.MoleculeFile = MoleculeFile
        
        self.AtomTypes = AtomTypes
        self.AnchorAtom = AnchorAtom

        self.ConvertOnly = ConvertOnly
        self.ProcessOnly = ProcessOnly
        self.Gen3D = Gen3D
        self.Target = Target
        
        self.Callback = Callback
        
        self.FlexAIDWRKInstall_Dir = os.path.join(self.FlexAID.FlexAIDInstall_Dir,'WRK')
        
        self.StartAtomIndex = StartAtomIndex
        
        #self.start()
        self.run()
        
    def run(self):
        
        self.FlexAID.ProcessRunning = True
        
        if not self.Copy_MoleculeFile():
            #self.FlexAID.DisplayMessage("  Processing %s molecule..." % molecule, 0)
            rv = self.process()
        
        self.Callback(False, '', 0, '', 0, False, False, 0, False, '')
        self.FlexAID.ProcessRunning = False
        
        return rv
        
    '''
    @summary: process: Processes the ligand (generates input files for FlexAID)
    '''  
    def process(self):

        # Set the command-line arguments to process ligand
        commandline = '"' + self.FlexAID.Process_LigandExecutable + '"'

        commandline += ' -f ' + '"' + self.TmpMoleculeFile + '"'
        
        if self.Target:
            # use default outputting for the target
            commandline += ' -target'
        else:
            #commandline += ' -o ' + '"' + os.path.join(self.FlexAID.FlexAIDTempProject_Dir,'LIGAND') + '"'
        
            if self.Gen3D:
                commandline += ' --gen3D'

            if self.ConvertOnly:
                commandline += ' -c'
            else:
                if self.ProcessOnly:
                    commandline += ' -p'
                
                commandline += ' --atom_index ' + str(self.StartAtomIndex)
                commandline += ' -ref'
                
                #if self.AtomTypes == 'Sobolev':
                #    commandline += ' --old_types'
                #elif self.AtomTypes == 'Gaudreault':
                #    commandline += ' --new_types'
                #elif self.AtomTypes == 'Sybyl':
                #    commandline += ' --babel_types'
                
                if self.AnchorAtom != -1:
                    commandline += ' --force_gpa ' + str(self.AnchorAtom)
        
        print(commandline)
        
        # Execute command-line
        try:
        
            if self.FlexAID.OSid == 'WIN':
                self.set_environment('BABEL_DATADIR', os.path.join(self.FlexAIDWRKInstall_Dir,'data'))
                self.FlexAID.Run = Popen(commandline, shell=False, stderr=PIPE, stdout=PIPE)
                                                
            elif self.FlexAID.OSid == 'LINUX':
                self.set_environment('LD_LIBRARY_PATH', os.path.join(self.FlexAIDWRKInstall_Dir,'libs'))
                self.set_environment('BABEL_LIBDIR', os.path.join(self.FlexAIDWRKInstall_Dir,'formats'))
                self.set_environment('BABEL_DATADIR', os.path.join(self.FlexAIDWRKInstall_Dir,'data'))
                self.FlexAID.Run = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)
                
            elif self.FlexAID.OSid == 'MAC':
                self.set_environment('DYLD_LIBRARY_PATH', os.path.join(self.FlexAIDWRKInstall_Dir,'libs'))
                self.set_environment('BABEL_LIBDIR', os.path.join(self.FlexAIDWRKInstall_Dir,'formats'))
                self.set_environment('BABEL_DATADIR', os.path.join(self.FlexAIDWRKInstall_Dir,'data'))
                self.FlexAID.Run = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)

            (out,err) = self.FlexAID.Run.communicate()
            #print "out", out
            #print "err", err
            
            if self.FlexAID.Run.returncode != 0:
                self.FlexAID.ProcessError = True
                return 1
            
        except:
            self.FlexAID.ProcessError = True
            return 1
        
        self.FlexAID.Run = None
        
        return 0

    def Copy_MoleculeFile(self):

        try:
            self.TmpMoleculeFile = os.path.join(self.FlexAID.FlexAIDTempProject_Dir,
                                                os.path.split(self.MoleculeFile)[1])

            copy(self.MoleculeFile, self.TmpMoleculeFile)
            
        except IOError:
            return 1
        except Error:
            pass
        
        
        return 0
    
    '''
    @summary: set_environment sets a environment variable
    '''  
    def set_environment(self, variable, value):

        os.putenv(variable, value)
        if os.getenv(variable) != value:
            os.environ[variable] = value
    
