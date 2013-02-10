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
import General

from subprocess import Popen, PIPE


class ProcLig:


    def __init__(self, top, StartAtomIndex, AtomTypes, AnchorAtom, ConvertOnly, Gen3D):
        
        self.top = top
        self.FlexAID = self.top.top

        self.LigandPath = self.FlexAID.IOFile.LigandPath.get()

        self.AtomTypes = AtomTypes
        self.AnchorAtom = AnchorAtom

        self.ConvertOnly = ConvertOnly
        self.Gen3D = Gen3D
        
        self.FlexAIDWRKInstall_Dir = os.path.join(self.FlexAID.FlexAIDInstall_Dir,'WRK')

        #Init the INP file path
        self.SimLigInpPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.inp')
        self.SimLigICPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.ic')

        self.StartAtomIndex = StartAtomIndex

        self.run()


    def run(self):
        
        self.FlexAID.ProcessRunning = True
        
        if not self.Copy_LigandFile():

            self.FlexAID.DisplayMessage("  Executing Process Ligand...", 0)
            
            if not self.process():
            
                if self.ConvertOnly:
                    self.FlexAID.DisplayMessage("  The ligand was converted successfully", 0)
                else:
                    self.FlexAID.DisplayMessage("  The ligand was processed successfully", 0)

            else:
                self.top.ProcessError = True

                if self.ConvertOnly:
                    self.FlexAID.DisplayMessage("  ERROR: The ligand could not be converted", 1)
                else:
                    self.FlexAID.DisplayMessage("  ERROR: The ligand could not be processed", 1)

        else:
            # could not copy file
            self.top.ProcessError = True

            self.FlexAID.DisplayMessage("  ERROR: The ligand could not be copied", 1)

        self.top.ProcessLigand(False, 0, '', 0, False, 0)
        self.FlexAID.ProcessRunning = False


    def Copy_LigandFile(self):
        
        try:
            copy(self.LigandPath, self.FlexAID.FlexAIDSimulationProject_Dir)            
        except IOError:
            return 1
        except Error:
            pass
        
        self.SimLigPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir, os.path.split(self.LigandPath)[1])

        return 0

       
    '''
    @summary: process: Processes the ligand (generates input files for FlexAID)
    '''  
    def process(self):

        # Set the command-line arguments to process ligand
        commandline = '"' + self.FlexAID.Process_LigandExecutable + '"'

        commandline += ' -f ' + '"' + self.SimLigPath + '"'
        commandline += ' -o ' + '"' + os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG') + '"'
        
        if self.Gen3D:
            commandline += ' --gen3D'
        
        if self.ConvertOnly:
            commandline += ' -c'
        else:
            commandline += ' --atom_index ' + str(self.StartAtomIndex)
            commandline += ' -ref'
            
            if self.AtomTypes == 'Sobolev':
                commandline += ' --old_types'
            elif self.AtomTypes == 'Gaudreault':
                commandline += ' --new_types'
            elif self.AtomTypes == 'Sybyl':
                commandline += ' --babel_types'
            
            if self.AnchorAtom != -1:
                commandline += ' --force_gpa ' + str(self.AnchorAtom)
        
        print(commandline)
        
        # Execute command-line
        try:
        
            if self.FlexAID.OSid == 'WIN':
                self.set_environment('BABEL_DATADIR', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'data') + '"')
                self.FlexAID.Run = Popen(commandline, shell=False, stderr=PIPE, stdout=PIPE)
                
            elif self.FlexAID.OSid == 'LINUX':
                self.set_environment('LD_LIBRARY_PATH', os.path.join(self.FlexAIDWRKInstall_Dir,'libs'))
                self.set_environment('BABEL_LIBDIR', os.path.join(self.FlexAIDWRKInstall_Dir,'formats'))
                self.FlexAID.Run = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)
                
            elif self.FlexAID.OSid == 'MAC':
                self.set_environment('DYLD_LIBRARY_PATH', os.path.join(self.FlexAIDWRKInstall_Dir,'libs'))
                self.set_environment('BABEL_LIBDIR', os.path.join(self.FlexAIDWRKInstall_Dir,'formats'))
                self.FlexAID.Run = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)

            (out,err) = self.FlexAID.Run.communicate()
            
            if self.FlexAID.Run.returncode != 0:
                self.top.ProcessError = True
            
        except:
            return 1

        self.FlexAID.Run = None
        
        if out.find('Done.') != -1:
            self.top.ReferencePath.set(os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG_ref.pdb'))
            return 0
        
        return 1
        
    
    '''
    @summary: set_environment sets a environment variable
    '''  
    def set_environment(self, variable, value):

        os.putenv(variable, value)
        if os.getenv(variable) == None:
            os.environ[variable] = value
        