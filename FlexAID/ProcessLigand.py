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

#import threading
import os
import shutil
import General

from subprocess import Popen, PIPE


class ProcLig():


    def __init__(self, top, StartAtomIndex, AtomTypes, AnchorAtom):
        
        #threading.Thread.__init__(self)

        self.top = top
        self.FlexAID = self.top.top

        self.ProtPath = self.FlexAID.IOFile.ProtPath
        self.LigandPath = self.FlexAID.IOFile.LigandPath

        self.AtomTypes = AtomTypes
        self.AnchorAtom = AnchorAtom

        self.FlexAIDWRKInstall_Dir = os.path.join(self.FlexAID.FlexAIDInstall_Dir,'WRK')

        #Init the INP file path
        self.SimLigInpPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.inp')
        self.SimLigICPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.ic')

        self.ResSeq = 0
        self.StartAtomIndex = StartAtomIndex

        self.FlexAID.DisplayMessage("  Executing Process Ligand...",0)

        #self.start()
        self.run()

 

    '''
    @summary: run: Side thread
    '''
    def run(self):
        
        #self.Set_LigandPDB()

        if not self.Copy_LigandFile():

            if not self.process():

                # ligand extraction success
                self.top.fProcessLigand = True
                self.top.ResSeq = 9999

                self.FlexAID.DisplayMessage("  The ligand '" + self.FlexAID.IOFile.LigandName.get() + "' was processed successfully",0)

            else:
                # ligand extraction failed
                self.top.ProcessError = True
                self.FlexAID.DisplayMessage("  ERROR: The ligand '" + self.FlexAID.IOFile.LigandName.get() + "' could not be processed",1)

        else:
            # could not copy file
            self.top.ProcessError = True
            self.FlexAID.DisplayMessage("  ERROR: The ligand '" + self.FlexAID.IOFile.LigandName.get() + "' could not be copied",1)


        self.top.ProcessLigand(False, 0)
        self.FlexAID.ProcessRunning = False

    '''
    @summary: Set_LigandPDB: Prepare the ligand PDB file for the ligand_extractor  
    '''                 
    def Copy_LigandFile(self):
        
        try:
            shutil.copy(self.LigandPath, self.FlexAID.FlexAIDSimulationProject_Dir)

            self.SimLigPath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir, os.path.split(self.LigandPath)[1])
        except:
            return 1

        return 0

    '''
    @summary: Set_LigandPDB: Prepare the ligand PDB file for the ligand_extractor  
    '''         
    def Set_LigandPDB(self):
        
        #Replace the index number in the ligand pdb file
        file = open(self.LigandPath, 'r')
        ligFile = file.readlines()
        file.close()
        
        modif_file = open(self.SimLigPath, 'w')
        
        for line in ligFile:
            
            if line.startswith('ATOM') or line.startswith('HETATM'):
                tmpLine = ""
                if line.startswith('ATOM'):
                    #Change the Atom for HetAtm
                    tmpLine = 'HETATM'
                else:
                    #Change the Atom Index no
                    tmpLine = line[0:6]
                    
                NoAtom = int(line[7:11]) + self.INDEXoffset                
                tmpLine += str(NoAtom).rjust(5, ' ')
                tmpLine += line[11:17]
                tmpLine += 'LIG A'
                tmpLine += str(self.ResSeq).rjust(4, ' ')
                tmpLine += line[26:]
                
                modif_file.write(tmpLine)
                
            elif line.startswith('CONECT'):
                #Change the Atom Index no
                tmpLine = line[0:6]
                NbConnect = ((len(line) - 7)/5) - 1
                AtomNo = int(line[7:11]) + self.INDEXoffset
                tmpLine += str(AtomNo).rjust(5, ' ')
                
                colNo = 12
                for i in range(0,NbConnect):
                    AtomNo = int(line[colNo:colNo+4]) + self.INDEXoffset
                    tmpLine += str(AtomNo).rjust(5, ' ')
                    colNo += 5
                
                tmpLine += '\n'   
                modif_file.write(tmpLine)
                
            else:                                       
                modif_file.write(line)                                    
        
        modif_file.close()      
        
       
    '''
    @summary: process: Processes the ligand (generates input files)
    '''  
    def process(self):

        # Set the command-line arguments to process ligand

        if self.FlexAID.OSid == 'WIN':
            commandline = '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'Process_Ligand.exe') + '"'
        else:
            commandline = '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'Process_Ligand') + '"'

        commandline += ' -f ' + '"' + self.SimLigPath + '"'
        #commandline += ' -rh '
        commandline += ' -o ' + '"' + os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG') + '"' 
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
            
        # Execute command-line
        try:
            if self.FlexAID.OSid == 'WIN':
                # Set environment variable here
                os.putenv('BABEL_DATADIR', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'data') + '"')
                pipe = Popen(commandline, shell=False, stderr=PIPE, stdout=PIPE)
            elif self.FlexAID.OSid == 'LINUX':
                os.putenv('LD_LIBRARY_PATH', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'libs') + '"')
                os.putenv('BABEL_LIBDIR', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'formats') + '"')
                pipe = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)
            elif self.FlexAID.OSid == 'MAC':
                os.putenv('DYLD_LIBRARY_PATH', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'libs') + '"')
                os.putenv('BABEL_LIBDIR', '"' + os.path.join(self.FlexAIDWRKInstall_Dir,'formats') + '"')
                pipe = Popen(commandline, shell=True, stderr=PIPE, stdout=PIPE)
            else:
                self.FlexAID.DisplayMessage("  ERROR: Your platform is not supported by NRGsuite",1)
                return 1
        except:
            return 1

        the_outerr = pipe.stderr.read()
        the_output = pipe.stdout.read()

        if the_output.find('Done.') != -1:
            self.top.ReferencePath = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG_ref.pdb')
            return 0

        else:
            return 1
        
    
