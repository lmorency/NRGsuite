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

from subprocess import Popen, PIPE
from datetime import datetime

import os
import re
import glob
import shutil
import hashlib

import Result
import General

if __debug__:
    import Constraint

class Manage:

    NUMBER_RESULTS = 10

    def __init__(self, top):
        
        #print("New instance of Manage Class")
        self.top = top
        
        self.FlexAID = self.top.top
        self.IOFile = self.FlexAID.IOFile
        self.Config1 = self.FlexAID.Config1
        self.Config2 = self.FlexAID.Config2
        self.Config3 = self.FlexAID.Config3
        self.GAParam = self.FlexAID.GAParam
        
        self.LOGFILE = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'log.txt')
        self.LOGFILETMP = self.LOGFILE + '.tmp'

        self.Now = str(datetime.now())
        self.Now = self.Now[0:int(self.Now.rfind('.'))]
        self.Now = self.Now.replace(':','-')
        self.Now = self.Now.replace(' ','-')

        self.Target = self.IOFile.TargetName.get()
        self.Ligand = self.IOFile.LigandName.get()
        self.COMPLEX = self.IOFile.Complex.get().upper()
        
        self.TmpFile = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'file.tmp')

        self.TargetFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'TARGET.inp.pdb')
        self.RefFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIGAND_ref.pdb')
        self.INPFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIGAND.inp')
        self.ICFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIGAND.ic')

        self.VarAtoms = list()
        self.listTmpPDB = list()
        self.RecAtom = dict()
        self.DisAngDih = dict()
        self.dictCoordRef = dict()
        
    ''' ==============================================================================
    @summary: Reference_Folders: Create folder references with the now timestamp
    ============================================================================== '''          
    def Reference_Folders(self):

        self.FlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,self.COMPLEX,self.Now)
        self.INPFlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAIDRunSimulationProject_Dir,self.Ligand + '.inp')
        self.ICFlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAIDRunSimulationProject_Dir,self.Ligand + '.ic')
        
        self.BINDINGSITE = os.path.join(self.FlexAIDRunSimulationProject_Dir,'binding_site.pdb')
        
        self.PAUSE = os.path.join(self.FlexAIDRunSimulationProject_Dir,'.pause')
        self.STOP = os.path.join(self.FlexAIDRunSimulationProject_Dir,'.stop')
        self.ABORT = os.path.join(self.FlexAIDRunSimulationProject_Dir,'.abort')
        self.UPDATE = os.path.join(self.FlexAIDRunSimulationProject_Dir,'.update')
        self.READ = os.path.join(self.FlexAIDRunSimulationProject_Dir,'.read')

        self.CONFIG = os.path.join(self.FlexAIDRunSimulationProject_Dir,'CONFIG.inp')
        self.ga_inp = os.path.join(self.FlexAIDRunSimulationProject_Dir,'ga_inp.dat')

    ''' ==============================================================================
    @summary: Create_Folders: Creation AND/OR copy of the required files  
    ============================================================================== '''          
    def Create_Folders(self):
       
        try:
            if not os.path.isdir(self.FlexAIDRunSimulationProject_Dir):
                os.makedirs(self.FlexAIDRunSimulationProject_Dir)
        except:
            return False

        return True
 
    ''' ==============================================================================
    @summary: Executable_Exists: Check if FlexAID is present
    ============================================================================== '''          
    def Executable_Exists(self):
        
        try:
            if not os.path.isfile(self.FlexAID.FlexAIDExecutable):
                return False
        except:
            return False

        return True

    ''' ==============================================================================
    @summary: Move_Files: Moves the required input files
    ============================================================================== '''          
    def Move_Files(self):

        try:
            shutil.copy(self.INPFlexAIDSimulationProject_Dir, self.FlexAIDRunSimulationProject_Dir)
            shutil.copy(self.ICFlexAIDSimulationProject_Dir, self.FlexAIDRunSimulationProject_Dir)
            shutil.copy(self.RefFlexAIDSimulationProject_Dir, self.FlexAIDRunSimulationProject_Dir)
            shutil.copy(self.TargetFlexAIDSimulationProject_Dir, self.FlexAIDRunSimulationProject_Dir)
        except:
            return False

        return True

    ''' ==============================================================================
    @summary: Clean: Clean unecessary files and unset state files
    ============================================================================== '''          
    def Clean(self):

        try:
            if os.path.isfile(self.PAUSE):
                os.remove(self.PAUSE)

            if os.path.isfile(self.ABORT):
                os.remove(self.ABORT)

            if os.path.isfile(self.STOP):
                os.remove(self.STOP)

            if os.path.isfile(self.UPDATE):
                os.remove(self.UPDATE)

            if os.path.isfile(self.READ):
                os.remove(self.READ)

            if os.path.isfile(self.LOGFILE):
                os.remove(self.LOGFILE)

            if os.path.isfile(self.LOGFILETMP):
                os.remove(self.LOGFILETMP)

        except:
            return False

        return True

    ''' ==============================================================================
    @summary: : Load_ResultFiles: Loads the results files after a simulation terminated
    ============================================================================== '''          
    def Load_ResultFiles(self):
        
        pattern = os.path.join(self.FlexAIDRunSimulationProject_Dir,'RESULT_*')
        for file in glob.glob(pattern):
            
            m = re.search("RESULT_(\d+)\.pdb$", file)
            if m:
                TOP = int(m.group(1)) + 1
                Res = Result.Result(file, TOP)
                
                self.top.Vars.ResultsContainer.Results.append(Res)
                continue
            
            m = re.search("RESULT_INI\.pdb$", file)
            if m:
                Res = Result.Result(file, 'REF')
                
                if self.Config2.UseReference.get():
                    self.top.Vars.ResultsContainer.Results.append(Res)
                    
                continue                

            m = re.search("RESULT_par\.res$", file)
            if m:
                self.top.Vars.ResultsContainer.ResultParams = file
                continue
    
    ''' ==================================================================================
    @summary: Hash_CONFIG: Hashes the contents of the CONFIG recursively
    ================================================================================== '''
    def Hash_CONFIG(self):

        hasher = hashlib.md5()
        
        hasher = General.hashfile_update(self.IOFile.TargetPath.get(), hasher)
        hasher = General.hashfile_update(self.INPFlexAIDSimulationProject_Dir, hasher)

        #hasher.update('SPACER 0.375\n')
        
        rngOpt = self.Config1.RngOpt.get()
        
        if rngOpt == 'LOCCEN':
            line = 'RNGOPT LOCCEN'
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[0]
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[1]
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[2]
            line += ' %.3f\n' % self.Config1.Vars.BindingSite.Sphere.Radius
            
            hasher.update(line)

        elif rngOpt == 'LOCCLF':
            self.Config1.Generate_CleftBindingSite()
            
            hasher = General.hashfile_update(self.Config1.CleftTmpPath, hasher)
        
        if self.Config1.Vars.TargetFlex.Count_SideChain() > 0:
            self.Create_FlexFile(self.TmpFile)

            hasher = General.hashfile_update(self.TmpFile, hasher)
                        
        if len(self.Config2.Vars.dictConstraints):
            self.Create_ConsFile(self.TmpFile)

            hasher = General.hashfile_update(self.TmpFile, hasher)

        if self.Config2.IntTranslation.get():
            hasher.update('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - -1\n')
            
        if self.Config2.IntRotation.get():
            hasher.update('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - 0\n')
        
        if self.Config2.UseReference.get():
            hasher.update('RMSDST ' + self.IOFile.ReferencePath.get() + '\n')
        
        # Ligand flexibility
        order = sorted(self.IOFile.Vars.dictFlexBonds.keys())
        hasher.update(self.Add_FlexBonds(order))

        if self.IOFile.AtomTypes.get() == 'Sobolev': # 8 atom types only
            hasher.update('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','scr_bin.dat') + '\n')

        elif self.IOFile.AtomTypes.get() == 'Gaudreault': # 12 atom types
            hasher.update('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','M6_cons_3.dat') + '\n')

        elif self.IOFile.AtomTypes.get() == 'Sybyl': # 26 atom types
            hasher.update('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','SYBYL_emat.dat') + '\n')

        # permeability of atoms
        Permea = 1.00 - float(self.Config3.Permeability.get())
        hasher.update('PERMEA ' + str(Permea) + '\n')
        
        # heterogroups consideration
        if self.Config3.IncludeHET.get():
            hasher.update('INCHET' + '\n')
            
            if self.Config3.ExcludeHOH.get():
                hasher.update('RMVHOH' + '\n')
    
        #config_file.write('VARDIS ' + self.Config3.DeltaDistance.get() + '\n')
        hasher.update('VARANG ' + self.Config3.DeltaAngle.get() + '\n')
        hasher.update('VARDIH ' + self.Config3.DeltaDihedral.get() + '\n')
        hasher.update('VARFLX ' + self.Config3.DeltaDihedralFlex.get() + '\n')

        hasher.update('SLVTYP ' + str(self.Config3.SolventTypeIndex.get()) + '\n')
        if self.Config3.SolventTypeIndex.get() == 0:
            hasher.update('SLVPEN ' + self.Config3.SolventTerm.get() + '\n')
        
        hasher.update('MAXRES ' + str(self.NUMBER_RESULTS) + '\n')
                
        return hasher.digest()

    ''' ==================================================================================
    @summary: Create_CONFIG: Creation of the CONFIG.inp
    ================================================================================== '''
    def Create_CONFIG(self):

        # Save the data to the configuration file
        config_file = open(self.CONFIG, 'w')

        config_file.write('PDBNAM ' + self.IOFile.TargetPath.get() + '\n')
        
        config_file.write('INPLIG ' + self.INPFlexAIDSimulationProject_Dir + '\n')
        config_file.write('METOPT GA\n')

        #config_file.write('BPKENM ' + self.Config1.BPE.get() + '\n')
        config_file.write('COMPLF ' + self.Config3.CompFct.get() + '\n')
        
        rngOpt = self.Config1.RngOpt.get()

        config_file.write('SPACER 0.375\n')
        
        if rngOpt == 'LOCCEN':
            line = 'RNGOPT LOCCEN'
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[0]
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[1]
            line += ' %.3f' % self.Config1.Vars.BindingSite.Sphere.Center[2]
            line += ' %.3f\n' % self.Config1.Vars.BindingSite.Sphere.Radius
            config_file.write(line)

        elif rngOpt == 'LOCCLF':
            self.Config1.Generate_CleftBindingSite()
            self.Copy_BindingSite()
            
            line = 'RNGOPT LOCCLF '
            line += self.BINDINGSITE + '\n'
            config_file.write(line)

        if self.Config1.Vars.TargetFlex.Count_SideChain() > 0:
            FlexSCFile = os.path.join(self.FlexAIDRunSimulationProject_Dir,'flexSC.lst')
            self.Create_FlexFile(FlexSCFile)
            config_file.write('FLEXSC ' + FlexSCFile + '\n')

        if len(self.Config2.Vars.dictConstraints):
            ConsFile = os.path.join(self.FlexAIDRunSimulationProject_Dir,'cons.lst')
            self.Create_ConsFile(ConsFile)
            config_file.write('CONSTR ' + ConsFile + '\n')

        if self.Config2.IntTranslation.get():
            config_file.write('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - -1\n')
            
        if self.Config2.IntRotation.get():
            config_file.write('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - 0\n')
        
        if self.Config2.UseReference.get():
            config_file.write('RMSDST ' + self.IOFile.ReferencePath.get() + '\n')
        
        # Ligand flexibility
        order = sorted(self.IOFile.Vars.dictFlexBonds.keys())
        config_file.write(self.Add_FlexBonds(order))

        if self.IOFile.AtomTypes.get() == 'Sobolev': # 8 atom types only
            config_file.write('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','scr_bin.dat') + '\n')

        elif self.IOFile.AtomTypes.get() == 'Gaudreault': # 12 atom types
            config_file.write('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','M6_cons_3.dat') + '\n')

        elif self.IOFile.AtomTypes.get() == 'Sybyl': # 26 atom types
            config_file.write('IMATRX ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','SYBYL_emat.dat') + '\n')

        # permeability of atoms
        Permea = 1.00 - float(self.Config3.Permeability.get())
        config_file.write('PERMEA ' + str(Permea) + '\n')
        
        # heterogroups consideration
        if self.Config3.IncludeHET.get():
            config_file.write('INCHET' + '\n')
            
            if self.Config3.ExcludeHOH.get():
                config_file.write('RMVHOH' + '\n')

        #config_file.write('VARDIS ' + self.Config3.DeltaDistance.get() + '\n')
        config_file.write('VARANG ' + self.Config3.DeltaAngle.get() + '\n')
        config_file.write('VARDIH ' + self.Config3.DeltaDihedral.get() + '\n')
        config_file.write('VARFLX ' + self.Config3.DeltaDihedralFlex.get() + '\n')


        config_file.write('SLVTYP ' + str(self.Config3.SolventTypeIndex.get()) + '\n')
        if self.Config3.SolventTypeIndex.get() == 0:
            config_file.write('SLVPEN ' + self.Config3.SolventTerm.get() + '\n')

        config_file.write('STATEP ' + self.FlexAIDRunSimulationProject_Dir  + '\n')
        config_file.write('TEMPOP ' + self.FlexAID.FlexAIDSimulationProject_Dir  + '\n')
                
        config_file.write('DEPSPA ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps') + '\n')

        config_file.write('MAXRES ' + str(self.NUMBER_RESULTS) + '\n')
        
        config_file.write('NRGSUI\n')
        config_file.write('NRGOUT 60\n')
        
        config_file.write('GRDBUF 1000\n')

        config_file.write('ENDINP\n')

        config_file.close()
        
        return
        
    ''' ==================================================================================
    @summary: Copy_BindingSite: Copies from the temp folder to results folder
    ================================================================================== '''
    def Copy_BindingSite(self):
        
        try:
            shutil.copy(self.Config1.CleftTmpPath, self.BINDINGSITE)
        except:
            return False

        return True

    ''' ==================================================================================
    @summary: Create_ga_inp: Creation of the ga_inp.dat file  
    ================================================================================== '''
    def Create_ga_inp(self, bContinue=False):
        
        # Save the data to the configuration file
        gaInp_file = open(self.ga_inp, 'w')

        gaInp_file.write('NUMCHROM ' + str(self.GAParam.NbChrom.get()) + '\n')
        gaInp_file.write('NUMGENER ' + str(self.GAParam.NbGen.get()) + '\n')
        
        gaInp_file.write('ADAPTVGA ' + str(self.GAParam.UseAGA.get()) + '\n')
        if self.GAParam.UseAGA.get():
            gaInp_file.write('ADAPTKCO')
            gaInp_file.write(' %.2f' % float(self.GAParam.AGAk1.get()))
            gaInp_file.write(' %.2f' % float(self.GAParam.AGAk2.get()))
            gaInp_file.write(' %.2f' % float(self.GAParam.AGAk3.get()))
            gaInp_file.write(' %.2f\n' % float(self.GAParam.AGAk4.get()))

        gaInp_file.write('CROSRATE %.3f\n' % float(self.GAParam.CrossRate.get()))
        gaInp_file.write('MUTARATE %.3f\n' % float(self.GAParam.MutaRate.get()))

        if bContinue:
            gaInp_file.write('POPINIMT IPFILE ' + self.top.Vars.ResultsContainer.ResultParams + '\n')
        else:
            gaInp_file.write('POPINIMT RANDOM\n')

        gaInp_file.write('FITMODEL ' + self.GAParam.FitModel.get() + '\n')
        if self.GAParam.FitModel.get() == 'PSHARE':
            gaInp_file.write('SHAREALF %.2f\n' % float(self.GAParam.FitAlpha.get()))
            gaInp_file.write('SHAREPEK %.2f\n' % float(self.GAParam.FitPeak.get()))
            gaInp_file.write('SHARESCL %.2f\n' % float(self.GAParam.FitScale.get()))

        gaInp_file.write('REPMODEL ' + self.GAParam.RepModel.get() + '\n')
        if self.GAParam.RepModel.get() == 'BOOM':
            gaInp_file.write('BOOMFRAC %.2f\n' % float(self.GAParam.RepB.get()))
        elif self.GAParam.RepModel.get() == 'STEADY':
            gaInp_file.write('STEADNUM ' + str(self.GAParam.RepSS.get()) + '\n')

        if self.GAParam.RepDup.get():
            gaInp_file.write('DUPLICAT\n')

        gaInp_file.write('PRINTCHR ' + str(self.GAParam.NbTopChrom.get()) + '\n')

        gaInp_file.write('OUTGENER 1\n')

        gaInp_file.close()
    
    ''' ==================================================================================
    FUNCTION Create_FlexFile: Create the Flex file list that contain the residue selected 
                              for the Flexible Side Chain.
    ==================================================================================  '''          
    def Create_FlexFile(self,outfile):
        
        FlexFile = open(outfile, 'w')
                  
        for item in self.Config1.Vars.TargetFlex.listSideChain:
            lastch = len(item) - 1
            ResName = item[0:3]
            ChainID = item[lastch]
            ResSeq = item[3:(lastch)]

            Line = 'RESIDU ' + str(ResSeq).rjust(4) + ' ' + ChainID + ' ' + ResName + '\n'
            FlexFile.write(Line)
            
        if os.path.isfile(os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps','rotobs.lst')):
            FlexFile.write('ROTOBS\n')

        FlexFile.close()

    ''' ==================================================================================
    FUNCTION Create_ConsFile: Creates the constraints file
    ==================================================================================  '''          
    def Create_ConsFile(self,outfile):
                    
        ConsFile = open(outfile, 'w')
                  
        for key in sorted(self.Config2.Vars.dictConstraints, key=str.lower):
            
            atom1 = self.Config2.parse_cons(self.Config2.Vars.dictConstraints[key][0])
            atom2 = self.Config2.parse_cons(self.Config2.Vars.dictConstraints[key][1])
            
            ConsFile.write('COVALENT ')
            self.Print_Constraint(ConsFile,atom1)
            ConsFile.write(':')
            self.Print_Constraint(ConsFile,atom2)
            ConsFile.write(':')
            ConsFile.write('%.2f' % float(self.Config2.Vars.dictConstraints[key][5]))
            ConsFile.write('\n')
            
        ConsFile.close()            

    ''' ==================================================================================
    FUNCTION Print_Constraint: Prints a constraint atom
    ==================================================================================  '''          
    def Print_Constraint(self,outfile,atom):

        outfile.write('%4d' % int(atom[2]) + ' ')
        outfile.write(atom[3] + ' ')
        outfile.write('%3s' % atom[1] + ' ')
        outfile.write('%5d' % int(atom[0]))

    ''' ==================================================================================
    FUNCTION Add_FlexBonds: Adds optimization lines in CONFIG
    ==================================================================================  '''          
    def Add_FlexBonds(self, list):
        
        Content = ''
        ForcedLines = ''

        # Append extra OPTIMZ lines in CONFIG file
        for k in list:
            strk = str(k)

            # If bond is flexible
            if self.IOFile.Vars.dictFlexBonds[k][0]:
                Content += 'OPTIMZ ' + str(self.IOFile.ResSeq.get()) + ' - ' + strk + '\n'
        
        return Content
        
    ''' ==================================================================================
    FUNCTION Modify_Input: Changes the ICDATA line for the new path and changes the new types
    ==================================================================================  '''          
    def Modify_Input(self):

        # Store file content
        try:
            inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'r')
            lines = inpFile.readlines()
            inpFile.close()

            # Re-write file
            inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'w')
            for line in lines:
                if line.startswith('HETTYP'):
                    index = line[6:11].strip()
                    newline = line[:11] + self.IOFile.Vars.dictAtomTypes[index][1].rjust(2, ' ') +  line[13:]
                    inpFile.write(newline)
                else:
                    inpFile.write(line)

            inpFile.close()
            
        except:
            return 1
            
        return 0
        
    ''' ==================================================================================
    FUNCTION Get_CoordRef: Get the list of initial coordinates of the reference
    ==================================================================================  '''          
    def Get_CoordRef(self):

        self.dictCoordRef.clear()

        try:
            file = open(self.IOFile.ReferencePath.get(),'r')
            self.ReferenceLines = file.readlines()
            file.close()

            for Line in self.ReferenceLines:
                
                if Line.startswith('HETATM'):
                    index = int(Line[6:11].strip())
                    CoordX = float(Line[30:38].strip())
                    CoordY = float(Line[38:46].strip())
                    CoordZ = float(Line[46:54].strip())

                    self.dictCoordRef[index] = [ CoordX, CoordY, CoordZ ]
                    
        except:
            return 1
            
        return 0
        
    ''' ==================================================================================
    FUNCTION CreateTempPDB: Create X copy of the ligand PDB file
    ==================================================================================  '''          
    def CreateTempPDB(self):
        
        # Create the custom path for the temporary pdb files
        for i in range (int(self.GAParam.NbTopChrom.get()) + 1):
            self.listTmpPDB.append(os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIGAND' + str(i) + '.pdb'))
    
    ''' ==================================================================================
    FUNCTION Get_RecAtom: Create a dictionary containing the atoms neighbours
    ==================================================================================  '''
    def Get_RecAtom(self):
    
        self.RecAtom.clear()

        try:
            file = open(self.INPFlexAIDSimulationProject_Dir)
            inpLines = file.readlines()
            file.close()
            
            # Creation of a dictionary containing the 3 neighbours of each atom
            for line in inpLines:
                if line.startswith('HETTYP'):
                    noLine = int(line[6:11])
                    self.RecAtom[noLine] = [int(line[21:26]), int(line[26:31]), int(line[31:36])]
        except:
            return 1
                
        return 0
    
    ''' ==================================================================================
    FUNCTION Get_DisAngDih: Create a dictionary containing the internal coordinates
    ==================================================================================  '''
    def Get_DisAngDih(self):

        self.DisAngDih.clear()
                
        try:
            file = open(self.ICFlexAIDSimulationProject_Dir)
            icLines = file.readlines()
            file.close()
            
            # Creation of a dictionary containing the 3 neighbors of each atoms
            for line in icLines:
                if line[0:6] != 'REFPCG':
                    noLine = int(line[0:5])
                    self.DisAngDih[noLine] = [float(line[7:15]), float(line[16:24]), float(line[25:33])]
        except:
            return 1
            
        return 0

    ''' ==================================================================================
    FUNCTION Get_VarAtoms: Get the atoms that need their values to be modified
    ==================================================================================  '''
    def Get_VarAtoms(self):
        
        del self.VarAtoms[:]
        
        NoAtom3 = 0
        NoAtom2 = 0
        NoAtom1 = 0

        file = open(self.INPFlexAIDSimulationProject_Dir)
        inpLines = file.readlines()
        file.close()
        
        # Creation of a dictionary containing the 3 neighbors of each atoms
        for line in inpLines:
            if line.startswith('HETTYP'):
                noAtom = int(line[6:11])
                if int(line[31:36]) == 0:
                    if int(line[26:31]) == 0:
                        if int(line[21:26]) == 0:
                            NoAtom3 = noAtom
                        else:
                            NoAtom2 = noAtom
                    else:
                        NoAtom1 = noAtom
        
        self.VarAtoms.extend( [ NoAtom3, NoAtom2, NoAtom1 ] )
