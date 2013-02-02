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

if __debug__:
    import Constraint

class Manage():

    def __init__(self,top):
        
        #print("New instance of Manage Class")
        self.top = top
        
        self.FlexAID = self.top.top
        self.IOFile = self.FlexAID.IOFile
        self.Config1 = self.FlexAID.Config1
        self.Config2 = self.FlexAID.Config2
        self.Config3 = self.FlexAID.Config3
        self.GAParam = self.FlexAID.GAParam

        self.PAUSE = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'.pause')
        self.STOP = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'.stop')
        self.ABORT = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'.abort')
        
        self.LOGFILE = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'sim.log')
        self.LOGFILETMP = self.LOGFILE + '.tmp'

        self.Now = str(datetime.now())
        self.Now = self.Now[0:int(self.Now.rfind('.'))]
        self.Now = self.Now.replace(':','-')
        self.Now = self.Now.replace(' ','-')

        self.Protein = self.IOFile.ProtName.get()
        self.Ligand = self.IOFile.LigandName.get()
        
        self.RefFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG_ref.pdb')
        self.INPFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.inp')
        self.ICFlexAIDSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,'LIG.ic')

        self.FlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAID.FlexAIDSimulationProject_Dir,self.Protein + '-' + self.Ligand,self.Now)
        self.INPFlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAIDRunSimulationProject_Dir,self.Ligand + '.inp')
        self.ICFlexAIDRunSimulationProject_Dir = os.path.join(self.FlexAIDRunSimulationProject_Dir,self.Ligand + '.ic')
        
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
            if self.FlexAID.OSid == 'WIN':
                p = Popen('copy "' + self.INPFlexAIDSimulationProject_Dir + '" "' + self.INPFlexAIDRunSimulationProject_Dir + '"' , shell=True)
                p.wait()
                p = Popen('copy "' + self.ICFlexAIDSimulationProject_Dir + '" "' + self.ICFlexAIDRunSimulationProject_Dir + '"' , shell=True)
                p.wait()
            else:
                p = Popen('cp "' + self.INPFlexAIDSimulationProject_Dir + '" "' + self.INPFlexAIDRunSimulationProject_Dir + '"' , shell=True)
                p.wait()
                p = Popen('cp "' + self.ICFlexAIDSimulationProject_Dir + '" "' + self.ICFlexAIDRunSimulationProject_Dir + '"' , shell=True)
                p.wait()

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

            if os.path.isfile(self.LOGFILE):
                os.remove(self.LOGFILE)

            if os.path.isfile(self.LOGFILETMP):
                os.remove(self.LOGFILETMP)

        except:
            return False

        return True

    ''' ==================================================================================
    @summary: Create_CONFIG: Creation of the CONFIG.inp
    ================================================================================== '''
    def Create_CONFIG(self):        
        
        # Save the data to the configuration file
        config_file = open(self.CONFIG, 'w')

        config_file.write('PDBNAM ' + self.IOFile.ProtPath.get() + '\n')
        
        config_file.write('INPLIG ' + os.path.join(self.FlexAIDRunSimulationProject_Dir,self.Ligand + '.inp') + '\n')
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
            line = 'RNGOPT LOCCLF '
            line += self.Config1.CleftTmpPath + '\n'
            config_file.write(line)        

        
        if self.Config1.Vars.TargetFlex.Count_SideChain() > 0:
            FlexSCFile = os.path.join(self.FlexAIDRunSimulationProject_Dir,'flexSC.lst')
            self.Create_FlexFile(FlexSCFile)
            config_file.write('FLEXSC ' + FlexSCFile + '\n')

        if len(self.Config2.Vars.dictCovConstraints):
            ConsFile = os.path.join(self.self.FlexAIDRunSimulationProject_Dir,'cons.lst')
            self.Create_ConsFile(ConsFile)
            config_file.write('CONSTR ' + ConsFile + '\n')

        if self.Config2.IntTranslation.get():
            config_file.write('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - -1\n')
            
        if self.Config2.IntRotation.get():
            config_file.write('OPTIMZ ' + str(self.IOFile.ResSeq.get())  + ' - 0\n')
        
        # Ligand flexibility
        order = sorted(self.Config2.Vars.dictFlexBonds.keys())
        self.Add_FlexBonds(config_file,order)

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

        config_file.write('STATEP ' + self.FlexAID.FlexAIDSimulationProject_Dir  + '\n')

        config_file.write('DEPSPA ' + os.path.join(self.FlexAID.FlexAIDInstall_Dir,'deps') + '\n')

        config_file.write('MAXRES 30'  + '\n')

        config_file.write('NRGSUI\n')
        config_file.write('ENDINP\n')

        config_file.close()
        

    ''' ==================================================================================
    @summary: Create_ga_inp: Creation of the ga_inp.dat file  
    ================================================================================== '''
    def Create_ga_inp(self):
        
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
    #def Create_ConsFile(self,outfile):
    #                
    #    ConsFile = open(outfile, 'w')
    #              
    #    for key in iter(self.Config2.Vars.dictCovConstraints):
    #        
    #        atom1 = Constraint.parse_cons(self.Config2.Vars.dictCovConstraints[key][0])
    #        atom2 = Constraint.parse_cons(self.Config2.Vars.dictCovConstraints[key][1])
    #        
    #        ConsFile.write('COVALENT ')
    #        self.Print_Constraint(ConsFile,atom1)
    #        ConsFile.write(':')
    #        self.Print_Constraint(ConsFile,atom2)
    #        ConsFile.write(':')
    #        ConsFile.write('%.2f' % float(self.Config2.Vars.dictCovConstraints[key][5]))
    #        ConsFile.write('\n')
    #        
    #    ConsFile.close()            

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
    def Add_FlexBonds(self, FilePtr,list):
                    
        ForcedLines = ''

        # Append extra OPTIMZ lines in CONFIG file
        for k in list:
            strk = str(k)

            # If bond is flexible
            if self.Config2.Vars.dictFlexBonds[k][0]:
                FilePtr.write('OPTIMZ ' + str(self.IOFile.ResSeq.get()) + ' - ' + strk + '\n')
                            
                # If the bond is forced
                if self.Config2.Vars.dictFlexBonds[k][1]:

                    ForcedLines += 'FLEDIH' + strk.rjust(3, ' ') + ' '
                    n = self.Config2.Vars.dictFlexBonds[k][2]            
                    for i in range(0, n):
                        ForcedLines += str(self.Config2.Vars.dictFlexBonds[k][i+3]).rjust(5, ' ')
                    ForcedLines += '\n'                
        
        self.Add_Forced(ForcedLines)

    ''' ==================================================================================
    FUNCTION Add_Forced: Adds FLEDIH lines in ligand inp
    ==================================================================================  '''          
    def Add_Forced(self, ForcedLines):

        # Append extra FLEDIH lines (Forced flexible bonds)
        if ForcedLines != '':

            # Store file content
            inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'r')
            lines = inpFile.readlines()
            inpFile.close()

            # Re-write file
            inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'w')
            for line in lines:
                if line.startswith('GPATOM'):
                    inpFile.write(ForcedLines)
                inpFile.write(line)
            inpFile.close()


    ''' ==================================================================================
    FUNCTION Modify_Input: Changes the ICDATA line for the new path and changes the new types
    ==================================================================================  '''          
    def Modify_Input(self):

        # Store file content
        inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'r')
        lines = inpFile.readlines()
        inpFile.close()

        # Re-write file
        inpFile = open(self.INPFlexAIDSimulationProject_Dir, 'w')
        for line in lines:
            if line.startswith('HETTYP'):
                index = line[6:11].strip()
                newline = line[:11] + self.Config2.Vars.dictAtomTypes[index][1].rjust(2, ' ') +  line[13:]
                inpFile.write(newline)
            else:
                inpFile.write(line)

        inpFile.close()
