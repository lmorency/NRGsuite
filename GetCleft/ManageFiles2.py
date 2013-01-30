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
@title: GetCleft - Interface - Default tab

@summary: This is the interface of GetCleft application accessible via the
          PyMOL menu (NRGsuite).

@organization: Najmanovich Research Group
@creation date:  Oct. 19, 2010
'''

import os
import shutil
import glob
import CleftObj

#=========================================================================================
'''                        ---  DIRECTORY MANAGEMENT  ---                          '''
#=========================================================================================     
class Manage:

    def __init__(self,top):
        
        self.top = top

    # Delete the clf files of the GetCleft output
    def delete_Temp(self):
    
        for CleftFile in glob.glob(os.path.join(self.top.GetCleftTempProject_Dir,'*_clf_*')):
            try:
                os.remove(CleftFile)
            except:
                pass
    
    # Store the temporary clefts
    def store_Temp(self, Target):
        
        for CleftFile in glob.glob(os.path.join(self.top.GetCleftTempProject_Dir,'*_sph_*')):
            CleftName = os.path.splitext(os.path.basename(CleftFile))[0]
            
            Cleft = CleftObj.CleftObj()
            Cleft.CleftFile = CleftFile
            Cleft.CleftName = CleftName
            Cleft.UTarget = Target.upper()
            Cleft.Set_CleftMD5()

            self.top.Default.TempBindingSite.Add_Cleft(Cleft)

    # When creating a new partition object, copy the partition in the temp directory (do not save)
    def copy_TempPartition(self, FromFile, DestFile):
        
        try:
            shutil.copy(FromFile, DestFile)
        except:
            self.top.DisplayMessage("  ERROR: Could not copy temporary partition file.", 1)
    
    # Copy files from the Temp to the Save directory
    def save_Temp(self):
            
        nCopy = 0
        for Cleft in iter(self.top.Default.TempBindingSite.listClefts):

            SavePath = os.path.join(self.top.GetCleftSaveProject_Dir, Cleft.UTarget)
            
            if not os.path.isdir(SavePath):
                os.makedirs(SavePath)
                
            DestFile = os.path.join(SavePath,os.path.basename(Cleft.CleftFile))
            try:
                shutil.copy(Cleft.CleftFile, DestFile)
                nCopy = nCopy + 1
            except:
                pass
                #self.top.DisplayMessage("  ERROR: Could not copy file '" + Cleft.CleftFile + "'", 2)

            Cleft.CleftFile = DestFile

        if nCopy:
            self.top.DisplayMessage("  Successfully copied (" + str(nCopy) + ") cleft object file(s)", 0)
        else:
            self.top.DisplayMessage("  No cleft object file(s) were copied", 0)

        self.top.CopySession = True
        
        return
        

    # Copies the temporary Partition file
    def copy_TempPartition(self, PartitionFile, OutputFile):
        
        try:
            shutil.copy(PartitionFile, OutputFile)            
            #self.top.DisplayMessage("Partition file saved to " + OutputFile, 0)
        except:
            self.top.DisplayMessage("Could not copy temporary partition file", 2)

    # Clean temporary directories and bindingsite
    def Clean(self):

        for file in glob.glob(os.path.join(self.top.GetCleftTempProject_Dir,'*.pdb')):
            try:
                os.remove(file)
            except OSError:
                self.top.DisplayMessage("  ERROR: Could not remove temporary file(s)", 1)
                pass

