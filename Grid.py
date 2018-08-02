'''
    NRGsuite: PyMOL molecular tools interface
    Copyright (C) 2011 Gaudreault, F., Morency, LP. & Najmanovich, R.

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

import threading
import Geometry

class Grid(threading.Thread):
    
    def __init__(self, top, CleftFile, OutputFile, Spacer, Estimate):
        
        threading.Thread.__init__(self)

        self.top = top

        self.top.ProcessError = False

        self.CleftFile = CleftFile
        self.Spacer = float(Spacer)

        self.OutputFile = OutputFile
        self.Estimate = Estimate

        self.dictGridPoints = dict()
        self.dictSpheres = dict()

        # Starts thread
        self.start()

    #=======================================================================
    """ Function Run: Generate a grid Based on the Cleft File """
    #=======================================================================         
    def run(self):

        if self.read_Cleft():
            self.top.DisplayMessage("Could not read Cleft Sphere file for generating Grid", 1)
            self.top.ProcessError = True
            self.top.GenGridRunning(False)
            return

        if self.build_Grid():
            self.top.DisplayMessage("Error while generating the grid from the Spheres", 1)
            self.top.ProcessError = True
            self.top.GenGridRunning(False)
            return

        if len(self.dictSpheres) == 0:
            self.top.DisplayMessage("Could not find any Sphere(s) in the Cleft file provided", 1)
            self.top.ProcessError = True
            self.top.GenGridRunning(False)
            return
        
        '''
        if self.OutputFile != '' and self.write_Grid():
            self.top.DisplayMessage("Error while outputting the grid to a file", 1)
            self.top.ProcessError = True
            self.top.GenGridRunning(False)
            return
        
        if self.Estimate and self.estimate_Volume():
            self.top.DisplayMessage("Error while estimating the size of the grid", 1)
            self.top.ProcessError = True
            self.top.GenGridRunning(False)
            return
        '''

        self.top.GenGridRunning(False)
        
    #=======================================================================
    """ read_Cleft: reads the coordinates of the spheres of the Cleft Sphere file """
    #=======================================================================         
    def read_Cleft(self):
        
        try:
            file = open(self.CleftFile, 'r')
            PDBLines = file.readlines()
            file.close()

            #0         1         2         3         4         5         6         7                  
            #01234567890123456789012345678901234567890123456789012345678901234567890123456789
            #ATOM    247  C   SPH Z   1       7.931   2.550 -14.373  1.00  2.02 
            for Line in PDBLines:
                if Line.startswith('ATOM  '):

                    Index = Line[6:11].strip()
                    CoordX = float(Line[30:38].strip())
                    CoordY = float(Line[38:46].strip())
                    CoordZ = float(Line[46:54].strip())
                    Radius = float(Line[60:66].strip())

                    self.dictSpheres[Index] = [ Radius, [ CoordX, CoordY, CoordZ ] ] 

        except:
            return 1
        
        return 0

    #=======================================================================
    """ build_Grid: builds the grid using the Spheres """
    #=======================================================================         
    def build_Grid(self):
        
        try:

            for sph in self.dictSpheres.keys():
                
                if float(1.0 / self.Spacer) - float(int(1.0 / self.Spacer)) > 0.001:
                    xmin = float(int((self.dictSpheres[sph][1][0] - self.dictSpheres[sph][0]) / self.Spacer)) * self.Spacer;
                    ymin = float(int((self.dictSpheres[sph][1][1] - self.dictSpheres[sph][0]) / self.Spacer)) * self.Spacer;
                    zmin = float(int((self.dictSpheres[sph][1][2] - self.dictSpheres[sph][0]) / self.Spacer)) * self.Spacer;
                    xmax = float(int((self.dictSpheres[sph][1][0] + self.dictSpheres[sph][0]) / self.Spacer) + 1.0) * self.Spacer;
                    ymax = float(int((self.dictSpheres[sph][1][1] + self.dictSpheres[sph][0]) / self.Spacer) + 1.0) * self.Spacer;
                    zmax = float(int((self.dictSpheres[sph][1][2] + self.dictSpheres[sph][0]) / self.Spacer) + 1.0) * self.Spacer;

                else:
                    xmin = float(int(self.dictSpheres[sph][1][0] - self.dictSpheres[sph][0] - self.Spacer));
                    ymin = float(int(self.dictSpheres[sph][1][1] - self.dictSpheres[sph][0] - self.Spacer));   
                    zmin = float(int(self.dictSpheres[sph][1][2] - self.dictSpheres[sph][0] - self.Spacer));
                    xmax = float(int(self.dictSpheres[sph][1][0] + self.dictSpheres[sph][0] + self.Spacer + 1.0))
                    ymax = float(int(self.dictSpheres[sph][1][1] + self.dictSpheres[sph][0] + self.Spacer + 1.0))
                    zmax = float(int(self.dictSpheres[sph][1][2] + self.dictSpheres[sph][0] + self.Spacer + 1.0))


                x = xmin
                y = ymin
                z = zmin
                
                sqrrad = sph.Radius
                
                while z < zmax:
                    while y < ymax:
                        while x < xmax:

                            key  = '%8.3f' % x
                            key += '%8.3f' % y
                            key += '%8.3f' % z

                            if key not in self.dictGridPoints:

                                sqrrad = self.dictSpheres[sph][0] * self.dictSpheres[sph][0]
                                if Geometry.sqrdistance( self.dictSpheres[sph][1], [ x, y, z ] ) < sqrrad:

                                    self.dictGridPoints[key] = ''

                            x += self.Spacer

                        x = xmin
                        y += self.Spacer

                    x = xmin
                    y = ymin
                    z += self.Spacer

        except:
            return 1

        return 0

    #=======================================================================
    """ write_Grid: outputs the grid in PDB format """
    #=======================================================================         
    def write_Grid(self):

        try:
            outfile = open(self.OutputFile, 'w')
            outfile.write('REMARK ParentCleft ' + self.CleftFile + '\n') 

            i = 1
            for vertex in sorted(self.dictGridPoints.keys()):
                outfile.write('ATOM  ')
                outfile.write('%5d' % i)
                outfile.write('  C   GRD A   1    ')
                outfile.write('%8.3f' % float(vertex[0:8]))
                outfile.write('%8.3f' % float(vertex[8:16]))
                outfile.write('%8.3f' % float(vertex[16:24]))
                outfile.write('  1.00  1.00           C  ')
                outfile.write('\n')
                i += 1

            outfile.close()

        except:
            return 1

        return 0

    #=======================================================================
    """ estimate_Volume: estimates the Grid size in Angstroms (call build_Grid before) """
    #=======================================================================         
    def estimate_Volume(self):
        
        try:

            n = len(self.dictGridPoints)
            s = self.Spacer
            s3 = s * s * s

            if n > 8:
                a = 1
                b = (n-8)/4
            else:
                a = n/8
                b = 0

            #print "a", a
            #print "b", b
            #print "n", n
            #print "s", s

            v = a * s3 + b * s3

            self.top.CleftVolume.set( v )

        except:
            return 1

        return 0
