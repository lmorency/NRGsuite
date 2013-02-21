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

# coding: utf-8 
'''
@title: Geometry.py

@summary: Module that define some tools used by FlexAID.

@contain: distance, angle, dihedralAngle

@organization: Najmanovich Research Group
@creation date:  oct. 13, 2010
'''

import math
    
'''******************************************************************************
  SUBROUTINE middle: Calculates the center of geometry between 2 points
******************************************************************************'''
def middle(pointA, pointB):
  
    return [ ( pointA[0] + pointB[0] ) / 2.0 ,
             ( pointA[1] + pointB[1] ) / 2.0 ,
             ( pointA[2] + pointB[2] ) / 2.0 
           ]
    
'''******************************************************************************
  SUBROUTINE distance: Calculates the cartesian distance between two point in 
                       3 dimensions.
******************************************************************************'''
def distance(pointA, pointB):
  
    d = 0.0
  
    for i in range(0, 3):
        d += (pointA[i]-pointB[i])**2
    
    d = math.sqrt(d)
    
    return d

'''******************************************************************************
  SUBROUTINE distance: Calculates the cartesian squared distance between two point in 
                       3 dimensions.
******************************************************************************'''
def sqrdistance(pointA, pointB):
  
    d = 0.0
  
    for i in range(0, 3):
        d += (pointA[i] - pointB[i])**2
    
    return d

'''*******************************************************************************
  SUBROUTINE angle: Calculates the valence angle between three consecutives atoms.
*******************************************************************************'''
def angle(pointA, pointB, pointC):

    cosa = 0.0
    absu = 0.0
    absv = 0.0
    
    for i in range(0, 3):
        cosa += (pointA[i]-pointB[i])*(pointC[i]-pointB[i])
        absu += (pointA[i]-pointB[i])*(pointA[i]-pointB[i])
        absv += (pointC[i]-pointB[i])*(pointC[i]-pointB[i])
    
    cosa = cosa / math.sqrt(absu*absv)
    cosa = math.acos(cosa)*(180.0/math.pi)
    
    return cosa
 
'''*******************************************************************************
  SUBROUTINE dihedralAngle: Calculates the torsional angle between 4 consecutives atoms
*******************************************************************************'''
def dihedralAngle(pointA, pointB, pointC, pointD):
  
    # Init the arrays of 3 floats
    t = [0.0, 0.0, 0.0]
    u = [0.0, 0.0, 0.0]
    w = [0.0, 0.0, 0.0]
    m = [0.0, 0.0, 0.0]
    n = [0.0, 0.0, 0.0]
    v = [0.0, 0.0, 0.0]
    
    # Init the float values
    absm = 0.0
    absn = 0.0
    absv = 0.0
    absu = 0.0
    costheta = 0.0
    theta = 0.0
    q = 0.0
    
    for i in range(0, 3):
        t[i] = pointA[i]-pointB[i]
        u[i] = pointC[i]-pointB[i]
        w[i] = pointD[i]-pointB[i]
    
    
    m[0] = (t[1]*u[2]) - (t[2]*u[1])
    m[1] = (t[2]*u[0]) - (t[0]*u[2])
    m[2] = (t[0]*u[1]) - (t[1]*u[0])
    
    n[0] = (w[1]*u[2]) - (w[2]*u[1])
    n[1] = (w[2]*u[0]) - (w[0]*u[2])
    n[2] = (w[0]*u[1]) - (w[1]*u[0])
    
    absm = math.sqrt((m[0]*m[0])+(m[1]*m[1])+(m[2]*m[2]))
    absn = math.sqrt((n[0]*n[0])+(n[1]*n[1])+(n[2]*n[2]))
    
    costheta = ((m[0]*n[0]) + (m[1]*n[1]) + (m[2]*n[2]))/(absm*absn)
    
    theta = math.acos(costheta)
    
    v[0] = (m[1]*n[2]) - (m[2]*n[1])
    v[1] = (m[2]*n[0]) - (m[0]*n[2])
    v[2] = (m[0]*n[1]) - (m[1]*n[0])
    
    absv = math.sqrt((v[0]*v[0])+(v[1]*v[1])+(v[2]*v[2]))
    absu = math.sqrt((u[0]*u[0])+(u[1]*u[1])+(u[2]*u[2]))
    
    q = ((v[0]*u[0]) + (v[1]*u[1]) + (v[2]*u[2]))/(absv*absu)
    
    theta = q*theta*(180.0/math.pi)
    
    return theta

'''
@summary: SUBROUTINE buildcc: builds the cartesianlen( coordinates of the tot atoms
          present in array list according to the reconstruction data.

@return: PDBCoord for each atom of the ligand (dictionary)
'''
def buildcc(ListAtom,RecAtom,DisAngDih,Ori):

    tot = len(ListAtom)
    PDBCoord = {}       

    # init the list (x, y, z)
    x = [0.0, 0.0, 0.0, 0.0]
    y = [0.0, 0.0, 0.0, 0.0]
    z = [0.0, 0.0, 0.0, 0.0]

    for an in range(0, tot):
        NoAtom = ListAtom[an] #Number of the atom treated

        for i in range(1, 4):
            #print('No Atom: ' + str(NoAtom) + ' -> ' + str(self.RecAtom[NoAtom][i-1]))
            j = RecAtom[NoAtom][i - 1]
            #print('j: ' + str(j))

            if(j != 0):
                x[i] = float(PDBCoord[j][0])
                y[i] = float(PDBCoord[j][1])
                z[i] = float(PDBCoord[j][2])
            elif(i == 1):
                x[i] = 1.0 + float(Ori[0])
                y[i] = 0.0 + float(Ori[1])
                z[i] = 0.0 + float(Ori[2])
            elif(i == 3):
                x[i] = 0.0 + float(Ori[0])
                y[i] = 1.0 + float(Ori[1])
                z[i] = 0.0 + float(Ori[2])
            else:
                x[i] = 0.0 + float(Ori[0])
                y[i] = 0.0 + float(Ori[1])
                z[i] = 0.0 + float(Ori[2])
        # END of FOR(i)         

        a = y[1] * (z[2] - z[3]) + y[2] * (z[3] - z[1]) + y[3] * (z[1] - z[2])
        b = z[1] * (x[2] - x[3]) + z[2] * (x[3] - x[1]) + z[3] * (x[1] - x[2])
        c = x[1] * (y[2] - y[3]) + x[2] * (y[3] - y[1]) + x[3] * (y[1] - y[2])
        op = math.sqrt((a * a) + (b * b) + (c * c))

        cx = float(a) / op
        cy = float(b) / op
        cz = float(c) / op
        #print('cx= ' + str(cx) + ' cy= ' + str(cy) + ' cz= ' + str(cz))

        a = x[2] - x[1]
        b = y[2] - y[1]
        c = z[2] - z[1]

        d = float(1.0) / (math.sqrt((a * a) + (b * b) + (c * c)))

        #print('d : ' + str(d))

        op = float(DisAngDih[NoAtom][0]) * d
        xn = a * op
        yn = b * op
        zn = c * op
        #print('d= ' + str(d) + ' op= ' + str(op) + ' xn= ' + str(xn) + ' yn= ' + str(yn) + ' zn= ' + str(zn))

        a = cx * cx
        b = cy * cy
        c = cz * cz
        #print('ang= ' + str(self.DisAngDih[NoAtom][1]))

        angPI = float(DisAngDih[NoAtom][1]) * math.pi / 180.0
        ct = math.cos(angPI)
        st = -1.0 * (math.sin(angPI))

        op = 1.0 - ct
       # print('ct= ' + str(ct) + ' st= ' + str(st) + ' op= ' + str(op))

        xk = (cx * cz * op - cy * st) * zn + ((1.0 - a) * ct + a) * xn + (cx * cy * op + cz * st) * yn
        yk = (cy * cx * op - cz * st) * xn + ((1.0 - b) * ct + b) * yn + (cy * cz * op + cx * st) * zn
        zk = (cz * cy * op - cx * st) * yn + ((1.0 - c) * ct + c) * zn + (cz * cx * op + cy * st) * xn
        #print('xk=' + str(xk) + ' yk=' + str(yk) + ' zk=' + str(zk))
        #print('dih= ' + str(self.DisAngDih[NoAtom][2]))

        dihPI = float(DisAngDih[NoAtom][2]) * math.pi / 180.0
        ct = math.cos(dihPI)
        st = math.sin(dihPI)

        op = 1.0 - ct

        cx = (x[2] - x[1]) * d
        cy = (y[2] - y[1]) * d
        cz = (z[2] - z[1]) * d  

        a = cx * cx
        b = cy * cy
        c = cz * cz
        #print('a= ' + str(a) + ' b= ' + str(b) + ' c= ' + str(c))

        x[0] = (((cx * cz * op) - (cy * st)) * zk) + ((((1.0 - a) * ct) + a) * xk) + (((cx * cy * op) + (cz * st)) * yk) + x[1]
        y[0] = (((cy * cx * op) - (cz * st)) * xk) + ((((1.0 - b) * ct) + b) * yk) + (((cy * cz * op) + (cx * st)) * zk) + y[1]
        z[0] = (((cz * cy * op) - (cx * st)) * yk) + ((((1.0 - c) * ct) + c) * zk) + (((cz * cx * op) + (cy * st)) * xk) + z[1]

        #print('NO Atom: ' + str(NoAtom) + ' X: ' + str(x[0]) + ' Y: ' + str(y[0]) + ' Z: ' + str(z[0]))

        #3 floating numbers        
        PDBCoord[NoAtom] = [x[0], y[0], z[0]]              


    #END of FOR(an)
    return PDBCoord    

'''
@summary: SUBROUTINE rmsd: calculates RMSD between predicted and reference
'''
def rmsd(dictCoord, dictCoordRef):
    
    tot = 0
    sum = 0.0

    try:
        for index in dictCoord.keys():
            sum += sqrdistance( dictCoord[index], dictCoordRef[index] )
            tot += 1

        sum /= float(tot)
        
        return ( math.sqrt(sum) )

    except:

        return 'N/A'
