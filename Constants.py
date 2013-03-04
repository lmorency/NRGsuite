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

from collections import defaultdict


MAX_LIGAND_ATOMS = 100

nFlexBonds = dict(zip(('VAL','LEU','ILE','MET', 
                       'PHE','ASN','ASP','GLN',  
                       'GLU','HIS','LYS','ARG', 
                       'SER','THR','TYR','TRP', 
                       'CYS'),                  
                      (1,2,2,3,2,2,2,3,3,       
                       2,4,4,1,1,2,2,1)))

# Side-chain dihedrals build list
setDihedrals = defaultdict(list)

setDihedrals = ({ 'VAL':['N','CA','CB','CG1'],
                  'LEU':['N','CA','CB','CG',
                         'CA','CB','CG','CD1'],
                  'ILE':['N','CA','CB','CG1',
                         'CA','CB','CG1','CD1'],
                  'MET':['N','CA','CB','CG',
                         'CA','CB','CG','SD',
                         'CB','CG','SD','CE'],
                  'PHE':['N','CA','CB','CG',
                         'CA','CB','CG','CD1'],
                  'ASN':['N','CA','CB','CG',
                         'CA','CB','CG','OD1'],
                  'ASP':['N','CA','CB','CG',
                         'CA','CB','CG','OD1'],
                  'GLN':['N','CA','CB','CG',
                         'CA','CB','CG','CD',
                         'CB','CG','CD','OE1'],
                  'GLU':['N','CA','CB','CG',
                         'CA','CB','CG','CD',
                         'CB','CG','CD','OE1'],
                  'HIS':['N','CA','CB','CG',
                         'CA','CB','CG','ND1'],
                  'LYS':['N','CA','CB','CG',
                         'CA','CB','CG','CD',
                         'CB','CG','CD','CE',
                         'CG','CD','CE','NZ',],
                  'ARG':['N','CA','CB','CG',
                         'CA','CB','CG','CD',
                         'CB','CG','CD','NE',
                         'CG','CD','NE','CZ',],
                  'SER':['N','CA','CB','OG'],
                  'THR':['N','CA','CB','OG1'],
                  'TYR':['N','CA','CB','CG',
                         'CA','CB','CG','CD1'],
                  'TRP':['N','CA','CB','CG',
                         'CA','CB','CG','CD1'],
                  'CYS':['N','CA','CB','SG'] })
