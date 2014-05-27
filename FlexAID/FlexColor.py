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
@title: Color.py

@summary: Module that define all colors used from Pymol simulation.

@contain: Get_Pymol_HeatColor, Get_RGB_HeatColor, GetHeatColorList, CreateNoList

@organization: Najmanovich Research Group
@creation date:  oct. 13, 2010
'''

from pymol import cmd

NBCOLOR = 20

'''
@summary: Get_Pymol_HeatColor: Based on a heat scale color chart, return the color
                        related on the colorId number pass has argument.
                        Related link: http://www.pymolwiki.org/index.php/Color_Values
@param pColorId: Integer
@return: String - Color name used by Pymol
'''
def Get_Pymol_HeatColor(pColorId):

    #From RED=0 to BLUE=19
    ColorList = ['red',
                 'br8',
                 'tv_red',
                 'oxygen',
                 'iron',
                 'tv_orange',
                 'sulfur',
                 'gold',
                 'yelloworange',
                 'neodymium',
                 'limon',
                 'chartreuse',
                 'tv_green',
                 'limegreen',
                 'teal',
                 'rhodium',
                 'slate',
                 'tv_blue',
                 'blue',
                 'density']

    if (pColorId<0) or (pColorId>19):
        return 'grey50'
    else:
        return ColorList[pColorId]


def Get_RGB_HeatColor(pColorId):

    #From RED=0 to BLUE=19
    ColorList = ['#FF0000',
                 '#E61A33',
                 '#FF1515',
                 '#FF4D4D',
                 '#E06633',
                 '#FF8C26',
                 '#E6C740',
                 '#FFD124',
                 '#FFDE5E',
                 '#C7FFC7',
                 '#BFFF40',
                 '#80FF00',
                 '#33FF33',
                 '#00FF80',
                 '#00BFBF',
                 '#0A7DAB',
                 '#8080FF',
                 '#4D4DFF',
                 '#0000FF',
                 '#1A1A99']

    if (pColorId<0) or (pColorId>19):
        return '#7D7D7D'
    else:
        return ColorList[pColorId]


'''
@summary: GetHeatColorList: Based on a heat scale color chart, return a color
                            list depending on the total of colors.
              Related link: http://www.pymolwiki.org/index.php/Color_Values
@param pTotColor: Integer
@param boolRGB: Boolean - If RGB or NOT
@return: List - List of colors (Name or RGB)
'''
def GetHeatColorList(pTotColor, boolRGB):

    ColorList = []
    noList = list()
    TotalColorList = 20     # Total of colors available

    if (pTotColor == 1):
        if boolRGB:
            ColorList.append(Get_RGB_HeatColor(0))
        else:
            ColorList.append(Get_Pymol_HeatColor(0))

    elif (pTotColor > 1) and (pTotColor < 21):

        noList = CreateNoList(pTotColor, TotalColorList)

        if boolRGB:
            for no in noList:
                ColorList.append(Get_RGB_HeatColor(no))

        else:
            for no in noList:
                ColorList.append(Get_Pymol_HeatColor(no))

    else:       # Number of color displayed higher then 20

        if boolRGB:

            for i in range(0, pTotColor):
                ColorList.append(Get_RGB_HeatColor(i))

        else:

            for i in range(0, pTotColor):
                ColorList.append(Get_Pymol_HeatColor(i))

    return ColorList


'''
@summary: CreateNoList: Create a list of number to get well repartitioned color list.
@param pTotColor: Integer
@param pTotalColorList: Integer
@return: noList - List of numbers
'''
def CreateNoList(pTotColor, pTotalColorList):

    noList = list()

    Modulo = (pTotalColorList-1) % (pTotColor-1)
    Partition = (pTotalColorList-Modulo-1)/ (pTotColor-1)

    stepStart = 0
    stepEnd = pTotalColorList - 1

    for i in range(0, pTotColor):

        if ((i % 2) == 0):
            noList.append(stepStart)
            stepStart = stepStart + Partition
        else:
            noList.append(stepEnd)
            stepEnd = stepEnd - Partition

    noList.sort()

    return noList
