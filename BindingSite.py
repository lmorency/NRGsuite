import re
import copy
import CleftObj
import SphereObj

class BindingSite(object):

    def __init__(self):

        self.Type = 0

        self.Sphere = None
        self.listClefts = list()

    ''' ==================================================================================
    FUNCTION Unset: Sets the binding-site as undefined
    ================================================================================== '''
    def Unset(self):

        self.Type = 0

    ''' ==================================================================================
    FUNCTION Set_Sphere: Sets the binding-site defined by a sphere
    ================================================================================== '''
    def Set_Sphere(self):
        
        self.Type = 1
        
    ''' ==================================================================================
    FUNCTION Set_Cleft: Sets the binding-site defined by one or more cleft(s)
    ================================================================================== '''
    def Set_Cleft(self):
        
        self.Type = 2

    ''' ==================================================================================
    FUNCTION Index_Cleft: Assign an index to each cleft in the list
    ================================================================================== '''
    def Index_Cleft(self):
        
        i = 1
        for Cleft in self.listClefts:
            Cleft.Index = i
            i = i + 1
        
    ''' ==================================================================================
    FUNCTION Add_Cleft: Adds a cleft only if it doesnt exist in the list
    ================================================================================== '''
    def Add_Cleft(self, NewCleft):
        
        for Cleft in self.listClefts:
            if Cleft.CleftMD5 == NewCleft.CleftMD5:
                return
        
        self.listClefts.append(NewCleft)
        self.Index_Cleft()

    ''' ==================================================================================
    FUNCTION Get_CleftName: Gets the cleft object matching the name
    ================================================================================== '''
    def Get_CleftName(self, CleftName):
        
        for Cleft in self.listClefts:
            if Cleft.CleftName == CleftName:
                return Cleft

        return None

    ''' ==================================================================================
    FUNCTION Get_SortedCleftNames: Returns a list of clefts sorted by CleftName
    ================================================================================== '''
    def Get_SortedCleftNames(self):
        
        convert = lambda text: int(text) if text.isdigit() else text.lower()
        alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
        
        CleftNames = list()
        for Cleft in self.listClefts:
            CleftNames.append(str(Cleft.CleftName))

        return sorted(CleftNames, key = alphanum_key)

    ''' ==================================================================================
    FUNCTION Remove_Cleft: Remove a cleft only if it exists in the list
    ================================================================================== '''
    def Remove_Cleft(self, RemCleft):
        
        for Cleft in self.listClefts:
            if Cleft.CleftMD5 == RemCleft.CleftMD5:
                self.listClefts.remove(Cleft)
                break

        self.Index_Cleft()

    ''' ==================================================================================
    FUNCTION Remove_CleftName: Remove a cleft only if it exists in the list
    ================================================================================== '''
    def Remove_CleftName(self, RemCleftName):
        
        for Cleft in self.listClefts:
            if Cleft.CleftName == RemCleftName:
                self.listClefts.remove(Cleft)
                break

        self.Index_Cleft()

    ''' ==================================================================================
    FUNCTION Clear_Cleft: Clears the clefts from the list
    ================================================================================== '''
    def Clear_Cleft(self):
        
        del self.listClefts[:]

    ''' ==================================================================================
    FUNCTION Count_Cleft: Count the number of clefts in the bindingsite
    ================================================================================== '''
    def Count_Cleft(self):
        
        return len(self.listClefts)

    ''' ==================================================================================
    FUNCTION Clear_Sphere: Clears the Spheres object
    ================================================================================== '''
    def Clear_Sphere(self):
        
        self.Sphere = None

    ''' ==================================================================================
    FUNCTION Clear: Clears the binding-site
    ================================================================================== '''
    def Clear(self):
        
        self.Clear_Cleft()
        self.Clear_Sphere()
        self.Unset()

    ''' ==================================================================================
    FUNCTION Copy: Copies an instance of a class
    =================================================================================  '''    
    def Copy(self):
        
        return copy.deepcopy(self)
