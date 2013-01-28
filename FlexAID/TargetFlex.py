import copy

class TargetFlex:

    def __init__(self):

        self.listSideChain = list()

    ''' ==================================================================================
    FUNCTION Add_SideChain: Adds a flexible side chain to array list
    ==================================================================================  '''           
    def Add_SideChain(self, Residue):

        if self.listSideChain.count(Residue) == 0:
            self.listSideChain.append(Residue)

    ''' ==================================================================================
    FUNCTION Remove_SideChain: Removes a flexible side chain to array list
    ==================================================================================  '''           
    def Remove_SideChain(self, Residue):

        if self.listSideChain.count(Residue) != 0:
            self.listSideChain.remove(Residue)

    ''' ==================================================================================
    FUNCTION Clear_SideChain: Empties the list of flexible side chain
    ==================================================================================  '''           
    def Clear_SideChain(self):

        del self.listSideChain[:]

    ''' ==================================================================================
    FUNCTION Count_SideChain: Returns the number of flexible side chain
    ==================================================================================  '''           
    def Count_SideChain(self):

        return len(self.listSideChain)

    ''' ==================================================================================
    FUNCTION Copy: Copies an instance of a class
    =================================================================================  '''    
    def Copy(self):
        
        return copy.deepcopy(self)
