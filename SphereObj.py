class SphereObj:

    ''' ==================================================================================
    FUNCTION __init__: Default constructor
    =================================================================================  '''    
    def __init__(self, Radius=0.0, MaxRadius=0.0, Center=[0.0,0.0,0.0]):

        self.Radius = Radius
        self.MaxRadius = MaxRadius
        self.Center = Center

    ''' ==================================================================================
    FUNCTION Set_Center: Sets the center of the sphere
    =================================================================================  '''    
    def Set_Center(self, List):
        
        self.Center = list(List)

    ''' ==================================================================================
    FUNCTION Set_Radius: Sets the radius of the sphere
    =================================================================================  '''    
    def Set_Radius(self, val):
        
        self.Radius = val

    ''' ==================================================================================
    FUNCTION Set_MaxRadius: Sets the maximum radius of the sphere
    =================================================================================  '''    
    def Set_MaxRadius(self, val):
        
        self.MaxRadius = val

    ''' ==================================================================================
    FUNCTION Reset: Reset the sphere
    =================================================================================  '''    
    def Reset(self):
        
        self.Radius = 0.0
        self.Center = [ 0.0, 0.0, 0.0 ]

    ''' ==================================================================================
    FUNCTION Copy: Copies an instance of a class
    =================================================================================  '''    
    def Copy(self):
        
        return SphereObj(self.Radius,self.MaxRadius,self.Center)

    ''' ==================================================================================
    FUNCTION Print: Prints the instance of a class
    =================================================================================  '''    
    def Print(self):
        
        print self
        print "Radius:", self.Radius
        print "MaxRadius:", self.MaxRadius
        print "Center:", self.Center
