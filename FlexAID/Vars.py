import copy

# Vars is a base class to override pickling

class Vars:
    
    dict_vars = {}
    
    # PICKLING
    # All class objects are pickled normally except for
    # TkinterVars. The object methods are skipped.
    def __getstate__(self):
    
        dict_data = self.__class__.__dict__.copy()

        for k in dict_data.keys():
            v = dict_data[k]
            # object is not a method
            if not k.startswith('__'):
                # object is a StringVar, IntVar, etc.
                name = v.__class__.__name__
                if name == 'StringVar' or name == 'IntVar' or \
                   name == 'BooleanVar' or name == 'DoubleVar':
                   
                    del dict_data[k]
                    
                    try:
                        dict_data['_' + k] = v.get()
                    except:
                        pass
                        
                elif k == 'dict_vars':
                    del dict_data[k]
                    
            else:
                del dict_data[k]
        
        # Copy instance variables
        dict_data = dict(dict_data.items() + self.__dict__.items())
        
        print "__getstate__"
        print "dict_data", dict_data
        
        return dict_data
    
        
    # UNPICKLING
    # First unpickle objects, list and dict
    # Then Vars because they are bound to controls
    def __setstate__(self, dict_data):

        Vars.dict_vars.clear()

        for k in dict_data.keys():
            v = dict_data[k]

            # object is a StringVar, IntVar, etc.            
            if k.startswith('_'):
                svar = k[1:]
                Vars.dict_vars[svar] = v
                
                del dict_data[k]

        self.__dict__.update(dict_data)
        
        return
    
    # refresh the values of Vars to trigger Tracers
    def refresh(self):
    
        for k, v in Vars.dict_vars.iteritems():
            var = self.__class__.__dict__.get(k)
            if var is not None:
                var.set(v)

        self.dict_vars.clear()

        
