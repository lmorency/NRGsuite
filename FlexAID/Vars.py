# Vars is a base class to override pickling

class Vars:

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
                
                #elif name == 'BindingSite':
                #    print v.listClefts
                
            else:
                del dict_data[k]
                
        #print dict_data
        return dict_data
    
        
    # UNPICKLING
    def __setstate__(self, dict_data):
    
        #print dict_data
        for k in dict_data.keys():
            v = dict_data[k]
            # object is a StringVar, IntVar, etc.
            if k.startswith('_'):
                svar = k[1:]
                
                try:
                    var = getattr(self.__class__,svar)
                    var.set(v)
                except:
                    pass
                                
            else:
                name = v.__class__.__name__
                #if name == 'BindingSite':
                #    print v.listClefts
                
                try:
                    var = getattr(self.__class__,k)
                    var = v
                except:
                    pass
                        
        return
        
