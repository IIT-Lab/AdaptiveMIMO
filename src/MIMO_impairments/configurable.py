class Configurable:
    def load_config(self, ConfigClass):
        if ConfigClass:
            conf_attr = [k for k in dir(ConfigClass) if not k.startswith('_')]
            obj_attr = [k for k in dir(self ) if not k.startswith('_')]

            for obj_k in obj_attr:
                if(obj_k in conf_attr):
                    if obj_k in self.__class__.__dict__ and isinstance(getattr(self.__class__, obj_k), property):
                        getattr(self.__class__, obj_k).fset(self, getattr(ConfigClass, obj_k))
                    else:
                        setattr(self, obj_k, getattr(ConfigClass, obj_k))
        else:
            return


