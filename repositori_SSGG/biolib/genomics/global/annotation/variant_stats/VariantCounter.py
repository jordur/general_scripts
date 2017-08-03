"""Created on Aug 13, 2013
@author: SGNET\gmarco"""


class VariantCounter(object):
    
    def __init__(self):
            self.__homo_var = 0
            self.__hetero_var = 0
            self.__total_variants = 0

    def get_total_variants(self):
        return self.__total_variants
          
    def get_homo_var(self):
        return self.__homo_var

    def get_hetero_var(self):
        return self.__hetero_var

    def del_total_variants(self):
        del self.__total_variants

    def del_homo_var(self):
        del self.__homo_var

    def del_hetero_var(self):
        del self.__hetero_var
        
    def increment_homo_var(self):
        self.__homo_var += 1
        self.__total_variants += 1
        
    def increment_hetero_var(self):
        self.__hetero_var += 1
        self.__total_variants += 1
        
    homo_var = property(get_homo_var, del_homo_var, "homo_var's docstring")
    hetero_var = property(get_hetero_var, del_hetero_var, "hetero_var's docstring")
    total_variants = property(get_total_variants, del_total_variants, "total_variants's docstring")