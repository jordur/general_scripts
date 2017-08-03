"""
Created on Aug 13, 2013

@author: SGNET\gmarco
"""


class Sample(object):
    """
    classdocs
    """

    def __init__(self, name, genotype_col, depth_col, ratio_col):
        self.name = name
        self.genotype_col = genotype_col
        self.depth_col = depth_col
        self.ratio_col = ratio_col

    def get_name(self):
        return self.__name

    def get_depth_col(self):
        return self.__depth_col

    def get_ratio_col(self):
        return self.__ratio_col

    def get_genotype_col(self):
        return self.__genotype_col

    def get_variant_counter(self):
        return self.__variant_counter

    def get_consequence_counter(self):
        return self.__consequence_counter

    def get_consequence_col(self):
        return self.__consequence_col

    def set_name(self, value):
        self.__name = value

    def set_depth_col(self, value):
        self.__depth_col = value

    def set_ratio_col(self, value):
        self.__ratio_col = value

    def set_genotype_col(self, value):
        self.__genotype_col = value

    def set_variant_counter(self, value):
        self.__variant_counter = value

    def set_consequence_counter(self, value):
        self.__consequence_counter = value

    def set_consequence_col(self, value):
        self.__consequence_col = value

    def del_name(self):
        del self.__name

    def del_depth_col(self):
        del self.__depth_col

    def del_ratio_col(self):
        del self.__ratio_col

    def del_genotype_col(self):
        del self.__genotype_col

    def del_variant_counter(self):
        del self.__variant_counter

    def del_consequence_counter(self):
        del self.__consequence_counter

    def del_consequence_col(self):
        del self.__consequence_col

    name = property(get_name, set_name, del_name, "name's docstring")
    depth_col = property(get_depth_col, set_depth_col, del_depth_col, "depth_col's docstring")
    ratio_col = property(get_ratio_col, set_ratio_col, del_ratio_col, "ratio_col's docstring")
    genotype_col = property(get_genotype_col, set_genotype_col, del_genotype_col, "genotype_col's docstring")