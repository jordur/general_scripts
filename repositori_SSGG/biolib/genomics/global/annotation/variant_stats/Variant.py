"""
Created on Aug 26, 2013

@author: SGNET\gmarco
"""


class Variant(object):
    """
    classdocs
    """

    def __init__(self, chr, pos, ref_allele, var_allele):
        self.chr = chr
        self.pos = int(pos)
        self.ref_allele = ref_allele
        self.var_allele = var_allele
        self.__forward_in_three_range = 0
        self.__forward_in_five_range = 0
        self.__reverse_in_three_range = 0
        self.__reverse_in_five_range = 0
        self.__total = 0

    def get_chr(self):
        return self.__chr

    def get_pos(self):
        return self.__pos

    def get_ref_allele(self):
        return self.__ref_allele

    def get_var_allele(self):
        return self.__var_allele

    def get_total(self):
        return self.__total

    def get_forward_in_three_range(self):
        return self.__forward_in_three_range

    def get_forward_in_five_range(self):
        return self.__forward_in_five_range

    def get_reverse_in_three_range(self):
        return self.__reverse_in_three_range

    def get_reverse_in_five_range(self):
        return self.__reverse_in_five_range

    def set_chr(self, value):
        self.__chr = value

    def set_pos(self, value):
        self.__pos = value

    def set_ref_allele(self, value):
        self.__ref_allele = value

    def set_var_allele(self, value):
        self.__var_allele = value

    def set_forward_in_three_range(self, value):
        self.__forward_in_three_range = value

    def set_forward_in_five_range(self, value):
        self.__forward_in_five_range = value

    def set_reverse_in_three_range(self, value):
        self.__reverse_in_three_range = value

    def set_reverse_in_five_range(self, value):
        self.__reverse_in_five_range = value

    def del_chr(self):
        del self.__chr

    def del_pos(self):
        del self.__pos

    def del_ref_allele(self):
        del self.__ref_allele

    def del_var_allele(self):
        del self.__var_allele

    def del_total(self):
        del self.__total

    def del_forward_in_three_range(self):
        del self.__forward_in_three_range

    def del_forward_in_five_range(self):
        del self.__forward_in_five_range

    def del_reverse_in_three_range(self):
        del self.__reverse_in_three_range

    def del_reverse_in_five_range(self):
        del self.__reverse_in_five_range

    def increment_total(self):
        self.__total += 1

    def increment_forward_in_three_range(self):
        self.__forward_in_three_range += 1
        self.__total += 1

    def increment_forward_in_five_range(self):
        self.__forward_in_five_range += 1
        self.__total += 1

    def increment_reverse_in_three_range(self):
        self.__reverse_in_three_range += 1
        self.__total += 1

    def increment_reverse_in_five_range(self):
        self.__reverse_in_five_range += 1
        self.__total += 1

    total = property(get_total, del_total, "total's docstring")
    chr = property(get_chr, set_chr, del_chr, "chr's docstring")
    pos = property(get_pos, set_pos, del_pos, "pos's docstring")
    ref_allele = property(get_ref_allele, set_ref_allele, del_ref_allele, "ref_allele's docstring")
    var_allele = property(get_var_allele, set_var_allele, del_var_allele, "var_allele's docstring")
    forward_in_three_range = property(get_forward_in_three_range, set_forward_in_three_range,
                                      del_forward_in_three_range, "forward_in_three_range's docstring")
    forward_in_five_range = property(get_forward_in_five_range, set_forward_in_five_range, del_forward_in_five_range,
                                     "forward_in_five_range's docstring")
    reverse_in_three_range = property(get_reverse_in_three_range, set_reverse_in_three_range,
                                      del_reverse_in_three_range, "reverse_in_three_range's docstring")
    reverse_in_five_range = property(get_reverse_in_five_range, set_reverse_in_five_range, del_reverse_in_five_range,
                                     "reverse_in_five_range's docstring")
        

    