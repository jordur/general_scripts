#!/usr/bin/env python

"""
Created on Aug 9, 2013

@author: SGNET\gmarco
"""


class Pathway(object):
    """Pathway class object"""

    def __init__(self, pw_id, name, revision, species):
        """Pathway class constructor"""
        self.pw_id = pw_id
        self.name = name
        self.revision = revision
        self.species = species
        self.ontology_report = 'N/A'
        self.ensembl_gene_list = []
        self.common_ensembl_gene_list = []
        self.common_hgnc_list = []
        self.ontology_list = []
        self.html_report = ''

    #GETs
    def get_pw_id(self):
        return self.pw_id

    def get_name(self):
        return self.name

    def get_revision(self):
        return self.revision

    def get_species(self):
        return self.species

    def get_ensembl_list(self):
        return self.ensembl_gene_list

    def get_common_ensembl_gene_list(self):
        return self.common_ensembl_gene_list

    def get_common_hgnc_list(self):
        return self.common_hgnc_list

    def get_ontology_list(self):
        return self.ontology_list

    def get_html_report(self):
        return self.html_report

    def get_ontology_report(self):
        return self.ontology_report

    #SETs
    def set_name(self, value):
        self.name = value

    def set_revision(self, value):
        if not isinstance(value, int):
            raise ValueError("Revision is not an integer !")
        self.revision = value

    def set_species(self, value):
        self.species = value

    def set_pw_id(self, value):
        self.pw_id = value

    def set_ensembl_list(self, value):
        self.ensembl_gene_list = value

    def set_common_ensembl_gene_list(self, value):
        self.common_ensembl_gene_list = value

    def set_common_hgnc_list(self, value):
        self.common_hgnc_list = value

    def set_ontology_list(self, value):
        self.ontology_list = value

    def set_html_report(self, value):
        self.html_report = value

    def set_ontology_report(self, value):
        self.ontology_report = value

    #DELs
    def del_pw_id(self):
        del self.pw_id

    def del_name(self):
        del self.name

    def del_revision(self):
        del self.revision

    def del_species(self):
        del self.species

    def del_ensembl_list(self):
        del self.ensembl_gene_list

    def del_common_ensembl_gene_list(self):
        del self.common_ensembl_gene_list

    def del_common_hgnc_list(self):
        del self.common_hgnc_list

    def del_ontology_list(self):
        del self.ontology_list

    def del_html_report(self):
        del self.html_report

    def del_ontology_report(self):
        del self.ontology_report

    def __repr__(self):
        #print self.pwId+"\t"+ self.name + "\t" + str(self.revision) + "\t" + self.species + "\n" 
        return "%s, %s, %i, %s" % (self.pw_id, self.name, self.revision, self.species)

    def __eq__(self, other):
        return self.pw_id == other.pw_id

    def __hash__(self):
        return hash(('pwid', self.pw_id))