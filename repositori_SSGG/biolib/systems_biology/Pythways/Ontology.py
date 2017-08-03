#!/usr/bin/env python

"""
Created on Aug 21, 2013

@author: SGNET\gmarco
"""


class Ontology(object):

    def __init__(self, ontology_id, name, ontology):
        """Ontology class constructor"""
        self.ontology_id = ontology_id
        self.name = name
        self.ontology = ontology

    def get_id(self):
        return self.ontology_id

    def get_name(self):
        return self.name

    def get_ontology(self):
        return self.ontology

    def set_ontology_id(self, value):
        self.ontology_id = value

    def set_name(self, value):
        self.name = value

    def set_ontology(self, value):
        self.ontology = value

    def del_ontology_id(self):
        del self.ontology_id

    def del_name(self):
        del self.name

    def del_ontology(self):
        del self.ontology