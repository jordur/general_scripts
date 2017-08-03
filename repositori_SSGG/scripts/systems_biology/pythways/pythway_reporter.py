#!/usr/bin/env python

'''
Created on Aug 13, 2013
@author: gmarco
'''

#Python imports
from jinja2 import Environment, FileSystemLoader
import argparse
import os.path
import sys
import time
import datetime
import cPickle as pickle

#Self imports
from Pythway_utils import *
from Pathway import Pathway
from Ontology import Ontology

def menu ():
    
    """Input menu for user."""
    
    parser = argparse.ArgumentParser(description='Pythway reporter v0.1',
                                     epilog="################################")

    required = parser.add_argument_group('Required options')
        
    required.add_argument('-i', dest='pathway_object_file', action = 'store', 
                        help='Input Pathway pickle object file: -i /home/my_pathway.p',
                        required=True)
    
    return parser.parse_args()

def file_exists (input_file):
    
    """Checks if input_file exists"""
    
    if not os.path.exists(input_file): raise IOError("The file %s does not exist!" %input_file)

def working_message (args):
    """Returns message with the options that the
    user has chosen"""
    
    message = 'Pythway reporter v0.1'
    message += '\n=============================='
    message += '\nInput file: %s' %args.pathway_object_file
    message += '\n=============================='
    
    return message

if __name__ == '__main__':
     
    #Parse user menu options.
    args = menu()
    
    #Check input file exists
    file_exists(args.pathway_object_file)
    
    #Print working message.
    print working_message (args)
    
    #Load related pathway list from pickle file.
    related_pathway_list = pickle.load( open(args.pathway_object_file, "rb" ))
    
    now = datetime.datetime.now()
    date = now.strftime("%d-%m-%Y %H:%M")
    
    #Jinja2 test
    env = Environment(loader=FileSystemLoader('.'))
    template = env.get_template('pathway_table_template.html')
    res = template.render(pathway_list = related_pathway_list, date = date)
    
    f = open('/home/likewise-open/SGNET/gmarco/Desktop/pathway_report.html','w')
    f.write(res)
    
    print 'DONE'