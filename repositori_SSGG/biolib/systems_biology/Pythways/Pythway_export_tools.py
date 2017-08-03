#!/usr/bin/env python

"""
Created on Aug 30, 2013
@author: gmarco
"""

import os
import datetime
import xmlrpclib
import shutil

from jinja2 import Environment, FileSystemLoader


def createdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

    return path


def removedir(path):
    if os.path.exists(path):
        shutil.rmtree(path)


def list_gpml_files(pathway_export_path, pathway_export_format):
    gpml_file_list = []
    for file_item in os.listdir(pathway_export_path):
        if file_item.endswith('.' + pathway_export_format):
            gpml_file_list.append(file_item)

    return gpml_file_list


def pathvisio_rpc_export(pathway_export_path, gpml_file_list, pathway_html_path):
    server = xmlrpclib.ServerProxy("http://10.0.0.30:7777")
    pathways_export_html_path = pathway_html_path + '/pathways'
    createdir(pathways_export_html_path)
    for gpml_file in gpml_file_list:
        gpml_file_path = pathway_export_path + '/' + gpml_file
        gpml_file_abs_path = os.path.abspath(gpml_file_path)
        server.PathVisio.exportPathwayHtml(gpml_file_abs_path, "", pathways_export_html_path)
    print 'Generated HTML Pathway graphics in HTML.'


def ontology_pathway_html_report(pathway, ontology_list, ontology_html_path, script_path):
    #Set report date time format
    now = datetime.datetime.now()
    date = now.strftime("%d-%m-%Y %H:%M")

    #Open html and set ontology_report
    html_output = ontology_html_path + '/{0}_ontology_report.html'.format(pathway.get_pw_id())
    pathway.set_ontology_report('./ontology/{0}_ontology_report.html'.format(pathway.get_pw_id()))

    #Jinja2 load path, template and variables, render template
    env = Environment(loader=FileSystemLoader(script_path))
    template = env.get_template('ontology_template.html')
    res = template.render(pathway=pathway, ontology_list=ontology_list, date=date).encode('utf-8')

    #Write rendered template to html, close file
    f = open(html_output, 'w')
    f.write(res)
    f.close()


def jinja_html_report(related_ensembl_pathways, pathway_html_path, script_path):
    #Set report date time format
    now = datetime.datetime.now()
    date = now.strftime("%d-%m-%Y %H:%M")

    #Jinja2 load path, template and variables
    env = Environment(loader=FileSystemLoader(script_path))
    template = env.get_template('pathway_table_template.html')
    res = template.render(pathway_list=related_ensembl_pathways, date=date).encode('utf-8')

    #Open, write and close html
    html_output = pathway_html_path + '/pathway_report.html'
    f = open(html_output, 'w')
    f.write(res)
    f.close()

    print 'Generated HTML report: {0}'.format(html_output)
