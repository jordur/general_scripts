#!/usr/bin/env python

"""
Created on Aug 12, 2013

@author: SGNET\gmarco
"""

import os
import base64
from cogent.db.ensembl import HostAccount
from cogent.db.ensembl import Genome
from SOAPpy import SOAPProxy
from Pathway import Pathway
from Ontology import Ontology

if __name__ == '__main__':
    pass


#Wikipathways WS API
url = 'http://www.wikipathways.org/wpi/webservice/webservice.php'
namespace = 'http://www.wikipathways.org/webservice'
server = SOAPProxy(url, namespace)


def is_subset(list1, list2):
    value = set(list1) <= set(list2)
    return value


def is_contained(list1, list2):
    if [i for i in list1 if i in list2]:
        return True
    else:
        return False


def common_genes(list1, list2):

    common_gene_list = list(set(list1).intersection(list2))
    common_gene_list = set(common_gene_list)
    
    return common_gene_list


def ensembl_to_hgnc(gene_list):
    account = HostAccount('blackrussian', 'bioinfo', 'A29bcd1234#', port=3306)
    human = Genome('human', Release=73, account=account)

    hgnc_list = []

    for gene in gene_list:
        hgnc_list.append(human.getGeneByStableId(StableId=gene).Symbol)

    hgnc_list = set(hgnc_list)

    return hgnc_list


def hgnc_to_ensembl_id(hgnc_list):
    account = HostAccount('blackrussian', 'bioinfo', 'A29bcd1234#', port=3306)
    human = Genome('human', Release=73, account=account)

    ensembl_stable_id_list = []

    for gene in hgnc_list:
        gene_query = human.getGenesMatching(Symbol=gene)
        for gene_obj in gene_query:
            if gene_obj.Symbol == gene:
                ensembl_stable_id_list.append(gene_obj.StableId)

    #Remove duplicates
    ensembl_stable_id_list = set(ensembl_stable_id_list)

    #Keep list elements starting with 'ENSG'
    ensembl_stable_id_list = [x for x in ensembl_stable_id_list if x.startswith('ENSG')]

    return ensembl_stable_id_list


def print_output(ws_output):
    #Loops through a list of dictionary items
    index = 1
    for item in ws_output:
        #calls select dictionary keys to print out values
        output = '   ' + str(index) + ')' + 'species:' + item['species'] + '\t '
        output += 'id:' + item['id'] + '\t ' + 'name:' + item['name']
        print output
        index += 1


def get_pathway_info(pwid):
    pathway_info = server.getPathwayInfo(pwId=pwid)
    pathway_id = pathway_info['id']
    pathway_name = pathway_info['name']
    if not pathway_id:
        print 'no'
    if not pathway_name:
        print 'no'

    return pathway_info


def get_pathway_xref_list(pathway, code):
    pathway_id = pathway.get_pw_id()
    ext_ref_code = code
    server.config.argsOrdering = {'getXrefList': ('pwId', 'code')}
    xref_list = server.getXrefList(pwId=pathway_id, code=ext_ref_code)
    #server_res = server.getXrefList(pwId, code) 

    return xref_list


def get_pathway_as(pathway, export_type, path):
    pathway_id = pathway.get_pw_id()
    pathway_revision = pathway.get_revision()
    try:
        server.config.argsOrdering = {'getPathwayAs': ('fileType', 'pwId', 'revision')}
        server_res = server.getPathwayAs(fileType=export_type, pwId=pathway_id, revision=pathway_revision)

        decoded_content = base64.b64decode(server_res)

        pathway_file_name = pathway_id + '.' + export_type
        gpml_output = os.path.join(path, pathway_file_name)

        with open(gpml_output, 'w') as f:
            f.write(decoded_content)

    except Exception, e:
        print e

    return 1


def create_ontology_list_by_pathway(pathway):
    ontology_list = []
    pathway_id = pathway.get_pw_id()

    server_res = server.getOntologyTermsByPathway(pwId=pathway_id)

    for element in server_res:
        try:
            ontology_id = element['id']
            name = element['name']
            ontology = element['ontology']

            #Casting Str type
            ontology_id = str(ontology_id)
            name = str(name)
            ontology = str(ontology)

            ontology_list.append(Ontology(ontology_id, name, ontology))

        except TypeError:
            continue

    return ontology_list


def find_pathways_by_xref(ids, codes):
    ensembl_id = ids
    ext_ref_code = codes

    server.config.argsOrdering = {'findPathwaysByXref': ('ids', 'codes')}
    server_res = server.findPathwaysByXref(ids=ensembl_id, codes=ext_ref_code)

    return server_res


def create_pathway_from_ensembl_gene_id(ensembl_gene_id):
    sc = 'En'
    pathway_list = []

    # Define the order of args: needed for this service
    server.config.argsOrdering = {'findPathwaysByXref': ('ids', 'codes')}

    # findPathwaysByXref (multiple args; returns list of dictionary references)    
    server_res = server.findPathwaysByXref(codes=sc, ids=ensembl_gene_id)

    #print type(server_res)

    if isinstance(server_res, list):
        for element in server_res:
            pw_id = str(element['id'])
            name = str(element['name'])
            revision = int(float(element['revision']))
            species = str(element['species'])

            #Create pathway object
            pathway_object = Pathway(pw_id, name, revision, species)

            if pathway_object not in pathway_list:
                pathway_list.append(pathway_object)

    else:
        pw_id = server_res.id
        name = server_res.name
        revision = int(server_res.revision)
        species = server_res.species

        #Create pathway object
        pathway_object = Pathway(pw_id, name, revision, species)

        if pathway_object not in pathway_list:
            pathway_list.append(pathway_object)
            #     for element in server_res:
            #         if not isinstance(element, primitiveTypes):
            #             pw_id = server_res.id
            #             name = server_res.name
            #             revision = int(server_res.revision)
            #             species = server_res.species
            #         else:
            #             pw_id = str(element['id'])
            #             name = str(element['name'])
            #             revision  = int(float(element['revision']))
            #             species =  str(element['species'])

            #Add it to list if it doesn't exist.

    return pathway_list