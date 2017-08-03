#!/usr/bin/env python

"""
Created on Aug 13, 2013
@author: gmarco
"""

#Python imports
import logging
import pydevd
import argparse
import time
import cPickle as Pickle
from termcolor import colored


#Self imports
from Pythway_utils import *
from Pythway_export_tools import *


def menu():
    """Input menu for user."""

    parser = argparse.ArgumentParser(description='Pythways finder v0.1',
                                     epilog="###############################")

    required = parser.add_argument_group('Required options')

    required.add_argument('-i', dest='input_ensembl_genes', action='store',
                          help='Input ENSEMBL gene ID list. One gene per line. ie: -i /home/my_genes.txt',
                          required=True)

    optional = parser.add_argument_group('Non-mandatory options')

    optional.add_argument('-e', dest='pathway_export_format', action='store',
                          help='Format you want to export the result Pathways: gpml, png, svg. ie: -e gpml',
                          required=False, default='gpml')

    optional.add_argument('-hgnc', dest='hgnc', action='store_true',
                          help='Use HGNC Symbol as input and not Ensembl Gene ID.',
                          required=False)

    optional.add_argument('-o', dest='output_file', action='store', default='./pythway_results.txt',
                          help="Output results file. If path to file doesn't exists it will be automatically created. "
                               "ie: -output /home/user/pathway/results.txt",
                          required=False)

    optional.add_argument('-d', dest='remote_debug', action='store_true',
                          help="Remote debug action.", required=False)

    return parser.parse_args()


def check_pathway_format(option):
    if option not in ['gpml', 'png', 'svg']:
        raise ValueError('The Pathway export format {0} is not supported.'.format(option))


def file_exists(input_file):
    """Checks if input_file exists"""

    if not os.path.exists(input_file):
        raise IOError("The file %s does not exist!" % input_file)


def ensembl_gene_input_filter(gene_input_file):
    #Load gene list into input_gene_list and remove empty lines and also \n\r
    gene_input_file = filter(None, [line.rstrip() for line in open(gene_input_file, 'r')])

    #Remove duplicates in list
    gene_input_file = set(gene_input_file)

    #Remove '-' values for lanes with no gene from Ensembl annotation.
    gene_input_file = filter(lambda a: a != '-', gene_input_file)

    #Remove rest of items not starting with 'EN' string.
    gene_input_file = filter(lambda a: a.startswith('EN'), gene_input_file)

    return gene_input_file


def hgnc_gene_input_filter(gene_input_file):

    #Load gene list into input_gene_list and remove empty lines and also \n\r
    gene_input_list = filter(None, [line.rstrip() for line in open(gene_input_file, 'r')])

    #Remove duplicates in list
    gene_input_list = set(gene_input_list)

    #Remove '-' values for lanes with no gene from Ensembl annotation.
    gene_input_list = filter(lambda a: a != '-', gene_input_list)

    #Convert HGNC to Ensembl Gene ID.
    ens_gene_list = hgnc_to_ensembl_id(gene_input_list)

    #Remove duplicates in list
    ens_gene_list = set(ens_gene_list)

    return ens_gene_list


def candidate_genes_input_filter(candidate_genes_file):
    #Load gene file list into input_gene_list and remove empty lines and also \n\r
    candidate_genes_file = filter(None, [line.rstrip() for line in open(candidate_genes_file, 'r')])

    #Remove header starting with '#'
    candidate_genes_file = [x for x in candidate_genes_file if not x.startswith('#')]

    #Split each element of list by '|' store first value
    candidate_genes_file = [i.split('|', 1)[0] for i in candidate_genes_file]

    #Remove duplicates in list
    candidate_genes_file = set(candidate_genes_file)

    return candidate_genes_file


def output_mgmt(output_file_path):
    """Creates output directories, even if nested in
    case that they don't exist."""

    if output_file_path == '.':
        raise ValueError(
            'The path {0} is not a valid option. Leave output option blank or specify a valid output file name.'.format(
                output_file_path))

    base_path = os.path.dirname(output_file_path)

    if base_path == '':
        base_path = '.'

    if not os.path.exists(base_path):
        print "Creating output directory..."
        os.makedirs(base_path)

    base_path = os.path.abspath(base_path)
    return base_path


def working_message(arg_list, gene_count):
    """Returns message with the options that the
    user has chosen"""

    message = 'Pythways v0.1'
    message += '\n=============================='
    message += '\nInput file: %s' % arg_list.input_ensembl_genes
    message += '\nProcessing %i genes.' % gene_count
    message += '\nPathway export format: %s' % arg_list.pathway_export_format
    message += '\nOutput: %s' % arg_list.output_file
    message += '\n=============================='

    return message


if __name__ == '__main__':

    #Enable logging file for exceptions.
    logging.basicConfig(format='%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.ERROR,
                        filename='genes_error.log')
    with open('genes_error.log', 'w'):
        pass

    #Parse user menu options.
    args = menu()

    #Check for remote debug option
    if args.remote_debug:
        user_ip = os.environ['USERIP']
        pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)

    related_pathway_list = []
    related_ensembl_pathways = []

    #Start timer counter.
    start_time = time.time()

    #Get current date
    now = datetime.datetime.now()
    date = now.strftime("%d_%m_%Y")

    #Check Pathway export format is supported
    pathway_export_format = args.pathway_export_format.lower()
    check_pathway_format(pathway_export_format)

    #Check input file exists
    file_exists(args.input_ensembl_genes)

    #Filter gene input list, remove duplicates and suppress '-'
    if args.hgnc:
        ensembl_filtered_input_gene_list = hgnc_gene_input_filter(args.input_ensembl_genes)
    else:
        ensembl_filtered_input_gene_list = ensembl_gene_input_filter(args.input_ensembl_genes)

    #Print working message.
    print working_message(args, len(ensembl_filtered_input_gene_list))

    #Create output dirs if needed.
    output_base_path = output_mgmt(args.output_file)

    #Create output dir for pathway export
    pathway_export_path = createdir(output_base_path + '/pathway_export')

    #Create output dir for pathway html
    pathway_html_path = createdir(output_base_path + '/html_report')

    #Create output dir for ontology html
    ontology_html_path = createdir(pathway_html_path + '/ontology')

    #Create Pickle dir
    # pickle_output_path = createdir(output_base_path + '/pickle')
    
    # #Pickle output
    # pickle_output = pickle_output_path + '/' + date + '_related_ensembl_pathways.p'

    #Open output file.
    output = open(args.output_file, 'w')

    #Print output file header
    output.write('#PW_ID\tPW_Name\tRevision\tRelated_HGNC_list\n')

    #Script path
    script_path = os.path.dirname(__file__)

    print 'Script is working. This may take a while since Wikipathways WS are slow as hell !'

    for ensembl_gene_id in ensembl_filtered_input_gene_list:
        #gene_iteration_counter += 1
        try:
            pathway_list = create_pathway_from_ensembl_gene_id(ensembl_gene_id)
        except AttributeError:
            logging.error("Wrong return type from %s" % ensembl_gene_id)

        #Skip to next pathway_list if no pathway_list is returned.
        try:
            if not pathway_list:
                continue

        except NameError:
            continue

        for pathway in pathway_list:
            #We remove duplicates and add ensembl gene list for each pathway.
            if pathway not in related_pathway_list:
                related_pathway_list.append(pathway)
                ensembl_id_list = get_pathway_xref_list(pathway, 'En')
                pathway.set_ensembl_list(ensembl_id_list)

                #Get common genes between input gene list and pathway
                common_gene_list = common_genes(ensembl_filtered_input_gene_list, ensembl_id_list)
                pathway.set_common_ensembl_gene_list(common_gene_list)

                #Translate common Ensembl Gene IDs to HGNC
                hgnc_common_list = ensembl_to_hgnc(common_gene_list)
                pathway.set_common_hgnc_list(hgnc_common_list)

    for pathway in related_pathway_list:
        #We check if input genes are contained in any of related obtained pathways found in related_pathway_list
        ensembl_gene_list = pathway.get_ensembl_list()

        if is_contained(ensembl_filtered_input_gene_list, ensembl_gene_list):
            #Get Ontology terms for pathway, create a list of Ontology objects and set it to it's corresponding Pathway
            ontology_list = create_ontology_list_by_pathway(pathway)
            pathway.set_ontology_list(ontology_list)
            related_ensembl_pathways.append(pathway)
            pathway_print_string = "\t".join([pathway.get_pw_id(), pathway.get_name(), str(pathway.get_revision()),
                                              ','.join(pathway.get_common_hgnc_list())])
            output.write(pathway_print_string + '\n')

    for pathway in related_ensembl_pathways:

        res = get_pathway_as(pathway, pathway_export_format, pathway_export_path)
        ontology_list = pathway.get_ontology_list()

        if len(ontology_list) > 0:
            ontology_pathway_html_report(pathway, ontology_list, ontology_html_path, script_path)


    #Save Pathway objects to pickle
    # print 'Generated pickle object file: {0}'.format(pickle_output)
    # Pickle.dump(related_ensembl_pathways, open(pickle_output, "wb"))

    gpml_file_list = list_gpml_files(pathway_export_path, pathway_export_format)
    pathvisio_rpc_export(pathway_export_path, gpml_file_list, pathway_html_path)
    removedir(pathway_export_path)
    jinja_html_report(related_ensembl_pathways, pathway_html_path, script_path)

    time = int(time.time() - start_time)
    print colored('DONE. Time elapsed: %s', 'green', attrs=['bold']) % str(datetime.timedelta(seconds=time))