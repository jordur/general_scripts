#!/usr/bin/env python
"""
SYNOPSIS
    vcf2DBNLVar.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    vcf2DBNLVar: Script that parses variant-calling annotated files and generates the corresponding GeneSys files
EXAMPLES
    vcf2DBNLVar.py -i variants_annotated.vcf -o variants_GeneSys.psv
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    vcf2DBNLVar.py v0.1
"""
import json
import httplib2
import pydevd
import sys
import os
import traceback
import optparse
import time
import requests
import vcf
from genotype import chromosome_number2str, get_chromosome_number
from variant import Variant


# Global parameters
#uri = 'http://genesys-w2p-test.sistemasgenomicos.com/dbnlvar/api/rest/'
uri = 'http://genesys-w2p-test.sistemasgenomicos.com/dbnlvar/api/'
auth = ('oscar.rodriguez@sistemasgenomicos.com', 'prueba')
separator = '|'


# Functions definitions
def read_table(connection, uri, data, method="GET"):
    resp, content = connection.request(uri + "chromosome/3.json", "GET", body="",
                                       headers={'content-type': 'application/json'})
    print content


def kaviar_test():
    uri = 'http://db.systemsbiology.net/kaviar/cgi-pub/Kaviar.pl'

    # using separate chromosome and coordinate fields
    # get variant information for the following positions on chromosome 1  for hg18
    freeze = "hg18"
    coordinates = "1"  # 0 for 0-based, 1 for 1-based
    chr = "chr1"
    positions = "77119, 555911, 714019"

    format = "json"

    args = {'frz': freeze, 'onebased': coordinates, 'chr': chr, 'pos': positions, 'format': format}
    header = requests.get(uri)
    print header.headers
    print header.request.headers()

    resp = requests.get(uri, args=args)
    print resp


def check_record(table=None, value=None):
    resp = requests.get(uri + table + '/' + str(value) + '.json', auth=auth)
    return True if resp.status_code == 200 else False


def create_record(payload=None):
    headers = {'Content-Type': 'application/json'}
    resp = requests.post(uri + 'rest/INSERT/chromosome/X.json', headers=headers, auth=auth, data=json.dumps(payload))
    return True if resp.status_code == 200 else False


def load_record(payload=None, template=None):
    headers = {'Content-Type': 'application/json'}
    record_json = {}
    record_json['ALT'] = []
    for substitution in payload.ALT:
        record_json['ALT'].append(str(substitution))
    record_json['CHROM'] = payload.CHROM
    record_json['FILTER'] = payload.FILTER
    record_json['FORMAT'] = payload.FORMAT
    record_json['INFO'] = payload.INFO
    record_json['POS'] = payload.POS
    record_json['QUAL'] = payload.QUAL
    record_json['REF'] = payload.REF
    #record_json[''] = payload.

    template.write_record(payload)
    resp = requests.post(uri + 'custom/LOAD', headers=headers, auth=auth, data=json.dumps(record_json))
    print json.dumps(record_json)
    return True if resp.status_code == requests.codes.ok else False


def load_consequence(consequence=None):
    headers = {'Content-Type': 'application/json'}
    #resp = requests.post(uri + 'custom/LOAD', headers=headers, auth=auth, data=json.dumps(payload))
    # Data definitions
    chrom = get_chromosome_number(consequence['Chr']['value'])
    chrom_name = chromosome_number2str(chrom)
    hgnc = consequence['SYMBOL']['value']
    ensembl_id = consequence['Gene']['value']
    transcript_ensembl_id = consequence['Feature_ID']['value']
    hgmd_name = consequence['RefSeq_ID']['value']

    # chromosome
    payload = {"id": chrom, "name": chrom_name}
    resp = requests.post(uri + 'rest/INSERT/chromosome/X.json', headers=headers, auth=auth, data=json.dumps(payload))
    # gene
    payload = {"chromosome_id": chrom, "HGNC": hgnc, "ensembl_id": ensembl_id}
    resp = requests.post(uri + 'rest/INSERT/gene/X.json', headers=headers, auth=auth, data=json.dumps(payload))
    gene_id = resp.id
    # transcript
    payload = {"genes_id": gene_id, "HGMD_name": hgmd_name, "ensembl_id": transcript_ensembl_id}
    resp = requests.post(uri + 'rest/INSERT/transcript/X.json', headers=headers, auth=auth, data=json.dumps(payload))
    return True if resp.status_code == 200 else False


# ****************************** main body *********************************
def main():
    global options, args
    # **********************
    # store in DBNLVar
    # **********************

    # Define connection parameters andd perform connection:
    connection = httplib2.Http(".cache")
    connection.add_credentials('oscar.rodriguez@sistemasgenomicos.com', 'prueba')

    # Open annotation file and parse each line in it
    annotation_vcf = vcf.Reader(open(options.input_vcf, 'r'))

    # Load metadata in variant object
    variant = Variant(samples=annotation_vcf.samples)

    for record in annotation_vcf:
        # Load variant information in DBNLVar, from consequences
        variant.get_from_record(record=record)
        for consequence in variant.consequences:
            resp = load_consequence(consequence=consequence)

        quit()
        # Store consequence non-relating data in DBNLVar
        if not check_record(table='chromosome', value=chrom_to_number(record.CHROM)):
            #payload = {'id': chrom_to_number(record.CHROM), 'name': number_to_chrom(chrom_to_number(record.CHROM))}
            load_record(payload=record)

        # Store consequence relating data in DBNLVar
        for consequence in record.INFO['CSQ']:
            for index, annotation in enumerate(consequence.split(separator)):
                payload = {}

        resp = requests.get(uri + 'chromosome/id/3/24.json', auth=auth)
        print resp.json()
        if resp.status_code == 200:
            content = resp.json()['content']
        else:
            print "ERROR: Problem in query"
            raise
        print content
        # Read variant rows


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-i', '--input_vcf', action='store', type='string', dest='input_vcf', default=None,
                          help='Path to input vcf variant calling file')
        #parser.add_option('-o', '--output_vcf', action='store', type='string', dest='output_vcf',
        #                  default='output.vcf', help='Output name')

        (options, args) = parser.parse_args()
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if not options.input_vcf:
            parser.error('ERROR: missing arguments!! Please check usage')
        if options.verbose:
            print time.asctime()

        # MAIN and global
        sender = 'genetonic@sistemasgenomicos.com'
        receivers = [os.environ['USERMAIL']]
        main()

        if options.verbose:
            print time.asctime()
        if options.verbose:
            print 'INFO: Duration of script run:',
        if options.verbose:
            print str(time.time() - start_time) + ' seconds'
        sys.exit(0)
    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        sys.exit(1)