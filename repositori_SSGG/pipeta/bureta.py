#!/usr/bin/env python
"""
SYNOPSIS
    TODO bureta [-h,--help] [--version]
DESCRIPTION
    This script replaces matraz.pl perl script. They will coexist
    until bureta is tested and works bug free.

    This script generates JSON config files for PipeTA.
    This docstring will be printed by the script if there
    is an error or if the user requests help (-h or --help).
EXAMPLES
    TODO: Show some examples of how to use this script.
EXIT STATUS
    TODO: List exit codes
AUTHOR
    TODO: Guillermo Marco <guillermo.marco@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    bureta.py v0.0.1
"""

import sys
import os
import re
import traceback
import argparse
import time
import pydevd
import json

matraz_config_file = os.environ['MATRAZ_CONFIG']
gluster = os.environ['PROJECTS']
gluster_2 = os.environ['PROJECTS2']


def convert_unicode_to_str(json_input):
    if isinstance(json_input, dict):
        return {convert_unicode_to_str(key): convert_unicode_to_str(value) for key, value in json_input.iteritems()}
    elif isinstance(json_input, list):
        return [convert_unicode_to_str(element) for element in json_input]
    elif isinstance(json_input, unicode):
        return json_input.encode('utf-8')
    else:
        return json_input


def config_reader(json_config, option):
    #Load JSON data, unicode
    json_data = json.load(open(json_config))

    #Conver unicode data to string
    json_data = convert_unicode_to_str(json_data)

    three_prime_primers = five_prime_primers = ''

    try:
        target = json_data['target_reference'][option]['target']
        capture = json_data['target_reference'][option]['capture']
        target_chrs = json_data['target_reference'][option]['target_chrs']
        chr_split = json_data['target_reference'][option]['chr_split']
        prefix = json_data['target_reference'][option]['prefix']
        capture_system = json_data['target_reference'][option]['capture_system']
        pipeline = json_data['target_reference'][option]['pipeline']
        aligner = json_data['target_reference'][option]['aligner']

        #Check if three and five primmers keys are defined in JSON config
        if 'three_prime_primers' in json_data['target_reference'][option] and 'five_prime_primers' in \
                json_data['target_reference'][option]:
            three_prime_primers = json_data['target_reference'][option]['three_prime_primers']
            five_prime_primers = json_data['target_reference'][option]['five_prime_primers']

    except KeyError:
        print 'ERROR: {0} not found in Bureta analysis database config file.'.format(option)

    print 'INFO: Analysis configuration loaded successfully.'
    return target, capture, target_chrs, chr_split, prefix, capture_system, pipeline, aligner, \
        three_prime_primers, five_prime_primers


def project_exists(project_code):
    gluster_project_path = os.path.join(gluster, project_code, 'rawdata')
    gluster_2_project_path = os.path.join(gluster_2, project_code, 'rawdata')

    if os.path.isdir(gluster_project_path):
        return gluster_project_path

    elif os.path.isdir(gluster_2_project_path):
        return gluster_2_project_path

    else:
        print 'ERROR: {0} project not found in Gluster/Gluster2.'.format(project_code)
        sys.exit(1)


def get_samples(rawdata_path):
    for sample_dir in os.listdir(rawdata_path):
        sample_dir_prefix = re.search('Sample_(.*)', sample_dir).group(1)
        print '\n', sample_dir
        print sample_dir_prefix

        subdir_path = os.path.join(rawdata_path, sample_dir)

        for sample_file in os.listdir(subdir_path):
            if sample_file.endswith('.fastq.gz' or 'fastq'):
                lane_number = re.search('_(L[\d]*)_', sample_file).group(1)
                read = re.search('_(R\d)_', sample_file).group(1)
                segment = re.search('([\d]*)\.fastq[\.gz]??', sample_file).group(1)

                print lane_number, read, segment

                sample_file_path = os.path.join(subdir_path, sample_file)
                print sample_file_path


def print_json_header():
    return 1
    # $$json{"name"} = "$bf";
	# $$json{"path"} = "$dir/$bf";
	# $$json{"references"} = {};
	# $$json{"references"}{"genomes"} = {};
	# $$json{"references"}{"genomes"}{"human"} = {};
	# $$json{"references"}{"genomes"}{"human"}{"fasta"} = "/share/references/genomes/human/hg19/reference/human_hg19.fa";
	# $$json{"references"}{"genomes"}{"human"}{"bwa_fasta"} = "/share/references/genomes/human/hg19/reference/bwa/hg19.fa";
	# $$json{"references"}{"genomes"}{"human"}{"bowtie2_reference"} = "/share/references/genomes/human/hg19/reference/bowtie2/hg19";
	# $$json{"references"}{"genomes"}{"human"}{"novoalign_index"} = "/share/references/genomes/human/hg19/reference/novoindex/hg19.ndx";
	# $$json{"references"}{"genomes"}{"human"}{"vcf"} = "/share/references/realign_recalibrate_hapmap/common_all.vcf";



def main():
    global args
    # TODO: Do something more interesting here...
    print 'Main'

    #Check if BF directory project exists in Gluster/Gluster2
    project_rawdata_path = project_exists(args.bf_project_code)

    #Load configuration from JSON config
    target, capture, target_chrs, chr_split, prefix, capture_system, pipeline, aligner, \
    three_prime_primers, five_prime_primers = config_reader(matraz_config_file, args.analysis_option)

    get_samples(project_rawdata_path)

    print 'DONE'


if __name__ == '__main__':
    try:
        start_time = time.time()

        parser = argparse.ArgumentParser(description='bureta.py', version='alpha version',
                                         epilog='This script replaces matraz.pl perl script.\
                                         They will coexist until bureta is tested and works bug free.\
                                         This script generates JSON config files for PipeTA.')

        parser.add_argument('-bf', dest='bf_project_code', action='store',
                            help='Bioinformatics project code',
                            required=True)

        parser.add_argument('-a', dest='analysis_option', action='store',
                            help="Analysis type to create JSON config.",
                            required=True)

        parser.add_argument('-d', dest='debug', action='store_true',
                            help="Remote debug action.", required=False)

        args = parser.parse_args()

        #parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
        #                               version='$Id$')
        #parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        #parser.add_option('-bf', '--bionformatics', action='store_true', default=False, help='BF project')
        #parser.add_option('-a', '--analysis', action='store_true', default=False, help='analysis option')
        #parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        #(options, args) = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
        if args.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        main()

    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        sys.exit(1)