#!/usr/bin/env python
"""
SYNOPSIS
    collecting_vcf_variants_vcf.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    collecting_vcf_variants_vcf: Script that merges variant-calling files from different deployed variant-callers.
    Unless otherwise stated, the first given variant-calling file will be the template file for outputting information.
    IMPORTANT: It's assumed that all .vcf input files are position-sorted!!
EXAMPLES
    collecting_vcf_variants_vcf.py -i variants1.vcf -i variants2.vcf -o merged.vcf
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    collecting_vcf_variants_vcf.py v0.1
"""
import shlex
import subprocess
import pydevd
import sys
import os
import traceback
import optparse
import time
import vcf
from vcf.parser import _Format as VcfFormat
from vcf.parser import _Info as VcfInfo
from genotype import get_chromosome_number, chromosome_number2str

# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args

    # Be sure to get files bgzipped and tabix indexed
    for vcf_file in options.input_vcf:
        #if not os.path.isfile(vcf_file + '.gz'):
        command_line = "bgzip -c " + vcf_file + " > " + vcf_file + ".gz"
        shlex.split(command_line)
        retcode = subprocess.check_output(command_line, stderr=subprocess.STDOUT, shell=True)
        print retcode
        #if not os.path.isfile(vcf_file + '.gz.tbi'):
        command_line = "tabix -f -p vcf " + vcf_file + ".gz"
        shlex.split(command_line)
        retcode = subprocess.check_output(command_line, stderr=subprocess.STDOUT, shell=True)
        print retcode

    # First vcf file will be the template file for outputting parameters
    template_vcf = vcf.Reader(open(options.input_vcf[0], 'r'))

    # Add essential fields in both formats and infos (header information)
    template_vcf.formats['FREQ'] = VcfFormat('FREQ', 1, 'String', 'Variant allele frequency')
    template_vcf.infos['SFREQ'] = VcfInfo('SFREQ', 1, 'Float', 'Maximum variant allele frequency of all samples')
    template_vcf.infos['SDP'] = VcfInfo('SDP', 1, 'Integer', 'Maximum sequencing depth of all samples')

    # Create a list of sorted variant-sites containing chr and position
    variant_sites = []
    for vcf_file in options.input_vcf:
        tmp_vcf = vcf.Reader(open(vcf_file, 'r'))
        for record in tmp_vcf:
            new_variant_site = (get_chromosome_number(record.CHROM), record.POS)
            if not new_variant_site in variant_sites:
                variant_sites.append(new_variant_site)
    variant_sites.sort(key=lambda variant: (variant[0], variant[1]))

    # Open all files for random access
    input_vcf = []
    for index, vcf_file in enumerate(options.input_vcf):
        input_vcf.append(vcf.Reader(open(vcf_file + '.gz', 'r')))
        # Perform tests and checks
        if index > 0 and input_vcf[index].samples != template_vcf.samples:
            print "INFO: not same sample list in", vcf_file

        # Add necessary FORMAT or INFO fields definitions in template
        for info in input_vcf[index].infos:
            if not info in template_vcf.infos:
                template_vcf.infos[info] = input_vcf[index].infos[info]
        for myformat in input_vcf[index].formats:
            if not myformat in template_vcf.formats:
                template_vcf.formats[myformat] = input_vcf[index].formats[myformat]

    # Open output handles
    output_vcf = vcf.Writer(open(options.output_vcf, 'w'), template_vcf, lineterminator='\n')
    output_indels_vcf = vcf.Writer(open(options.output_vcf+'_indels.vcf', 'w'), template_vcf, lineterminator='\n')
    output_snps_vcf = vcf.Writer(open(options.output_vcf+'_snps.vcf', 'w'), template_vcf, lineterminator='\n')

    # Now parse each variant-site and fetch information from vcfs:
    for my_variant_site in variant_sites:
        records = []
        for my_vcf in input_vcf:
            try:
                for record in my_vcf.fetch(chromosome_number2str(my_variant_site[0]), my_variant_site[1],
                                           my_variant_site[1]):
                    # vcf.fetch returns also next position if described, must be therefore removed
                    if record.POS == my_variant_site[1]:
                        records.append(record)
            except KeyError:
                # This exception is raised when the primary key is not found in one of the files. No actions required
                pass
        # master_records are those records for being output to merged vcf. A master record will be created for each
        # group of variants from a same variant site that can be merged
        master_records = [records[0]]
        for record in records[1:]:
            add_to_master = False
            already_added = False
            for master_record in master_records:
                if master_record != record:
                    add_to_master = True
                else:
                    if not master_record.merge(record):
                        add_to_master = True
                    else:
                        already_added = True
            if add_to_master and not already_added:
                master_records.append(record)

        for master_record in master_records:
            output_vcf.write_record(master_record)
            if master_record.is_snp:
                output_snps_vcf.write_record(master_record)
            elif master_record.is_indel:
                output_indels_vcf.write_record(master_record)


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-i', '--input_vcf', action='append', type='string', dest='input_vcf', default=None,
                          help='Path to input vcf variant calling files (it can be repeated as many times as desired)')
        parser.add_option('-o', '--output_vcf', action='store', type='string', dest='output_vcf',
                          default='merged.vcf', help='Output name')

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