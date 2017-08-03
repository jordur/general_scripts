#!/usr/bin/env python
"""
SYNOPSIS
    repairing_vcf.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    repairing_vcf.py: Script that parses VarScan variant-calling files searching multiple variants alleles, since
    VarScan current version (v2.3.6) seems to have trouble with these sites
EXAMPLES
    repairing_vcf.py -i variants_varscan.vcf -o variants_varscan_repaired.vcf
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    repairing_vcf.py v0.1
"""
import pydevd
import sys
import os
import traceback
import optparse
import time
import vcf


# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args

    # Open and parse each line of the vcf file
    input_vcf = vcf.Reader(open(options.input_vcf, 'r'))

    # Create output folder if not existent
    if not os.path.exists(options.output_folder):
        os.makedirs(options.output_folder)

    # Open output files for writing
    output_files = {}
    for sample in input_vcf.samples:
        output_files[sample] = open(os.path.join(options.output_folder, sample), 'w')

    for record in input_vcf:
        for sample in record.samples:
            if sample.aaf is not None:
                for aaf in sample.aaf:
                    output_files[sample.sample].write(str(aaf) + '\n')
            else:
                output_files[sample.sample].write('.\n')


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-i', '--input_vcf', action='store', type='string', dest="input_vcf", default=None,
                          help='Path to input vcf variant calling file')
        parser.add_option('-o', '--output_folder', action='store', type='string', dest='output_folder',
                          default='./sample_freqs', help='Output folder')

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