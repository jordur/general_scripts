#!/usr/bin/env python
"""
SYNOPSIS
    vcf2GeneSys.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    vcf2GeneSys: Script that parses variant-calling annotated files and generates the corresponding GeneSys files
EXAMPLES
    vcf2GeneSys.py -i variants_annotated.vcf -o variants_GeneSys.psv
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    vcf2GeneSys.py v0.1
"""
import pydevd
import sys
import os
import traceback
import optparse
import time
import vcf
from variant import Variant


# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args
    separator = '|'

    # Open and parse each line of the vcf file
    input_vcf = vcf.Reader(open(options.input_vcf, 'r'))
    if options.non_model:
        variant = Variant(samples=input_vcf.samples, organism_type='non_model', ploidy=options.ploidy)
    else:
        variant = Variant(samples=input_vcf.samples, ploidy=options.ploidy)

    # Open output file
    with open(options.output_vcf, 'w') as output_psv:
        # Generate output file header
        #variant = ConsequenceType(input_vcf.samples)
        output_psv.write(variant.create_psv_header(separator=separator))

        # Now parse lines in .vcf and output with new format:
        for record in input_vcf:
            # Only output sites that hasn't been filtered out
            if len(record.FILTER) == 0:
                #for consequence in range(0, len(record.INFO['CSQ'])):
                variant.get_from_record(record=record)
                output_psv.write(variant.put_to_psv(separator=separator))


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-n', '--non-model', action='store_true', default=False, help='use in case of non-model '
                                                                                        'organisms')
        parser.add_option('-p', '--ploidy', action='store', type='int', dest='ploidy', default=2,
                          help='Organism ploidy. Enter number of copies of chromosome sets. Default is diploid (2)')
        parser.add_option('-i', '--input_vcf', action='store', type='string', dest='input_vcf', default=None,
                          help='Path to input vcf variant calling file')
        parser.add_option('-o', '--output_vcf', action='store', type='string', dest='output_vcf',
                          default='output.vcf', help='Output name')

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