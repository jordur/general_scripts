#!/usr/bin/env python
"""
SYNOPSIS
    normalize_variants_vcf.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    normalize_variants_vcf: Script that parses variant-calling files from different deployed variant-callers, and
    generates or normalizes the FREQ information for each sample
EXAMPLES
    normalize_variants_vcf.py -i variants.vcf -o normalized_variants.vcf
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    normalize_variants_vcf.py v0.1
"""
import pydevd
import sys
import os
import traceback
import optparse
import time
import vcf
from vcf.parser import _Format as VcfFormat
from vcf.parser import _Info as VcfInfo
from variant import norm_freq

# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args

    # Open and parse each line of the vcf file
    input_vcf = vcf.Reader(open(options.input_vcf, 'r'))
    # If an FREQ field already exists in FORMAT or INFO, it has to be stored and be used when importing from input
    former_vcfformat_freq = input_vcf.formats['FREQ'] if 'FREQ' in input_vcf.formats else None
    former_vcfinfo_sfreq = input_vcf.infos['FREQ'] if 'FREQ' in input_vcf.infos else None
    former_vcfinfo_sdp = input_vcf.infos['DPS'] if 'DPS' in input_vcf.infos else None
    input_vcf.formats['FREQ'] = VcfFormat('FREQ', None, 'String', 'Variant allele frequency')
    input_vcf.infos['SFREQ'] = VcfInfo('SFREQ', 1, 'Float', 'Maximum variant allele frequency of all samples')
    input_vcf.infos['SDP'] = VcfInfo('SDP', 1, 'Integer', 'Maximum sequencing depth of all samples')
    output_vcf = vcf.Writer(open(options.output_vcf, 'w'), input_vcf, lineterminator='\n')
    if former_vcfformat_freq is not None:
        input_vcf.formats['FREQ'] = former_vcfformat_freq
    if former_vcfinfo_sfreq is not None:
        input_vcf.infos['SFREQ'] = former_vcfinfo_sfreq
    if former_vcfinfo_sdp is not None:
        input_vcf.infos['SDP'] = former_vcfinfo_sdp

    for record in input_vcf:
        if not 'FREQ' in record.FORMAT.split(':'):
            record.add_format('FREQ')

        # Default values for added INFO fields
        site_freq = None
        site_depth = 0

        # iterate over all call objects of record
        for call in record.samples:
            # Allele frequency and Depth evaluation among samples
            try:
                site_freq = max(site_freq, max(call.aaf)) if call.aaf is not None else site_freq
                site_depth = max(call.depth, site_depth) if call.depth is not None else site_depth
            except Exception:
                print "ERROR: unforeseen exception when normalizing record:", record
                raise
            call.add_format('FREQ', norm_freq(call.aaf))
            # TODO: unfortunately GATK filtering doesn't yet deal correctly with "None" (.) values
            if site_freq is None or site_freq == '.':
                site_freq = 0
            record.add_info('SFREQ', site_freq)
            record.add_info('SDP', site_depth)
        output_vcf.write_record(record)


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
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