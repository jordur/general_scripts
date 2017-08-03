#!/usr/bin/env python
"""
SYNOPSIS
    freq_plots.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    freq_plots: Script that plots an histogram of aaf (alternative allelic frequencies)
EXAMPLES
    freq_plots.py -i freqs.tsv -o plot.jpg
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    freq_plots.py v0.1
"""
import pydevd
import csv
import sys
import os
import traceback
import optparse
import time
import matplotlib.pyplot as plt
#import numpy as np


# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args

    # Open input file
    with open(options.aaf_file, 'rb') as tsvfile:
        dialect = csv.Sniffer().sniff(tsvfile.read(1024))
        tsvfile.seek(0)
        sample_freqs = []
        sample_names = []
        for row in csv.reader(tsvfile, dialect):
            if sample_names == []:
                for sample in row:
                    sample_names.append(sample)
            else:
                for sample in row:
                    sample_freqs.append(sample)

        n, bins, patches = plt.hist(sample_freqs[:], 50, normed=1, facecolor='g', alpha=0.75)
        plt.show()


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-i', '--aaf_file', action='store', type='string', dest='aaf_file', default=None,
                          help='Path to .tsv file containing samples alternative allele frequencies')
        parser.add_option('-o', '--output_file', action='store', type='string', dest='output_file',
                          default='output.jpg', help='Output name')

        (options, args) = parser.parse_args()
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if not options.aaf_file:
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