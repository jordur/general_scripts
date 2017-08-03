#!/usr/bin/env python
"""
SYNOPSIS
    test_hgvs.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    test_hgvs: Script that parses variant-calling annotated files and generates the corresponding GeneSys files
EXAMPLES
    test_hgvs.py -i variants_annotated.vcf -o variants_GeneSys.psv
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    test_hgvs.py v0.1
"""
import pydevd
import sys
import os
import traceback
import optparse
import time
import vcf
from variant import Variant
import hgvs
import hgvs.utils
from pygr.seqdb import SequenceFileDB


# Functions definitions
# Read genome sequence using pygr.
genome = SequenceFileDB('/share/references/genomes/human/hg19/reference/human_hg19.fa')

# Read RefSeq transcripts into a python dict.
with open('/share/references/genes/human/Homo_sapiens.GRCh37.73.modified_genePred') as infile:
    transcripts = hgvs.utils.read_transcripts(infile)


# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    # Following code line is for compatibility with Ensembl GenePred files:
    # Ensembl transcripts don't include transcript version in .gtf files (gene set files). To assure version
    # compatibility, Ensembl annotation and gtf files must have the same version!!
    name = name.split('.')[0] if name[:4] == 'ENST' else name
    return transcripts.get(name)


# ****************************** main body *********************************
def main():
    global options, args
    separator = '|'

    # Parse the HGVS name into genomic coordinates and alleles.
    #chrom, offset, ref, alt = hgvs.parse_hgvs_name('ENST00000515609.1:c.30G>T', genome, get_transcript=get_transcript)
    #print chrom, offset, ref, alt

    # Format an HGVS name.
    chrom, offset, ref, alt = ('chr2', 179616770, 'GAA', 'G')
    transcript = get_transcript('ENST00000359218.5')
    hgvs_name = hgvs.format_hgvs_name(chrom, offset, ref, alt, genome, transcript)
    print hgvs_name
    chrom, offset, ref, alt = ('chr2', 179616770, 'GAA', 'GA')
    transcript = get_transcript('ENST00000359218.5')
    hgvs_name = hgvs.format_hgvs_name(chrom, offset, ref, alt, genome, transcript)
    hgvs_var = hgvs.HGVSName(hgvs_name)
    hgvs_str = 'ENST00000359218.5:c.10597+1079_10597+1080delTTinsT'
    hgvs_var2 = hgvs.HGVSName(hgvs_str)

    print hgvs_name
    quit()

    # Open and parse each line of the vcf file
    input_vcf = vcf.Reader(open(options.input_vcf, 'r'))
    variant = Variant(samples=input_vcf.samples)

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