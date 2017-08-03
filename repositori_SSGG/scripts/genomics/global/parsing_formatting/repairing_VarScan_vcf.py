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
from vcf.parser import _Format as VcfFormat
from vcf.model import _Record as Record
from vcf.model import _Substitution as Substitution
from vcf.model import _Call as Call
from variant import norm_freq


# Functions definitions


# ****************************** main body *********************************
def main():
    global options, args

    # Create template for outputting from input file
    template_vcf = vcf.Reader(open(options.input_vcf, 'r'))
    template_vcf.formats['FREQ'] = VcfFormat('FREQ', None, 'String', 'Variant allele frequency')

    # Open and parse each line of the vcf file
    input_vcf = vcf.Reader(open(options.input_vcf, 'r'))
    repaired_vcf = vcf.Writer(open(options.output_vcf, 'w'), template_vcf, lineterminator='\n')
    for record in input_vcf:
        if len(record.ALT) > 1:
            # Problem 1 in VarScan files: alternative allele repeated. Output unique alternative alleles
            new_alt = []
            for index, alternative in enumerate(record.ALT):
                if not alternative in new_alt:
                    new_alt.append(alternative)
            for call in record.samples:
                if call.gt_nums is not None:
                    gt_alleles = []
                    for allele in call.gt_bases.split(call.gt_phase_char()):
                        try:
                            gt_alleles.append(str(new_alt.index(allele)+1))
                        except ValueError:
                            gt_alleles.append('0')
                    call.add_format('GT', call.gt_phase_char().join(gt_alleles))
                    call.gt_nums = call.gt_phase_char().join(gt_alleles)
            record.ALT = new_alt
            record.alleles = [record.REF] + new_alt

            # Problem 2 in VarScan files: only a FREQ value is stored, even in case of multiple alternative alleles
            for call in record.samples:
                # Assign allele frequency to those alleles predicted by variant caller
                freqs = []
                for index in range(1, len(record.ALT) + 1):
                    if call.gt_nums == '0/0' or call.gt_nums == '0|0':
                        freqs.append(norm_freq(call.aaf[0]))
                    elif call.gt_alleles is None:
                        freqs.append(None)
                    elif str(index) in call.gt_alleles and ('0' in call.gt_alleles or call.gt_type == 2):
                        freqs.append(norm_freq(call.aaf[0]))
                    elif str(index) not in call.gt_alleles:
                        freqs.append(None)
                    else:
                        print "ERROR: unforeseen case when obtaining freqs!!"
                        print "record:", record
                        raise
                call.add_format('FREQ', freqs)
        else:
            for call in record.samples:
                call.add_format('FREQ', norm_freq(call.aaf))
        # Problem 3: hybrid records (both SNP and Indel in a single record)
        snv_record = None
        indel_record = None
        if record.is_indel and len(record.ALT) > 1:
            # Creation of lists containing indels and snvs
            indel_alternatives = []
            snv_alternatives = []
            for index, alt in enumerate(record.ALT):
                if len(alt) == len(record.REF) and alt.sequence[1:] == record.REF[1:]:
                    if not alt in snv_alternatives:
                        snv_alternatives.append(alt)
                else:
                    if not alt in indel_alternatives:
                        indel_alternatives.append(alt)
            for index, alt in enumerate(record.ALT):
                if len(alt) == len(record.REF):
                    # Check if alternative allele could be expressed as SNV
                    if alt.sequence[1:] == record.REF[1:]:
                        print "INFO: splitting hybrid record (containing both SNP and indel)", record
                        my_alternative = Substitution(alt.sequence[0])
                        my_samples = []
                        for call in record.samples:
                            my_samples.append(Call(call.site, call.sample, call.data.copy()))
                        my_record = Record(record.CHROM, record.POS, record.ID, record.REF[0], [my_alternative],
                                           record.QUAL, record.FILTER, record.INFO, record.FORMAT,
                                           record._sample_indexes, my_samples)
                        for call in my_record.samples:
                            call.add_format('FREQ', call.get_format('FREQ')[index])
                            if call.gt_nums is not None:
                                gt_alleles = []
                                for allele in call.gt_bases.split(call.gt_phase_char()):
                                    try:
                                        gt_alleles.append(str(snv_alternatives.index(allele)+1))
                                    except ValueError:
                                        gt_alleles.append('0')
                                call.add_format('GT', call.gt_phase_char().join(gt_alleles))
                                call.gt_nums = call.gt_phase_char().join(gt_alleles)
                        if snv_record is None:
                            snv_record = my_record
                        else:
                            if not snv_record.merge(my_record):
                                print "ERROR: Impossible to split record", record
                                raise
                else:
                    # Indel
                    my_samples = []
                    for call in record.samples:
                        my_samples.append(Call(call.site, call.sample, call.data.copy()))
                    my_record = Record(record.CHROM, record.POS, record.ID, record.REF, [alt],
                                       record.QUAL, record.FILTER, record.INFO, record.FORMAT,
                                       record._sample_indexes, my_samples)
                    for call in my_record.samples:
                        call.add_format('FREQ', call.get_format('FREQ')[index])
                        if call.gt_nums is not None:
                            gt_alleles = []
                            for allele in call.gt_bases.split(call.gt_phase_char()):
                                try:
                                    gt_alleles.append(str(indel_alternatives.index(allele)+1))
                                except ValueError:
                                    gt_alleles.append('0')
                            call.add_format('GT', call.gt_phase_char().join(gt_alleles))
                            call.gt_nums = call.gt_phase_char().join(gt_alleles)
                    if indel_record is None:
                        indel_record = my_record
                    else:
                        if not indel_record.merge(my_record):
                            print "ERROR: Impossible to split record", record
                            raise
        if snv_record is not None:
            repaired_vcf.write_record(snv_record)
            repaired_vcf.write_record(indel_record)
        else:
            repaired_vcf.write_record(record)


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-i', '--input_vcf', action='store', type='string', dest="input_vcf", default=None,
                          help='Path to input vcf variant calling file')
        parser.add_option('-o', '--output_vcf', action='store', type='string', dest='output_vcf',
                          default='repaired_output.vcf', help='Output name')

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