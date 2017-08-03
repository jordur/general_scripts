#!/usr/bin/env python
"""
SYNOPSIS
    genotype.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    genotype: Module to get information from annotation files
EXAMPLES
    genotype.py
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    genotype.py v0.1
"""

from vcf.model import _Call


def chromosomes():
    """Function that returns the human chromosomes names"""
    return ' '.join(str(i) for i in xrange(1,24))," X Y M"


def get_chromosome_number(chromosome):
    """ Function that returns a sorted number for each chromosome (X=23, Y=24, M=25)"""
    if isinstance(chromosome, str):
        if chromosome.startswith('chr'):
            chromosome = chromosome[3:]
    if chromosome == "X" or chromosome == "x":
        return 23
    elif chromosome == "Y" or chromosome == "y":
        return 24
    elif chromosome == "M" or chromosome == "m":
        return 25
    else:
        try:
            chr = int(chromosome)
        except ValueError:
            chr = chromosome
        return chr


def chromosome_number2str(chromosome):
    """Returns the string chrXX corresponding to a int chromosome"""
    try:
        chrom = int(chromosome)
        chrom = 'X' if chrom == 23 else 'Y' if chrom == 24 else 'M' if chrom == 25 else chrom
        return 'chr'+str(chrom)
    except ValueError:
        return chromosome
    except Exception:
        print "ERROR: problem in chromosome_number2str for chromosome", chromosome
        raise


def get_genotype(call=None, var_allele=None, ploidy=2):
    """ Function that returns the predicted genotype inferred from allelic frequency """
    try:
        if call.aaf is None:
            freq = None
        else:
            if isinstance(var_allele, str):
                var_allele = call.site.ALT.index(var_allele)
            if var_allele is None:
                freq = max(af if af != '.' else 0 for af in call.aaf)
            else:
                if call.aaf[var_allele] == '.':
                    freq = None
                else:
                    freq = float(call.aaf[var_allele])
    except AttributeError:
        return 'Uncovered', float(0)
    except IndexError:
        print 'WARNING: Impossible to guess genotype for call in site', call.site
        return 'Unknown', float(0)
    except Exception:
        print "ERROR in call.site", call.site
        raise
    if freq is None:
        freq = 0
        if call.depth is None or call.depth == 0:
            genotype = 'Uncovered'
        else:
            genotype = 'Unknown'
    elif call.site.is_indel:
        if ploidy == 2:
            if freq == 0:
                genotype = 'Homo_ref'
            elif 0 < freq < 0.3:
                genotype = 'UNC_Hetero'
            elif 0.3 <= freq < 0.6:
                genotype = 'P_Hetero'
            elif 0.6 <= freq:
                genotype = 'P_Homo_var'
            else:
                print "ERROR: Unforeseen case in genotype.py, call.site", call.site
                raise Exception
        elif ploidy == 1:
            if freq == 0:
                genotype = 'Homo_ref'
            elif 0 < freq < 0.3:
                genotype = 'UNC_Homo'
            elif 0.3 <= freq:
                genotype = 'P_Homo_var'
            else:
                print "ERROR: Unforeseen case in genotype.py, call.site", call.site
                raise Exception
        else:
            print "ERROR: unforeseen ploidiy in genotype.py, call.site", call.site
            raise Exception
    else:
        if ploidy == 2:
            if freq == 0:
                genotype = 'Homo_ref'
            elif 0 < freq < 0.12:
                genotype = 'P_Homo_ref'
            elif 0.12 <= freq < 0.35:
                genotype = 'UNC_Hetero'
            elif 0.35 <= freq < 0.65:
                genotype = 'P_Hetero'
            elif 0.65 <= freq < 0.85:
                genotype = 'UNC_Homo'
            elif 0.85 <= freq:
                genotype = 'P_Homo_var'
            else:
                print "ERROR: Unforeseen case in genotype.py, call.site", call.site
                raise Exception
        elif ploidy == 1:
            if freq == 0:
                genotype = 'Homo_ref'
            elif 0 < freq < 0.12:
                genotype = 'P_Homo_ref'
            elif 0.12 <= freq < 0.85:
                genotype = 'UNC_Homo'
            elif 0.85 <= freq:
                genotype = 'P_Homo_var'
            else:
                print "ERROR: Unforeseen case in genotype.py, call.site", call.site
                raise Exception
        else:
            print "ERROR: unforeseen ploidiy in genotype.py, call.site", call.site
            raise Exception
    return genotype, freq

if __name__ == '__main__':
    True