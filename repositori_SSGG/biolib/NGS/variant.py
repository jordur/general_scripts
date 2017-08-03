#!/usr/bin/env python
"""
SYNOPSIS
    variant.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    Module to get information from annotation files
EXAMPLES
    variant.py
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    variant.py v0.1
"""
#import hgvs
#import hgvs.utils
from pygr.seqdb import SequenceFileDB
#from hgvs import InvalidHGVSName
from genotype import get_genotype
from collections import OrderedDict
from fisher import pvalue
from math import log


# Read genome sequence using pygr.
#genome = SequenceFileDB('/share/references/genomes/human/hg19/reference/human_hg19.fa')

# Read RefSeq transcripts into a python dict.
#with open('/share/references/genes/human/Homo_sapiens.GRCh37.73.modified_genePred') as infile:
#    transcripts = hgvs.utils.read_transcripts(infile)


# Provide a callback for fetching a transcript by its name.
def get_transcript(name):
    # Following code line is for compatibility with Ensembl GenePred files:
    # Ensembl transcripts don't include transcript version in .gtf files (gene set files). To assure version
    # compatibility, Ensembl annotation and gtf files must have the same version!!
    name = name.split('.')[0] if name[:4] == 'ENST' else name
    return transcripts.get(name)


def norm_freq(freq):
    if isinstance(freq, list):
        result = []
        for myfreq in freq:
            result.append(norm_freq(myfreq))
        return result
    else:
        if freq is not None and freq != '.':
            return '{0:.3f}'.format(
                float(freq.replace('%', '')) / 100 if isinstance(freq, str) and '%' in freq else float(freq))
        else:
            return None


class Call(OrderedDict):
    """Class for call information"""

    def __init__(self):
        super(Call, self).__init__()
        self['Genotype'] = None
        self['Depth'] = None
        self['Var/Depth'] = None

    @property
    def genotype(self):
        return self['Genotype']

    @genotype.setter
    def genotype(self, mygenotype):
        self['Genotype'] = mygenotype

    @property
    def depth(self):
        return self['Depth']

    @depth.setter
    def depth(self, mydepth):
        self['Depth'] = mydepth

    @property
    def freq(self):
        return self['Var/Depth']

    @freq.setter
    def freq(self, myfreq):
        self['Var/Depth'] = myfreq


class ConsequenceType(OrderedDict):
    """Class definition for consequence types (a variant may contain several consequence types)"""

    def __init__(self, samples=OrderedDict(), organism_type='model'):
        super(ConsequenceType, self).__init__()
        if organism_type == 'model':
            self.add_arg('SYMBOL', 'HGNC_symbol')
            self.add_arg('Chr', 'Chr')
            self.add_arg('Pos', 'Pos')
            self.add_arg('Ref_Allele', 'Ref_Allele')
            self.add_arg('Var_Allele', 'Var_Allele')
            self['Samples'] = OrderedDict()
            for sample in samples:
                self['Samples'][sample] = Call()
            self.add_arg('Strand_Bias', 'Strand_Bias')
            self.add_arg('Gene', 'Gene')
            self.add_arg('Gene_description', 'Gene_description')
            self.add_arg('HGVSc_name', 'HGVSc_name')
            self.add_arg('ALLELE_NUM', 'ALLELE_NUM')
            self.add_arg('Intron', 'Intron')
            self.add_arg('Exon', 'Exon')
            self.add_arg('HGVSp', 'HGVSp_name')
            self.add_arg('Variant_effect', 'Variant_effect')
            self.add_arg('ALL_MAF', 'ALL_MAF')
            self.add_arg('AFR_MAF', 'AFR_MAF')
            self.add_arg('AMR_MAF', 'AMR_MAF')
            self.add_arg('ASN_MAF', 'ASN_MAF')
            self.add_arg('EUR_MAF', 'EUR_MAF')
            self.add_arg('InterPro_IDs', 'InterPro_IDs')
            self.add_arg('InterPro_descriptions', 'InterPro_descriptions')
            self.add_arg('HGMD_info', 'HGMD_info')
            self.add_arg('Related_publication', 'Related_publication')
            self.add_arg('Existing_variation', 'Existing_variation')
            self.add_arg('Feature_type', 'Feature_type')
            self.add_arg('Feature_ID', 'Feature_ID')
            self.add_arg('RefSeq_ID', 'RefSeq_ID')
            self.add_arg('CCDS_ID', 'CCDS_ID')
            self.add_arg('Canonical_isoform', 'Canonical_isoform')
            self.add_arg('Conservation_score', 'Conservation_score')
            self.add_arg('Grantham_distance', 'Grantham_distance')
            self.add_arg('Condel_prediction', 'Condel_prediction')
            self.add_arg('SIFT_prediction', 'SIFT_prediction')
            self.add_arg('PolyPhen_prediction', 'PolyPhen_prediction')
            self.add_arg('Affected_prot_domains', 'Affected_prot_domains')
            self.add_arg('Regulatory_Motif_name', 'Regulatory_Motif_name')
            self.add_arg('Regulatory_Motif_position', 'Regulatory_Motif_position')
            self.add_arg('Regulatory_High_Inf_Pos', 'Regulatory_High_Inf_Pos')
            self.add_arg('Regulatory_Motif_Score_Change', 'Regulatory_Motif_Score_Change')
            self.add_arg('Flanking_sequence', 'Flanking_sequence')
        elif organism_type == 'non_model':
            self.add_arg('Chr', 'Chr')
            self.add_arg('Pos', 'Pos')
            self.add_arg('Ref_Allele', 'Ref_Allele')
            self.add_arg('Var_Allele', 'Var_Allele')
            self['Samples'] = OrderedDict()
            for sample in samples:
                self['Samples'][sample] = Call()
            self.add_arg('Strand_Bias', 'Strand_Bias')
            self.add_arg('Allele', 'Allele')
            self.add_arg('Gene', 'Gene')
            self.add_arg('Feature', 'Feature_ID')
            self.add_arg('Feature_type', 'Feature_type')
            self.add_arg('Variant_effect', 'Variant_effect')
            self.add_arg('cDNA_position', 'cDNA_position')
            self.add_arg('CDS_position', 'CDS_position')
            self.add_arg('Protein_position', 'Protein_position')
            self.add_arg('Amino_acids', 'Amino_acids')
            self.add_arg('Codons', 'Codons')
            self.add_arg('Existing_variation', 'Existing_variation')
            self.add_arg('ALLELE_NUM', 'Allele_num')
            self.add_arg('DISTANCE', 'Distance')

    def add_arg(self, name=None, psv_name=None):
        self[name] = {}
        self[name]['psv_name'] = psv_name
        self[name]['value'] = None

    def get_arg(self, name=None):
        return self[name]['value'] if name is not None else None

    def get_arg_psv_name(self, name=None):
        return self[name]['psv_name'] if name is not None else None

    def set_arg(self, name=None, value=None):
        self[name]['value'] = value

    def sample_call(self, name=None):
        return self['Samples'][name] if name is not None else None


def normalize_variant_site_notation(ref, alt, type='hgvs'):
    ref_normalized = ''
    alt_normalized = ''
    if alt == '' or ref == '':
        return ref, alt
    elif len(ref) >= len(alt):
        j = 0
        discrepant_zone = False
        for i in range(0, len(ref)):
            if j < len(alt) and ref[i] == alt[j] and (not discrepant_zone or len(ref) - i <= len(alt) - j):
                j += 1
            else:
                discrepant_zone = True
                ref_normalized += ref[i]
                if len(alt) == len(ref):
                    alt_normalized += alt[i]
                    # VCF format uses most left-aligned genomic coordinate
            discrepant_zone = True if type == 'VCF' else discrepant_zone
        if j != len(alt):
            return ref, alt
        else:
            return ref_normalized, alt_normalized
    else:
        alt_normalized, ref_normalized = normalize_variant_site_notation(alt, ref, type=type)
        return ref_normalized, alt_normalized


def compare_changes(ref1, alt1, ref2, alt2):
    return True if (ref1 == ref2 and alt1 == alt2) or (
        reverse_complement(ref1) == ref2 and reverse_complement(alt1) == alt2) else False


def obtain_vcf_alternative_allele(ref1, alt1, ref2, valt2):
    """Returns the vcf alternative allele relative to the hgvs consequence
    """
    r1, a1 = normalize_variant_site_notation(ref1, alt1, type='hgvs')
    for allele in valt2:
        r2, a2 = normalize_variant_site_notation(ref2, allele, type='VCF')
        if compare_changes(r1, a1, r2, a2):
            return allele
    print "WARNING: No suitable VCF-alternative allele (ref:", ref2, ", alt:", valt2, "found for hgvs mutation (ref:", \
        ref1, ", alt:", alt1, ")."
    return valt2


def reverse_complement(sequence):
    """Return the reverse complement of the dna string."""

    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}

    reverse_complement_sequence = ''

    sequence_list = list(sequence)
    sequence_list.reverse()

    for letter in sequence_list:
        reverse_complement_sequence += complement[letter.upper()]
    return reverse_complement_sequence


class Variant:
    """Class definition for variants"""

    def __init__(self, samples=OrderedDict(), organism_type='model', ploidy=2):
        self._organism_type = organism_type
        self._ploidy = ploidy
        self._taxonomy = {'genome': '', 'specie': ''}
        self._sequencing = {'NGS_id': '', 'seq_technology': '', 'platform_name': '', 'platform_version': '',
                            'capture_technology': ''}
        self._consequences = []
        self._consequences.append(ConsequenceType(organism_type=organism_type))
        self._samples = OrderedDict()
        for sample in samples:
            self._samples[sample] = Call()

    def get_from_record(self, record=None, separator='|'):
        if record is not None:
            self.consequences = []
            # Consequences in variant may be different to consequences in vcf, since different alleles will be
            # be considered as different consequences
            for consequence_number in range(0, len(record.INFO['CSQ'])):
                consequence = record.INFO['CSQ'][consequence_number]
                args = consequence.split(separator)
                try:
                    if self.organism_type == 'model':
                        alt = int(args[4])
                    else:
                        alt = int(args[11])
                    if alt > len(record.ALT):
                        print "WARNING: reported allele", alt, "not in alternative alleles list for record", record
                        alt = len(record.ALT) - 1
                        if alt > 0:
                            print "ERROR: Impossible to guess alternative allele!!"
                            quit()
                except IndexError:
                    print "ERROR: consequence not containing enough information even for getting alternative " \
                          "allele:", consequence
                    print "record:", record
                    raise
                except ValueError:
                    print "WARNING: wrong allele information in record", record
                    alt = []
                    for i in range(0, len(record.ALT)):
                        alt.append(i)
                if not isinstance(alt, list):
                    alt = [alt]
                for my_alt in alt:
                    self.consequences.append(ConsequenceType(self._samples, organism_type=self.organism_type))
                    index = 0
                    for arg_name in self.consequences[consequence_number].keys():
                        try:
                            arg = args[index]
                        except IndexError:
                            print "ERROR: No information for parameter", arg_name, "in CSQ. Index:", index
                            print "record:", record
                        if arg_name == 'Chr':
                            self.consequences[consequence_number].set_arg(arg_name, record.CHROM)
                        elif arg_name == 'Pos':
                            self.consequences[consequence_number].set_arg(arg_name, record.POS)
                        elif arg_name == 'Ref_Allele':
                            self.consequences[consequence_number].set_arg(arg_name, record.REF)
                        elif arg_name == 'Var_Allele':
                            self.consequences[consequence_number].set_arg(arg_name, str(record.ALT[my_alt-1]))
                        elif arg_name == 'Samples':
                            for sampindex, sample in enumerate(record.samples):
                                call = record.samples[sampindex]
                                self.sample_add(name=sample.sample, call=Call())
                                genotype, freq = get_genotype(call, ploidy=self.ploidy)
                                self.samples[sample.sample].genotype = genotype
                                self.samples[sample.sample].freq = freq

                                genotype, freq = get_genotype(call, var_allele=my_alt-1, ploidy=self.ploidy)
                                self.consequences[consequence_number]['Samples'][
                                    sample.sample].genotype = genotype
                                self.consequences[consequence_number]['Samples'][sample.sample].freq = freq
                                try:
                                    self.samples[sample.sample].depth = call.depth if call.depth is not None else '-'
                                    self.consequences[consequence_number]['Samples'][
                                        sample.sample].depth = call.depth if call.depth is not None else '-'
                                except AttributeError:
                                    self.samples[sample.sample].depth = '-'
                                    self.consequences[consequence_number]['Samples'][sample.sample].depth = '-'
                        elif arg_name == 'Strand_Bias':
                            # Obtain strand bias value and output it
                            try:
                                self.consequences[consequence_number].set_arg('Strand_Bias', '{0:.3f}'.format(
                                    record.INFO['FS']))
                            except KeyError:
                                # Assume that the variant site has been called by VarScan
                                ref_minus = 0
                                ref_plus = 0
                                alt_minus = 0
                                alt_plus = 0
                                for sampindex, sample in enumerate(record.samples):
                                    call = record.samples[sampindex]
                                    ref_plus += call.__getitem__('RDF') if call.__getitem__('RDF') is not None else 0
                                    ref_minus += call.__getitem__('RDR') if call.__getitem__('RDR') is not None else 0
                                    alt_plus += call.__getitem__('ADF') if call.__getitem__('ADF') is not None else 0
                                    alt_minus += call.__getitem__('ADR') if call.__getitem__('ADR') is not None else 0
                                strand_bias = '{0:.3f}'.format(
                                    -10 * log(pvalue(ref_plus, ref_minus, alt_plus, alt_minus).two_tail))
                                strand_bias = '0.000' if strand_bias == '-0.000' else strand_bias
                                self.consequences[consequence_number].set_arg('Strand_Bias', strand_bias)
                        else:
                            index += 1
                            if arg_name == 'ALL_MAF':
                                # ALL_MAF field includes by default the alternative allele that must be removed in psv
                                if ':' in arg:
                                    arg = arg[arg.index(':') + 1:]
                            elif arg_name == 'Flanking_sequence' and arg != '':
                                # Only include considered allele
                                begin = arg.index('[')
                                end = arg.index(']')
                                arg = arg[:begin] + '[' + record.REF + '/' + str(record.ALT[my_alt-1]) + arg[end:]
                            elif arg_name == 'Variant_effect' or arg_name == 'Existing_variation':
                                # Variant effect and Existing_variation: must be comma-separated
                                arg = arg.replace('&', ',')
                            self.consequences[consequence_number].set_arg(arg_name, arg)

    def create_psv_header(self, separator='|'):
        string = '#'
        for arg_name, value in self.consequences[0].items():
            if arg_name == 'Samples':
                for sample_name, sample in self.samples.items():
                    for call_name, call in sample.items():
                        string = string + sample_name + '_' + call_name + separator
            elif arg_name != 'ALLELE_NUM':
                string = string + self.consequences[0].get_arg_psv_name(arg_name) + separator
        length = len(separator)
        string = string[:-length]
        return string + '\n'

    def put_to_psv(self, separator='|'):
        lines = []
        for consequence in self.consequences:
            string = ''
            for arg_name, value in consequence.items():
                if arg_name == 'Samples':
                    for sample_name, sample in consequence['Samples'].items():
                        for call_name, call in sample.items():
                            string = string + str(call) + separator
                elif arg_name != 'ALLELE_NUM':
                    arg = str(consequence.get_arg(arg_name)) if consequence.get_arg(arg_name) != '' else '-'
                    string = string + arg + separator
            length = len(separator)
            string = string[:-length]
            lines.append(string)
        return '\n'.join(lines) + '\n'

    @property
    def taxonomy(self):
        return self._taxonomy

    @taxonomy.setter
    def taxonomy(self, **kwargs):
        for kwarg in kwargs:
            self._taxonomy[kwarg]

    @property
    def sequencing(self):
        return self._sequencing

    @property
    def consequences(self):
        return self._consequences

    @consequences.setter
    def consequences(self, consequences=None):
        self._consequences = consequences

    @property
    def samples(self):
        return self._samples

    def sample_add(self, name=None, call=None):
        if not name in self._samples:
            self._samples[name] = call

    @property
    def organism_type(self):
        return self._organism_type

    @property
    def ploidy(self):
        return self._ploidy