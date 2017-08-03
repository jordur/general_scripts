"""
Created on Aug 13, 2013

@author: SGNET\gmarco
"""


class ConsequenceCounter(object):
    """
    classdocs
    """
    #todo: GET, SET, DEL minusculas, increment mayus

    def __init__(self):
        self.__total_consequences = 0
        self.__transcript_ablation = 0
        self.__splice_donor_variant = 0
        self.__splice_acceptor_variant = 0
        self.__stop_gained = 0
        self.__frameshift_variant = 0
        self.__stop_lost = 0
        self.__initiator_codon_variant = 0
        self.__inframe_insertion = 0
        self.__inframe_deletion = 0
        self.__missense_variant = 0
        self.__transcript_amplification = 0
        self.__splice_region_variant = 0
        self.__incomplete_terminal_codon_variant = 0
        self.__synonymous_variant = 0
        self.__stop_retained_variant = 0
        self.__coding_sequence_variant = 0
        self.__mature_miRNA_variant = 0
        self.__5_prime_UTR_variant = 0
        self.__3_prime_UTR_variant = 0
        self.__non_coding_exon_variant = 0
        self.__nc_transcript_variant = 0
        self.__intron_variant = 0
        self.__NMD_transcript_variant = 0
        self.__upstream_gene_variant = 0
        self.__downstream_gene_variant = 0
        self.__TFBS_ablation = 0
        self.__TFBS_amplification = 0
        self.__TF_binding_site_variant = 0
        self.__regulatory_region_variant = 0
        self.__regulatory_region_ablation = 0
        self.__regulatory_region_amplification = 0
        self.__feature_elongation = 0
        self.__feature_truncation = 0
        self.__intergenic_variant = 0

    def get_total_consequences(self):
        return self.__total_consequences

    def get_transcript_ablation(self):
        return self.__transcript_ablation

    def get_splice_donor_variant(self):
        return self.__splice_donor_variant

    def get_splice_acceptor_variant(self):
        return self.__splice_acceptor_variant

    def get_stop_gained(self):
        return self.__stop_gained

    def get_frameshift_variant(self):
        return self.__frameshift_variant

    def get_stop_lost(self):
        return self.__stop_lost

    def get_initiator_codon_variant(self):
        return self.__initiator_codon_variant

    def get_inframe_insertion(self):
        return self.__inframe_insertion

    def get_inframe_deletion(self):
        return self.__inframe_deletion

    def get_missense_variant(self):
        return self.__missense_variant

    def get_transcript_amplification(self):
        return self.__transcript_amplification

    def get_splice_region_variant(self):
        return self.__splice_region_variant

    def get_incomplete_terminal_codon_variant(self):
        return self.__incomplete_terminal_codon_variant

    def get_synonymous_variant(self):
        return self.__synonymous_variant

    def get_stop_retained_variant(self):
        return self.__stop_retained_variant

    def get_coding_sequence_variant(self):
        return self.__coding_sequence_variant

    def get_mature_miRNA_variant(self):
        return self.__mature_miRNA_variant

    def get_5_prime_UTR_variant(self):
        return self.__5_prime_UTR_variant

    def get_3_prime_UTR_variant(self):
        return self.__3_prime_UTR_variant

    def get_non_coding_exon_variant(self):
        return self.__non_coding_exon_variant

    def get_nc_transcript_variant(self):
        return self.__nc_transcript_variant

    def get_intron_variant(self):
        return self.__intron_variant

    def get_NMD_transcript_variant(self):
        return self.__NMD_transcript_variant

    def get_upstream_gene_variant(self):
        return self.__upstream_gene_variant

    def get_downstream_gene_variant(self):
        return self.__downstream_gene_variant

    def get_TFBS_ablation(self):
        return self.__TFBS_ablation

    def get_TFBS_amplification(self):
        return self.__TFBS_amplification

    def get_TF_binding_site_variant(self):
        return self.__TF_binding_site_variant

    def get_regulatory_region_variant(self):
        return self.__regulatory_region_variant

    def get_regulatory_region_ablation(self):
        return self.__regulatory_region_ablation

    def get_regulatory_region_amplification(self):
        return self.__regulatory_region_amplification

    def get_feature_elongation(self):
        return self.__feature_elongation

    def get_feature_truncation(self):
        return self.__feature_truncation

    def get_intergenic_variant(self):
        return self.__intergenic_variant

    def set_total_consequences(self, value):
        self.__total_consequences = value

    def set_transcript_ablation(self, value):
        self.__transcript_ablation = value

    def set_splice_donor_variant(self, value):
        self.__splice_donor_variant = value

    def set_splice_acceptor_variant(self, value):
        self.__splice_acceptor_variant = value

    def set_stop_gained(self, value):
        self.__stop_gained = value

    def set_frameshift_variant(self, value):
        self.__frameshift_variant = value

    def set_stop_lost(self, value):
        self.__stop_lost = value

    def set_initiator_codon_variant(self, value):
        self.__initiator_codon_variant = value

    def set_inframe_insertion(self, value):
        self.__inframe_insertion = value

    def set_inframe_deletion(self, value):
        self.__inframe_deletion = value

    def set_missense_variant(self, value):
        self.__missense_variant = value

    def set_transcript_amplification(self, value):
        self.__transcript_amplification = value

    def set_splice_region_variant(self, value):
        self.__splice_region_variant = value

    def set_incomplete_terminal_codon_variant(self, value):
        self.__incomplete_terminal_codon_variant = value

    def set_synonymous_variant(self, value):
        self.__synonymous_variant = value

    def set_stop_retained_variant(self, value):
        self.__stop_retained_variant = value

    def set_coding_sequence_variant(self, value):
        self.__coding_sequence_variant = value

    def set_mature_mi_rna_variant(self, value):
        self.__mature_miRNA_variant = value

    def set_5_prime_UTR_variant(self, value):
        self.__5_prime_UTR_variant = value

    def set_3_prime_UTR_variant(self, value):
        self.__3_prime_UTR_variant = value

    def set_non_coding_exon_variant(self, value):
        self.__non_coding_exon_variant = value

    def set_nc_transcript_variant(self, value):
        self.__nc_transcript_variant = value

    def set_intron_variant(self, value):
        self.__intron_variant = value

    def set_NMD_transcript_variant(self, value):
        self.__NMD_transcript_variant = value

    def set_upstream_gene_variant(self, value):
        self.__upstream_gene_variant = value

    def set_downstream_gene_variant(self, value):
        self.__downstream_gene_variant = value

    def set_TFBS_ablation(self, value):
        self.__TFBS_ablation = value

    def set_TFBS_amplification(self, value):
        self.__TFBS_amplification = value

    def set_TF_binding_site_variant(self, value):
        self.__TF_binding_site_variant = value

    def set_regulatory_region_variant(self, value):
        self.__regulatory_region_variant = value

    def set_regulatory_region_ablation(self, value):
        self.__regulatory_region_ablation = value

    def set_regulatory_region_amplification(self, value):
        self.__regulatory_region_amplification = value

    def set_feature_elongation(self, value):
        self.__feature_elongation = value

    def set_feature_truncation(self, value):
        self.__feature_truncation = value

    def set_intergenic_variant(self, value):
        self.__intergenic_variant = value

    def del_total_consequences(self):
        del self.__total_consequences

    def del_transcript_ablation(self):
        del self.__transcript_ablation

    def del_splice_donor_variant(self):
        del self.__splice_donor_variant

    def del_splice_acceptor_variant(self):
        del self.__splice_acceptor_variant

    def del_stop_gained(self):
        del self.__stop_gained

    def del_frameshift_variant(self):
        del self.__frameshift_variant

    def del_stop_lost(self):
        del self.__stop_lost

    def del_initiator_codon_variant(self):
        del self.__initiator_codon_variant

    def del_inframe_insertion(self):
        del self.__inframe_insertion

    def del_inframe_deletion(self):
        del self.__inframe_deletion

    def del_missense_variant(self):
        del self.__missense_variant

    def del_transcript_amplification(self):
        del self.__transcript_amplification

    def del_splice_region_variant(self):
        del self.__splice_region_variant

    def del_incomplete_terminal_codon_variant(self):
        del self.__incomplete_terminal_codon_variant

    def del_synonymous_variant(self):
        del self.__synonymous_variant

    def del_stop_retained_variant(self):
        del self.__stop_retained_variant

    def del_coding_sequence_variant(self):
        del self.__coding_sequence_variant

    def del_mature_mi_rna_variant(self):
        del self.__mature_miRNA_variant

    def del_5_prime_UTR_variant(self):
        del self.__5_prime_UTR_variant

    def del_3_prime_UTR_variant(self):
        del self.__3_prime_UTR_variant

    def del_non_coding_exon_variant(self):
        del self.__non_coding_exon_variant

    def del_nc_transcript_variant(self):
        del self.__nc_transcript_variant

    def del_intron_variant(self):
        del self.__intron_variant

    def del_NMD_transcript_variant(self):
        del self.__NMD_transcript_variant

    def del_upstream_gene_variant(self):
        del self.__upstream_gene_variant

    def del_downstream_gene_variant(self):
        del self.__downstream_gene_variant

    def del_TFBS_ablation(self):
        del self.__TFBS_ablation

    def del_TFBS_amplification(self):
        del self.__TFBS_amplification

    def del_TF_binding_site_variant(self):
        del self.__TF_binding_site_variant

    def del_regulatory_region_variant(self):
        del self.__regulatory_region_variant

    def del_regulatory_region_ablation(self):
        del self.__regulatory_region_ablation

    def del_regulatory_region_amplification(self):
        del self.__regulatory_region_amplification

    def del_feature_elongation(self):
        del self.__feature_elongation

    def del_feature_truncation(self):
        del self.__feature_truncation

    def del_intergenic_variant(self):
        del self.__intergenic_variant

    def increment_transcript_ablation(self):
        self.__transcript_ablation += 1
        self.__total_consequences += 1

    def increment_splice_donor_variant(self):
        self.__splice_donor_variant += 1
        self.__total_consequences += 1

    def increment_splice_acceptor_variant(self):
        self.__splice_acceptor_variant += 1
        self.__total_consequences += 1

    def increment_stop_gained(self):
        self.__stop_gained += 1
        self.__total_consequences += 1

    def increment_frameshift_variant(self):
        self.__frameshift_variant += 1
        self.__total_consequences += 1

    def increment_stop_lost(self):
        self.__stop_lost += 1
        self.__total_consequences += 1

    def increment_initiator_codon_variant(self):
        self.__initiator_codon_variant += 1
        self.__total_consequences += 1

    def increment_inframe_insertion(self):
        self.__inframe_insertion += 1
        self.__total_consequences += 1

    def increment_inframe_deletion(self):
        self.__inframe_deletion += 1
        self.__total_consequences += 1

    def increment_missense_variant(self):
        self.__missense_variant += 1
        self.__total_consequences += 1

    def increment_transcript_amplification(self):
        self.__transcript_amplification += 1
        self.__total_consequences += 1

    def increment_splice_region_variant(self):
        self.__splice_region_variant += 1
        self.__total_consequences += 1

    def increment_incomplete_terminal_codon_variant(self):
        self.__incomplete_terminal_codon_variant += 1
        self.__total_consequences += 1

    def increment_synonymous_variant(self):
        self.__synonymous_variant += 1
        self.__total_consequences += 1

    def increment_stop_retained_variant(self):
        self.__stop_retained_variant += 1
        self.__total_consequences += 1

    def increment_coding_sequence_variant(self):
        self.__coding_sequence_variant += 1
        self.__total_consequences += 1

    def increment_mature_miRNA_variant(self):
        self.__mature_miRNA_variant += 1
        self.__total_consequences += 1

    def increment_5_prime_UTR_variant(self):
        self.__5_prime_UTR_variant += 1
        self.__total_consequences += 1

    def increment_3_prime_UTR_variant(self):
        self.__3_prime_UTR_variant += 1
        self.__total_consequences += 1

    def increment_non_coding_exon_variant(self):
        self.__non_coding_exon_variant += 1
        self.__total_consequences += 1

    def increment_nc_transcript_variant(self):
        self.__nc_transcript_variant += 1
        self.__total_consequences += 1

    def increment_intron_variant(self):
        self.__intron_variant += 1
        self.__total_consequences += 1

    def increment_NMD_transcript_variant(self):
        self.__NMD_transcript_variant += 1
        self.__total_consequences += 1

    def increment_upstream_gene_variant(self):
        self.__upstream_gene_variant += 1
        self.__total_consequences += 1

    def increment_downstream_gene_variant(self):
        self.__downstream_gene_variant += 1
        self.__total_consequences += 1

    def increment_TFBS_ablation(self):
        self.__TFBS_ablation += 1
        self.__total_consequences += 1

    def increment_TFBS_amplification(self):
        self.__TFBS_amplification += 1
        self.__total_consequences += 1

    def increment_TF_binding_site_variant(self):
        self.__TF_binding_site_variant += 1
        self.__total_consequences += 1

    def increment_regulatory_region_variant(self):
        self.__regulatory_region_variant += 1
        self.__total_consequences += 1

    def increment_regulatory_region_ablation(self):
        self.__regulatory_region_ablation += 1
        self.__total_consequences += 1

    def increment_regulatory_region_amplification(self):
        self.__regulatory_region_amplification += 1
        self.__total_consequences += 1

    def increment_feature_elongation(self):
        self.__feature_elongation += 1
        self.__total_consequences += 1

    def increment_feature_truncation(self):
        self.__feature_truncation += 1
        self.__total_consequences += 1

    def increment_intergenic_variant(self):
        self.__intergenic_variant += 1
        self.__total_consequences += 1