#!/usr/bin/env python

'''
Created on Aug 22, 2013
@author: Guillermo Marco Puche
@email: guillermo.marco@sistemasgenomicos.com
'''

#Python imports
import argparse
import os.path
import sys
import pydevd
import time
import datetime
import re

#Self imports
from Sample import Sample
from ConsequenceCounter import ConsequenceCounter

def menu ():
    
    """Input menu for user."""
    
    
    parser = argparse.ArgumentParser(description='Consequence Stats v0.1',
                                     epilog="################################")

    required = parser.add_argument_group("Required options")
        
    required.add_argument('-i', dest='input_psv', action = 'store', 
                        help='Input: *.psv annotation file. ie: /home/annotation.psv',
                        required=True)
    
    required.add_argument('-d', dest='depth_threshold', type=float, action = 'store', 
                        help="Depth threshold. ie: -d 10"
                        ,required=True)
    
    required.add_argument('-freq1', dest='freq_1', type=float, action = 'store', 
                        help="Freq ref threshold. ie: -fref 0.12",  
                        required=True)
    
    required.add_argument('-freq2', dest='freq_2', type=float, action = 'store', 
                        help="Freq var threshold. ie: -fvar 0.75",
                        required=True)
    
    required.add_argument('-o', dest='output_file', action = 'store', default = './consequence_stats.txt',
                        help="Output file. If path to file doesn't exists it will be automatically created. ie: -output /home/results/consequence_stats.txt",
                        required=False)

    required.add_argument('-t', dest='remote_debug', action='store_true',
                          help="Remote debug action.", required=False)
    
    return parser.parse_args()

def create_Samples_from_header(line):
    
    """Creates sample objects for input file.
    It uses CN and CP input lists from users.
    
    It creates a sample object for each one containing the name
    and the column index for genotype, depth amd var/depth (ratio) columns    
    """
    samples_name_list = []
    sample_return_list = []
    
    geno_sufix = '.replicate1_Genotype'
    depth_sufix = '.replicate1_Depth'
    ratio_sufix = '.replicate1_Var/Depth'
    consequence = 'Variant_effect'
    
    columns = line.split('|')

    for i in columns:
        res = re.search('(.*)\.replicate1_Genotype$', i)
        if res: samples_name_list.append(res.group(1))

    for sample_name in samples_name_list:
        genotype_col = columns.index(sample_name+geno_sufix) 
        depth_col = columns.index(sample_name+depth_sufix)
        ratio_col = columns.index(sample_name+ratio_sufix)
        
        #We get consequence column
        consequence_col = columns.index(consequence)
        
        #We create the sample object
        new_sample = Sample(sample_name, genotype_col, depth_col, ratio_col)
        
        #We set the consequence column
        new_sample.set_consequence_col(consequence_col)
        
        #We create and set a ConsequenceCounter object to the Sample.
        consequence_counter = ConsequenceCounter()
        
        #We set the ConsequenceCounter object to Sample object
        new_sample.set_consequence_counter(consequence_counter)
        
        sample_return_list.append(new_sample)
            
    return sample_return_list

def check_samples_depth (line, sample_list, depth_threshold):
    
    """Checks sample depths are above depth_threshold
    Returns true to skip line if there's at least one sample
    above depth threshold.
    
    Returns False (to continue program flow) and also:
    """
    
    sample_depth_true = []
    
    columns = line.split('|')

    for sample in sample_list:
        #We store samples that are greater or equal to depth_threshold
        if columns[sample.get_depth_col()] == "-": columns[sample.get_depth_col()] = 0.0
        if float(columns[sample.get_depth_col()]) >= depth_threshold: sample_depth_true.append(sample)
    
    #Return True to skip_line if no cp_depth_true samples
    if len(sample_depth_true) == 0: return True, None
    
    else:
        #Return False to not skip line and continue main execution and upcomming modules
        return False, sample_depth_true

def categorize_variants(line, samples_depth_true, freq1, freq2):
    
    #Convert to float just in case..
    freq1 = float(freq1)
    freq2 = float(freq2)
    
    #Split line
    columns = line.split('|')
    
    for sample in samples_depth_true:
        
        #Check if variant depth/var ratio freq is between freq1 and freq2 if so it's hetero.
        if freq1 <= float(columns[sample.get_ratio_col()]) < freq2:
            count_line_consequences(line, sample)
        
        #Else check if depth/var ratio freq is greater than freq2 it's homo.
        elif float(columns[sample.get_ratio_col()]) >= freq2:
            count_line_consequences(line, sample)

def count_line_consequences(line, sample):
  
    consequence_terms = ['transcript_ablation','splice_donor_variant','splice_acceptor_variant','stop_gained','frameshift_variant','stop_lost','initiator_codon_variant','inframe_insertion','inframe_deletion','missense_variant','transcript_amplification','splice_region_variant','incomplete_terminal_codon_variant','synonymous_variant','stop_retained_variant','coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant','non_coding_exon_variant','nc_transcript_variant','intron_variant','NMD_transcript_variant','upstream_gene_variant','downstream_gene_variant','TFBS_ablation','TFBS_amplification','TF_binding_site_variant','regulatory_region_variant','regulatory_region_ablation','regulatory_region_amplification','feature_elongation','feature_truncation','intergenic_variant']
    #Split line
    columns = line.split('|')
    
    #print columns[sample.get_consequence_col()]
    for consequence_term in consequence_terms:
        #consequence_term_lc = consequence_term.lower()
        consequence_method = 'increment_'+consequence_term
        res = re.search(consequence_term, columns[sample.get_consequence_col()])
        
        if res: getattr(sample.get_consequence_counter(), consequence_method)()

def file_exists (input_file):
    
    """Checks if input_file exists"""
    
    if not os.path.exists(input_file): raise IOError("The file %s does not exist!" %input_file)

def output_mgmt (output_file_path):
    
    """Creates output directories, even if nested in
    case that they don't exist."""
    
    if output_file_path == '.':
        raise ValueError('The path {0} is not a valid option. Leave output option blank or specify a valid output file name.'.format(output_file_path))
    
    output_base_path = os.path.dirname(output_file_path)
    
    if not os.path.exists(output_base_path):
        print "Creating output directory..."
        os.makedirs(output_base_path)
    
    return output_base_path

def working_message (args):
    """Returns message with the options that the
    user has chosen"""
    
    message = 'Consequence Stats v0.1 working...'
    message += '\n=============================='
    message += '\nInput file: %s' %args.input_psv
    message += '\nDepth Threshold: %s' %args.depth_threshold
    message += '\nFreq 1 threshold: %s' %args.freq_1
    message += '\nFreq 2 threshold: %s' %args.freq_2
    message += '\nOutput: %s' %args.output_file
    message += '\n=============================='
    
    return message


def calc_percentage(count, total):
    percentage = (1.0 * count / total) * 100
    rounded_percentage = round(percentage, 2)
    return rounded_percentage


def generate_report (sample_list, output_file):
    
    consequence_terms = ['transcript_ablation','splice_donor_variant','splice_acceptor_variant','stop_gained','frameshift_variant','stop_lost','initiator_codon_variant','inframe_insertion','inframe_deletion','missense_variant','transcript_amplification','splice_region_variant','incomplete_terminal_codon_variant','synonymous_variant','stop_retained_variant','coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant','3_prime_UTR_variant','non_coding_exon_variant','nc_transcript_variant','intron_variant','NMD_transcript_variant','upstream_gene_variant','downstream_gene_variant','TFBS_ablation','TFBS_amplification','TF_binding_site_variant','regulatory_region_variant','regulatory_region_ablation','regulatory_region_amplification','feature_elongation','feature_truncation','intergenic_variant']
    
    #Construct the header
    header = '#Variant_effect\t'
    
    for sample in sample_list:
        header+=sample.get_name()+'\t'
    
    header+='\n'
    
    #Write header to file
    output_file.write(header)
    
    #Write row per consequence
    total_cons = 0
    for sample in sample_list:
        total_cons += (sample.get_consequence_counter().get_total_consequences())

    total_percent = 0
    for consequence_term in consequence_terms:
        consequence_values = []
        output_file.write(consequence_term+'\t')
        #consequence_term_lc = consequence_term.lower()
        consequence_method = 'get_'+consequence_term

        for sample in sample_list:
            conseq_counter = (int(getattr(sample.get_consequence_counter(), consequence_method)()))
            conseq_percentage = calc_percentage(conseq_counter, total_cons)
            total_percent += conseq_percentage
            consequence_values.append(str(conseq_percentage))
            #consequence_values.append(str(getattr(sample.get_consequence_counter(), consequence_method)()))


            
        output_file.write('\t'.join(consequence_values))
        output_file.write('\n')

if __name__ == '__main__':
     
    #Start timer counter.
    start_time = time.time()
    
    #Parse user menu options.
    args = menu()
    
    #Check input file exists
    file_exists(args.input_psv)
    
    #Check if freq2 is > than freq1
    if args.freq_2 < args.freq_1:
        raise ValueError('ERROR: freq2 must be greater than freq1 !')

    #Check for remote debug option
    if args.remote_debug:
        user_ip = os.environ['USERIP']
        pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)

    #Print working message.
    print working_message (args)
    
    #Create output dirs if needed.
    output_base_path = output_mgmt(args.output_file)
    
    #Open output file.
    output_file = open(args.output_file, 'w')
    
    #We start to parse the file
    with open(args.input_psv, 'r') as input_file:
                
        for line in input_file: 

            #If header line.                        
            if line.startswith('#'):
                #We create our sample objects list
                sample_list = create_Samples_from_header(line)
                    
            else:
                split_line = line.split('|')
                    
                skip_variant, samples_depth_true = check_samples_depth(line, sample_list, args.depth_threshold)
                if skip_variant: continue
                
                categorized = categorize_variants(line, samples_depth_true, args.freq_1, args.freq_2)                    
                
        
    generate_report(sample_list, output_file)
          
    time = int(time.time() - start_time)
    print ('DONE. Time elapsed: %s') %str(datetime.timedelta(seconds=time))