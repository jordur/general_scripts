#!/usr/bin/env python

'''
Created on Aug 22, 2013
@author: Guillermo Marco Puche
@email: guillermo.marco@sistemasgenomicos.com
'''

import argparse
import os.path
import sys
import time
from Sample import Sample
from VariantCounter import VariantCounter
import re
#from variant_tools import *

def menu ():
    
    """Input menu for user."""
    
    
    parser = argparse.ArgumentParser(description='Variant Stats v0.1',
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
    
    required.add_argument('-o', dest='output_file', action = 'store', default = './variant_stats.txt',
                        help="Output file. If path to file doesn't exists it will be automatically created. ie: -output ./results/variant_stats.txt",
                        required=False)
    
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
    
    columns = line.split('|')

    for i in columns:
        res = re.search('(.*)\.replicate1_Genotype$', i)
        if res: samples_name_list.append(res.group(1))

    for sample_name in samples_name_list:
        genotype_col = columns.index(sample_name+geno_sufix) 
        depth_col = columns.index(sample_name+depth_sufix)
        ratio_col = columns.index(sample_name+ratio_sufix)
            
        new_sample = Sample(sample_name, genotype_col, depth_col, ratio_col)
        
        #We create and set a VariantCounter object to the Sample.      
        new_sample.set_variant_counter(VariantCounter())
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
            sample.get_variant_counter().increment_hetero_var()
        
        #Else check if depth/var ratio freq is greater than freq2 it's homo.
        elif float(columns[sample.get_ratio_col()]) >= freq2:
            sample.get_variant_counter().increment_homo_var()

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
    
    message = 'Variant Stats v0.1 working...'
    message += '\n=============================='
    message += '\nInput file: %s' %args.input_psv
    message += '\nDepth Threshold: %s' %args.depth_threshold
    message += '\nFreq 1 threshold: %s' %args.freq_1
    message += '\nFreq 2 threshold: %s' %args.freq_2
    message += '\nOutput: %s' %args.output_file
    message += '\n=============================='
    
    return message


def generate_report (sample_list, output_file):
    
    output_file.write('#Sample_Name\tTotal_Variants\tTotal_Hetero\tTotal_Homo')
    
    for sample in sample_list:
        output_file.write('\n{0}\t{1}\t{2}\t{3}'.format(sample.get_name(), sample.get_variant_counter().total_variants, sample.get_variant_counter().get_hetero_var(), sample.get_variant_counter().get_homo_var()))
        
    output_file.close()

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
    
    #Print working message.
    print working_message (args)
    
    #Create output dirs if needed.
    output_base_path = output_mgmt(args.output_file)
    
    #Open output file.
    output_file = open(args.output_file, 'w')
    
    #Initialize to empty list previous_variant_identifiers
    previous_variant_identifiers = [None,None,None,None]
    
    unique = []
    
    #We start to parse the file
    with open(args.input_psv, 'r') as input_file:
                
        for line in input_file: 

            #If header line.                        
            if line.startswith('#'):
                #We create our sample objects list
                sample_list = create_Samples_from_header(line)
                    
            else:
                split_line = line.split('|')
                variant_identifiers = [split_line[1],split_line[2],split_line[3],split_line[4]]
                
                #If line is ""new"" (we consider new a line with different Chr Pos Ref_Allele & Var_Allele values)
                if variant_identifiers != previous_variant_identifiers:
                    
                    skip_variant, samples_depth_true = check_samples_depth(line, sample_list, args.depth_threshold)
                    if skip_variant: continue
                    
                    categorize_variants(line, samples_depth_true, args.freq_1, args.freq_2)
                    
                    unique.append(line)                    
                    
                #Update current line identifier
                previous_variant_identifiers = variant_identifiers
                
        #Output variables
        total_unique_variants = len(unique)
        
        generate_report(sample_list, output_file)
        #'Total unique variants: {0}'.format(total_unique_variants)
        
        print 'DONE. Time elapsed: %.2f seconds.' %(time.time() - start_time) 