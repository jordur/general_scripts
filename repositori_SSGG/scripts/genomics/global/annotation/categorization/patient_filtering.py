#!/usr/bin/env python

'''
Created on Aug 13, 2013
@author: SGNET\gmarco
'''

import argparse
import os.path
import sys
import time
from Sample import Sample

def menu ():
    
    """Input menu for user."""
    
    
    
    parser = argparse.ArgumentParser(description='Patient filtering v0.1',
                                     epilog="################################")

    required = parser.add_argument_group("Required options")
        
    required.add_argument('-i', dest='input_psv', action = 'store', 
                        help='Input: *.psv annotation file. ie: /home/annotation.psv',
                        required=True)
    
    required.add_argument('-cn', dest='cn_samples', action = 'store', nargs='+',
                        help='CN Samples list. ie: -cn Sample_1 Sample_2 Sample3',
                        required=True)
    
    required.add_argument('-cp', dest='cp_samples', action = 'store', nargs='+',
                        help='CP Samples list. ie: -cp Sample_4 Sample_5 Sample6',
                        required=True)
    
    required.add_argument('-d', dest='depth_value', type=float, action = 'store', 
                        help="Depth threshold. ie: -d 10"
                        ,required=True)
    
    required.add_argument('-fref', dest='freq_ref', type=float, action = 'store', 
                        help="Freq ref threshold. ie: -fref 0.12",  
                        required=True)
    
    required.add_argument('-fvar', dest='freq_var', type=float, action = 'store', 
                        help="Freq var threshold. ie: -fvar 0.75",
                        required=True)
    
    required.add_argument('-o', dest='output_file', action = 'store', default = './patient_filtered.psv',
                        help="Output file. If path to file doesn't exists it will be automatically created. ie: -output ./results/filtered.psv",
                        required=False)
    
    return parser.parse_args()

def create_Samples_from_header(line, sample_list):
    
    """Creates sample objects for input file.
    It uses CN and CP input lists from users.
    
    It creates a sample object for each one containing the name
    and the column index for genotype, depth amd var/depth (ratio) columns    
    """

    sample_object_list = []
    geno_sufix = '.replicate1_Genotype'
    depth_sufix = '.replicate1_Depth'
    ratio_sufix = '.replicate1_Var/Depth'
    
    columns = line.split('|')
    #if not set(sample_list).issubset(set(columns)): quit("ERROR: Specified samples: %s couldn't be found in input annotation file." %sample_list)
    for sample_name in sample_list:
        try:
            genotype_col = columns.index(sample_name+geno_sufix) 
            depth_col = columns.index(sample_name+depth_sufix)
            ratio_col = columns.index(sample_name+ratio_sufix)
            
            new_sample = Sample(sample_name, genotype_col, depth_col, ratio_col)
            sample_object_list.append(new_sample)
            
        except ValueError:
            print "ERROR: Specified samples: %s couldn't be found in input annotation file." %sample_list
            quit()
            
    return sample_object_list

def check_samples_depth (line, cp_samples, cn_samples, depth_threshold):
    
    """Checks sample depths are above depth_threshold
    Returns true to skip line if there's at least one sample
    above depth threshold.
    
    Returns False (to continue program flow) and also:
    
    - cp_depth_true list with cp sample objects that are above or equal depth threshold.
    - cn_informative list with cn sample objects that are above or equal depth threshold.
    - cn_non_informative list with cn sample objects that are below to depth threshold.
    """
    cp_depth_true = []
    cn_informative = []
    cn_non_informative = []
    
    columns = line.split('|')

    for cp_sample in cp_samples:
        #We store cp_samples that are greater or equal to depth_threshold
        if columns[cp_sample.get_depth_col()] == "-": columns[cp_sample.get_depth_col()] = 0.0
        if float(columns[cp_sample.get_depth_col()]) >= depth_threshold: cp_depth_true.append(cp_sample) 
                
    for cn_sample in cn_samples:
        #We categorize cn_samples that are greater or equal to depth_threshold as informative
        if columns[cn_sample.get_depth_col()] == "-": columns[cn_sample.get_depth_col()] = 0.0
        if float(columns[cn_sample.get_depth_col()]) >= depth_threshold: cn_informative.append(cn_sample)
        #We categorize cn_samples that are strictly minor to depth_threshold as non informative
        else: cn_non_informative.append(cn_sample)
    
    #Return True to skip_line if no cp_depth_true samples
    if len(cp_depth_true) == 0: return True, None, None, None
    
    else:
        #Return False to not skip line and continue main execution and upcomming modules
        return False, cp_depth_true, cn_informative, cn_non_informative

def check_cp_depth_true_var_depth_ratio(line, cp_depth_true, freq_ref_threshold):
    
    """Check for cp samples with enough depth that the freq_ref is above threshold
    or return True to skip the line.
    """
    columns = line.split('|')
    
    for cp_sample in cp_depth_true:
        #Check for each sample in cp_depth_true if var/depth ratio is lesser than freq_ref_threshold
        #Return true to skip_line to skip line
        if float(columns[cp_sample.get_ratio_col()]) < freq_ref_threshold: return True 

    #Do not skip line        
    return False

def check_cn_informative_var_depth_ratio(line, cn_informative, freq_ref_threshold):
    
    """Check for cn_informative samples with enough depth that the freq_ref is above threshold
    or return True to skip the line.
    """
    
    columns = line.split('|')
    
    for cn_informative_sample in cn_informative:
        #Check for each sample in cn_informative if var/depth ratio is greater equal than freq_ref_threshold
        #Return true to skip_line to skip line
        if float(columns[cn_informative_sample.get_ratio_col()]) >= freq_ref_threshold: return True 

    #Do not skip line
    return False

def file_exists (input_file):
    
    """Checks if input_file exists"""
    
    if not os.path.exists(input_file): raise IOError("The file %s does not exist!" %input_file)

def output_mgmt_old (output_file_path):
    
    """Creates output directories, even if nested in
    case that they don't exist."""
    
    output_base_path = os.path.dirname (output_file_path)
    
    if not os.path.exists(output_base_path):
        print "Creating output directory..."
        os.makedirs(output_base_path)
        
    return 1

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
    
    message = 'Patient filtering v0.1 working...'
    message += '\n=============================='
    message += '\nInput file: %s' %args.input_psv
    message += '\nCN_Samples: %s' %args.cn_samples
    message += '\nCP_Samples: %s' %args.cp_samples
    message += '\nDepth Threshold: %s' %args.depth_value
    message += '\nFreq ref threshold: %s' %args.freq_ref
    message += '\nFreq var threshold: %s' %args.freq_var
    message += '\nOutput: %s' %args.output_file
    message += '\n=============================='
    
    return message

def percentage_calc(read, total):
    
    """Changes input types to float and returns percentage."""
    
    return ((float(input_file.tell())/float(total_file_size))*100)    

def status(completed, completed_ant):
    
    """Prints percentage completed and flushes stdout.
    Returns completed_ant to update counter in main."""
    
    if completed - completed_ant > 10:
        completed_ant = float(completed)
        sys.stdout.write("Working:%3d%% completed.\r" % completed)
        sys.stdout.flush()
        return completed_ant
    return completed_ant


if __name__ == '__main__':
    
    #Start timer counter.
    start_time = time.time()
    
    #Parse user menu options.
    args = menu()
    
    #Check input file exists
    file_exists(args.input_psv)
    
    #Print working message.
    print working_message (args)
    
    #Create output dirs if needed.
    output_mgmt(args.output_file)
    
    #Open output file.
    salida = open(args.output_file, 'w')

    #Get total size of input file.
    total_file_size = int(os.path.getsize(args.input_psv))
    
    with open(args.input_psv, 'r') as input_file:
        
        #Completed counters to 0.
        completed_ant = 0.0
        completed = 0.0
        
        for line in input_file: 

            #If header line.                        
            if line.startswith('#'):
                completed = percentage_calc(input_file.tell(),total_file_size)
                 
                cn_samples = create_Samples_from_header(line, args.cn_samples)
                cp_samples = create_Samples_from_header(line, args.cp_samples)
                salida.write(line)
            
            else:
                #Update counters for percentage completed display.
                completed = percentage_calc(input_file.tell(),total_file_size)
                completed_ant = status(completed,completed_ant)
                
                #Check depths or skip line.
                skip_line, cp_depth_true, cn_informative, cn_non_informative = check_samples_depth(line, cp_samples, cn_samples, args.depth_value)
                if skip_line: continue
                
                #Check var/depth ratio or skip line.
                skip_line = check_cp_depth_true_var_depth_ratio(line, cp_depth_true, args.freq_ref)
                if skip_line: continue
                
                #Check that number cn_informative samples is not zero.
                if len(cn_informative)>0:
                    skip_line = check_cn_informative_var_depth_ratio(line, cn_informative, args.freq_ref)
                    if skip_line: continue
                    
                    salida.write(line)
                
                #If all are non informative print line, we can't discard it.
                else: salida.write(line)

        print 'DONE. Time elapsed: %.2f seconds.' %(time.time() - start_time) 