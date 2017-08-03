#!/usr/bin/env python

'''
Created on Aug 22, 2013
@author: Guillermo Marco Puche
@email: guillermo.marco@sistemasgenomicos.com
'''

#Python imports
import pysam
import argparse
import os.path
import sys
import datetime
import time

#self classes
from Variant import Variant

def menu ():
    
    """Input menu for user."""
    
    
    
    parser = argparse.ArgumentParser(description='Haloplex v0.1',
                                     epilog="################################")

    required = parser.add_argument_group("Required options")
        
    required.add_argument('-i', dest='input_psv', action = 'store', 
                        help='Input: *.psv annotation file. ie: /home/annotation.psv',
                        required=True)
    
    required.add_argument('-bam', dest='input_bam', action = 'store',
                        help='BAM File requieres that .bam.bai (bam index) is also in BAM directory.',
                        required=True)
    
    required.add_argument('-r', dest='position_range', type=int, action = 'store', 
                        help="Number of positions to get variants in bam respect the variant position in annotation. ie: -r 100"
                        ,required=True)
    
    required.add_argument('-3', dest='three_range', type=int, action = 'store', 
                        help="3' number of bases to look in range. ie -3 10",
                        required=True)
    
    required.add_argument('-5', dest='five_range', type=int, action = 'store', 
                        help="5' number of bases to look in range ie: -5 10",
                        required=True)
    
    required.add_argument('-o', dest='output_file', action = 'store', default = './variants_report.psv',
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


def strand_check(strand):
    if strand == 'forward' or strand == 'f': return False
    elif strand == 'reverse' or strand == 'r': return True
    else: raise ValueError('Incorrect strand argument.')
    
def working_message (args):
    """Returns message with the options that the
    user has chosen"""
    
    message = 'Haloplex v0.1 working...'
    message += '\n=============================='
    message += '\nInput file: %s' %args.input_psv
    message += '\nOutput: %s' %args.output_file
    message += '\n=============================='
    
    return message

def sam_stuff(variant, position_range, three_range, five_range, input_bam, output_file):
    samfile = pysam.Samfile(input_bam, "rb")
    
    chr = variant.get_chr()
    pos = variant.get_pos()
    start_pos = pos - position_range
    end_pos = pos + position_range


    forward_three = []
    forward_five = []
    forward_out = []
    
    reverse_three = []
    reverse_five = []
    reverse_out = []
    
    output_file.write(('{0}\t{1}\n').format(variant.get_chr(), variant.get_pos()))
    
    for alignedread in samfile.fetch(chr, start_pos, end_pos):
        
        if alignedread.pos <= pos <= alignedread.aend:
            
            #If reverse strand 3' <---- 5'
            if alignedread.is_reverse:
                three_start = alignedread.pos
                three_end = alignedread.pos+three_range
                three_interval = range(three_start, three_end+1)
                
                five_start = alignedread.aend-five_range
                five_end = alignedread.aend
                five_interval = range(five_start, five_end+1)
                
                if pos in three_interval:
                    variant.increment_reverse_in_three_range()
                    reverse_three.append(alignedread)
                    
                elif pos in five_interval: 
                    variant.increment_reverse_in_five_range()
                    reverse_five.append(alignedread)
                else:
                    variant.increment_total()
                    reverse_out.append(alignedread)
                
            #If forward strand 5' ----> 3'
            else:
                three_start = alignedread.aend-three_range
                three_end = alignedread.aend
                three_interval = range(three_start, three_end+1)
                
                five_start = alignedread.pos
                five_end = alignedread.pos+five_range
                five_interval = range(five_start, five_end+1)
                
                if pos in three_interval:
                    variant.increment_forward_in_three_range()
                    forward_three.append(alignedread)
                    
                elif pos in five_interval: 
                    variant.increment_forward_in_five_range()
                    forward_five.append(alignedread)
                else:
                    variant.increment_total()
                    forward_out.append(alignedread)

    output_file.write('FORWARD\n')
    output_file.write('Total 3\' 5\'\n')                    
    output_file.write(('{0}\t{1}\t{2}\n').format(variant.get_total(), variant.get_forward_in_three_range(), variant.get_forward_in_five_range()))
        
    output_file.write('3\':\n')
        
    for elemento in forward_three:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))
        
    output_file.write('5\':\n')
    
    for elemento in forward_five:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))

    output_file.write('OUT FORWARD:\n')
    for elemento in forward_out:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))
    

    output_file.write('REVERSE\n')
    output_file.write('Total 3\' 5\'\n')                    
    output_file.write(('{0}\t{1}\t{2}\n').format(variant.get_total(), variant.get_reverse_in_three_range(), variant.get_reverse_in_five_range()))
        
    output_file.write('3\':\n')
        
    for elemento in reverse_three:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))
        
    output_file.write('5\':\n')
    
    for elemento in reverse_five:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))

    output_file.write('OUT FORWARD:\n')
    for elemento in reverse_out:
        output_file.write(('{0}\t{1}\t{2}\n').format(elemento.qname, elemento.pos, elemento.aend))

                 
    return None
            
    #samfile.close()

def variant_reader (input_psv):
    
    variant_list = []
    
    #Initialize to empty list previous_variant_identifiers
    previous_variant_identifiers = [None,None,None,None]
    
    #We start to parse the file
    with open(input_psv, 'r') as input_file:
                
        for line in input_file: 

            #If header line.                        
            if line.startswith('#'):
                continue
                    
            else:
                split_line = line.split('|')
                variant_identifiers = [split_line[1],split_line[2],split_line[3],split_line[4]]
                
                #If line is ""new"" (we consider new a line with different Chr Pos Ref_Allele & Var_Allele values)
                if variant_identifiers != previous_variant_identifiers:
                    #add 1 to position since Ensembl uses 1 based position and pysam uses 0 based position
                    new_variant = Variant(split_line[1],int(split_line[2])+1,split_line[3],split_line[4])
                    variant_list.append(new_variant)                           
                    
                #Update current line identifier
                previous_variant_identifiers = variant_identifiers
    
    return variant_list

if __name__ == '__main__':
    
    #Start timer counter.
    start_time = time.time()
    
    #Parse user menu options.
    args = menu()
    
    #Check input psv annotation exists.
    file_exists(args.input_psv)
    
    #Check that input bam exists.
    file_exists(args.input_bam)
    
    #Print working message.
    print working_message (args)
    
    #Create output dirs if needed.
    output_mgmt(args.output_file)
    
    #Open output file.
    output_file = open(args.output_file, 'w')
    
    #Create Variant Objects
    variant_list = variant_reader(args.input_psv)

    #output_file.write('#Chr\tPosition\tTotal_Variants\t3_prime_range\t5_prime_range\n')
    for variant in variant_list:
        
        if variant.get_chr() != 13 and variant.get_pos() != 32913485:
            print variant.get_chr(), variant.get_pos()
            continue

        else:
            print 'ENTRO', variant.get_chr(), variant.get_pos()
        #print variant.get_chr(), variant.get_pos(), variant.get_total(), variant.get_in_three_range(), variant.get_in_five_range()
        sam_stuff(variant, args.position_range, args.three_range, args.five_range, args.input_bam, output_file)
        #print variant.get_chr(), variant.get_pos(), variant.get_total(), variant.get_in_three_range(), variant.get_in_five_range()
        #output_file.write(('{0}\t{1}\t{2}\t{3}\t{4}\n').format(variant.get_chr(), variant.get_pos(), variant.get_total(), variant.get_in_three_range(), variant.get_in_five_range()))
            
    time = int(time.time() - start_time)
    print ('DONE. Time elapsed: %s') %str(datetime.timedelta(seconds=time))