#!/usr/bin/env python

'''
Created on Sep 05, 2013
@author: gmarco
'''

#Python imports
import argparse
import os.path
import time
import datetime
import string
import random
from subprocess import call
import subprocess

def menu ():
    
    """Input menu for user."""
    
    parser = argparse.ArgumentParser(description='Reads on target bedtools v0.1',
                                     epilog="################################")

    required = parser.add_argument_group('Required options')
        
    required.add_argument('-i', dest='input_bam', action = 'store', 
                        help='BAM input file. ie: -i Sample_1234_Panel.bam',
                        required=True)
    
    required.add_argument('-p', dest='target_region', action = 'store', 
                        help='Target regions BED file. ie: -p target_regions.bed',
                        required=True)
    required.add_argument('-f', dest='trash_dir', action = 'store', 
                        help='Trash directory. ie: -f ./trash/',
                        required=True)
    
    required.add_argument('-o', dest='output', action = 'store',
                        help="Output results file. ie: -o ./trash/stats/on_target.stat.all",
                        required=True)
    
    return parser.parse_args()

def working_message():
    print 'Reads on target bedtools v0.1'
    print '=============================='
    print 'Working...'

if __name__ == '__main__':
      
    #Start timer counter.
    start_time = time.time()
    
    #Store menu options namespace.
    args = menu()
    
    #Parse menu options and normalize them.
    input_bam = args.input_bam
    target = args.target_region
    input_bam_basename = os.path.splitext(os.path.basename(input_bam))[0]
    trash_dir = os.path.normpath(args.trash_dir)
    
    #Open output file to write results.
    output = open(args.output, 'w')
    
    #Trash intersect bam.
    bam_on_target = os.path.join(trash_dir, input_bam_basename+'_on_target.bam')
    
    #Print working message
    working_message()
    
    #Open trash file and generate intersect BAM.
    with open(bam_on_target, 'w') as intersect_output:
        call(['intersectBed', '-abam', '{0}'.format(input_bam), '-b', '{0}'.format(target)],stdout=intersect_output)

    #Count total reads on intersect file.
    samtools_view_count = subprocess.Popen(['samtools', 'view', '-c', bam_on_target], stdout=subprocess.PIPE)
    stdout, stderr = samtools_view_count.communicate()
    reads_on_target = stdout.strip()
    
    #Write to output file.
    output.write('reads_on_target\t{0}'.format(reads_on_target))
    
    #Get elapsed time.
    time = int(time.time() - start_time)
    print ('DONE. Time elapsed: {0}'.format(datetime.timedelta(seconds=time)))
     