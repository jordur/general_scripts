#!/usr/bin/env python
"""
SYNOPSIS
	sg_mirmap_test [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
	Script for testing the miRmap library
EXAMPLES: 
	TODO Show some examples of how to use this script.
EXIT STATUS
	TODO: List exit codes
AUTHOR
	TODO: Arbol Rodriguez <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
	This script belongs to Sistemas Genomicos S.L.
VERSION
	sg_mirmap_test.py v0.1.0
"""

import sys, os, traceback, optparse
import time
import re
import mirmap
#from pexpect import run, spawn

def main ():

	global options, args

	# ****************************** main body *********************************
	print 'Script for testing the miRmap library'

	# Open files
	target = open(options.target, 'r')
	mirna = open(options.mirna, 'r')

	# Show and save files contents:
	print 'Target file:'
	seq_target = target.read()
	print seq_target
	print 'miRNA sequence:'
	seq_mirna = mirna.read()
	print seq_mirna

	# Perform different actions:
	mim = mirmap.mm(seq_target, seq_mirna)
	mim.find_potential_targets_with_seed(allowed_lengths=[6,7], allowed_gu_wobbles={6:0,7:0}, allowed_mismatches={6:0,7:0}, take_best=True)
	print 'mim.end_sites: {}' . format(mim.end_sites)                                    # Coordinate(s) (3' end) of the target site on the target sequence
	mim.eval_tgs_au(with_correction=False)			# TargetScan features manually evaluated with
	mim.eval_tgs_pairing3p(with_correction=False)	# a non-default parameter.
	mim.eval_tgs_position(with_correction=False)
	print 'mim.prob_binomial: {}' . format(mim.prob_binomial)		# mim's attribute: the feature is automatically computed
	print 'mim.report: {}' . format(mim.report())

	
if __name__ == '__main__':
	try:
		start_time = time.time()
		parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='sg_mirmap_test.py v0.1.0')
		parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
		parser.add_option ('-t', '--target', action='store', dest="target", help='Specifies the target sequence')
		parser.add_option ('-m', '--mirna', action='store', dest="mirna", help='Specifies the miRNA sequence')
		(options, args) = parser.parse_args()
		#if len(args) < 1:
		#    parser.error ('missing argument')
		if options.verbose: print time.asctime()
		main()
		if options.verbose: print time.asctime()
		if options.verbose: print 'Duration of script run:',
		if options.verbose: print (time.time() - start_time) / 60.0
		sys.exit(0)
	except KeyboardInterrupt, e: # Ctrl-C
		raise e
	except SystemExit, e: # sys.exit()
		raise e
	except Exception, e:
		print 'ERROR, UNEXPECTED EXCEPTION'
		print str(e)
		traceback.print_exc()
		os._exit(1)
