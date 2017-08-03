#!/usr/bin/env python
"""
SYNOPSIS
    PipeTA.py [-h,--help] [-v,--verbose] [--version] [-c,--config]
DESCRIPTION
    PipeTA: Pipeline To Analyze, a script for retrieving general configuration parameters for launching pipelines
EXAMPLES
    PipeTA.py
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    PipeTA.py v0.1.0
"""

import sys, os, traceback, optparse
import time
import re
import json
import os
#from pexpect import run, spawn

def main ():

    global options, args

    # ****************************** main body *********************************
    print 'PipeTA: Pipeline To Analyze, a script for retrieving general configuration parameters for launching pipelines'
    
    config = json.load(open(options.config))
    
    #print(os.environ)
    os.environ["capture"] = config["references"]["panels"][analysis]["capture"]
    os.environ["target"] = config["references"]["panels"][analysis]["target"]
    os.environ["involved_chrs"] = config["references"]["panels"][analysis]["involved_chrs"]
    
    
    print(os.environ)
    
if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='0.1.0')
        parser.add_option ('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option ('-c', '--config', action='store_true', dest="config", default="/share/apps/etc/pipelines/config.json", help='Json configuration file')
        (options, args) = parser.parse_args()
        if len(args) < 1:
            parser.error ('missing argument')
        else: analysis = args[0]
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
