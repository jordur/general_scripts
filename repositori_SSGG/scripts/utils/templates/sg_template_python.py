#!/usr/bin/env python
"""
SYNOPSIS
    TODO sg_template_python [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    TODO This describes how to use this script. This docstring
    will be printed by the script if there is an error or
    if the user requests help (-h or --help).
EXAMPLES
    TODO: Show some examples of how to use this script.
EXIT STATUS
    TODO: List exit codes
AUTHOR
    TODO: Oscar Rodriguez <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    sg_template_python.py v0.1.0
"""

import sys
import os
import traceback
import optparse
import time
import pydevd


def main():
    global options, args
    # TODO: Do something more interesting here...
    print 'This is the python template script. Please use it for developing of new python scripts'


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='$Id$')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        (options, args) = parser.parse_args()
        #if len(args) < 1:
        #    parser.error ('missing argument')
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if options.verbose:
            print time.asctime()
        main()
        if options.verbose:
            print time.asctime()
            print 'Duration of script run:',
            print (time.time() - start_time) / 60.0
        sys.exit(0)
    except KeyboardInterrupt, e:  # Ctrl-C
        raise e
    except SystemExit, e:  # sys.exit()
        raise e
    except Exception, e:
        print 'ERROR, UNEXPECTED EXCEPTION'
        print str(e)
        traceback.print_exc()
        sys.exit(1)