#!/usr/bin/env python
"""
SYNOPSIS
    prepare4genesys.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    prepare4genesys: Module that processes PipeTA analysis, so that they can be directy imported in GeneSysTM
EXAMPLES
    prepare4genesys.py -c project.json
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    prepare4genesys.py v0.1
"""
import pydevd
import sys
import os
import traceback
import optparse
import time
from essay import Essay
from essay import FileData

# Global parameters
module = "prepare4genesys"


# ****************************** main body *********************************
def main():
    global options, args

    # Read essay state json
    essay = Essay()
    essay.read_essay(options.essay_path)

    print 'Launching prepare4genesys jobs for essay', essay.get_name()

    # Create GeneSys compatible file
    annotation = FileData(name=essay.get_name() + '_annotation.vcf', data_type='analysis', modifier='-i ',
                          add_args='annotation', path=essay.get_path())
    genesys_psv = FileData(name=essay.get_name() + '_GeneSys.psv', data_type='trash', modifier='-o ',
                           add_args='prepare4genesys', path=essay.get_path())
    hold_jobs = essay.submit_from_essay(job_name='genesys', command='vcf2GeneSys.py', input_data=annotation,
                                        output_data=genesys_psv, module=module)

    # Export result files via jenkins task (cURL) (only if not exome)
    if essay.get_target_name() is not None and not (
            essay.get_target_name().startswith('exome') or essay.get_target_name().startswith('exoma')):
        if 'gluster2' in essay.get_path():
            param = '\'http://Bioinfo:411977dc3264743d7b8b4cb380978e14@10.0.0.82:8080/view/Bioinfo/job/Pipeta_Baming_G2/buildWithParameters?token=SOCELTOKENG2&path=' + essay.get_path() + '&project=' + essay.get_project_name() + '&analysis=' + essay.get_name() + '\''
        else:
            param = '\'http://Bioinfo:411977dc3264743d7b8b4cb380978e14@10.0.0.82:8080/view/Bioinfo/job/Pipeta_Baming_G1/buildWithParameters?token=SOCELTOKENG1&path=' + essay.get_path() + '&project=' + essay.get_project_name() + '&analysis=' + essay.get_name() + '\''
        hold_jobs = essay.submit_from_essay(job_name='curl2UGM', command='curl', add_args=param, module=module,
                                            hold_jobs=hold_jobs)

        # Call the vcf2DBNLVar script
        #hold_jobs.append(
        #    essay.submit_from_essay(job_name='vcf2DBNLVar', command='vcf2DBNLVar', input_data=annotation,
        #                            module=module))
    return hold_jobs


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        parser.add_option('-e', '--essay_path', action='store', type='string', dest="essay_path", default="",
                          help='Path to essay folder containing state.json file')
        parser.add_option('-q', '--queue', action='store', type='string', dest="queue", default="low",
                          help='SGE queue name (currently low, med, mem or high available)')
        parser.add_option('-p', '--priority', action='store', type='int', dest="priority", default=0,
                          help='SGE priority (value from 0 to 1024)')
        (options, args) = parser.parse_args()
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if not options.essay_path:
            parser.error('ERROR: missing arguments!! Please check usage')
        if options.verbose:
            print time.asctime()

        # MAIN and global
        sender = 'genetonic@sistemasgenomicos.com'
        receivers = [os.environ['USERMAIL']]
        main()

        if options.verbose:
            print time.asctime()
        if options.verbose:
            print 'INFO: Duration of script run:',
        if options.verbose:
            print str(time.time() - start_time) + ' seconds'
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