#!/usr/bin/env python
"""
SYNOPSIS
    probeta.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    probeta: Probe To Analyze, a launcher, status-checker, controller and reports creator for in silico analysis
    This script follows the whole workflow of a given PipeTA pipeline, performing, among others, the following actions:
    - checks jobs states
    - updates these states in the essay state.json (together with errors information)
    - checks errors in both SGE and logs files
    - creates a performance report
    - sends e-mail notifications when problems in analysis appear, when analysis finish, etc.
    - allows relaunching whole modules (simply set module state to 'relaunch')
    Please check its new information in the essay state.json file and in:
    - error.log: error information from both SGE and logs files
    - performance.tsv: file that gathers computing information (duration, CPU use, slots, memory use, etc.) for all the jobs launched
    Probeta can be used even on already finished projects with minor (or none) changes, depending on the version of PipeTA used in the analysis.
    Probeta is able to call PipeTA when needed, to relaunch analysis 
    
EXAMPLES
    probeta.py -c project.json
    probeta.py -r -v -c project.json # This option relaunches all jobs that had errors
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    probeta.py v0.1.0
"""

import sys
import os
import traceback
import optparse
import time
import json
import shlex
import subprocess
import drmaa
import smtplib
import pydevd


def DecodeStatus(state):
    """ Function for decoding SGE status """
    return {
        drmaa.JobState.UNDETERMINED: 'process status cannot be determined',
        drmaa.JobState.QUEUED_ACTIVE: 'job is queued and active',
        drmaa.JobState.SYSTEM_ON_HOLD: 'job is queued and in system hold',
        drmaa.JobState.USER_ON_HOLD: 'job is queued and in user hold',
        drmaa.JobState.USER_SYSTEM_ON_HOLD: 'job is queued and in user and system hold',
        drmaa.JobState.RUNNING: 'job is running',
        drmaa.JobState.SYSTEM_SUSPENDED: 'job is system suspended',
        drmaa.JobState.USER_SUSPENDED: 'job is user suspended',
        drmaa.JobState.DONE: 'job finished normally',
        drmaa.JobState.FAILED: 'job finished, but failed',
    }


def HoldState(job):
    """ Function that returns the hold state (if waiting for other jobs, for example) of a given job
    This function considers that if errors in previous jobs occur, current job shouldn't be launched"""


def get_dependent_jobs(jobid, essay, s):
    jobList = []
    for module in essay['modules'].keys():
        for job in essay['modules'][module]['jobs'].keys():
            if str(jobid) != str(essay['modules'][module]['jobs'][job]):
                for myjobid in essay['modules'][module]['jobs'][job]['hold_jobs']:
                    if str(myjobid) == str(jobid):
                        jobList.append(str(essay['modules'][module]['jobs'][job]['jobid']))
                        if essay['modules'][module]['jobs'][job]['state'] != "process":
                            if essay['modules'][module]['jobs'][job]['SGE_id'] != 'unknown' and \
                                            essay['modules'][module]['jobs'][job]['SGE_id'] != "":
                                if options.verbose:
                                    print "INFO: holding job " + job + " (" + str(essay['modules'][module]['jobs'][job][
                                        'SGE_id']) + "), dependent of job " + jobid
                                holded = hold_job(essay['modules'][module]['jobs'][job]['jobid'], essay, s)
                                if holded != 0:
                                    if options.verbose:
                                        print "INFO: hold not possible, relaunching " + job
                                    set_job_for_relaunch(essay['modules'][module]['jobs'][job]['jobid'], essay, s)
                                else:
                                    if options.verbose:
                                        print "INFO: holded job" + job
                        else:
                            print "INFO: Impossible to hold process " + job + ". This is a process, NOT a SGE job!"
    return jobList


def set_job_for_relaunch(jobid, essay, s):
    try:
        module, job = get_module_and_job_from_job_id(jobid, essay)
        # Set new state to pending
        essay['modules'][module]['jobs'][job]['state'] = 'pending'
        if not essay['modules'][module]['jobs'][job].has_key('error_history'):
            essay['modules'][module]['jobs'][job]['error_history'] = {}
        SGE_id = essay['modules'][module]['jobs'][job]['SGE_id']
        essay['modules'][module]['jobs'][job]['SGE_id'] = ''
        essay['modules'][module]['jobs'][job]['error_history'][SGE_id] = {}
        # SGE errors
        if essay['modules'][module]['jobs'][job]['performance'].has_key('exit_status'):
            essay['modules'][module]['jobs'][job]['error_history'][SGE_id]['exit_status'] = \
                essay['modules'][module]['jobs'][job]['performance']['exit_status']
        if essay['modules'][module]['jobs'][job]['performance'].has_key('failed'):
            essay['modules'][module]['jobs'][job]['error_history'][SGE_id]['failed'] = \
                essay['modules'][module]['jobs'][job]['performance']['failed']
        if essay['modules'][module]['jobs'][job]['performance'].has_key('hostname'):
            essay['modules'][module]['jobs'][job]['error_history'][SGE_id]['hostname'] = \
                essay['modules'][module]['jobs'][job]['performance']['hostname']
        essay['modules'][module]['jobs'][job]['performance'] = {}
        # Log errors
        essay['modules'][module]['jobs'][job]['error_history'][SGE_id]['errors'] = \
            essay['modules'][module]['jobs'][job]['errors']
        essay['modules'][module]['jobs'][job]['errors'] = ""
        try:
            # Hold/relaunch dependent Jobs
            jobList = get_dependent_jobs(jobid, essay, s)
            return True
        except Exception, e:
            print "WARNING: dependent jobs of " + str(jobid) + " could not be relaunched: " + str(e)
            return False
    except Exception, e:
        print "WARNING: job with jobid " + str(jobid) + " couldn't be set to be relaunched"
        return False


def get_module_and_job_from_job_id(jobid, essay):
    for module in essay['modules'].keys():
        for job in essay['modules'][module]['jobs'].keys():
            if jobid == essay['modules'][module]['jobs'][job]['jobid']:
                return (module, job)
    return ("", "")


def hold_job(jobid, essay, s):
    module, job = get_module_and_job_from_job_id(jobid, essay)
    try:
        s.control(essay['modules'][module]['jobs'][job]['SGE_id'], drmaa.JobControlAction.HOLD)
        return 0
    except Exception, e:
        print 'ERROR: Impossible to hold process ' + str(essay['modules'][module]['jobs'][job]['SGE_id'])
        print 'ERROR: ' + str(e)
        return 1


def unhold_job(jobid, essay, s):
    module, job = get_module_and_job_from_job_id(jobid, essay)
    try:
        s.control(essay['modules'][module]['jobs'][job]['SGE_id'], drmaa.JobControlAction.RELEASE)
        return 0
    except Exception, e:
        print 'ERROR: Impossible to unhold process ' + essay['modules'][module]['jobs'][job]['SGE_id']
        print 'ERROR: ' + str(e)
        return 1


def pipeta(config, queue, priority, additional_args):
    """Launch pipeta.pl with arguments"""
    try:
        command_line = "pipeta.pl -m " + config + " -q " + queue + " -p " + priority + " " + additional_args
        print "INFO: Launching " + command_line
        pipeta_args = shlex.split(command_line)
        retcode = subprocess.check_output(command_line, stderr=subprocess.STDOUT, shell=True)
        print retcode
    except Exception, e:
        print "ERROR: Problems launching " + command_line
        quit()


def send_mail(subject, body):
    global sender, receivers
    message = "From: your lovely GeneTonic computation cluster <genetonic@sistemasgenomicos.com>\n"
    message += "To: " + str(receivers) + "\n"
    message += "Subject: PROBETA REPORT - " + subject + "\n"
    message += body
    try:
        smtpObj = smtplib.SMTP('localhost')
        smtpObj.sendmail(sender, receivers, message)
        print "INFO: Successfully sent email"
    except SMTPException:
        print "ERROR: unable to send email"


def get_log_folder(essay, module, job):
    # For compatibility with older versions of pipeta, the log folder can be also obtained from job name
    if essay['modules'][module]['jobs'][job].has_key('log'):
        log_folder = essay['modules'][module]['jobs'][job]['log']
    else:
        log_folder = job
        log_folder.replace('jobs', 'logs', 1)
        tmp_str = log_folder.split('/')
        log_folder = ""
        for name in tmp_str[:-1]:
            log_folder = log_folder + name + '/'
    return log_folder


def check_error_in_line(line):
    check_words = [{'in': 'ERROR', 'not_in': []},
                   {'in': 'FAULT', 'not_in': ['DEFAULT']},
                   {'in': 'FATAL', 'not_in': []},
                   {'in': 'UNINITIALIZED', 'not_in': []},
                   {'in': 'FAIL', 'not_in': ['FAILING BADMATEFILTER', 'FAILING BADCIGARFILTER',
                                             'WERE FAILED BY THE STRAND-FILTER', 'FAILING DUPLICATEREADFILTER',
                                             'FAILING FAILSVENDORQUALITYCHECKFILTER', 'FAILING MALFORMEDREADFILTER',
                                             'FAILING MAPPINGQUALITYUNAVAILABLEFILTER',
                                             'FAILING MAPPINGQUALITYZEROFILTER', 'FAILING NOTPRIMARYALIGNMENTFILTER',
                                             'FAILING PLATFORM454FILTER', 'FAILING UNMAPPEDREADFILTER']}]
    uppercase_line = line.upper()
    error_in_line = False
    for word_check in check_words:
        if word_check['in'] in uppercase_line:
            not_in_check = True
            for not_in in word_check['not_in']:
                if not_in in uppercase_line:
                    not_in_check = False
            if not_in_check:
                error_in_line = True
    return error_in_line


def main():
    global options, args

    # ****************************** main body *********************************
    print '''ProbeTA: Probe To Analyze, a launcher and status-checker for in silico analysis... v1.0'''

    # Read configuration json
    try:
        config = json.load(open(options.config))
    except ValueError:
        print "ERROR: Impossible to read json configuration file " + options.config

    # Launch pipeta.pl if it hasn't already been launched (state.json file of any of the analysis doesn't exist)
    stateFound = True
    for essay in config['essays'].keys():
        if not os.path.isfile(config['path'] + "/" + essay + "/state.json"):
            stateFound = False
    if not stateFound and not options.getdata:
        pipeta(options.config, options.queue, str(options.priority), "")
    elif options.getdata:
        # pipeta is launched in create_only_mode!! (just files & folders will be created)
        pipeta(options.config, options.queue, str(options.priority), "-c")

    # Start SGE session
    s = drmaa.Session()
    s.initialize()

    # Check which SGE jobs are running or queued_active
    finished = False
    relaunch = False
    subject = "project: " + config['name']
    body = subject + "\n"
    while not finished:
        finished = True
        changes = False
        # Read state.json files from all essays
        print "INFO: reading state.json"
        subject = "project: " + config['name']
        for essay in config['essays'].keys():
            subject = subject + ", essay: " + essay
            body = subject + "\n"
            with open(config['path'] + "/" + essay + "/state.json", 'r') as input:
                state = json.load(input)

            # Get all jobs_ids from essay
            if options.verbose:
                print "INFO: getting job ids"
            for module in state['modules'].keys():
                if state['modules'][module]['state'] == 'relaunch':
                    moduleRelaunch = True
                    state['modules'][module]['state'] = 'pending'
                else:
                    moduleRelaunch = False
                for job in state['modules'][module]['jobs'].keys():

                    # Relaunch job if the whole module must be relaunched
                    if moduleRelaunch:
                        relaunchJob = True
                    else:
                        relaunchJob = False
                    if options.verbose:
                        print "INFO: starting with job " + job
                        # For compatibility with older pipeta versions, if SGE_id doesn't exit, this identificator is got from state
                    if state['modules'][module]['jobs'][job].has_key('SGE_id'):
                        id = state['modules'][module]['jobs'][job]['SGE_id']
                    else:
                        id = state['modules'][module]['jobs'][job]['state'].split('=')[1]

                    # If "getdata" mode is set, look for existent jobIDs from log files and get the biggest one:
                    if options.getdata:
                        logFolder = get_log_folder(state, module, job)
                        id = 'unknown'
                        jobname = job.split('/')[-1]
                        if options.verbose:
                            print "INFO: processing job name: " + jobname
                        for file in os.listdir(logFolder):
                            if jobname in file.split('/')[-1] and '.o' in file:
                                if id == 'unknown':
                                    id = file.split('.o')[-1]
                                elif file.split('.o')[-1] > id:
                                    id = file.split('.o')[-1]
                        if options.verbose:
                            print "INFO: for module " + module + ", job " + job + ", following id is found: " + id
                        state['modules'][module]['jobs'][job]['SGE_id'] = id

                    # Process only jobs submitted directly to SGE
                    if state['modules'][module]['jobs'][job]['state'][:7] != "process" and \
                                    state['modules'][module]['jobs'][job]['state'][:8] != "proccess":

                        # Check job status. If no status is returned, it's assumed that job is finished
                        try:
                            jobState = s.jobStatus(id)
                        except Exception, e:
                            jobState = "finished"

                        # Actions for running and hold jobs:
                        if jobState == "queued_active" or jobState == "running" or jobState == "system_on_hold" or jobState == "user_system_on_hold" or jobState == "user_on_hold" or jobState == "system_suspended":
                            finished = False
                            if options.verbose:
                                print "INFO: actions for running and holding job: " + job + " (" + id + ") with state: " + jobState
                                # If not running or active, set a user hold, for controlling purposes
                            if jobState != "queued_active" and jobState != "running" and jobState != "finished" and jobState != "system_suspended":
                                s.control(id, drmaa.JobControlAction.HOLD)
                                # If currently only hold by user, the job is released
                            if jobState == "user_on_hold":
                                s.control(id, drmaa.JobControlAction.RELEASE)

                        # Initialize json structure parameters in case that are not yet existent
                        # (for compatibility with older pipeta releases)
                        if not state['modules'][module]['jobs'][job].has_key('errors'):
                            state['modules'][module]['jobs'][job]['errors'] = ''
                        if not state['modules'][module]['jobs'][job].has_key('performance'):
                            state['modules'][module]['jobs'][job]['performance'] = {}

                        # Actions for finished jobs:
                        if jobState == "finished":
                            if options.verbose:
                                print "INFO: actions for finished job " + job + " (" + id + ")"

                            if state['modules'][module]['jobs'][job]['performance'] == {}:
                                changes = True
                                if options.verbose:
                                    print "INFO: getting SGE info for finished job " + job + " (" + id + ")"
                                    # Get SGE job info
                                tries = 0
                                info_available = False
                                if id != '' and id != 'unknown':
                                    while not info_available and tries < 3 and id != 'unknown':
                                        try:
                                            command_line = "qacct -j " + id
                                            # s.jobStatus(id) #This approach cannot be used since the job hasn't been
                                            # created in the drmaa session
                                            retcode = subprocess.check_output(command_line, stderr=subprocess.STDOUT,
                                                                              shell=True)
                                            info_available = True
                                        except Exception, e:
                                            if options.verbose:
                                                print "WARNING: No SGE information (via qacct) for job " + job + " (" + \
                                                      id + "). Exception: " + str(e)
                                            tries += 1
                                            time.sleep(10)
                                            pass

                                if info_available:
                                    # Transform the qacct information to a dictionary that can later be stored
                                    l = retcode.split('\n')
                                    l.pop(0)
                                    for param in l:
                                        l2 = param.split(' ', 1)
                                        try:
                                            l2[1]
                                        except Exception, e:
                                            break
                                        else:
                                            key = l2[0]
                                            value = l2[1].strip()
                                            state['modules'][module]['jobs'][job]['performance'][key] = value
                                else:
                                    # If no information from qacct (SGE) could be retrieved, it's assumed that the job
                                    # must be relaunched
                                    finished = False
                                    changes = True
                                    relaunchJob = True
                                    print 'WARNING:', job, ', with ID=', id, \
                                        ' has not SGE information available. It will be relaunched'

                            # SGE bad exit status
                            if options.verbose:
                                print "INFO: checking SGE bad exit status for job " + job + " (" + id + ")"
                                # There may exist jobs killed before running, and thus not having any qacct information
                            if state['modules'][module]['jobs'][job]['performance'].has_key('exit_status'):
                                exitstatus = state['modules'][module]['jobs'][job]['performance']['exit_status']
                            else:
                                exitstatus = "666: killed before running"
                                state['modules'][module]['jobs'][job]['performance']['exit_status'] = exitstatus
                                body = body + "\nProblem with job " + job + " (" + id + "). Exit status: " + exitstatus
                            if state['modules'][module]['jobs'][job]['performance'].has_key('failed'):
                                failed = state['modules'][module]['jobs'][job]['performance']['failed']
                            else:
                                failed = "666: killed before running"
                                state['modules'][module]['jobs'][job]['performance']['failed'] = failed
                                body = body + "\nProblem with job " + job + " (" + id + "). Failed status: " + failed
                            if exitstatus != "0" or failed != "0":
                                changes = True
                                finished = False

                                # 'relaunchJob' variable is set to true, so that pipeta gets launched again
                                relaunchJob = True
                                print 'WARNING: ' + job + ', with ID=' + id + ' did not finished properly. Exit_status=' + exitstatus + ', failed=' + failed

                            # Errors from LOG FILES are gathered
                            if options.verbose:
                                print "INFO: checking log file for job " + job + " (" + id + ")"
                                # Open Log-File
                            logFolder = get_log_folder(state, module, job)
                            for file in os.listdir(logFolder):
                                if file.endswith('.o' + id):
                                    logfile = open(state['modules'][module]['jobs'][job]['log'] + '/' + file, 'r')
                                    for line in logfile.readlines():
                                        if check_error_in_line(line):
                                            changes = True
                                            state['modules'][module]['jobs'][job]['errors'] = line
                                            body = body + "\nProblem with job " + job + " (" + id + "). ERROR: " + line

                                            # In case that jobs with ERRORS must be relaunched
                                            if options.relaunch:
                                                finished = False
                                                relaunchJob = True
                                            pass
                                        else:
                                            # For compatibility with previous probeta versions
                                            state['modules'][module]['jobs'][job]['errors'] = ""
                                    logfile.close()

                        # Actions for jobs whose state changed
                        if state['modules'][module]['jobs'][job]['state'] != jobState and not relaunch:
                            changes = True
                            if options.verbose:
                                print "INFO: job " + job + " (" + id + ") changed its stored state = " + \
                                      state['modules'][module]['jobs'][job][
                                          'state'] + ", to current state = " + jobState
                                # Update job state only for those jobs not going to be relaunched
                            if not relaunch:
                                if options.verbose:
                                    print "INFO: updating job state for job " + job + " (" + id + ")"
                                state['modules'][module]['jobs'][job]['state'] = jobState

                        # Relaunch job if necessary (states, SGE_id, etc is modified)
                        if relaunchJob and not options.getdata:
                            if options.verbose:
                                print "INFO: relaunching job " + job + " (" + id + ") and its dependent jobs"
                                # set_job_for_relaunch function relaunches current job and sets hold states or relaunches
                                # dependent jobs
                            relaunch = set_job_for_relaunch(state['modules'][module]['jobs'][job]['jobid'], state, s)

            # If there are changes in states, they get stored in state.json:
            errorsfound = ""
            if changes:
                # Update performance and ERRORS information
                if options.verbose:
                    print "INFO: changes found and therefore updating performance.tsv"
                with open(config['path'] + "/" + essay + "/performance.tsv", 'w') as performance:

                    # Open error logs file
                    try:
                        errorlog = open(config['path'] + "/" + essay + "/error.log", 'w')
                    except Exception, e:
                        print "ERROR: problem opening error logs file"
                        quit()

                    performance.write(
                        "module\tfile\tjob\telapsed time (secs)\tcpu time (secs*cpu)\t mem time (secs*GB)\tmax virtual memory (B)\tslots (un.)\tmean cpu usage (cpu)\tmean memory usage (GB)\n")
                    for module in state['modules'].keys():
                        for job in state['modules'][module]['jobs'].keys():
                            if state['modules'][module]['jobs'][job].has_key('performance') and \
                                            state['modules'][module]['jobs'][job]['performance'] != {}:
                                if state['modules'][module]['jobs'][job]['performance'][
                                    'failed'] != '666: killed before running':
                                    jobname = state['modules'][module]['jobs'][job]['performance']['jobname']
                                    jobnumber = state['modules'][module]['jobs'][job]['performance']['jobnumber']
                                    wallclock = state['modules'][module]['jobs'][job]['performance']['ru_wallclock']
                                    cputime = state['modules'][module]['jobs'][job]['performance']['cpu']
                                    memtime = state['modules'][module]['jobs'][job]['performance']['mem']
                                    maxvmem = state['modules'][module]['jobs'][job]['performance']['maxvmem']
                                    slots = state['modules'][module]['jobs'][job]['performance']['slots']
                                    if float(wallclock) != 0:
                                        meancpu = str(float(cputime) / float(wallclock))
                                        meanmem = str(float(memtime) / float(wallclock))
                                    else:
                                        meancpu = '0'
                                        meanmem = '0'
                                    performance.write(
                                        module + '\t' + jobname + '\t' + jobnumber + '\t' + wallclock + '\t' + cputime + '\t' + memtime + '\t' + maxvmem + '\t' + slots + '\t' + meancpu + '\t' + meanmem + '\n')
                                    if options.verbose:
                                        print "INFO: Outputting errors for job " + job + " (" + jobnumber + ")"
                                        # Output ERRORS:
                                    logserrors = state['modules'][module]['jobs'][job]['errors']
                                    if logserrors != "":
                                        errorsfound = 'ERROR: in job ' + job + ' (SGE ID ' + jobnumber + '): ' + logserrors
                                        errorlog.write(
                                            'ERROR: in job ' + job + ' (SGE ID ' + jobnumber + '): ' + logserrors)

                    # Close error log file
                    errorlog.close()

                # Update state.json
                if options.verbose:
                    print "INFO: changes found and therefore updating state.json"
                with open(config['path'] + "/" + essay + "/state.json", 'w') as outfile:
                    json.dump(state, outfile, sort_keys=True, indent=4, separators=(',', ': '))
            else:
                print "INFO: No changes..."

            # Relaunch PipeTA if necessary (SGE errors detected):
            if relaunch:
                if options.verbose:
                    print "WARNING: problems with finished jobs require PipeTA to be relaunched"
                if options.mail:
                    send_mail('WARNING: problems with analysis ' + config['name'],
                             'There were problems with your analysis ' + config[
                                 'name'] + '. Some jobs have been relaunched. Check it in ' + config[
                                 'path'] + '/' + essay + '\n' + body)
                pipeta(options.config, options.queue, str(options.priority), "")
                relaunch = False

            # Send mail in case of errors:
            if options.mail:
                if errorsfound != "":
                    send_mail('ERROR: errors appeared in analysis ' + config['name'],
                             'There were problems with your analysis ' + config[
                                 'name'] + '. Some jobs have been relaunched. Check it in ' + config[
                                 'path'] + '/' + essay + ':\n' + errorsfound + '\n' + body)

        if not finished:
            time.sleep(60)

    # Close SGE session                        
    s.exit()

    print "INFO: probeta successfully finished. Please check error.log and performance.tsv of each essay included in the project"
    if options.mail:
        send_mail(subject + ' :: Your analysis finished successfully!!',
                 body + '\n\nPlease feed me with additional computational problems. I am bored...\n')


if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'],
                                       version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-c', '--config', action='store', type='string', dest="config", default="",
                          help='Json configuration file')
        parser.add_option('-q', '--queue', action='store', type='string', dest="queue", default="low",
                          help='SGE queue name (currently low, med, mem or high available)')
        parser.add_option('-p', '--priority', action='store', type='int', dest="priority", default=0,
                          help='SGE priority (value from 0 to 1024)')
        parser.add_option('-r', '--relaunch', action='store_true', default=False,
                          help='option for re-launching jobs with errors in logs')
        parser.add_option('-m', '--mail', action='store_true', default=False,
                          help='option for sending mail when finished or errors appeared')
        parser.add_option('-g', '--getdata', action='store_true', default=False,
                          help='option for only getting information (jobs IDs, etc.) from a previously launched analysis')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        (options, args) = parser.parse_args()
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if not options.config:
            parser.error('ERROR: missing arguments!! Please check usage')
        if options.verbose: print time.asctime()

        # MAIN and global
        sender = 'genetonic@sistemasgenomicos.com'
        receivers = [os.environ['USERMAIL']]
        main()

        if options.verbose: print time.asctime()
        if options.verbose: print 'INFO: Duration of script run:',
        if options.verbose: print str(time.time() - start_time) + ' seconds'
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