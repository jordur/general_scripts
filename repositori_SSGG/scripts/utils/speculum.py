#!/usr/bin/env python

"""
SYNOPSIS
    TODO speculum [-h,--help] [--version]
DESCRIPTION
    Copies project to corresponding DONE directory.
    After, launches Jenkins backup task.

EXAMPLES
    TODO: Show some examples of how to use this script.
EXIT STATUS
    TODO: List exit codes
AUTHOR
    TODO: Guillermo Marco <guillermo.marco@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    speculum.py v0.0.1
"""

import sys
import os
import traceback
import argparse
import time
import pydevd
import shutil
import jenkins
from termcolor import cprint

matraz_config_file = os.environ['MATRAZ_CONFIG']
gluster = os.environ['PROJECTS']
gluster_2 = os.environ['PROJECTS2']
gluster_done = os.environ['PROJECTS_DONE']
gluster_2_done = os.environ['PROJECTS2_DONE']


def convert_unicode_to_str(json_input):
    if isinstance(json_input, dict):
        return {convert_unicode_to_str(key): convert_unicode_to_str(value) for key, value in json_input.iteritems()}
    elif isinstance(json_input, list):
        return [convert_unicode_to_str(element) for element in json_input]
    elif isinstance(json_input, unicode):
        return json_input.encode('utf-8')
    else:
        return json_input


def project_exists(project_code):
    gluster_project_path = os.path.join(gluster, project_code)
    gluster_done_project_path = os.path.join(gluster_done, project_code)

    gluster_2_project_path = os.path.join(gluster_2, project_code)
    gluster_2_done_project_path = os.path.join(gluster_2_done, project_code)

    if os.path.isdir(gluster_project_path):
        if os.path.exists(gluster_done_project_path):
            sys.exit('WARN: {0} already exists, aborting backup.'.format(gluster_done_project_path))
        return 1, gluster_project_path, gluster_done_project_path

    elif os.path.isdir(gluster_2_project_path):
        if os.path.exists(gluster_2_done_project_path):
            sys.exit('WARN: {0} already exists, aborting backup.'.format(gluster_2_done_project_path))
        return 2, gluster_2_project_path, gluster_2_done_project_path

    else:
        sys.exit('ERROR: {0} project not found in Gluster/Gluster2.'.format(project_code))


def copy_work_to_done(project_name, project_path, done_project_path):
    print 'INFO: Copying project: {0}'.format(project_name)
    print'INFO: Source: {0}'.format(project_path)
    print'INFO: Destination: {0}'.format(done_project_path)
    print 'INFO: Please be patient this can take a while..'
    shutil.copytree(project_path, done_project_path, symlinks=False, ignore=shutil.ignore_patterns('trash'))
    cprint('INFO: Successfully copied clean project from work to DONE.', 'green')
    return 1


def remove_trash_from_done(done_project_path):
    print 'INFO: Removing trash folders from DONE..'
    for root, dirs, files in os.walk(done_project_path):
        for subdir in dirs:
            if subdir == 'trash':
                shutil.rmtree(os.path.join(root, subdir))
    print 'DONE: Trash cleaned from project DONE copy.'


def remove_dir(path):
    try:
        shutil.rmtree(path)

    except (OSError, 13):
        cprint('WARNING: Couldn\'t delete source directory. Check your permissions.', 'red')
        cprint('WARNING: You must delete source {0} manually.'.format(path), 'red')


def jenkins_copy_to_backup(environment, project_name):
    print 'INFO: Launched Jenkins backup task to Mojito server.'
    job_name = 'Backup_de_DONE_Gluster{0}'.format(environment)
    jenkins_handler = jenkins.Jenkins('http://jenkins.sistemasgenomicos.com:8080/', 'Bioinfo', 'bi0inf02013')
    jenkins_handler.build_job(job_name, {'project': project_name},
                              'G{0}Tokendone'.format(environment))

    #job = jenkins_handler.get_job_name(job_name)
    #last_build = job.get_last_build()
    #
    #print job
    #print last_build

    #job_info = jenkins_handler.get_job_info(job_name)
    #print "a"


def print_logo():
    #http://patorjk.com/software/taag/#p=display&h=2&v=2&f=Small%20Slant&t=Speculum%20v0%20.1
    cprint('''   ____                  __
  / __/__  ___ ______ __/ /_ ____ _
 _\ \/ _ \/ -_) __/ // / / // /  ' \\
/___/ .__/\__/\__/\_,_/_/\_,_/_/_/_/
   /_/
   You can press CTRL+C to abort process.
   ''', 'cyan')


def main():
    global args
    try:
        #Print ASCII logo
        print_logo()

        #Check if BF directory project exists in Gluster/Gluster2
        project_env, project_path, done_project_path = project_exists(args.bf_project_name)

        #Copy files from work (Ongoing) to DONE
        copy_work_to_done(args.bf_project_name, project_path, done_project_path)

    except KeyboardInterrupt:  # Ctrl-C:
        cprint('WARNING: User keyboard interrupt detected. Aborting backup.', 'red')
        print 'INFO: Removing backup destination dir: {0}'.format(done_project_path)
        remove_dir(done_project_path)
        sys.exit(2)

    #Copy files to backup server with Jenkins
    jenkins_copy_to_backup(project_env, args.bf_project_name)

    print 'INFO: Deleting source dir: {0}'.format(project_path)
    remove_dir(project_path)

    cprint('DONE: Enjoy backup ;)\n', 'green')


if __name__ == '__main__':
    try:
        start_time = time.time()

        parser = argparse.ArgumentParser(
            description='Copies project to corresponding DONE directory,\
            removes source then launches Jenkins backup task.',
            version='alpha version')

        parser.add_argument('-bf', dest='bf_project_name', action='store',
                            help='Bioinformatics project directory name',
                            required=True)

        parser.add_argument('-d', dest='debug', action='store_true',
                            help="Remote debug action.", required=False)

        args = parser.parse_args()
        if args.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)

        main()

    except SystemExit, e:  # sys.exit()
        #print 'System unexpected exit detected.'
        raise e

    except Exception, e:
        print 'ERROR: Aborting backup !'
        print str(e)
        traceback.print_exc()
        sys.exit(1)