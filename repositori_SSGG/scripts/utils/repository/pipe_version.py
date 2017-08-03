#!/usr/bin/env python
"""
SYNOPSIS
    pipe_version.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    pipe_version: script for obtaining pipeline version and repositories commit_IDs
EXAMPLES
    pipe_version.py -p DNA-seq
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    pipe_version v0.1
"""

import subprocess
import sys
import os
import traceback
import optparse
import time
import git
from datetime import datetime
import pydevd


class RepoData:
    """Common base class for repositories"""
    def __init__(self, name, path):
        self.name = name
        self.path = path
        # Get last update date
        try:
            self.repo = git.Repo(path)
            lastpull = self.repo.git.log('--simplify-by-decoration', pretty="%ai %d", max_count=1)
        except Exception, e:
            # The path is not a git repo:
            print 'INFO: Not a git-repository or not possible to get git information', name, '. Exception: '+str(e)
            print 'INFO: Getting version information from path name or by launching application without parameters'
            print "INFO: Launching " + self.name
            proc = subprocess.Popen(self.name, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
            output = proc.stdout.read()
            self.date = None
            self.lasttag = None
            self.version = 0
            for line in output.split('\n'):
                uppercaseline = line.upper()
                if 'VERSION' in uppercaseline:
                    # check separator:
                    separator = 'VERSION:' if 'VERSION:' in uppercaseline else 'VERSION'
                    tmp = uppercaseline.split(separator)[1].split('.')
                    version = tmp[0]+'.'+tmp[1].split('-')[0]
                    self.version = float(version)
                    return
                if self.name.upper() in uppercaseline:
                    if line.split()[1].startswith('v'):
                        try:
                            self.version = float(line.split()[1][1:])
                            return
                        except Exception, e:
                            # Not a version number
                            self.version = 0
                    else:
                        try:
                            self.version = float(line.split()[1])
                            return
                        except Exception, e:
                            # Not a version number
                            self.version = 0
            return

        # If the path is a git repo:
        data = lastpull.split()
        self.date = datetime.strptime(data[0]+" "+data[1], '%Y-%m-%d %H:%M:%S')
        # Get last tag version
        try:
            self.lasttag = self.repo.git.describe()
        except Exception, e:
            print "WARNING: No tag-information found for repository", name, ". Exception: "+str(e)
        data = self.lasttag.split('-')
        if data[0].startswith('v'):
            data[0] = float(data[0][1:])
        else:
            data[0] = float(data[0])
        self.version = data[0]

    def get_path(self):
        return self.path

    def display_repository(self):
        print "Repository name:", self.name, ", path:", self.path, ", version:", self.version, ", tag: ", self.lasttag,\
            ", date:", self.date

    def set_version(self, version):
        self.version = version

    def get_version(self):
        return self.version

    def set_date(self, date):
        self.date = date

    def get_date(self):
        return self.date


def main():
    global options, args

    # ****************************** main body *********************************
    if options.verbose:
        print '''pipe_version: script for checking pipeline versions...'''
    
    # Define existent repostiories
    scripts_path = os.environ['SCRIPTS']
    repos = {}
    repos['scripts'] = RepoData("scripts", scripts_path)
    repos['biolib'] = RepoData("biolib", scripts_path+"/biolib")
    repos['pipeta'] = RepoData("pipeta", scripts_path+"/utils/PipeTA")
    repos['panels'] = RepoData("panels", scripts_path+"/genomics/panels/panel_design/")

    # And applications
    apps_path = os.environ['APPS']
    repos['gatk'] = RepoData("gatk_toolkit", apps_path+'/GATK')
    repos['bwa'] = RepoData("bwa", apps_path+'/bwa')
    repos['VarScan'] = RepoData("VarScan", apps_path+'/VarScan')
    repos['samtools'] = RepoData("samtools", apps_path+'/samtools')

    pipe_version = 0
    date = datetime.fromordinal(1)
    for repo in repos.keys():
        # Get each repository last tag
        repos[repo].display_repository()
        if repos[repo].get_date() is not None and repos[repo].get_date() > date:
            date = repos[repo].get_date()
        pipe_version += repos[repo].get_version()
        
    print "Pipeline version: ", pipe_version, ", date: ", date
    return pipe_version

if __name__ == '__main__':
    try:
        start_time = time.time()
        parser = optparse.OptionParser(formatter=optparse.TitledHelpFormatter(), usage=globals()['__doc__'], version='0.1.0')
        parser.add_option('-v', '--verbose', action='store_true', default=False, help='verbose output')
        parser.add_option('-p', '--pipe', action='store', type='string', dest="pipe", default="", help='Pipeline name (DNA-reseq, ...)')
        parser.add_option('-d', '--debug', action='store_true', default=False, help='allows remote debugging')
        (options, args) = parser.parse_args()
        if options.debug:
            user_ip = os.environ['USERIP']
            pydevd.settrace(user_ip, port=58484, stdoutToServer=True, stderrToServer=True)
        if not options.pipe:
            parser.error('ERROR: missing arguments!! Please check usage')
        if options.verbose:
            print time.asctime()
        
        # MAIN and global
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