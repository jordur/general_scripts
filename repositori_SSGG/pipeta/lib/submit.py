#!/usr/bin/env python
"""
SYNOPSIS
    submit.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    submit.py: Module for submitting jobs to SGE (and in the future, for whatever other job engines)
EXAMPLES
    TODO
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    submit.py v0.1
"""
import drmaa


# Functions definitions
def submit_job(job_name=None, command=None, add_args=None, cwd=None, queue='shudra.q', hold_jobs=None):
    # Start SGE session
    s = drmaa.Session()
    s.initialize()
    jt = s.createJobTemplate()
    jt.jobName = job_name
    jt.remoteCommand = command
    if add_args is not None:
        jt.args = add_args
    jt.joinFiles = True
    jt.nativeSpecification = '-shell y -V -j y -S /bin/bash -q ' + queue
    jt.workingDirectory = cwd
    if hold_jobs is not None:
        jt.nativeSpecification += ' -hold_jid ' + hold_jobs
    jobid = s.runJob(jt)
    print 'INFO: Job ' + job_name + ' has been submitted with id ' + jobid
    s.deleteJobTemplate(jt)

    # Close SGE session
    s.exit()
    # Return SGE job id
    return jobid