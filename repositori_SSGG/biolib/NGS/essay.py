#!/usr/bin/env python
"""
SYNOPSIS
    essay.py [-h,--help] [-v,--verbose] [--version]
DESCRIPTION
    Module for essay control in python
EXAMPLES
    TODO:
EXIT STATUS
    TODO: List exit codes
AUTHOR
    Arbol <oscar.rodriguez@sistemasgenomicos.com>
LICENSE
    This script belongs to Sistemas Genomicos S.L.
VERSION
    essay.py v0.1
"""
import subprocess
import sys
import os
import traceback
import json
from submit import submit_job

data_types = ['analysis', 'trash']


# Class definitions
class FileData():
    """ Stores relevant information for any file involved in analysis
    """

    def __init__(self, name=None, data_type=None, modifier=None, path=None, add_args=None):
        self._name = name
        self._type = data_type if data_type in data_types else None
        self._modifier = modifier
        self._add_args = add_args
        add_args = add_args if isinstance(add_args, (list, tuple)) else [add_args]
        self._path = path + '/' + data_type + '/' + '/'.join(add_args) + '/' + name if path is not None else None

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, name):
        self._name = name

    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, data_type):
        self._type = data_type

    @property
    def modifier(self):
        if self._modifier is None:
            return self._name
        else:
            return self._modifier + self.path

    @modifier.setter
    def modifier(self, modifier):
        self._modifier = modifier

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        self._path = path

    def obtain_path(self, essay_path):
        add_args = self.add_args if isinstance(self.add_args, (list, tuple)) else [self.add_args]
        self._path = essay_path + '/' + self.type + '/' + '/'.join(add_args) + '/' + self.name

    @property
    def add_args(self):
        return self._add_args


class Essay():
    """Essay class definition, for dealing with bioinformatics essays"""

    def __init__(self, state=None):
        self._state = state

    def get_name(self):
        return self._state['name']

    def get_project_name(self):
        return self._state['project_name']

    def get_path(self):
        return self._state['path']

    def get_target_name(self):
        try:
            target_name = self._state['target_reference']['name']
        except KeyError:
            target_name = None
        return target_name

    def get_samples(self):
        try:
            sample_list = self._state['samples_order']
        except KeyError:
            sample_list = None
        return sample_list

    def read_essay(self, essay_path):
        try:
            with open(essay_path + '/state.json', 'r') as input_state_file:
                self._state = json.load(input_state_file)
        except Exception, e:
            print 'ERROR: state.json file could not be opened for essay in ', essay_path, '!!'
            print str(e)
            traceback.print_exc()
            sys.exit(1)

    def write_state(self):
        try:
            with open(self._state['path'] + '/state.json', 'w') as state_file:
                json.dump(self._state, state_file, sort_keys=True, indent=4, separators=(',', ': '))
        except Exception, e:
            print 'ERROR: cannot write state to ' + self._state['path'] + '/state.json!!'
            print str(e)
            traceback.print_exc()
            sys.exit(1)

    def add_arg(self, info_type=None, folder_type=None, keys=None, value=None):
        """ function for adding new arguments in state.json. "args" contains all json keys, and the last argument
        contains the value that is going to be assigned to the key"""
        if info_type is None:
            print 'ERROR: No info type has been given to store information in state.json'
            raise
        mydict = self._state[info_type] if folder_type is None else self._state[info_type][folder_type]

        if len(keys) > 1:
            for arg in keys[:-1]:
                if arg not in mydict.keys():
                    mydict[arg] = {}
                mydict = mydict[arg]
        if value is None:
            if keys[-1] not in mydict.keys():
                mydict[keys[-1]] = {}
        else:
            mydict[keys[-1]] = value

    def get_arg(self, info_type=None, folder_type=None, keys=None):
        """ function for getting arguments in state.json. "args" contains all json keys, and the last argument
        contains the value"""
        try:
            if info_type is None:
                print 'ERROR: No info type has been given to retrieve information from state.json'
                raise
            mydict = self._state[info_type] if folder_type is None else self._state[info_type][folder_type]
            if len(keys) > 1:
                for arg in keys[:-1]:
                    if arg not in mydict.keys():
                        mydict[arg] = {}
                    mydict = mydict[arg]
                return mydict[keys[-1]]
        except KeyError:
            print 'ERROR: Key arguments can\'t be found in state.json'

    def add_job(self, module=None, job=None, sge_id=None, log=None, hold_jobs=None):
        """ Gets a new job ID, by looking all existent jobs numbers and assigning next number"""
        job_id = 0
        try:
            assert isinstance(self._state['modules'][module]['jobs'][job]['jobid'], int)
            job_id = self._state['modules'][module]['jobs'][job]['jobid']
        except KeyError or AssertionError:
            if not module in self._state['modules'].keys():
                self._state['modules'][module] = {}
            if not 'jobs' in self._state['modules'][module].keys():
                self._state['modules'][module]['jobs'] = {}
            for mymodule in self._state['modules']:
                if 'jobs' in self._state['modules'][mymodule].keys():
                    for myjob in self._state['modules'][mymodule]['jobs']:
                        if myjob != job and 'jobid' in self._state['modules'][mymodule]['jobs'][myjob].keys():
                            job_id = max(job_id, int(self._state['modules'][mymodule]['jobs'][myjob]['jobid']))
            if not job in self._state['modules'][module]['jobs'].keys():
                self._state['modules'][module]['jobs'][job] = {}
            job_id += 1
        self._state['modules'][module]['jobs'][job]['jobid'] = job_id
        self._state['modules'][module]['jobs'][job]['state'] = 'running'
        hold_jobs = [] if hold_jobs is None else hold_jobs
        self._state['modules'][module]['jobs'][job]['hold_jobs'] = hold_jobs
        self._state['modules'][module]['jobs'][job]['SGE_id'] = sge_id
        self._state['modules'][module]['jobs'][job]['log'] = log
        return self._state['modules'][module]['jobs'][job]['jobid']

    def add_file(self, folder_type=None, path=None, name=None):
        """ Adds file (or folder) to tree structure of essay and creates the folder where the file is included
        """
        if folder_type is None:
            print 'ERROR: No type of information for the folder creation has been given (trash, analysis, logs, jobs,' \
                  ' etc...)'
            raise
        folder_path = self._state['path'] + '/' + folder_type
        path = [path] if isinstance(path, str) else path
        for arg in path:
            folder_path += '/' + arg
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
        if name is not None:
            folder_path += '/' + name
            self.add_arg(info_type='tree', folder_type=folder_type, keys=path, value=name)
        else:
            self.add_arg(info_type='tree', folder_type=folder_type, keys=path)
        return folder_path

    def get_queue(self):
        return self._state['queue']

    def get_jobs_for_module(self, module):
        job_list = []
        for mymodule in self._state['flow'][module]['previous'].split(','):
            for myjob in self._state['modules'][mymodule]['jobs']:
                if self._state['modules'][mymodule]['jobs'][myjob]['jobid'] not in job_list:
                    job_list.append(self._state['modules'][mymodule]['jobs'][myjob]['jobid'])
        return job_list

    def get_sge(self, job):
        for mymodule in self._state['modules']:
            for myjob in self._state['modules'][mymodule]['jobs']:
                if self._state['modules'][mymodule]['jobs'][myjob]['jobid'] == job:
                    return self._state['modules'][mymodule]['jobs'][myjob]['SGE_id']
        return None

    def get_sge_ids(self, hold_jobs):
        sge_ids = []
        hold_jobs = str(hold_jobs) if isinstance(hold_jobs, int) else hold_jobs
        for job in hold_jobs.split(','):
            sge_ids.append(self.get_sge(int(job)))
        if sge_ids is [None]:
            return None
        else:
            return ','.join(sge_ids)

    def submit_from_essay(self, job_name=None, command=None, input_data=None, output_data=None, add_args=None,
                          add_path=None, module=None, run_in_headnode=False, hold_jobs=None):
        """ Submits jobs for essays, with default conventions:
        analysis: definitive results
        trash: temporary results
        logs: log information coming from jobs being launched
        jobs: for testing purposes, jobs being launched are also stored
        """
        try:
            assert isinstance(self._state['modules'][module]['jobs'][job_name]['state'], unicode)
            if self._state['modules'][module]['jobs'][job_name]['state'] == 'pending':
                raise AssertionError
            return self._state['modules'][module]['jobs'][job_name]['jobid']
        except KeyError or AssertionError:
            # Define arguments
            job_args = []
            if input_data is not None:
                input_data.obtain_path(self.get_path())
                job_args.extend([input_data.modifier])
            if output_data is not None:
                output_data.obtain_path(self.get_path())
                job_args.extend([output_data.modifier])
                self.add_file(folder_type=output_data.type, path=output_data.add_args)
            if add_args is not None:
                job_args.extend([add_args])

            path = [module]
            if add_path is not None:
                if isinstance(add_path, str):
                    path.extend([add_path])
                else:
                    path.extend(add_path)

            workdir = self.add_file(folder_type='logs', path=path)
            job = self.add_file(folder_type='jobs', path=path, name=job_name)

            if hold_jobs is None:
                hold_jobs = self.get_jobs_for_module(module)
            elif isinstance(hold_jobs, int):
                hold_jobs = [hold_jobs]

            tmp_command = command + ' ' + ' '.join(job_args)
            with open(job, 'w') as job_file:
                job_file.write(tmp_command)
            if run_in_headnode:
                proc = subprocess.Popen(tmp_command)
                sge_id = proc.pid
            else:
                sge_id = submit_job(job_name=job_name, command=command, add_args=job_args, cwd=workdir,
                                    queue=self.get_queue(), hold_jobs=','.join(str(self.get_sge(x)) for x in hold_jobs))
            # Create job in state.json and create needed folders
            job_id = self.add_job(module=module, job=job_name, sge_id=sge_id, log=workdir,
                                  hold_jobs=hold_jobs)
            self.write_state()
            return job_id