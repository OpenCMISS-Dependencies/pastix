#!/usr/local/bin/python
# -*- coding: utf-8 -*-

""" This scripts launches regression test on several projects """
__autor__ = u"Xavier Lacoste "
__license__ = "CeCILL-C"

from Machines.machines import machines, remote_machines
from Projects.PaStiX import PaStiX
from Projects.Hips import Hips
from Projects.Scotch import Scotch
from optparse import OptionParser
from copy import deepcopy
import os
import re
import datetime
import subprocess
import socket
import subprocess
import pickle
import time


#try: 
from Conf.conf import _liste_run, _database_data, _projects
#except ImportError:
#    from Conf.default import _liste_run, _database_data, _pastix_dir

###############################
_tmp_dir         = u'/tmp/%s'
_log_plafrim     = u'%s.o*'
_regression_dir  = u'/tmp/%s/regression'
_regression_arch = u'/tmp/%s/regression.tar.gz'
_regression_remote_dir  = u'/tmp/%s-regression/regression'
_regression_remote_arch = u'/tmp/%s-regression.tar.gz'

from Database.database import fillDatabase

#####################################################
# This is the entry point of the Regression testing #
#####################################################
if __name__ == "__main__":
    # Defining the options of the program
    usage="""usage: %prog [options]

Launch compilations and run defined in Conf/conf.py.
"""
    username    = os.getenv('USER')

    parser = OptionParser(usage=usage)
    # Adding the options

    parser.add_option("-o", "--only_run", dest="only_run",
                      help="Only launch runs, suppose build has been done",
                      action="store_true", default=False)

    parser.add_option('-e', '--export', dest='export',
                      help='Create a release and test it',
                      action='store_true', default=False)

    parser.add_option('-u', '--username', dest='username',
                      help='username to use for export command',
                      metavar='string', type=str, default=[])

    # Parse the command line
    (options, args) = parser.parse_args()
    
            
    if __debug__:
        print "Remote regression"
    # generate an id (unique) for the tests : 14 numbers based on the time
    # (we assume we don't launch 2 tests in less than 1sec)
    ident = "%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d"% (
        datetime.datetime.now().timetuple()[0:6] )

    username    = os.getenv('USER')

    if (options.username):
        username = options.username

    threads  = []
    remotes  = []
    
    branches = {}
    for project in _projects.keys():
        branches[project] = []

    for run in _liste_run:
        for remote_machine_name in run[u'remote_machines']:
            if not remote_machine_name in remotes:
                remotes.append(remote_machine_name)

    # buildict = {
    #     remote_machine_name : {
    #         machinename : { 
    #             build_opt_set_name : {  
    #                 u'project'  : projectName, 
    #                 u'build_opt': build_set, 
    #                 u'binaries' :[binarylist]}}}}
    if options.only_run == False:
        build_dict = {}
        #project by machine
        pbm = {}
        for run in _liste_run:
            for machineconf in run[u'machines']:
                remote_machine_name = machines[machineconf][u'machinename']
                if not build_dict.has_key(remote_machine_name):
                    build_dict[remote_machine_name] = {}
                    pbm[remote_machine_name] = []
                if not run[u'project'] in pbm[remote_machine_name]:
                    pbm[remote_machine_name].append(run[u'project'])
                if not build_dict[remote_machine_name].has_key(machineconf):
                    build_dict[remote_machine_name][machineconf] = {}

                for build in run[u'compilations']:
                    if not  build_dict[remote_machine_name][machineconf].has_key(build[u'NAME']):
                        
                        conf = deepcopy(machines[machineconf])
                    
                        conf[u'ROOT'] = _projects[run[u"project"]][u"export_dir"] % (u'%s-regression' % ident)
                        
                        build_dict[remote_machine_name][machineconf][build[u'NAME']] = {
                            u'conf'     : conf,
                            u'project'  : run[u'project'],
                            u'build_opt': build,
                            u'binaries' : []}
                        branch = build[u'BRANCH']
                        if not branch in branches[run[u'project']]:
                            branches[run[u'project']].append(branch)
                    
                    for binary in run[u'binaries']: 
                        if (not binary in build_dict[remote_machine_name][machineconf][build[u'NAME']][u'binaries']):
                            build_dict[remote_machine_name][machineconf][build[u'NAME']][u'binaries'].append(binary)

    # Creating archive of the regression project 
    path =  os.path.dirname(os.path.realpath(__file__))
    cmd  = [u'mkdir', u'-p', os.path.dirname(_regression_arch % ident)]
    if __debug__:
        print "%s" % (str(cmd))
    subprocess.Popen(cmd).wait()
    cmd  = [u'tar', u'-c', os.path.basename(path), u'--exclude', '.git']
    cmd2 = [u'gzip']
    if __debug__:
        print "%s | %s" % (str(cmd), str(cmd2))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                         cwd=os.path.dirname(path))
    
    subprocess.Popen(cmd2, stdin=p.stdout,
                     stdout=open(_regression_arch % ident, u'w')).wait()
    
    
    for project in _projects.keys():
        if (options.export):
            p = _projects[project][u'func']
            if (not p == None) and (not len(branches[project]) == 0):
                p = p()
                p.export(_projects[project][u"export_name"] % ident, 
                         _projects[project][u'svnurl'] % username, 
                         branches[project], ident)
        else:
            p = _projects[project][u'func']
            d = _projects[project][u'dir']
            if (not p == None) and (not d == None):
                p = p()
                p.export(_projects[project][u"export_name"] % ident, d, 
                         branches[project], ident)


    # run_dict = { 
    #     remote_machine_name : {
    #         project_name : {
    #             machineconf : runs[]}}}
    run_dict = {}
    for run in _liste_run:
        for machineconf in run[u'machines']:
            remote_machine_name = machines[machineconf][u'machinename']
            
            if not run_dict.has_key(remote_machine_name):
                run_dict[remote_machine_name] = {}
            project = run[u'project']
            if not run_dict[remote_machine_name].has_key(project):
                run_dict[remote_machine_name][project] = {}
            if not run_dict[remote_machine_name][project].has_key(machineconf):
                run_dict[remote_machine_name][project][machineconf] = []

            run_dict[remote_machine_name][project][machineconf].append(deepcopy(run))


    for remote_machine_name in remotes:
        print u'Building and creating runs on %s' % remote_machine_name
        cmd = [u'mkdir', '-p', _tmp_dir % ident]
        rm = remote_machines[remote_machine_name]
        rm.remoteCallCommand(cmd).wait()
        filename  = u'/tmp/%s/pickle' % ident
        filename2 = u'/tmp/%s/pickle2' %ident
        pickle.dump(build_dict[remote_machine_name], open(filename, u'w'))
        rm .fileCopyTo(filename, filename).wait()
        pickle.dump(run_dict[remote_machine_name], open(filename2, u'w'))
        rm.fileCopyTo(filename2, filename2).wait()
        rm.fileCopyTo(_regression_arch % ident,
                      _regression_remote_arch % ident).wait()
        cmd = [u'mkdir', '-p', os.path.dirname(_regression_remote_dir % ident)]
        rm.remoteCallCommand(cmd).wait()
        cmd = [u'tar', u'-xzf', _regression_remote_arch %ident, u'-C', 
               os.path.dirname(_regression_remote_dir %ident)]
        rm.remoteCallCommand(cmd).wait()

        for pname in pbm[remote_machine_name]:
            p = rm.fileCopyTo(_projects[pname][u"export_name"] % ident, 
                              _projects[pname][u"export_name"] % (
                    u'%s-regression' % ident))
            p.wait()
            print p.communicate()
            cmd = [u'tar', u'-xzf', 
                   _projects[pname][u"export_name"] % (
                    u'%s-regression' % ident), 
                   u'-C', 
                   os.path.dirname(_projects[pname][u"export_dir"] % (
                        u'%s-regression'  % ident))]
            rm.remoteCallCommand(cmd).wait()

        cmd = [u"python"]
        if not __debug__:
            cmd.append(u'-0')
        cmd.append("%s/remote-regression.py" % (_regression_remote_dir % ident))
        cmd.append("--id")
        cmd.append(ident)
        cmd.append("--build_dict")
        cmd.append(filename)
        cmd.append("--run_dict")
        cmd.append(filename2)
        threads.append(remote_machines[remote_machine_name].remoteCallCommand(cmd))

    for t in threads:
        t.wait()
        print t.communicate()[0]

    threads = []
    print u"Now runs!"
    for machineconf in remotes:
        cmd = ["python"]
        if not __debug__:
            cmd.append(u'-0')
        cmd.append("%s/remote-regression.py" % (_regression_remote_dir % ident))
        cmd.append("--run")
        cmd.append("--id")
        cmd.append(ident)
        threads.append(remote_machines[machineconf].remoteCallCommand(cmd))

    for t in threads:
        t.wait()
        print t.communicate()[0]

    fillDatabase(ident, _liste_run, _database_data)
    
    for remote_machine_name in remotes:
        print "Cleaning"
        cmd = [u'rm', '-rf', _regression_remote_dir % ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()
        cmd = [u'rm', '-rf', _regression_remote_arch % ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()
        cmd = [u'rm', '-rf', _regression_dir % ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()
        cmd = [u'rm', '-rf', _regression_arch % ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()
        cmd = [u'rm', '-rf', _tmp_dir % ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()
        cmd = [u'rm', '-rf', _log_plafrim %ident]
        remote_machines[remote_machine_name].remoteCallCommand(cmd).wait()


                                         
