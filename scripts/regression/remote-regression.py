#!/usr/local/bin/python
# -*- coding: utf-8 -*-

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

from Conf.conf import _projects
for p in _projects.keys():
    _projects[p][u'dir'] = None

def builds_and_runs(ident, build_dict, run_dict):
    """ Launch builds using build_dict and runs using run_dict.
    
    build_dict must be like :
      build_dict = {
        machineconfname : {
          build_name : {
            conf     : build_configuration
            project : project_name
         }
       }
    }

    run_dict is given in parameter to project.Run()
    """
    myMachines = {}
    for machineconf in build_dict.keys():
        for build_name in build_dict[machineconf].keys():
            build = build_dict[machineconf][build_name]
            conf = build[u'conf']
            if not myMachines.has_key(build[u'project']):
                myMachines[build[u'project']] = {}
            myMachines[build[u'project']][machineconf] = conf
            print "[%s] Building option %s on %s [%s]" % (build[u'project'],
                                                          build_name, 
                                                          machineconf, 
                                                          conf[u'ROOT'])

            p = _projects[build[u'project']][u'func']()
            p.Configure(conf,
                        build[u'build_opt'],
                        build[u'binaries'],
                        ident)

            p.Build(conf)

    for project in run_dict.keys():
        p = _projects[project][u'func']()
        p.Run(myMachines[project],
              run_dict[project],
              ident)

def launch(file, machine):
    """ 
    launch the commands contained in the given file 
    of the given machine on the given machine.
    """
    fp = open(file, 'r')
    timeout = 0
    if machine.mpdclean != None :
        process = subprocess.Popen(machine.mpdclean, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        output  = process.communicate()[0]

    if machine.mpdboot != None:
        process = subprocess.Popen(machine.mpdboot, shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        output  = process.communicate()[0]
        print output
    while(True):
        line = fp.readline()
        if line == "": break
        while (True):
            p = subprocess.Popen(line, shell=True, 
                                 stdout=subprocess.PIPE, 
                                 stderr=subprocess.STDOUT)
            t_beginning = time.time()
            seconds_passed = 0
            while True:
                if p.poll() is not None:
                    break
                seconds_passed = time.time() - t_beginning
                if timeout != 0 and seconds_passed > timeout:
                    p.kill()
                    break
                time.sleep(0.1)

            output  = p.communicate()[0]
            retcode = p.poll()
            #print u'%s : %d' % (line, retcode)
            if machine.qstat == None or retcode == 0 : break
            time.sleep(machine.sleeptime)

    fp.close()

#####################################################
# This is the entry point of the Regression testing #
#####################################################
if __name__ == "__main__":
    # Defining the options of the program
    usage="""usage: %prog [options]

Launch compilations and run defined in Conf/conf.py.
"""
    machinename = socket.gethostbyaddr(socket.gethostname())[0]
    username    = os.getenv('USER')
    machinename = machinename+u'-'+username

    parser = OptionParser(usage=usage)
    # Adding the options
    parser.add_option("-i", "--id", dest="id",
                      help="regression identifiaction string",
                      metavar="string", type=str, default=[])

    parser.add_option("-b", "--build_dict", dest="build_dict",
                      help="regression build dictionary file path",
                      metavar="string", type=str, default=[])

    parser.add_option("-n", "--run_dict", dest="run_dict",
                      help="regression build dictionary file path",
                      metavar="string", type=str, default=[])

    parser.add_option("-r", "--run", dest="run",
                      help="run",
                      action="store_true", default=False)

    # Parse the command line
    (options, args) = parser.parse_args()
    
    if options.run:
        for project in _projects.keys():
            ident = options.id
            resultdir = remote_machines[machinename].resultdir
            filename = os.path.join(resultdir,
                                    ident,
                                    u'%s-launch.sh' % (project))
            if os.path.isfile(filename) : 
                launch(filename, remote_machines[machinename]) 
                remote_machines[machinename].waitJobs(ident)
                filename = os.path.join(resultdir, 
                                        ident,
                                        u'%s-execution.xml' % (project))
                p = _projects[project][u'func']()
                p.Parse(filename, ident, remote_machines[machinename])
    else:
        ident = options.id
        build_dict = pickle.load(open(options.build_dict))
        run_dict   = pickle.load(open(options.run_dict))
        builds_and_runs(ident, build_dict, run_dict)

    


