#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement # This isn't required in Python 2.6
import os
import subprocess
import sys
import hashlib
import re

from copy import deepcopy
import xml.etree.cElementTree as xml

def getText(nodelist):
    """Get text from an XML node"""
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)

    return ''.join(rc)



class Project:
    """ Class project :

    The class project must contains
      _liste_default_yes : compilation flag which are defined by default.
      _liste_default_not : compilation flag which are not defined by default.
      _optionsComp       : available compilation option
      _binaries          : available binaries
      _config            : configuration file template
    """
    _name = 'Default'
    _liste_default_yes = []
    _liste_default_not = []
    _optionsComp = {}
    _binaries = {}
    _config   = u""
    _config_name = u"config.in"
    _env = os.environ
    _env["LC_CTYPE"] = "en_US.UTF-8" 
    _env["LANG"]     = "en_US.UTF-8"
    _configFlags = {}
    _parse_util = {}
    _parse_count = {}
    def __init__(self):
        """ Initialisation of the class

        Add entries to _optionsComp using _liste_default_yes and
        _liste_default_not

        Replace __OPT__ in config file by a list of __keywords__ using
        _liste_default_yes and _liste_default_not keys as keywords.
        """
        for key in self._liste_default_yes:
            self._optionsComp[key] = {
                u'database_id' : key,
                u'default'     : {
                    u'replace'       : {u'__%s__' % (key) : u'-D%s' % (key)},
                    u'database_val'  : u'DEFINED'},
                False          : {
                    u'replace'       : {u'__%s__' % (key) : u''},
                    u'database_val'  : u'UNDEFINED',
                    u'searchInName'  : re.compile(u'_NO_%s' % key)}
                }

        for key in self._liste_default_not:
            self._optionsComp[key] = {
                u'database_id' : key,
                u'default'     : {
                    u'replace'       : {u'__%s__' % (key) : u''},
                    u'database_val'  : u'UNDEFINED'},
                True           : {
                    u'replace'       : {u'__%s__' % (key) : u'-D%s' % (key)},
                    u'database_val'  : u'DEFINED',
                    u'searchInName'  : re.compile(u'_%s' % key)}
                }


        for key in self._liste_default_yes:
            self._config = self._config.replace(u'__OPT__',
                                                u'__%s__ __OPT__' % key)
        for key in self._liste_default_not:
            self._config = self._config.replace(u'__OPT__',
                                                u'__%s__ __OPT__' % key)
        self._config = self._config.replace(u'__OPT__', u'')

        return

    def _addTextChild(self, node, key, value, attributes={}):
        """  Add a text child to an XML node. 
        
        The resulting XML will be :
        <node_key>
        ....
        <key attributes_key=attributes_val ...> value </key>
        </node_key>

        """
        child = xml.Element(key)
        for key2 in attributes.keys():
            child.attrib[key2] = attributes[key2]
        child.text  = value
        node.append(child)



    def getRevision(self):
        """ Get revision from build directory 

        Check if __REVISION__ file exists.
        If it exists the contained integer will be used as revision.
        Else, will get SVN revision consulting the repository. 

        If repository contains .svn, we use svn log -l 1 and
        search for the revision number in it.
        Else, we use git svn log -n 1. Thus, **it only works with SVN
        central repository...** 

        """
        rev_search = re.compile('r([0-9]+)')
        revision = -1
        if (os.path.isfile(os.path.join(self.branchdir, u'__REVISION__'))):
            return open(os.path.join(self.branchdir,
                                     u'__REVISION__'), 'r').read()
        if (os.path.isdir(os.path.join(self.branchdir, u'.svn'))):
            try:
                process = subprocess.Popen([u'svn', u'log', u'-l', u'1'],
                                           env=self._env,
                                           cwd=self.branchdir,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT)
                svn_log = process.communicate()[0]
                svn_log = svn_log.decode(u'utf8',
                                         u'replace')
                res = rev_search.search(svn_log)
                if (res):
                    revision = int(res.group(1))
            except:
                revision = -1
                return revision
        else:
            try:
                process = subprocess.Popen([u'git', u'rev-parse', 'HEAD'],
                                           env=self._env,
                                           cwd=self.branchdir,
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.STDOUT)
                git_hash = process.communicate()[0]
                git_hash = git_hash.decode(u'utf8', u'replace')
                revision = git_hash.strip()
            except:
                revision = -1
                return revision
        return revision

    def getDiff(self):
        """ Get differences from build directory 

        Check for a __DIFF__ file in the directory and use it as diff if
        it exists.
        
        Else, if .svn directory exists we use svn diff and if not we use 
        git diff remotes/git-svn.

        If an error occurs with SVN we return a u"Error on getting
        diff", and if git doesn't works u"Could not get diff" is returned.

        """
        svn_diff = u""
        if (os.path.isfile(os.path.join(self.branchdir, u'__DIFF__'))):
            return open(os.path.join(self.branchdir, u'__DIFF__'), 'r').read()
        if (os.path.isdir(os.path.join(self.branchdir, u'.svn'))):
            try:
                process = subprocess.Popen([u'svn', u'diff'],
                                           env=self._env,
                                           cwd=self.branchdir,
                                           stdout=subprocess.PIPE)
                svn_diff = process.communicate()[0]
                svn_diff = svn_diff.decode(u'utf8', u'replace')
            except:
                svn_diff = u"Error on getting diff"
        else:
            try:
                process = subprocess.Popen([u'git',u'diff',
                                            u'remotes/git-svn', u'.'],
                                           env=self._env,
                                           cwd=self.branchdir,
                                           stdout=subprocess.PIPE)
                svn_diff = process.communicate()[0]
                svn_diff = svn_diff.decode(u'utf8', u'replace')
            except:
                try:
                    process = subprocess.Popen([u'git',u'diff',
                                                u'remotes/orgin/master', u'.'],
                                               env=self._env,
                                               cwd=self.branchdir,
                                               stdout=subprocess.PIPE)
                    svn_diff = process.communicate()[0]
                    print u"yes"
                    svn_diff = svn_diff.decode(u'utf8', u'replace')
                except:
                    svn_diff = u"Could not get diff"
        return svn_diff

    def Configure(self,  machine, options, compilations, ident):
        """Configure the project on the given machine

        Replaces __keywords__ in config and build an XML file
        representing the compilation."""

        self.branchdir    = os.path.join(machine[u'ROOT'],options[u"BRANCH"])
        self.compilations = compilations
        self.options      = options
        self.user         = machine[u"username"]
        self.resultdir    = os.path.join(machine[u"result_dir"], ident)

        # Create the xml document
        self.doc_xml = xml.Element(u'TestCase')
        config_xml = xml.Element(u"config")
        self.doc_xml.append(config_xml)

        self._addTextChild(config_xml, u"PROJECT", self._name)
        self._addTextChild(config_xml, u"NAME",    options[u"NAME"])
        self._addTextChild(config_xml, u"MACHINE", machine[u"NAME"])
        date = self._idToDate(ident)
        self._addTextChild(config_xml, u"DATE",    date)
        self._addTextChild(config_xml, u"USER",    machine[u"username"])
        self.diff = self.getDiff()
        self._addTextChild(config_xml, "DIFF",     self.diff)
        try:
            self.diff = self.diff.decode(u'utf8', u'replace')
        except:
            self.diff = u"error in diff"
        self.md5diff = hashlib.md5(self.diff).hexdigest()
        self._addTextChild(config_xml, "MD5DIFF",  self.md5diff)

        self.install = os.path.join(self.resultdir, 
                                    self.user,
                                    machine["NAME"], 
                                    self._name,
                                    options["NAME"])

        revision = str(self.getRevision()).decode(u'utf8', u'replace')

        self._addTextChild(config_xml, u"REVISION", revision)

        for key in self._optionsComp.keys():
            if (options.has_key(key) and
                self._optionsComp[key].has_key(options[key])):
                opt_set = self._optionsComp[key][options[key]]
            else:
                opt_set = self._optionsComp[key][u'default']

            replace_dict = opt_set[u'replace']
            database_val = opt_set[u'database_val']
            database_id  = self._optionsComp[key][u'database_id']
            for key2 in replace_dict.keys():
                self._config   = self._config.replace(key2,
                                                      replace_dict[key2])
            self._addTextChild(config_xml, u"option",
                               database_val,
                               {u'key' : database_id})

        if (options.has_key("INTEGER")):
            if (options['INTEGER'] == "int32"):
                if (machine.has_key(u'SCOTCH_INT32')):
                    self._config = self._config.replace(u'__SCOTCH_DIR__',
                                                        machine[u'SCOTCH_INT32'])
                if (machine.has_key(u'METIS_INT32')):
                    self._config = self._config.replace(u'__METIS_DIR__',
                                                        machine[u'METIS_INT32'])
            elif (options['INTEGER'] == "int64"):
                if (machine.has_key(u'SCOTCH_INT64')):
                    self._config = self._config.replace(u'__SCOTCH_DIR__',
                                                        machine[u'SCOTCH_INT64'])
                if (machine.has_key(u'METIS_INT64')):
                    self._config = self._config.replace(u'__METIS_DIR__',
                                                        machine[u'METIS_INT64'])
            elif (options['INTEGER'] == "long"):
                if (machine.has_key(u'SCOTCH_LONG')):
                    self._config = self._config.replace(u'__SCOTCH_DIR__',
                                                        machine[u'SCOTCH_LONG'])
                if (machine.has_key(u'METIS_LONG')):
                    self._config = self._config.replace(u'__METIS_DIR__',
                                                        machine[u'METIS_LONG'])        
            else:
                if (machine.has_key(u'SCOTCH_INT')):
                    self._config = self._config.replace(u'__SCOTCH_DIR__',
                                                        machine[u'SCOTCH_INT'])
                if (machine.has_key(u'METIS_INT')):
                    self._config = self._config.replace(u'__METIS_DIR__',
                                                        machine[u'METIS_INT'])
        else:
            if (machine.has_key(u'SCOTCH_INT')):
                self._config = self._config.replace(u'__SCOTCH_DIR__',
                                                    machine[u'SCOTCH_INT'])
            if (machine.has_key(u'METIS_INT')):
                self._config = self._config.replace(u'__METIS_DIR__',
                                                    machine[u'METIS_INT'])

        for key in self._configFlags.keys():
            machine_key = self._configFlags[key][u'machine_key']
            if machine.has_key(machine_key) : tmp = machine[machine_key]
            else : tmp = self._configFlags[key][u'default']
            self._config = self._config.replace(key,  tmp)

        self.md5conf = hashlib.md5(self._config).hexdigest()

        self._addTextChild(config_xml, u"MD5CONF",    self.md5conf)
        self._addTextChild(config_xml, u"CONFIG_FILE", self._config)

        if (machine.has_key(u"HWLOC_HOME")):
            defined = u"DEFINED"
            self._config = self._config.replace(u"__HWLOC_HOME__",
                                                machine[u'HWLOC_HOME'])
        else:
            defined  = u"UNDEFINED"

        self._addTextChild(config_xml, u"option", defined, {u'key' : u"HWLOC"})

        return config_xml


    def Build(self, machine):
        """Build the project on the given machine"""
        # Copy config.in
        f=open(os.path.join(self.branchdir, self._config_name), 'w')
        f.write(self._config)
        f.close()

        search_warning = re.compile(u"warning\s*:",
                                    flags=re.IGNORECASE|re.MULTILINE)
        search_remark  = re.compile(u"remark\s*:",
                                    flags=re.IGNORECASE|re.MULTILINE)
        search_error   = re.compile(u"error\s*:",
                                    flags=re.IGNORECASE|re.MULTILINE)

        subprocess.Popen([u"mkdir", u"-p", self.install], env=self._env).wait()
        subprocess.Popen([u"cp",
                          os.path.join(self.branchdir, self._config_name),
                          os.path.join(self.install,   self._config_name)],
                         env=self._env).wait()

        retcode = 0
        env=self._env
        if not machine[u"MODULES" ] == None:
            if not machine[u"MODULES" ] == u"":
                for module in machine[u"MODULES"]:
                    print "env[u\"module\"] '%s'" % env[u"module"]
                    pipe = subprocess.Popen(u"module load %s; env" % (module),
                                            stdout=subprocess.PIPE, shell=True,
                                            env=env)
                    output = pipe.communicate()[0]
                    env={}
                    for line in output.splitlines():
                        words = line.split('=',1)
                        if len(words) == 2:
                            key = words[0]
                            env[key]=words[1]
                        else:
                            print line
                            env[key]+=u"\n"+line
        pipe = subprocess.Popen(u"pkg-config --libs-only-L | sed -e 's/-L//'",
                                stdout=subprocess.PIPE, shell=True,
                                env=env)
        output = pipe.communicate()[0]
        env[u'HWLOC_HOME']=output


        for libname in self._libraries.keys():
            make_xml = xml.Element(u"Build")
            self.doc_xml.append(make_xml)
            print u"building %s on %s with %s [%s]" %(libname, machine[u'NAME'],
                                                      self.options[u'NAME'],
                                                      machine[u'ROOT'])
            self._addTextChild(make_xml, u"OBJECT_TYPE", u"library")
            self._addTextChild(make_xml, u"OBJECT_NAME", libname)

            commandes = self._libraries[libname]

            log = ""

            for commande in commandes:
                commande = commande.replace(u'__MAKE__',
                                            machine['MAKE'])

                commande = commande.replace(u'__MAKE_J__',
                                            machine['MAKE_J'])
                if (retcode == 0):
                    log += u'## Running '+commande+u'\n'
                    make = subprocess.Popen(commande,
                                            stdout=subprocess.PIPE,
                                            stderr=subprocess.STDOUT,
                                            cwd=self.branchdir,
                                            shell=True, env=env)
                    log += make.communicate()[0].decode(u'utf8', u'replace')
                    retcode = make.returncode

                    if retcode < 0:
                        string = u"\"%s\" was terminated by signal %d"
                        strin = string % (commande, -retcode)
                        print >>sys.stderr, string
                        log += string
                    elif not retcode == 0:
                        string = u"\"%s\" returned %d" % (commande, retcode)
                        print >>sys.stderr, string
                        log += string

            warnings = 0
            warnings += len(search_warning.findall(log))
            warnings += len(search_remark.findall(log))
            errors = 0
            errors += len(search_error.findall(log))
            if not retcode == 0:
                errors = errors + 1
                
            self._addTextChild(make_xml, u"WARNINGS",
                               unicode(str(warnings)))
            self._addTextChild(make_xml, u"ERRORS",   unicode(str(errors)))
            self._addTextChild(make_xml, u"LOG",      log)
            self._addTextChild(make_xml, u"RETURN",   unicode(str(retcode)))

        for compilation in self.compilations:
            print u"building %s on %s with %s [%s]" %(compilation,
                                                      machine[u'NAME'],
                                                      self.options[u'NAME'],
                                                      machine[u'ROOT'])
            make_xml = xml.Element(u"Build")
            self.doc_xml.append(make_xml)
            self._addTextChild(make_xml, u"OBJECT_TYPE", u"binary")
            self._addTextChild(make_xml, u"OBJECT_NAME", compilation)
            cmd = [machine[u"MAKE"],
                   self._binaries[compilation][u'make_cmd']]
            log = "## %s\n" % (" ".join(cmd))
            build_dir = self._binaries[compilation]['build_dir']

            try:
                make = subprocess.Popen(cmd,
                                        env=env,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.STDOUT,
                                        cwd=os.path.join(self.branchdir,
                                                         build_dir))
            except OSError:
                print "Couldn't run %s in %s" % (cmd, os.path.join(self.branchdir,
                                                                   build_dir))
                raise

            log += make.communicate()[0].decode(u'utf8', u'replace')

            warnings = 0
            warnings += len(search_warning.findall(log))
            warnings += len(search_remark.findall(log))
            errors = 0
            errors += len(search_error.findall(log))

            retcode = make.returncode
            self._addTextChild(make_xml, u"WARNINGS",
                               unicode(str(warnings)))

            if not retcode == 0:
                errors = errors + 1

            self._addTextChild(make_xml, u"ERRORS",   unicode(str(errors)))
            self._addTextChild(make_xml, u"LOG",      log)
            self._addTextChild(make_xml, u"RETURN",   unicode(str(retcode)))

            if retcode < 0:
                string = u"\"%s\" was terminated by signal %d"
                string = string % (u" ".join(cmd), -retcode)
                print >>sys.stderr, string
            elif retcode == 0:
                binary_name = self._binaries[compilation][u'binary']
                subprocess.Popen([u"cp",
                                  os.path.join(self.branchdir,
                                               binary_name),
                                  self.install] ,
                                 env= env).wait()
            else:
                string = u"\"%s\" returned %d"
                string = string % (u' '.join(cmd), retcode)
                print >>sys.stderr, string

        filename = os.path.join(self.install, u"compilations.xml")
        xml.ElementTree(self.doc_xml).write(filename)

    def _runOne(self, binary, case, parameters,
                machine, options, revision, ident):
        raise NotImplementedError


    def Run(self, machines, run_dict, ident):
        """Run tests for the project"""
        launch =u""

        self.doc_xml = xml.Element(u"Runs")

        for machineconf in run_dict.keys():
            machine         = machines[machineconf]
            print u"Running on %s [%s]" % (machine[u'NAME'], machine[u'ROOT'])
            for run in run_dict[machineconf]:
                binaries        = run['binaries']
                parameters_list = run['parameters']
                cases           = run['cases']
                for options in run['compilations']:


                    self.user      = machine["username"]
                    self.resultdir = os.path.join(machine["result_dir"],
                                                  ident)



                    self.branchdir    = machine[u'ROOT']+options[u"BRANCH"]


                    self.diff = self.getDiff()
                    self.md5diff = hashlib.md5(self.diff).hexdigest()


                    revision = str(self.getRevision()).decode(u'utf8',
                                                              u'replace')

                    self.install      = os.path.join(self.resultdir,
                                                     self.user,
                                                     machine["NAME"],
                                                     self._name,
                                                     options["NAME"])
                    m = hashlib.md5()

                    for line in open(os.path.join(self.install,
                                                  self._config_name), "r"):
                        m.update(line)

                    self.md5conf = m.hexdigest()

                    for parameters in parameters_list:
                        for case in cases:
                            for binary in binaries:
                                launch += self._runOne(binary, case, parameters,
                                                       machine, options,
                                                       revision, ident)


        f=open(os.path.join(self.resultdir,
                            u'%s-launch.sh' % (self._name)), u'w')
        f.write(launch)
        f.close()
        filename = os.path.join(self.resultdir,
                                u'%s-execution.xml' % (self._name))
        xml.ElementTree(self.doc_xml).write(filename)



    def Parse(self, xmlfile, ident, machine):
        self.doc_xml = xml.parse(xmlfile).getroot()
        self.resultdir = os.path.join(machine.resultdir, ident)
        assert self.doc_xml.tag == u"Runs"
        for run_xml in self.doc_xml.getiterator(u'Run'):
            logfile = run_xml.find(u'logfile').text
            for count_key in self._parse_count.keys():
                self._parse_count[count_key][u"count"] = 0
            try:
                log = u""
                for line in open(logfile, u'r'):
                    log += line
                    for count_key in self._parse_count.keys():
                        count_re = self._parse_count[count_key][u"parse_key"]
                        if (count_re.search(line)):
                            self._parse_count[count_key][u"count"] += 1
                    for parse_key in self._parse_util.keys():
                        parse_re = self._parse_util[parse_key][u"parse_key"]
                        parse_type = self._parse_util[parse_key][u"type"]
                        resultat = parse_re.search(line)
                        if (resultat):
                            child = xml.Element(u"OUTPUT")
                            child.attrib[u"key"] = parse_key
                            if (parse_type == u"memory"):
                                child.attrib[u"type"] = u"BIGINT"
                                mem = float(resultat.group(1))
                                if resultat.group(2) == u"Ko":
                                    mem  = mem*1024
                                if resultat.group(2) == u"Mo":
                                    mem  = mem*1024*1024
                                if resultat.group(2) == u"Go":
                                    mem  = mem*1024*1024*1024
                                text = str(int(mem))
                            else:
                                child.attrib[u"type"] = parse_type
                                text  = resultat.group(1)
                            child.text = unicode(text)
                            run_xml.append(child)

                for count_key in self._parse_count.keys():
                    child = xml.Element(u"OUTPUT")
                    child.attrib[u"key"]  = count_key
                    child.attrib[u"type"] = u"INT"
                    count = str(self._parse_count[count_key][u"count"])
                    count = count.decode(u'utf8', u'ignore')
                    child.text = count
                    run_xml.append(child)

                child   = xml.Element(u"LOG")
                log = log.decode(u'utf8', u'ignore')
                child.text = log
                run_xml.append(child)
            except IOError:
                child   = xml.Element(u"LOG")
                log     = u"No log file"
                child.text = log
                run_xml.append(child)

        self._checkError()

        xml.ElementTree(self.doc_xml).write(xmlfile)

    def _checkError(self):
        """Count number of errors or warnings in output log"""
        raise NotImplementedError

    def listBuildFlags(self):
        ''' Return the list of available flags'''
        liste = {}
        for flag in self._optionsComp.keys():
            liste[flag] = []
            for value in self._optionsComp[flag].keys():
                if not (value == u'default' or
                        value == u'database_id'):
                    liste[flag].append(value)

        return liste


    def buildCompilationFromName(self, name, branch = u'/trunk'):
        ''' build a compilation dictionary from name '''
        dictionary = {u'NAME'   : name,
                      u'BRANCH' : branch}
        for flag in self._optionsComp.keys():
            opt = self._optionsComp[flag]
            for value in opt.keys():
                if not (value in [u'default', u'database_id']):
                    if opt[value].has_key(u'searchInName'):
                        if opt[value][u'searchInName'].search(name):
                            dictionary[flag] = value

        return dictionary


    def listCompilationNames(self):
        ''' Create a list of available keywords for compilation name '''
        liste = []
        for flag in self._optionsComp.keys():
            liste2 = [flag]
            opt = self._optionsComp[flag]
            for value in opt.keys():
                if not (value in [u'default', u'database_id']):
                    if opt[value].has_key(u'searchInName'):
                        liste2.append(opt[value][u'searchInName'].pattern)
            liste.append(liste2)

        return liste

    def _idToDate(self, ident):
        """Transform ident string into date string.

        String format : 
           ident : yyyymmddhhmmss
           Date  : yyyy-mm-dd hh:mm:ss
           """
        year  = ident[0:4]
        month = ident[4:6]
        day   = ident[6:8]
        hour  = ident[8:10]
        minut = ident[10:12]
        sec   = ident[12:14]
        yearmonday = u"-".join([year,
                                month,
                                day])
        hourminsec = u':'.join([hour,
                                minut,
                                sec])
        return u' '.join([ yearmonday,
                           hourminsec ])

