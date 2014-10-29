import xml.etree.cElementTree as xml

class FileNotSupported(Exception):
    def __init__(self,raison):
        self.raison = raison
    
    def __str__(self):
        return self.raison

try:
    import MySQLdb
except ImportError:
    def fillDatabase(ident):
        return True
else:
    import os
    import sys
    sys.path.append(u'..')
    from Projects.PaStiX import PaStiX

    from Machines.machines import remote_machines


#print Machines.machines.machines
    def fillDatabase(ident, liste_runs, database_data):
        database = MySQLdb.connect (host    = database_data[u"host"  ],
                                    user    = database_data[u"user"  ],
                                    passwd  = database_data[u"passwd"],
                                    db      = database_data[u"db"    ],
                                    charset = u'utf8')
        create_database(database)
        for run in liste_runs:
            for machineconf in run[u'remote_machines']:
                cmd = [u"find"]
                cmd.append(os.path.join(remote_machines[machineconf].resultdir,
                                        ident))
                cmd.append(u"-regex")
                cmd.append(u"\".*\.xml\"")
                machine = remote_machines[machineconf]
                files = machine.remoteCallCommand(cmd).communicate()[0]
                if (not files == None):
                    for f in files.splitlines():
                        print f
                        cp = machine.fileCopyFrom(f,u"/tmp/%s-tmp.xml" % ident)
                        cp.wait()
                        try:
                            xmlToSQL(database, u"/tmp/%s-tmp.xml" % ident, 
                                     database_data[u"db"])
                        except FileNotSupported:
                            print "File not supported : "+ f
                        except:
                            raise

        database.close()
        cmd = ["rm", "-rf"]
        cmd.append(os.path.join(machine.resultdir, ident))
        remote_machines[machineconf].remoteCallCommand(cmd).wait()

    def addIfNotExist(database, tableName, columnName, sqltype, db):
        """
        Add a new column, using the given type, in the given table if
        it doesn't exist.
        """

        cmd = u"""SELECT * FROM `information_schema`.`COLUMNS`
                  WHERE `COLUMN_NAME`='%s'
                  AND `TABLE_NAME`='%s'
                  AND `TABLE_SCHEMA`='%s'""" % (columnName, tableName, db)

        cursor = database.cursor ()
        cursor.execute(cmd)
        if (None ==  cursor.fetchone()):
            cursor2 = database.cursor ()
            cmd = u"""ALTER TABLE `%s`.`%s`
                  ADD COLUMN `%s` %s""" % (db, tableName, columnName, sqltype)
            cursor2.execute(cmd)
            cursor2.close()
        cursor.close ()

    def create_database(database):
        """
        Create database to store regression results if it
        doesn't exist.
        """
        create_build = """
        CREATE TABLE BUILD (
            PROJECT       VARCHAR(10),
            MACHINE       VARCHAR(50),
            USER          VARCHAR(20),
            REVISION      VARCHAR(40),
            DATE          DATETIME,
            MD5DIFF       VARCHAR(32),
            MD5CONF       VARCHAR(32),
            OBJECT_NAME   VARCHAR(50),
            OBJECT_TYPE   VARCHAR(10),
            NAME          VARCHAR(100),
            `RETURN`      INT,
            WARNINGS      INT,
            ERRORS        INT,
            CONFIG_FILE   TEXT,
            DIFF          TEXT,
            LOG           LONGTEXT,
            PRIMARY  KEY (
               `PROJECT`, `MACHINE` ,  `USER` ,  `REVISION` ,  `DATE` ,  `MD5DIFF` ,
               `MD5CONF`, `OBJECT_NAME`)
        ) CHARACTER SET utf8;"""



        create_run = """
        CREATE TABLE RUN (
            PROJECT       VARCHAR(10),
            BINARY_NAME   VARCHAR(50),
            CASE_NAME     VARCHAR(40),
            BUILD         VARCHAR(40),
            USER          VARCHAR(20),
            DATE          DATETIME,
            MACHINE       VARCHAR(50),
            MD5DIFF       VARCHAR(32),
            MD5CONF       VARCHAR(32),
            OPT_SET_NAME  VARCHAR(50),
            MPI_PROC_NBR  SMALLINT,
            THREAD_NBR    SMALLINT,
            REVISION      VARCHAR(40),
            LOG           LONGTEXT,
            PRIMARY KEY (
               `PROJECT`, `BINARY_NAME`, `CASE_NAME`, `USER`, `DATE`,
               `MACHINE`, `MD5CONF`, `MD5DIFF`, `OPT_SET_NAME`,
               `MPI_PROC_NBR`, `THREAD_NBR`, `REVISION`)
        ) CHARACTER SET utf8;"""

        cursor = database.cursor ()
        cursor.execute(u"SHOW TABLES")
        row = cursor.fetchone ()
        rows = []
        while (not row == None):
            rows.append(row[0])
            row = cursor.fetchone ()
        if (not rows or not "BUILD" in rows):
            cursor2 = database.cursor()
            cursor2.execute(create_build)
            cursor2.close ()
        if (not rows or not "RUN" in rows):
            cursor2 = database.cursor()
            cursor2.execute(create_run)
            cursor2.close ()
        cursor.close ()

    def xmlToSQL(database, filename, db):
        """File database from the given xml file"""
        doc_xml = xml.parse(filename).getroot()
        if (doc_xml.tag == u"TestCase"):

            config = doc_xml.find(u"config")

            keys = [u'PROJECT', u'NAME', u"MACHINE", u"USER",
                    u"REVISION", u"DATE",
                    u"MD5CONF", u"MD5DIFF", u'DIFF', u"CONFIG_FILE"]
            config_dict = {}
            for key in keys:
                config_dict[u"`%s`" % key]  = u"'%s'" % (config.find(key).text)

            for option in config.getiterator(u'option'):
                key=option.attrib["key"]
                config_dict[u"`%s`" % key] = u"'%s'" % option.text
                addIfNotExist(database, u"BUILD",
                              key,
                              u"VARCHAR(30)", db)

            for build in doc_xml.getiterator(u'Build'):
                build_dict= config_dict
                keys = [u"OBJECT_TYPE", u"OBJECT_NAME", u"LOG",
                        u"RETURN", u"WARNINGS", u"ERRORS"]
                for key in keys:
                    build_dict[u"`%s`" % (key)]  = u"'%s'" % \
                        (unicode(build.find(key).text).replace("'","\\'") )

                cmd  = u"REPLACE INTO `BUILD` ( "
                cmd += u", ".join(build_dict.keys())
                cmd += u") VALUES ("
                cmd += u", ".join(build_dict.values())
                cmd += u")"

                cursor = database.cursor ()
                try:
                    cursor.execute(cmd)
                except:
                    print "Error running : SQL command"

                    
                cursor.close ()

        elif (doc_xml.tag == "Runs"):
            for run in doc_xml.getiterator(u'Run'):
                run_dict = {}

                for option in run.getiterator(u'option'):
                    run_dict[u'`%s`' % (option.attrib[u"id"])] = \
                        u"'%s'" % (unicode(option.text).replace("'","\\'") )
                    addIfNotExist(database, u'RUN',
                                  option.attrib[u"id"],
                                  option.attrib[u"type"], db)


                keys = [u'PROJECT',   u"BINARY_NAME", u"BUILD", u"OPT_SET_NAME",
                        u"CASE_NAME", u"MPI_PROC_NBR",
                        u"THREAD_NBR", u"USER", u"MACHINE",
                        u"REVISION", u'LOG', u"DATE", u"MD5CONF", u"MD5DIFF"]
                for key in keys:
                    node = run.find(key)
                    if node == None:
                        print "ERROR: No tag %s" % key
                    else:
                        text = unicode(node.text).replace("'","\\'")
                        run_dict[u"`%s`" % (key)]  = u"'%s'" %(text)


                for output in run.getiterator(u'OUTPUT'):
                    key = output.attrib["key"]
                    addIfNotExist(database, u"RUN",
                                  key,
                                  output.attrib["type"], db)
                    run_dict[u"`%s`" % key] = u"'%s'" % unicode(output.text).replace("'","\\'")

                cmd  = "REPLACE INTO `RUN` ( "
                cmd += ", ".join(run_dict.keys())
                cmd += ") VALUES ("
                cmd += ", ".join(run_dict.values())
                cmd += ")"

                cursor = database.cursor ()
                cursor.execute(cmd)
                cursor.close ()

        else:

            raise FileNotSupported("XML file unrocognized : " + filename)
