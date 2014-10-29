import MySQLdb
import MySQLdb.cursors
import mod_python
import re

from mod_python import apache
from mod_python import Session
from cgi import escape
from urllib import unquote

from copy import deepcopy

_build_columns = [u'PROJECT', u'REVISION', u'DATE', u'MACHINE', u'USER', u'NAME',
                  u'OBJECT_NAME', u'WARNINGS', u'ERRORS', u'RETURN']
_same_build_columns = [u'REVISION', u'DATE', 
                       u'WARNINGS', u'ERRORS', u'RETURN']
_run_columns   = [u'PROJECT', u'BINARY_NAME', u'CASE_NAME', u'USER', u'DATE', u'MACHINE',
                  u'OPT_SET_NAME', 
                  u'MPI_PROC_NBR', u'THREAD_NBR', u'MATRIX_SIZE', u'FACT_TIME', u'SOLVE_TIME',
                  u'WARNINGS', u'ERRORS']
_same_run_columns   = [u'REVISION',  u'DATE',
                       u'FACT_TIME', u'SOLVE_TIME',
                       u'WARNINGS', u'ERRORS']

_maj_session_index = u"""javascript:lancerRequete('index.py/majSession?KEY=%s&VALUE=%s','GET');"""
_maj_session       = u"""javascript:lancerRequete('majSession?KEY=%s&VALUE=%s','GET');"""
_rm_column         = u"""javascript:lancerRequete('rmColumn?KEY=%s&INDEX=%s','GET');"""
_add_column        = u"""javascript:lancerRequete('addColumn?KEY=%s&COLUMN=%s&INDEX=%s','GET');"""

_database_data = {
    u"host"   : u"localhost",
    u"user"   : u"pastix",
    u"passwd" : u"PaStiX",
    u"db"     : u"pastix-python"}

_header = u"""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
   <STYLE>
#trier { background-color:white; color:black; border-collapse:collapse; BORDER:white 1px solid; FONT:12 Arial; TEXT-ALIGN:center }
#trier TR { background-color:#ffefd5 }
#trier .title { background-color:#bf2b2f; FONT:14 Arial; color:#ffffff; font-weight:bold }
   SPAN { FONT:bold 12 Arial; CURSOR:pointer }
#BODY { background-color:#FFF5E5 }
#trier TD { BORDER:white 1px solid; }
   </STYLE>
   <head>
   <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
   <title align=center>Results sumary</title> 
   <link rel="stylesheet" type="text/css" href="../pastix.css" />
   <link rel="stylesheet" type="text/css" href="../example.css" />
   <link href="../prettify/prettify.css" type="text/css" rel="stylesheet" />
   <script type="text/javascript" src="../prettify/prettify.js"></script>
   </head>
   <body onload="prettyPrint()" bgcolor="white">
   <script src="../sortable.js" type="text/javascript"
 onerror="alert('Error: failed to load ' + this.src)"></script>

   <script type="text/javascript">
   function lancerRequete(requete,methode)
     {
       if (window.XMLHttpRequest)
         {
            xhr_object = new XMLHttpRequest();
            xhr_object.open(methode, requete, false);
            xhr_object.send(null);
            //xhr_object.onreadystatechange = function() 
            //  { 
                 if(xhr_object.readyState == 4) 
                 {
                   //alert(xhr_object.responseText);
                 }
            //  }
         }
       else if(window.ActiveXObject)
         {
           xhr_object = new ActiveXObject("Microsoft.XMLHTTP");
           xhr_object.open(methode, requete, false);
           xhr_object.send(null);
           if(xhr_object.readyState == 4) 
             {
               //alert(xhr_object.responseText);
             }
         }
       else
         {
           alert('Votre navigateur ne supporte pas les objets XMLHTTPRequest...');
           return(false);
         }
       return(true);
     }
  </script>


"""
_header_no_sort = u"""
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
   <!--STYLE>
#trier { background-color:white; color:black; border-collapse:collapse; BORDER:white 1px solid; FONT:12 Arial; TEXT-ALIGN:center }
#trier TR { background-color:#ffefd5 }
#trier .title { background-color:#bf2b2f; FONT:14 Arial; color:#ffffff; font-weight:bold }
   SPAN { FONT:bold 12 Arial; CURSOR:pointer }
#BODY { background-color:#FFF5E5 }
#trier TD { BORDER:white 1px solid; }
   </STYLE-->
   <head>
   <meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
   <title align=center>Results sumary</title> 
   <link rel="stylesheet" type="text/css" href="../pastix.css" />
   <link rel="stylesheet" type="text/css" href="../example.css" />
   <link href="../prettify/prettify.css" type="text/css" rel="stylesheet" />
   <script type="text/javascript" src="../prettify/prettify.js"></script>
   <script type="text/javascript" src="../prettify/lang-logs-compile.js"></script>
   </head>
   <body onload="prettyPrint()">


"""
_footer = u"""</body></html>"""

_table_header = u"""
<table class="sortable">
  <thead>
"""
_table_header_end = u"""
  </thead>
  <tbody> 
"""

_table_footer = u"""
</tbody>
</table>
"""

_hideorshow = u"""
  <script>
	 function hideOrShow(id)
      {
        var obj;
        obj  = document.getElementById(id);
        if (obj.style.visibility=="collapse")
           obj.style.visibility="";
        else
           obj.style.visibility="collapse";
      }
   </script>
"""


def hello():
   s = u"""\
<html>
<body>
<h2>Hello World!</h2>
</body>
</html>
"""
   return s



def grep(string,list):
   expr = re.compile(string)
   line = 0
   output = u''
   for text in list:
      line = line + 1
      match = expr.search(text)
      if match != None:
         output += u'%s\n' % (match.string)
   return output

def _index(numberOfresults, project=None):
   if project:
      page_html = _header
      page_html += u"""
<h1>Last tests for project %s</h1>
<h2>Last builds</h2>
""" % project
   else:
      page_html = _header.replace(u"../", u"")
      page_html += u"""
<h1>Last tests</h1>
<h2>Last builds</h2>
"""
   database = MySQLdb.connect (host   = _database_data[u"host"  ],
                               user   = _database_data[u"user"  ],
                               passwd = _database_data[u"passwd"],
                               db     = _database_data[u"db"    ])
   if project:
      cmd = u"SELECT DISTINCT REVISION, MACHINE, `DATE` FROM `BUILD` WHERE `PROJECT` like '%s' ORDER BY `DATE` DESC LIMIT 0 , %d;" % (project, numberOfresults)
   else:
      cmd = u"SELECT DISTINCT PROJECT, REVISION, MACHINE, `DATE` FROM BUILD ORDER BY `DATE` DESC LIMIT 0 , %d;" % numberOfresults
   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd)
   row = cursor.fetchone()
   page_html += _table_header

   if project:
      page_html += u"""
    <tr>
      <th>Revision</th>
      <th>Machine</th>
      <th>Date</th>
      <th>Number of builds</th>
      <th>Number of warnings</th>
      <th>Number of errors</th>
    </tr>
"""
      link = u"""href='./print_builds' onClick="%s%s%s%s%s" """
      my_maj_session = _maj_session
   else:
      page_html += u"""
    <tr>
      <th>Project</th>
      <th>Revision</th>
      <th>Machine</th>
      <th>Date</th>
      <th>Number of builds</th>
      <th>Number of warnings</th>
      <th>Number of errors</th>
    </tr>
"""
      link = u"""href='./index.py/print_builds' onClick="%s%s%s%s%s" """
      my_maj_session = _maj_session_index

   page_html += _table_header_end
   if project:
      row_html = u"""<tr class="%s"><td>%s</td><td>%s</td><td>%s</td><td><a %s>%i</a></td><td><a %s>%i</a></td><td><a %s>%i</a></td></tr>"""
   else:
      row_html = u"""<tr class="%s"><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><a %s>%i</a></td><td><a %s>%i</a></td><td><a %s>%i</a></td></tr>"""
   while row :
      if project: 
         my_project = project
      else:
         my_project  = row[u'PROJECT']
      revision    = row[u'REVISION']
      machine     = row[u'MACHINE']
      date     = str(row[u'DATE'])
      
      cmd = """SELECT COUNT(*)      AS NB, 
                 SUM(BUILD.WARNINGS)  AS NB_WARNINGS,
                 SUM(BUILD.ERRORS)    AS NB_ERRORS
                 FROM `BUILD` WHERE `REVISION` LIKE '%s'
                 AND `MACHINE` LIKE '%s'
                 AND `DATE` LIKE '%s'"""
      cmd = cmd % (revision, machine, date)
      if project:
         cmd += """ AND `PROJECT` LIKE '%s'""" % project
      cursor2 = database.cursor()
      cursor2.execute(cmd)
      row2 = cursor2.fetchone()
      cursor2.close()
      links = []
      
      maj_rev    = my_maj_session %(u'REVISION', str(revision))
      maj_proj   = my_maj_session %(u'PROJECT', str(my_project))
      maj_mach   = my_maj_session %(u'MACHINE', machine)
      maj_date   = my_maj_session %(u'DATE', date)
      maj_filter = my_maj_session %(u'FILTER', u'all')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)

      maj_filter = my_maj_session %(u'FILTER', u'warning')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)

      maj_filter = my_maj_session %(u'FILTER', u'error')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)
      if project:
         liste = [revision, machine, date]
      else:
         liste = [my_project, revision, machine, date]
      for r in row2:
         liste.append(links.pop(0))
         liste.append(r)   

      if row2[2] > 0:
         liste.insert(0, u"error")
      elif row2[1] > 0:
         liste.insert(0, u"warning")
      else:
         liste.insert(0, u"ok")     
      page_html += row_html % tuple(liste)
      row = cursor.fetchone()
   cursor.close ()
   page_html += _table_footer


   page_html += u"""
<h2>Last runs</h2>
"""


   if project:
      cmd = u"SELECT DISTINCT REVISION, MACHINE, `DATE` FROM RUN WHERE `PROJECT` LIKE '%s' ORDER BY `DATE` DESC LIMIT 0 , %d;" % (project, numberOfresults)
   else:
      cmd = u"SELECT DISTINCT PROJECT, REVISION, MACHINE, `DATE` FROM RUN ORDER BY `DATE` DESC LIMIT 0 , %d;" % numberOfresults

   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd)
   row = cursor.fetchone()
   page_html += _table_header
   if project:
      page_html += u"""
    <tr>
      <th>Revision</th>
      <th>Machine</th>
      <th>Date</th>
      <th>Number of runs</th>
      <th>Number of warnings</th>
      <th>Number of errors</th>
    </tr>
"""
      link = u""" href='./print_runs' onclick="%s%s%s%s%s" """
   else:
      page_html += u"""
    <tr>
      <th>Project</th>
      <th>Revision</th>
      <th>Machine</th>
      <th>Date</th>
      <th>Number of runs</th>
      <th>Number of warnings</th>
      <th>Number of errors</th>
    </tr>
"""
      link = u""" href='./index.py/print_runs' onclick="%s%s%s%s%s" """

   page_html += _table_header_end

   if project:
      row_html = u"""<tr class="%s"><td>%s</td><td>%s</td><td>%s</td><td><a %s>%i</a></td><td><a %s>%i</a></td><td><a %s>%i</a></td></tr>"""
   else:
      row_html = u"""<tr class="%s"><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><a %s>%i</a></td><td><a %s>%i</a></td><td><a %s>%i</a></td></tr>"""
   while row :
      links = []
      if project: 
         my_project = project
      else:
         my_project  = row[u'PROJECT']
      revision = row[u'REVISION']
      machine  = row[u'MACHINE']
      date     = str(row[u'DATE'])
      cmd = u"""SELECT COUNT(*)      AS NB, 
                 SUM(RUN.WARNINGS)  AS NB_WARNINGS,
                 SUM(RUN.ERRORS)    AS NB_ERRORS
                 FROM `RUN` WHERE `REVISION` LIKE '%s'
                 AND `MACHINE` LIKE '%s'
                 AND `DATE` LIKE '%s'""" 
        
      cmd = cmd % (revision, machine, date)
      page_html += u"<!-- %s -->" % cmd
      cursor2 = database.cursor()
      cursor2.execute(cmd)
      row2 = cursor2.fetchone()
      cursor2.close()
      maj_rev    = my_maj_session %(u'REVISION', str(revision))
      maj_proj   = my_maj_session %(u'PROJECT', str(my_project))
      maj_mach   = my_maj_session %(u'MACHINE', machine)
      maj_date   = my_maj_session %(u'DATE', date)
      maj_filter = my_maj_session %(u'FILTER', u'all')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)

      maj_filter = my_maj_session %(u'FILTER', u'warning')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)

      maj_filter = my_maj_session %(u'FILTER', u'error')
      my_link    = link % (maj_proj, maj_rev, maj_mach, maj_date, maj_filter)
      links.append(my_link)

      if project:
         liste = [revision, machine, date]
      else:
         liste = [my_project, revision, machine, date]
      for r in row2:
         liste.append(links.pop(0))
         liste.append(r)

        
      if row2[2] > 0:
         liste.insert(0, u"error")
      elif row2[1] > 0:
         liste.insert(0, u"warning")
      else:
         liste.insert(0, u"ok")
      i = 0
      for l in liste:
         if l == None :
            liste[i] = -1
            liste[0] = u'error' 
         i = i+1
      page_html += row_html % tuple(liste)
        
      row = cursor.fetchone()
   cursor.close ()

   page_html += _table_footer
   page_html += _footer
   database.close()
   return page_html


def same_build(req):
   session = Session.Session(req)
   PROJECT     = session[u'PROJECT']
   MACHINE     = session[u'MACHINE']
   USER        = session[u'USER']
   MD5DIFF     = session[u'MD5DIFF']
   MD5CONF     = session[u'MD5CONF']
   OBJECT_NAME = session[u'OBJECT_NAME']
   NAME        = session[u'NAME']
   myfilter = 'all'
   if session.has_key(u'SAME_BUILD_COLUMNS'):
      columns = session[u'SAME_BUILD_COLUMNS']
   else:
      columns = _same_build_columns
 
   search_arg = {
      u'PROJECT'     : PROJECT,
      u'OBJECT_NAME' : OBJECT_NAME,
      u'USER'        : USER,
      u'MACHINE'     : MACHINE,
      u'MD5CONF'     : MD5CONF,
      u'MD5DIFF'     : MD5DIFF,
      u'NAME'        : NAME
      }
   link_logs = [u"""<a href='./print_build' onclick="%s%s%s%s%s%s%s">more info</a>""",
                u"""<a href='./print_builds' onclick="%s%s%s%s%s%s%s">build session</a>"""]
   page_html = _header  

   page_html += u"""
<h1>Builds of %s, %s, %s by %s on %s </h1>
""" %(PROJECT, OBJECT_NAME, NAME, USER, MACHINE)
   page_html += _table_header
   page_html += _table_header_end
   page_html += u"""<tr><td><a href=../index.py>Go back to global page</a></td>"""
   page_html += u"""<td><a href=chooseSameBuildColumns>Choose columns</a></td></tr>"""
   page_html += _table_footer

   array = _get_array_from_columns(search_arg, u'BUILD', columns, link_logs, 
                                   [u'MACHINE', u'USER', u'REVISION', u'DATE',
                                    u'MD5DIFF', u'MD5CONF', u'OBJECT_NAME'], myfilter)

   page_html += _array_to_table(array)

   page_html += _footer
   return page_html

def same_run(req):
   session = Session.Session(req)
   PROJECT      = session[u'PROJECT']
   BINARY_NAME  = session[u'BINARY_NAME']
   CASE_NAME    = session[u'CASE_NAME']   
   USER         = session[u'USER']        
   MACHINE      = session[u'MACHINE']     
   MD5CONF      = session[u'MD5CONF']     
   MD5DIFF      = session[u'MD5DIFF']     
   OPT_SET_NAME = session[u'OPT_SET_NAME']
   MPI_PROC_NBR = session[u'MPI_PROC_NBR']
   THREAD_NBR   = session[u'THREAD_NBR']  
   myfilter     = "all"
   if session.has_key(u'SAME_RUN_COLUMNS'):
      columns = session[u'SAME_RUN_COLUMNS']
   else:
      columns = _same_run_columns

   search_arg = {
      u'PROJECT'      : session[u'PROJECT'],
      u'BINARY_NAME'  : session[u'BINARY_NAME'],
      u'CASE_NAME'    : session[u'CASE_NAME'],   
      u'USER'         : session[u'USER'],        
      u'MACHINE'      : session[u'MACHINE'],     
      u'MD5CONF'      : session[u'MD5CONF'],     
      u'MD5DIFF'      : session[u'MD5DIFF'],     
      u'OPT_SET_NAME' : session[u'OPT_SET_NAME'],
      u'MPI_PROC_NBR' : session[u'MPI_PROC_NBR'],
      u'THREAD_NBR'   : session[u'THREAD_NBR'],  
      }
   page_html = _header 
   page_html += u"""
<h1>Runs of %s,%s,%s,%s by %s on %s ( %s threads, %s mpi nodes)</h1>
""" %(PROJECT, BINARY_NAME, CASE_NAME, OPT_SET_NAME, USER, MACHINE, THREAD_NBR, MPI_PROC_NBR)
   page_html += _table_header
   page_html += _table_header_end
   page_html += u"""<tr><td><a href=../index.py>Go back to global page</a></td>"""
   page_html += u"""<td><a href=chooseSameRunColumns>Choose columns</a></td></tr>"""
   page_html += _table_footer

   link_logs = [u"""<a href='./print_log' onclick="%s%s%s%s%s%s %s%s%s%s%s %s">log</a>""",
                u"""<a href='./print_runs' onclick="%s%s%s%s%s%s %s%s%s%s%s %s">run session</a>"""]
   array = _get_array_from_columns(search_arg, u'RUN', columns, link_logs, 
                                   [u'PROJECT',      u'BINARY_NAME',
                                    u'CASE_NAME',    u'USER',
                                    u'DATE',         u'MACHINE',
                                    u'MD5CONF',      u'MD5DIFF',
                                    u'OPT_SET_NAME', u'MPI_PROC_NBR',
                                    u'THREAD_NBR',   u'REVISION'], myfilter)
   page_html += _array_to_table(array)
   
   page_html += _footer
   
   return page_html


def index(req):
   page_html = _index(20)
   return page_html

def pastix(req):
   page_html = _index(20, "PaStiX")
   return page_html
def scotch(req):
   page_html = _index(20, "Scotch")
   return page_html
def hips(req):
   page_html = _index(20, "Hips")
   return page_html

    

def union(a, b):
    """ return the union of two lists """
    return list(set(a) | set(b))

def _get_array_from_columns(search_arg, table, columns, link, link_data, myfilter):

    database = MySQLdb.connect (host   = _database_data[u"host"  ],
                                user   = _database_data[u"user"  ],
                                passwd = _database_data[u"passwd"],
                                db     = _database_data[u"db"    ])

    
    cmd        = u"SELECT "
    my_columns = union(columns, link_data)
    cmd       += u"`%s`" %(my_columns.pop())
    for c in my_columns:
       cmd    += u", `%s`" %(c) 
    cmd       += u" FROM `%s` WHERE "
    search_str = u"`%s` LIKE '%s'"
    my_search_arg = deepcopy(search_arg)
    key, value = my_search_arg.popitem()
    cmd += search_str % (key, str(value))
    for key in my_search_arg.keys():
        cmd += u' AND '
        cmd += search_str % (key, str(search_arg[key]))
    if (myfilter == u'error'):
        cmd += u' AND `ERRORS` > 0 '
    if (myfilter == u'warning'):
        cmd += u' AND `WARNINGS` > 0 '
    cmd += u" ORDER BY `DATE` DESC LIMIT 0 , 10000;"
    

    cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
    cursor.execute(cmd % table)
    
    array = []
    tmp = []
    for key in columns:
        tmp.append(key)    
    tmp.append(u'more')
    tmp.append(u'same runs')
    array.append(tmp)  

   
    row = cursor.fetchone()

    while (row):
        tmp = []
        for key in columns:
            tmp.append(str(row[key]))
        liste = []
        for data in link_data:
            liste.append(_maj_session %(data, str(row[data])))
        tmp.append(link[0] % tuple(liste))
        tmp.append(link[1] % tuple(liste))
        array.append(tmp)
        row = cursor.fetchone()
    
    database.close()

    return array

def _array_to_table(array):
    try: 
       error_col = array[0].index(u"ERRORS")
    except ValueError:
       error_col = -1

    try:
       warn_col = array[0].index(u"WARNINGS")
    except ValueError:
       warn_col = -1
    if (warn_col != -1 or error_col != -1):
       string = u'<tr class="%s">'
    else:
       string = u'<tr>'

    string1 = u'<tr>'
    for col in array[0]:
       string += u'<td>%s</td>'
       string1 += u'<th>%s</th>'
    string += u'</tr>'
    string1 += u'</tr>'
    output = ""
    output += _table_header
    output += string1 % tuple(array.pop(0))
    output += _table_header_end
    for row in array:
       if (warn_col != -1 or error_col != -1):
          if error_col != -1 and (row[error_col] == 'None' or 
                                  int(row[error_col]) != 0):
             row.insert(0, u"error")
          elif warn_col != -1 and int(row[warn_col]) != 0:
             row.insert(0, u"warning")
          else:
             row.insert(0, u"ok")

       output += string % tuple(row)
    output += _footer

    return output

def print_builds(req):
    #return '%s %s %s' % (revision, machine, date)
    session = Session.Session(req)
    project  = session[u'PROJECT']
    revision = session[u'REVISION']
    machine  = session[u'MACHINE']
    date     = session[u'DATE']
    myfilter = 'all'
    if session.has_key(u'BUILD_COLUMNS'):
       columns = session[u'BUILD_COLUMNS']
    else:
       columns = _build_columns
 
    search_arg = {
        u'REVISION' : revision,
        u'MACHINE'  : machine,
        u'DATE'     : date
        }
    link_logs = [u"""<a href='./print_build' onclick="%s%s%s%s%s%s%s%s">more info</a>""",
                 u"""<a href='./same_build' onclick="%s%s%s%s%s%s%s%s">same builds</a>"""]
    page_html = _header  

    page_html += u"""
<h1>Builds of %s (r%s) on %s the %s </h1>
""" %(project, revision, machine, date)
    page_html += _table_header
    page_html += _table_header_end
    page_html += u"""<tr><td><a href=../index.py>Go back to global page</a></td>"""
    page_html += u"""<td><a href=chooseBuildColumns>Choose columns</a></td></tr>"""
    page_html += _table_footer

    array = _get_array_from_columns(search_arg, u'BUILD', columns, link_logs, 
                                    [u'MACHINE', u'USER', u'REVISION', u'DATE',  u'MD5DIFF',
                                     u'MD5CONF', u'OBJECT_NAME', 'NAME'], myfilter)

    page_html += _array_to_table(array)
    #page_html += array
    page_html += _footer
    return page_html
    

def print_runs(req):
   session = Session.Session(req)
   project  = session[u'PROJECT']
   revision = session[u'REVISION']
   machine  = session[u'MACHINE']
   date     = session[u'DATE']
   myfilter = session[u'FILTER']
   if session.has_key(u'RUN_COLUMNS'):
      columns = session[u'RUN_COLUMNS']
   else:
      columns = _run_columns

   search_arg = {
      u'REVISION' : revision,
      u'MACHINE'  : machine,
      u'DATE'     : date
      }
   page_html = _header 
   page_html += u"""
<h1>Runs of %s (r%s) on %s the %s </h1>
""" %(project, revision, machine, date)
   page_html += _table_header
   page_html += _table_header_end
   page_html += u"""<tr><td><a href=../index.py>Go back to global page</a></td>"""
   page_html += u"""<td><a href=chooseRunColumns>Choose columns</a></td></tr>"""
   page_html += _table_footer

   link_logs = [u"""<a href='./print_log' onclick="%s%s%s%s%s%s %s%s%s%s%s %s">log</a>""",
                u"""<a href='./same_run' onclick="%s%s%s%s%s%s %s%s%s%s%s %s">same runs</a>"""]
   array = _get_array_from_columns(search_arg, u'RUN', columns, link_logs, 
                                   [u'PROJECT',      u'BINARY_NAME',
                                    u'CASE_NAME',    u'USER',
                                    u'DATE',         u'MACHINE',
                                    u'MD5CONF',      u'MD5DIFF',
                                    u'OPT_SET_NAME', u'MPI_PROC_NBR',
                                    u'THREAD_NBR',   u'REVISION'], myfilter)
   page_html += _array_to_table(array)
   
   page_html += _footer
   
   return page_html

def print_log(req):
   session = Session.Session(req)
   PROJECT      = session[u'PROJECT']
   BINARY_NAME  = session[u'BINARY_NAME']
   CASE_NAME    = session[u'CASE_NAME']   
   USER         = session[u'USER']        
   DATE         = session[u'DATE']        
   MACHINE      = session[u'MACHINE']     
   MD5CONF      = session[u'MD5CONF']     
   MD5DIFF      = session[u'MD5DIFF']     
   OPT_SET_NAME = session[u'OPT_SET_NAME']
   MPI_PROC_NBR = session[u'MPI_PROC_NBR']
   THREAD_NBR   = session[u'THREAD_NBR']  
   REVISION     = session[u'REVISION']    

   page_html  = _header_no_sort 
   cmd =  u"SELECT `LOG` from RUN WHERE PROJECT='%s' AND BINARY_NAME='%s' AND CASE_NAME='%s' "
   cmd += u"AND USER='%s' AND `DATE`='%s' AND MACHINE='%s' "
   cmd += u"AND MD5CONF='%s' AND MD5DIFF='%s' AND OPT_SET_NAME='%s'"
   cmd += u"AND MPI_PROC_NBR='%s' AND THREAD_NBR='%s' AND REVISION='%s'"

   database = MySQLdb.connect (host   = _database_data[u"host"  ],
                               user   = _database_data[u"user"  ],
                               passwd = _database_data[u"passwd"],
                               db     = _database_data[u"db"    ])
   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd % (PROJECT, BINARY_NAME, CASE_NAME, USER, DATE, MACHINE, MD5CONF, MD5DIFF, 
                         OPT_SET_NAME, MPI_PROC_NBR, THREAD_NBR, REVISION))
   row = cursor.fetchone()
   #page_html += u'<code class="prettyprint lang-sh">%s</code>' % row['LOG'].replace(u"\n", u"\n<br>").replace(u' ', u'&nbsp;')
   if (row['LOG']):
      page_html += u'<code>%s</code>' % row['LOG'].replace(u"\n", u"\n<br>").replace(u' ', u'&nbsp;')
   else:
      page_html += u"<code>Empty log...</code>"
   page_html += _footer
   database.close()
   return page_html

def print_build(req, filename=None, key_grep=None):
   session = Session.Session(req)

   project     = session[u'PROJECT']
   MACHINE     = session[u'MACHINE']
   USER        = session[u'USER']
   REVISION    = session[u'REVISION']
   DATE        = session[u'DATE']
   MD5DIFF     = session[u'MD5DIFF']
   MD5CONF     = session[u'MD5CONF']
   OBJECT_NAME = session[u'OBJECT_NAME']

   remote = [u'LOG',     u'CONFIG_FILE', u'DIFF']
   warn   = [u'LOG']
   hide   = [u'MD5DIFF', u'MD5CONF']
   database = MySQLdb.connect (host   = _database_data[u"host"  ],
                               user   = _database_data[u"user"  ],
                               passwd = _database_data[u"passwd"],
                               db     = _database_data[u"db"    ])

   page_html  = _header_no_sort

   page_html += """
<h1>Build of %s [%s (r%s)] from %s on %s the %s </h1>
""" %(OBJECT_NAME, project, REVISION, USER, MACHINE, DATE)

   cmd =  u"SELECT * from BUILD WHERE MACHINE='%s' AND USER='%s' "
   cmd += u"AND REVISION='%s' AND `DATE`='%s' AND MD5DIFF='%s' "
   cmd += u"AND MD5CONF='%s' AND OBJECT_NAME='%s'"
   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd % (MACHINE, USER, REVISION, DATE, MD5DIFF, MD5CONF, OBJECT_NAME))
   row = cursor.fetchone()
   
   if filename == None:
      page_html += _table_header
      page_html += _table_header_end
        
      row_html      = u'<tr><td>%s</td><td>%s</td></tr>'
      row_html_href = u'<tr><td>%s</td><td><a href=\'./print_build?filename=%s\'>show</a></td></tr>'
      row_html_href_warn = u'<tr><td>%s</td><td>[<a href=\'./print_build?filename=%s\'>log</a>] [<a href=\'./print_build?filename=%s&key_grep=warning\'>warnings</a>] [<a href=\'./print_build?filename=%s&key_grep=error\'>errors</a>]</td></tr>'
      for key in row.keys():
         if key not in hide:
            if key in remote:
               if key in warn:
                  page_html += row_html_href_warn % (key, key, key, key)
               else:
                  page_html += row_html_href % (key, key)
            else:
               if not 'None' == str(row[key]):
                  page_html += row_html % (key, row[key])
         cursor.close()
      page_html += _table_footer

   else:
      if key_grep == None:
         output = u''
         for line in row[filename].decode(u'utf8', u'replace').splitlines():
            output += '%s\n' % line.rstrip() 
      else:
         output = grep(u'%s:.*' %key_grep,  row[filename].decode(u'utf8', u'replace').splitlines())
         #page_html += u'<pre class="prettyprint lang-sh">%s</pre>' % output #.replace(u"\n", u"\n<br>").replace(u' ', u'&nbsp;')
      page_html += u'<pre>%s</pre>' % output #.replace(u"\n", u"\n<br>").replace(u' ', u'&nbsp;')
   page_html += _footer
   return page_html

def handler(req):

   session = Session.Session(req)

   try:
      session[u'HITS'] += 1
   except:
      session[u'HITS'] = 1

   session.save()
   
   req.content_type = u'text/plain'
   req.write(u'Hits: %d\n' % session[u'HITS'])
   return apache.OK

def _chooseColumn(req, session_key, reference_list, page, name, table, page_ref):
    session = Session.Session(req)
    if not session.has_key(session_key):
       session[session_key] = deepcopy(reference_list)
    session.save()
    page_html  = _header 
    page_html += u"""
<script language="Javascript" type="text/javascript" > 
   function choix(formulaire) 
     { 
       var j; 
       var i = formulaire.boite1.selectedIndex; 
       if (i != 0) {
          lancerRequete('majSession?KEY=SELECTED&VALUE='.concat(formulaire.boite1.options[i].text),'GET');
       }
     }
 
</script> 
"""
    page_html += u"""
<h1>Modify %ss page's columns</h1>
<h2>Remove columns</h2>
<a href='%s'>Go back to %ss page</a><br>""" % (name, page, name)
    database = MySQLdb.connect (host   = _database_data[u"host"  ],
                                user   = _database_data[u"user"  ],
                                passwd = _database_data[u"passwd"],
                                db     = _database_data[u"db"    ])

    cmd = u'DESCRIBE %s' % table
    cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
    cursor.execute(cmd) 
  
    page_html += u"""<form name="formulaire"> 
<select name="boite1" onChange="choix(this.form);location.reload()">\n"""
    if not session.has_key('SELECTED'):
       page_html += u"""<option selected>...........Choose a column name...........</option>\n"""
    else:
       page_html += u"""<option>...........Choose a column name...........</option>\n"""

    line = u"""<option%s>%s</option>"""
    while (True):
       row = cursor.fetchone()
       if row == None : break
       column = row[u'Field']
       selected = u''
       if session.has_key(u'SELECTED') and session[u'SELECTED'] == column:
          selected = u' selected'
       page_html += line % (selected, column)

    page_html += u"""</select></form>"""

    page_html += _table_header
    page_html += _table_header_end
    index = 0
    for column in session[session_key]:
       rm_column  = _rm_column % (session_key, index)
       if session.has_key(u'SELECTED'):
          add_before = _add_column %(session_key, session[u'SELECTED'], index)
          add_after  = _add_column %(session_key, session[u'SELECTED'], index+1)
       else:
          add_before = _add_column %(session_key, 0, index)
          add_after  = _add_column %(session_key, 0, index+1)
       page_html += u'''
<tr>
  <td>%s</td>
  <td><a href='%s' onclick=%s>remove</a></td>
  <td><a href='%s' onclick=%s>add before</a></td>
  <td><a href='%s' onclick=%s>add after</a></td>
</tr>''' % (column, page_ref, rm_column, page_ref, add_before, page_ref, add_after)
       index = index +1
    page_html += _table_footer

   
    page_html += _footer
    return page_html

def chooseBuildColumns(req):   
   return _chooseColumn(req, u'BUILD_COLUMNS', _build_columns, u'print_builds', u'build', u'BUILD', u'chooseBuildColumns')

def chooseSameBuildColumns(req):   
   return _chooseColumn(req, u'BUILD_COLUMNS', _same_build_columns, u'same_build', u'build', u'BUILD', u'chooseSameBuildColumns')

def chooseRunColumns(req):   
   return _chooseColumn(req, u'RUN_COLUMNS', _run_columns, u'print_runs', u'run', u'RUN', u'chooseRunColumns')

def chooseSameRunColumns(req):   
   return _chooseColumn(req, u'SAME_RUN_COLUMNS', _same_run_columns, u'same_run', u'run', u'RUN', u'chooseSameRunColumns')



def majSession(req):
   session = Session.Session(req)
   GET = mod_python.util.parse_qs(req.args)
   session[GET[u'KEY'][0]] = GET[u'VALUE'][0]
   session.save()
   print session[GET[u'KEY'][0]]
   return str(GET)

def rmColumn(req):
   session = Session.Session(req)
   GET = mod_python.util.parse_qs(req.args)
   session[GET[u'KEY'][0]].pop(int(GET[u'INDEX'][0]))
   session.save()
   print session[GET[u'KEY'][0]]
   return str(GET)

def addColumn(req):
   session = Session.Session(req)
   GET = mod_python.util.parse_qs(req.args)
   session[GET[u'KEY'][0]].insert(int(GET[u'INDEX'][0]), GET[u'COLUMN'][0])
   session.save()
   print session[GET[u'KEY'][0]]
   return str(GET)


def sendMail(req):
   # Import smtplib for the actual sending function
   import smtplib
   import MimeWriter  
   import mimetools
   import StringIO  
   encoding = "base64"
   charset = "utf8"
 
   sender  = u'regression@python.fr'
   to      = u'xavier.lacoste@inria.fr'
   #declaration des buffers
   out = StringIO.StringIO() 
   html = _index(10).replace("<STYLE>", 
                             "<STYLE>%s" % open(u'/home/pastix/pastix-user/ricar/Scripts/regression/www/pastix.css').read())
   htmlin = StringIO.StringIO(html)
   txtin = StringIO.StringIO("Ne fonctionne qu'en HTML")
   
   #declaration et initialisation du writer
   writer = MimeWriter.MimeWriter(out)
   writer.addheader("Subject", "Ici le sujet du message")
   writer.addheader("MIME-Version", "1.0")
   writer.startmultipartbody("alternative")
   writer.flushheaders()
 
   #ajout de la partie text
   textPart = writer.nextpart()
   textPart.addheader("Content-Transfer-Encoding", encoding)
   pout = textPart.startbody("text/plain", [("charset", charset)])
   mimetools.encode(txtin, pout, encoding)
   txtin.close()
 
   #On ajoute la partie html
   htmlPart = writer.nextpart()
   htmlPart.addheader("Content-Transfer-Encoding", encoding)
   pout = htmlPart.startbody("text/html", [("charset", charset)])
   mimetools.encode(htmlin, pout, encoding)
   htmlin.close()
   
   #on clot le mail
   writer.lastpart()
   mail = out.getvalue()
   out.close()
   smtp = smtplib.SMTP()
   smtp.connect()
   smtp.sendmail(sender, [to], mail)
   smtp.close()
   

