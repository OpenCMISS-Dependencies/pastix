# Import smtplib for the actual sending function
import smtplib
import Database.database
import MySQLdb
import MySQLdb.cursors
import re
_database_data = {
    u"host"   : u"localhost",
    u"user"   : u"pastix",
    u"passwd" : u"PaStiX",
    u"db"     : u"pastix-python"}
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

_footer = u"""</body></html>"""
def grep(string,list):
   expr = re.compile(string)
   line = 0
   output = u''
   for text in list:
      line = line + 1
      match = expr.search(text)
      if match != None:
         try:
            output += u'%s\n' % (match.string)
         except:
            print match.string
            raise
         
   return output


def _index(numberOfresults):
   page_html = _header.replace(u"../", u"")
   page_html += u"""
<h1>Last tests</h1>
<h2>Last builds</h2>
"""
   search_warning = re.compile(u"warning\s*:",
                               flags=re.IGNORECASE|re.MULTILINE)
   search_remark  = re.compile(u"remark\s*:",
                               flags=re.IGNORECASE|re.MULTILINE)
   search_error   = re.compile(u"error\s*:",
                               flags=re.IGNORECASE|re.MULTILINE)
   database = MySQLdb.connect (host   = _database_data[u"host"  ],
                               user   = _database_data[u"user"  ],
                               passwd = _database_data[u"passwd"],
                               db     = _database_data[u"db"    ])
   cmd = u"SELECT `DATE` FROM BUILD ORDER BY `DATE` DESC LIMIT 0 , 1;"
   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd)
   row = cursor.fetchone()
   date = str(row['DATE'])

   cmd = u"SELECT DISTINCT `PROJECT`, `REVISION`, `MACHINE`"
   cmd += u" FROM `BUILD` WHERE `DATE` like \"%s\";" % date

   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
   cursor.execute(cmd)
   row = cursor.fetchone()
   page_html += _table_header
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
   page_html += _table_header_end
   row_html = u"""<tr class="%s">
   <td>%s</td>
   <td>%s</td>
   <td>%s</td>
   <td>%s</td>
   <td>%i</td>
   <td>%i</td>
   <td>%i</td>
</tr>"""
   errors = []
   while row :
      project  = row[u'PROJECT']
      revision = row[u'REVISION']
      machine  = row[u'MACHINE']
      cmd = """SELECT COUNT(*)      AS NB, 
                 SUM(BUILD.WARNINGS)  AS NB_WARNINGS,
                 SUM(BUILD.ERRORS)    AS NB_ERRORS
                 FROM `BUILD` WHERE `REVISION` LIKE '%s'
                 AND `PROJECT` LIKE '%s' 
                 AND `MACHINE` LIKE '%s'
                 AND `DATE` LIKE '%s'"""
      cmd = cmd % (revision, project,  machine, date)
      cursor2 = database.cursor()
      cursor2.execute(cmd)
      row2 = cursor2.fetchone()
      cursor2.close()
      liste = [project, revision, machine, date]
      for r in row2:
         liste.append(r)   
      if row2[2] > 0:
         liste.insert(0, u"error")
      elif row2[1] > 0:
         liste.insert(0, u"warning")
      else:
         liste.insert(0, u"ok")     
      page_html += row_html % tuple(liste)
      row = cursor.fetchone()
      cmd = u"""select LOG, NAME FROM `BUILD` WHERE `REVISION` LIKE '%s'
                 AND `PROJECT` LIKE '%s' 
                AND `MACHINE` LIKE '%s'
                AND `DATE` LIKE '%s'
                AND ( `WARNINGS` NOT LIKE 0 
                      OR `ERRORS` NOT LIKE 0)"""
      cmd = cmd % (revision, project, machine, date)
      cursor2 = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)
      cursor2.execute(cmd)
      row2 = cursor2.fetchone()
      while row2:
         l = grep(u'error:', row2[u'LOG'].decode(u'utf8', u'ignore').splitlines())#//search_error.findall(row2[0])
         for s in l.splitlines(): 
            errors.append([project, machine, row2[u'NAME'], s])
         l = grep(u'warning:', row2[u'LOG'].decode(u'utf8', u'ignore').splitlines()) #search_warning.findall(row2[0])
         for s in l.splitlines(): 
            errors.append([project, machine, row2[u'NAME'], s])
         row2 = cursor2.fetchone()
         
   cursor.close ()
   page_html += _table_footer

      
   cmd = u"SELECT DISTINCT PROJECT, REVISION, MACHINE"
   cmd += u" FROM RUN WHERE `DATE` like \"%s\";" % (date)

   cursor = database.cursor(cursorclass=MySQLdb.cursors.DictCursor)

   cursor.execute(cmd)
   row = cursor.fetchone()
   page_html += _table_header
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
   page_html += _table_header_end
   row_html = u"""<tr class="%s">
   <td>%s</td>
   <td>%s</td>
   <td>%s</td>
   <td>%s</td>
   <td>%i</td>
   <td>%i</td>
   <td>%i</td>
</tr>"""
   while row :
      project  = row[u'PROJECT']
      revision = row[u'REVISION']
      machine  = row[u'MACHINE']
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
      liste = [project, revision, machine, date]
      for r in row2:
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
   page_html += _footer 

   save = [u'',u'',u'']
   page_html += u'<h3> Bulding errors </h3>\n'
   for e in errors:
      if (save[0] != e[0] or
          save[1] != e[1] or
          save[2] != e[2]):
         if e[0] != '':
            page_html +=u'</pre><h4>%s on %s with %s :</h4><pre>\n' % (e[0], e[1], e[2])
         else:
            page_html +=u'<h4>%s on %s with %s :</h4><pre>\n' % (e[0], e[1], e[2])
         save[0] = e[0]
         save[1] = e[1]
         save[2] = e[2]
         
      page_html += u'%s\n' % e[3]

   page_html+= u'</pre>'
   database.close()
   return page_html

def sendMail():
   # Import smtplib for the actual sending function
   import smtplib
   import MimeWriter  
   import mimetools
   import StringIO  
   encoding = "base64"
   charset = "utf8"
 
   sender  = u'regression@python.fr'
   to      = [u'pastix-log@localhost']
   #declaration des buffers
   out = StringIO.StringIO() 
   html = _index(10).replace("<STYLE>", 
                             "<STYLE>%s" % open(u'/home/pastix/pastix-user/ricar/Scripts/regression/www/pastix.css').read())
   htmlin = StringIO.StringIO(html)
   txtin = StringIO.StringIO("Ne fonctionne qu'en HTML")
   
   #declaration et initialisation du writer
   writer = MimeWriter.MimeWriter(out)
   writer.addheader("Subject", "Regression")
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

   smtp.sendmail(sender, to, mail)
   smtp.close()
   
sendMail()
