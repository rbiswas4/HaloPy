import MySQLdb as mdb
import sys


def createcMdb(tablename,mass, concentration, z, case):
   if z==0.0: strz='0'
   if z==1.0: strz='1'
   if z==2.0: strz='2'
   var=tablename+strz
   length=len(mass)
   try:
      con = mdb.connect('localhost', 'mysqlu','','testcm');
      cur = con.cursor()
      cur.execute('drop table if exists '+var)
      if case=='cMdata': cur.execute( 'create table '+var+' (Id INT PRIMARY KEY AUTO_INCREMENT, Mass float, conc float)')
      if case=='fitdata': cur.execute( 'create table '+var+' (Id INT PRIMARY KEY AUTO_INCREMENT, amplitude float, indx float)')     
      for i in range(0,length):
         var1= mass[i]
         var2= concentration[i]
         if case=='cMdata':cur.execute("INSERT INTO "+var+" (Mass, conc) VALUES('%f','%f')" %(var1,var2))
         if case=='fitdata':cur.execute("INSERT INTO "+var+" (Amplitude, indx) VALUES('%f','%f')" %(var1,var2))
   except mdb.Error, e:
       print "Error %d: %s" % (e.args[0],e.args[1])
       sys.exit(1)
   finally:
      if con:
         con.close()

#if __name__ == "__main__":
#   createcMdb(1.0, 1.0, 0)
