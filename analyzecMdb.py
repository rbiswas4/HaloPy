import MySQLdb as mdb
import sys
import numpy as np

import matplotlib.pylab
from pylab import *
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib.mlab as mlab
from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import AutoLocator
from matplotlib.widgets import Slider, Button, RadioButtons
import discslider  as ds

def findmasslimit(z):
   if z==0.0: strz='0'
   if z==1.0: strz='1'
   if z==2.0: strz='2'
   var='conc_'+strz
   try:
      con = mdb.connect('localhost', 'mysqlu','','testcm');
      cur = con.cursor()
      cur.execute("select min(Mass), max(Mass) from "+var)
      row= cur.fetchone()
   except mdb.Error, e:
       print "Error %d: %s" % (e.args[0],e.args[1])
       sys.exit(1)
   finally:
      if con:
         con.close() 
   return row[0], row[1]
   
def findconc(lowmass, highmass, z):
   if z==0.0: strz='0'
   if z==1.0: strz='1'
   if z==2.0: strz='2'
   var='conc_'+strz
   M,c=[],[]
   try:
      con = mdb.connect('localhost', 'mysqlu','','testcm');
      cur = con.cursor()
      cur.execute("select Mass, conc from "+var+" where Mass>%le and Mass <%le order by conc" %(lowmass,highmass))
      numrows= int(cur.rowcount)
      for i in range(numrows):
         row= cur.fetchone()
         M.append(row[0])
         c.append(row[1])
         #print "%le %lf" %(row[0], row[1])   
   except mdb.Error, e:
       print "Error %d: %s" % (e.args[0],e.args[1])
       sys.exit(1)
   finally:
      if con:
         con.close()
   return np.array(M), np.array(c)

def comp_hist(c):
   """rc('font', family='sans serif', style='normal', variant='normal', stretch='normal', size='10.0')
   lab_fontsize = 9
   axes_fontsize = 12
   lineM = ['-k','--r',':b','-m','-g','-y']
   lineL = [':k',':r',':b',':m',':g',':y']
   lw1 = 2
   lw2 = 1
   majorLocatorx   =  MaxNLocator(8)
   majorLocatory   = MaxNLocator(10)
   minorLocator   = AutoMinorLocator() #MaxNLocator(50) #MultipleLocator(2)
   minorLocatorx   = AutoMinorLocator() #MultipleLocator(2)
   subplots_adjust(hspace=0., wspace=0.)"""
   ax = subplot(    1,     1,    1)
   ax.clear()
   n, bins, patches = ax.hist(c, 10, normed=1, facecolor='green', alpha=0.75)
   bincenters=0.5*(bins[1:]+bins[:-1])  
   meanc= np.mean(c)
   stdc= np.var(c)**0.5
   nsamp=len(c)
   y = mlab.normpdf( bincenters, meanc, stdc)
   plot(bincenters, y, '--r', linewidth=2,label="# samples= %d \n mean=%lf \nstd=%lf"%(nsamp,meanc,stdc))
   xlabel(r'$\rm conc$', fontsize=20)
   ylabel(r'P($\rm conc$)', fontsize=20)
   ax.legend(loc="upper right", #bbox_to_anchor=[0.7, 1.0],
           ncol=2, shadow=False)
   return bincenters,n

def print_hist(bincenters,n):
   print " histogram bin c vs counts"
   for i in range(len(bincenters)):
          print bincenters[i], n[i]

def update(val):
    lowmass = slowmass.val
    highmass= sbinsize.val+lowmass
    z= floor(sz.val)
    M,c=findconc(lowmass,highmass,z)  
    bincenters,n=comp_hist(c)  
    print "num of data= ", len(M) 
    print "range of mass from ",lowmass, " to ", highmass," at redshift= ",z  
    print "mean mass= %le std= %le mean c= %lf std= %lf" %(np.mean(M), np.var(M)**0.5, np.mean(c), np.var(c)**0.5)
    draw()  

def reset(event):
    slowmass.reset()
    sbinsize.reset()
    sz.reset()
 
def updateval():
    slowmass.on_changed(update)
    sbinsize.on_changed(update)
    sz.on_changed(update)
 
  #button.on_clicked(reset)
  #radio.on_clicked(update)

if __name__ == "__main__":
    """if len(sys.argv) !=4:
      print "usage: python analyzecMdb.py low-mass high-mass redshift"
      exit()
    lowmass=float (sys.argv[1])
    highmass=float(sys.argv[2])
    z=float(sys.argv[3])
    M,c=findconc(lowmass,highmass,z)
    
    bincenters,n=comp_hist(c)
    print_hist(bincenters,n)
    """   
    minmass, maxmass= findmasslimit(0)
    print "min and max mass", minmass, maxmass
    print " showing histogram plot..."
    ax = subplot(111)
    subplots_adjust(left=0.1, bottom=0.25)
    axcolor = 'lightgoldenrodyellow'
    f0=(maxmass-minmass)*0.1
    fz=1
    valfmt='%le'
    valinit='%le'
    #resetax = axes([0.025, 0.5, 0.15, 0.15])
    #button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')
    #rax = axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
    #radio = RadioButtons(rax, ('z=0':0, 'z=1':1, 'z=2':2), active=0)
    axlowmass = axes([0.12, 0.03, 0.65, 0.03], axisbg=axcolor)
    axbinsize  = axes([0.12, 0.08, 0.65, 0.03], axisbg=axcolor)
    axz= axes([0.12, 0.13, 0.65, 0.03], axisbg=axcolor)
    slowmass = ds.DiscreteSlider(axlowmass, 'Masslow', minmass, 0.5*maxmass, increment=f0,valinit= f0)
    sbinsize = ds.DiscreteSlider(axbinsize, 'binsize', minmass, 0.5*maxmass, increment= f0, valinit= f0)
    sz= ds.DiscreteSlider(axz,'redshift',0,2.5,increment=fz,valinit=fz)   

    updateval()
    show()

