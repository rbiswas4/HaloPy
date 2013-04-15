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
from scipy.optimize import curve_fit
import createcMdb as cMdb

def meanconc(num,z):
   if z==0.0: strz='0'
   if z==1.0: strz='1'
   model= ['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30','31','32','33','34','35','36','37']
   ax = subplot(    1,     1,    1)
   ax.clear()
   M,c=[],[]
   state="f"
   path="/Users/Suman/backup/research/concentration/"

   file= path+"wcdm_output_180_0619/cM.180."+state+"."+model[num]+".z"+strz
   f=open(file,"r")
   data=np.loadtxt(f)
   close,f
   x,y,yerr=data[:,2],data[:,5], data[:,6]
   ax.errorbar(x,y,yerr)
   M.extend(data[:,2])
   c.extend(data[:,5])
 
   file= path+"wcdm_output_365_0619/cM.365."+state+"."+model[num]+".z"+strz
   f=open(file,"r")
   data=np.loadtxt(f)
   close,f
   x,y,yerr=data[:,2],data[:,5], data[:,6]
   ax.errorbar(x,y,yerr)
   M.extend(data[:,2])
   c.extend(data[:,5])
   
   file= path+"wcdm_output_1300_0619/cM.1300."+state+"."+model[num]+".z"+strz
   f=open(file,"r")
   data=np.loadtxt(f)
   close,f
   x,y,yerr=data[:,2],data[:,5], data[:,6]
   M.extend(data[:,2])
   c.extend(data[:,5])
   M=np.array(M)
   c=np.array(c)
   ampl,indx=cMfit(M,c)
   #print ampl,indx
   ax.errorbar(x,y,yerr)
   plot(M,cMfunc(M,ampl,indx),'--k', linewidth=2,label='%2.2lf[1.12$(M/(5.10^{13}))^{0.3}+0.53]^{-%0.2lf}$'%(ampl,indx))
   xscale('log')
   xlim(1e12,2e15)
   ylim(1,8)
   xlabel(r"M M$_\odot$/h", fontsize=18)
   ylabel(r"c",fontsize=18)
   ax.legend(loc="upper right", #bbox_to_anchor=[0.7, 1.0],
           ncol=2, shadow=False)
   #plt.title(r"model= "+num, fontsize=18)
   #show()
   return ampl, indx


def print_cM(M,c,err):
   print " mean c-M"
   for i in range(len(c)):
          print M[i], n[i]

def cMfunc(m,ampl,ind):
  return ampl*(1.12*(m/5e13)**0.3+0.53)**-ind

def cMfit(M,c):
    ampl0, ind0=5.0, 0.3
    pfinal, covar=curve_fit(cMfunc, M, c, p0=[ampl0, ind0], maxfev=10000)
    return pfinal[0],pfinal[1]

def update(val):
    num = int(floor(snum.val))
    z= floor(sz.val)
    meanconc(num,z)
    print "model= ",num,"redshift= ",z  
    """print "num of data= ", len(M) 
    print "range of mass from ",lowmass, " to ", highmass," at redshift= ",z  
    print "mean mass= %le std= %le mean c= %lf std= %lf" %(np.mean(M), np.var(M)**0.5, np.mean(c), np.var(c)**0.5)"""
    draw()  

def reset(event):
    smodel.reset()
    sz.reset()
 
def updateval():
    snum.on_changed(update)
    sz.on_changed(update)
 
def createfitdb():
   ampl0,indx0,ampl1,indx1=[],[],[],[]
   for i in range (1,38):
       var,var1= meanconc(i,0.0)
       ampl0.append(var)
       indx0.append(var1)
       var,var1= meanconc(i,1.0)
       ampl1.append(var)
       indx1.append(var1)
    #for i in range (0,37): print i,ampl0[i],indx0[i],ampl1[i],indx1[i]
    #exit()
   cMdb.createcMdb('cmfit_', ampl0, indx0, 0.0,'fitdata')
   cMdb.createcMdb('cmfit_', ampl1, indx1, 1.0, 'fitdata')

if __name__ == "__main__": 
    #createfitdb()
    print " showing mean c-M plot..."
    ax = subplot(111)
    subplots_adjust(left=0.1, bottom=0.25)
    axcolor = 'lightgoldenrodyellow'
    fz=1
    valinit='%d'
    axnum  = axes([0.12, 0.05, 0.65, 0.03], axisbg=axcolor)
    axz= axes([0.12, 0.1, 0.65, 0.03], axisbg=axcolor)
    snum = ds.DiscreteSlider(axnum, 'model #', 1, 37.5, increment=fz,valinit= fz)
    sz= ds.DiscreteSlider(axz,'redshift',0,1.5,increment=fz,valinit=fz)   

    updateval()
    show()
