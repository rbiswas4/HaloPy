from numpy import *
from sys import *
from scipy import *
from params import *

def read_input_params (infile):
   inputparams= tuple(open(infile,"r"))
   return inputparams
  #
  
def readhalofiles(halofileroot, filenum):
   FOFcount= []
   SOcount=  []
   for i in range (0,filenum):
      infile= halofileroot+".sodproperties."+str(i)
      for line in open(infile,'r'):
           data = line.split(' ')
           #print data
           #exit()
           FOFcount.append(int(data[fofcount_col]))
           SOcount.append(int(data[socount_col]))              
   FOFcount = array(FOFcount)
   SOcount = array(SOcount)
   return FOFcount, SOcount

def get_bin_mean(a,bin_edges, massbins):
   mean_mass= []
   for i in range (0, massbins-1):
      b_start= bin_edges[i]
      b_end= bin_edges[i+1]
      ind_upper = nonzero(a >= b_start)[0]
      a_upper = a[ind_upper]
      a_range = a_upper[nonzero(a_upper < b_end)[0]]
      if(len(a_range)!=0): meana= mean(a_range)
      else: meana=0.5*(b_start+b_end)
      mean_mass.append(meana)
   return array(mean_mass)
     

def logrange(start, end, binnum):
   b_edge=[]
   binsize=log(end/start)/(binnum-1)
   for i in range(0,binnum):
      b_edge.append(start*exp(i*binsize))     
   return array(b_edge), binsize

def MF_fit(x, z):
  A=0.333/(1+z)**0.11
  a1= 0.788/(1+z)**0.01
  p= 0.807
  q= 1.795
  delc=1.686
  return A*(2*a1/pi)**0.5*exp(-a1*delc**2/2/x**2)*(1+(x**2/a1/delc**2)**p)*(delc/x)**q

def calc_mf(b_start, b_end, massbins, mass):
    bins, binsize= logrange(b_start, b_end, massbins) 
    hist,bin_edges=histogram(mass,bins=bins)
    meanmass= get_bin_mean(mass, bin_edges, massbins)
    return meanmass, hist, binsize, bin_edges
