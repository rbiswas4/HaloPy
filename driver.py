#!/usr/bin/env python
#CHANGES:
#Introduced two extra arguments to be able to write to an outputdirectory 
#and not do num particle corrections.
#introduced cb transfer function
import sys
import os
from MF import *
from cM import *
from fsigma import *
from growth_ode import *

from growth_rahulv2 import *
from astropy.cosmology import set_current
from astropy.cosmology import get_current
from astropy.cosmology import w0waCDM 
import argparse

import time

##############################################################################
		#Parsing input
##############################################################################
parser = argparse.ArgumentParser(prog="driver.py" ,
	formatter_class=argparse.RawTextHelpFormatter,
	description = """
	Code to do FoF/SOD mass functions and c-M relations. 
	USAGE:
	driver.py inputfile --correctnumparticles [True] --outdir [$PWD] 
	OUTPUT:  
	txt output files are created in dirname specified at commandline, 
	default is $PWD""")

parser.add_argument("--minparticles", 
	type = int ,
	default  = 400 ,
	help = "Minimum number of particles in a halo required to be counted" )
parser.add_argument(dest = "inputfile" ,
	type = str ,
	help = "relative or absolute path to inputfile")
parser.add_argument("--outdir" ,
	type = str ,
	default = './' ,
	help = "relative path to output directory")	
parser.add_argument('--maphalomasses',dest='maphalomasses',action='store_true')
parser.add_argument('--no-maphalomasses',dest='maphalomasses',action='store_false')
parser.set_defaults(maphalomasses=True)

args = parser.parse_args()
print sys.argv
infile  = args.inputfile 
# args.correctnumparticles 
outdir  = args.outdir
if not os.path.exists(outdir) :
	os.makedirs(outdir) 
outdir =  outdir + "/"
rhoc= 2.77536627e11; #Msun.h^2/Mpc^3

tin= time.time()
#infile= sys.argv[1]
inputparams= read_input_params(infile) # input parameter file
print "read input params file ", infile ,"\n"

case= int(inputparams[0])
massbins= int(inputparams[1])
filenum= int(inputparams[2])
zseek= float(inputparams[3])
if int(inputparams[4])==1:
    model_no=int(inputparams[5])
    cosmoparams= read_input_params("model_param.100.txt")  # read cosmo param from model_param.txt, usually useful for suite of sims like coyote, supercoyote
    tmp=cosmoparams[model_no].split()
    hubble= float(tmp[4])
    Omegam=float(tmp[1])/hubble**2
    ns= float(tmp[5])
    sigma8= float(tmp[3])
    w0= float(tmp[6])
    wa= float(tmp[7])
    Omeganu=float(tmp[8])/hubble**2
    Omegab=float(tmp[2])/hubble**2  
else:
    model_no=0
    hubble= float(inputparams[6])
    Omegam= float(inputparams[7])
    ns= float(inputparams[8])
    sigma8= float(inputparams[9])
    w0= float(inputparams[10])
    wa= float(inputparams[11])
    Omeganu= float(inputparams[12])
    Omegab= float(inputparams[13])

print hubble, Omegam, ns, sigma8, w0, wa
Omegac = Omegam - Omegab - Omeganu  
print Omegac, Omegab, Omeganu

print "VALS", inputparams[15], inputparams[16]
boxsize= float(inputparams[14])
particlenum= int(inputparams[15])
minparticle= int(inputparams[16])
Delta= float(inputparams[17])
fileroot= inputparams[18].rstrip()
halofile= inputparams[19].rstrip()
Tkfile= inputparams[20].rstrip()

	#massres is particle mass
massres= Omegam*rhoc*(boxsize/particlenum)**3 #Msun/h
simvol= (boxsize)**3 #Mpc^3/h^3
print 'mass res= %le' % (massres), 'Msun/h'    
print fileroot, halofile, Tkfile    
#k,Tk= read_Tk(Tkfile)
from camb_utils import cambio as cio 
k, Tk = cio.cbtransfer(Tkfile, Omegac, Omegab ) 
N= Pk_norm(sigma8, ns, k, Tk, hubble)

if Omeganu>0: Neff=3.04
else: Neff=0
avals, D , info  = growth(Omegam =Omegam, w0 =w0, wa =wa, Tcmb =2.725, h =hubble, Neff  = Neff, omeganuh2val = Omeganu*hubble**2 ) 

#scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0, 0.0) #wa(last parameter)=0.0 here, if wa varies, then modify lines around 21-35 to read in the value either from model_param.txt or from input file
 
Ds, lDs= interp_D(zseek,avals,D[:,0], D[:,1])

print "D(z) and log D/log z(z)= ", Ds, lDs, 'at z= ',zseek
print "time takes= ",time.time()-tin

if case == 0 or case == 2 or case==4:
        print "\n################## FOF MF ##########################################\n"
        fofproperties= inputparams[21].rstrip()
        halofileroot=fileroot+'/'+halofile
        fofcount_col= readheaders(halofileroot, fofproperties,'fof_halo_count')
        FOFcount= readhalofiles(halofileroot, fofproperties, fofcount_col,filenum)
        
        print "\ntotal # of fof clusters= ",len(FOFcount)
        FOFbinstart=min(FOFcount)*massres
        FOFbinstart=400.*massres
        FOFbinend=max(FOFcount)*massres
	if args.maphalomasses :
		print " Doing halo mass maps \n"
        	FOFmass=FOFcount*(1.0-1.0/FOFcount**0.65)*massres
	else:
		print "Not Correcting for number of particles in halos\n"
        	FOFmass=FOFcount*massres
        print "minFOF maxFOF"
        print "# particles= ",min(FOFcount), max(FOFcount)
	print "mass Msun/h= ", '%le %le' %(FOFbinstart, FOFbinend)
        file= open(outdir + halofile+'.fof_mf','w')
        print >>file, "# FOF Mass Msun/h, # clusters, dn/dln M (h/Mpc)^3, frac err, 1/sigma(M), f(sigma), fsigma_fit \n"
        meanmass, hist, binsize, bin_edges= calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
        for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >0:  print >> file, '%le %5d %7.4le %7.4lf %7.4lf %7.4lf %7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek))
	file.close()

if case >= 1 and case<=5:
        print "\n##################### SO MF ########################################\n"
        tin= time.time()
        sodproperties= inputparams[22].rstrip()
        halofileroot=fileroot+'/'+halofile
        socount_col= readheaders(halofileroot, sodproperties,'sod_halo_count')
        SOcount= readhalofiles(halofileroot, sodproperties, socount_col, filenum)   
	print "\ntotal # of SO clusters= ",len(SOcount)
        SObinstart=max(min(SOcount), minparticle)*massres
        SObinend=max(SOcount)*massres
        SOmass=SOcount*massres
        SOradius= (3.0*SOmass/4.0/pi/rhoc/Delta)**(1.0/3.0)
	print "minSO maxSO"
        print "# particles= ",min(SOcount), max(SOcount)
	print "mass Msun/h= ", '%le %le' %(SObinstart, SObinend)
        file= open(outdir + halofile+'.sod_mf','w')
        print>> file, "# SO Mass Msun/h # clusters dn/dln M (h/Mpc)^3 frac err 1/sigma(M) f(sigma)\n"
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
	   sigmaM= Ds*sigmaM
           if hist[i] >2:  print>> file, '%le %5d %7.4le %7.4lf %7.4lf %7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM))
        file.close()
	print "time takes= ",time.time()-tin

if case >= 3:
        print "\n####################### c-M relation #################################\n"
        def fitz(z, var):
          if var ==0:
            if z<=1: y= 1.08*(1+z)**0.27
            if z>1: y= 1.3*z**-0.19
          if var>0 and z<=1 and z>0: y=1.08*(1+z)**0.35 
          if var>0 and z==0: y=0.95
          return y
            # z0=1.08, z1=1.3, z2=1.14
        tin= time.time()      
        sodprofile= inputparams[23].rstrip()
        halofileroot=fileroot+'/'+halofile
	cts_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_count')
	radius_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_radius')
        overden_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_rho_ratio')        
        print "\ndoing c-M calculation now.."
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, halofileroot, sodprofile, filenum, SObinstart, cts_col, radius_col, overden_col)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)

	file= open(outdir + halofile+'.cM','w')
        print>>file, "# SO Mass Msun/h, # clusters, cmean, frac err, variance\n"
        for i in range(0,massbins-1):  
            if count[i] > 3 and variance[i] !=0.0 and cmean[i] !=0:  print >>file, '%le %5d %3.4lf %3.4lf %3.4lf' %(meanmass[i], count[i],cmean[i]*fitz(zseek, model_no), (cmeanerr[i]**2+cmean[i]**2/count[i])**0.5, (variance[i]/cmean[i]**2-1.0)**0.5)
        #for i in np.arange (0,2.2,0.2): print i, fitz(i)*1.0/(1+i)**(0.71*0.89)
        file.close()
	print "time takes= ",time.time()-tin


"""
if case == 5:
     print "\n ######################### create mysql db########################\n"
     createcMdb('conc',massSO, conc, zseek,'cMdata')
  
  """ 


