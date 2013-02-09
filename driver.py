import sys
from MF import *
from cM import *
from fsigma import *
from growth_ode import *
from params import *
import time

from pylab import *

tin= time.time()
infile= sys.argv[1]
inputparams= read_input_params(infile) 

case= int(inputparams[0])
massbins= int(inputparams[1])
filenum= int(inputparams[2])
zseek= float(inputparams[3])
hubble= float(inputparams[4])
Omegam= float(inputparams[5])
ns= float(inputparams[6])
sigma8= float(inputparams[7])
w0= float(inputparams[8])
boxsize= float(inputparams[9])
particlenum= int(inputparams[10])
minparticle= int(inputparams[11])
Delta= float(inputparams[12])
halofileroot= inputparams[13].rstrip()
Tkfile= inputparams[14].rstrip()
massres= Omegam*rhoc*(boxsize/particlenum)**3 #Msun/h
simvol= (boxsize/hubble)**3

FOFcount, SOcount= readhalofiles(halofileroot, filenum)   
print "\n# clusters FOF SO mass res Msun/h\n"
print len(FOFcount), len(SOcount), massres    
    
FOFbinstart=min(FOFcount)*massres
FOFbinend=max(FOFcount)*massres
SObinstart=max(min(SOcount), minparticle)*massres
SObinend=max(SOcount)*massres
FOFmass=FOFcount*massres
SOmass=SOcount*massres
SOradius= (3.0*SOmass/4.0/pi/rhoc/Delta)**(1.0/3.0)
print SOmass[0], SOradius[0]
print "\nminFOF maxFOF minSO maxSO\n"
print FOFbinstart, FOFbinend, SObinstart, SObinend
   
k,Tk= read_Tk(Tkfile)
N= Pk_norm(sigma8, ns, k, Tk, hubble)
scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0)
Ds, lDs= interp_D(zseek,scalefac,D0, D1)
print "growth= ", Ds, lDs
print "time takes= ",time.time()-tin
################## FOF MF ###########################################
if case == 0 or case == 2 or case==4:
        print "\nFOF Mass Msun/h # clusters dn/dlog M Mpc^-3 frac err 1/sigma(M) f(sigma), fsigma_fit \n"
        meanmass, hist, binsize, bin_edges= calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
        for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >0:  print meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/hubble**3/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek)
         
##################### SO MF ########################################
if case >= 1 and case<=4:
        tin= time.time()
        print "\nSO Mass Msun/h # clusters dn/dlog M Mpc^-3 frac err 1/sigma(M) f(sigma)\n"
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
	   sigmaM= Ds*sigmaM
           if hist[i] >2:  print meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/hubble**3/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM)
	print "time takes= ",time.time()-tin
####################### c-M relation #######################################
if case >= 3:
        tin= time.time()
        #radius, overden, cts= readhaloprof(halofileroot, filenum, profnum) 
        print "\n profile read..and fitting for c"
        #conc, concerr, norm, massSO= conc_each_halo(SOmass, SOradius, radius, overden, cts, SObinstart)
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, halofileroot, filenum, SObinstart)
        print len(conc), len(massSO), min(conc), max(conc), mean(norm), min(norm), max(norm)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)
        print "\nSO Mass Msun/h # clusters cmean err cmean variance"
        for i in range(0,massbins-1):  
            if count[i] > 2:  print meanmass[i], count[i],cmean[i], (cmeanerr[i]**2+cmean[i]**2/count[i])**0.5, (variance[i]/cmean[i]**2-1.0)**0.5
	print "time takes= ",time.time()-tin
