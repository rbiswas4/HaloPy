import sys
from MF import *
from cM import *
from fsigma import *
from growth_ode import *
import time

rhoc= 2.77536627e11; #Msun.h^2/Mpc^3

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
simvol= (boxsize)**3 #Mpc^3/h^3
print 'mass res= %le' % (massres), 'Msun/h'    
    
k,Tk= read_Tk(Tkfile)
N= Pk_norm(sigma8, ns, k, Tk, hubble)
scalefac, D0, D1= Da_ode(Omegam, 1-Omegam, w0)
Ds, lDs= interp_D(zseek,scalefac,D0, D1)
print "D(z) and log D/log z(z)= ", Ds, lDs, 'at z= ',zseek
print "time takes= ",time.time()-tin
print "\n################## FOF MF ##########################################\n"
if case == 0 or case == 2 or case==4:
        fofproperties= inputparams[15].rstrip()
        fofcount_col= readheaders(halofileroot, fofproperties,'fof_halo_count')
        FOFcount= readhalofiles(halofileroot, fofproperties, fofcount_col,filenum)
        
        print "\ntotal # of fof clusters= ",len(FOFcount)
        FOFbinstart=min(FOFcount)*massres
        FOFbinend=max(FOFcount)*massres
        FOFmass=FOFcount*(1.0-1.0/FOFcount**0.65)*massres
        print "minFOF maxFOF"
        print "# particles= ",min(FOFcount), max(FOFcount)
	print "mass Msun/h= ", '%le %le' %(FOFbinstart, FOFbinend)
        file= open(halofileroot+'.fof_mf','w')
        print >>file, "# FOF Mass Msun/h, # clusters, dn/dln M (h/Mpc)^3, frac err, 1/sigma(M), f(sigma), fsigma_fit \n"
        meanmass, hist, binsize, bin_edges= calc_mf(FOFbinstart, FOFbinend, massbins, FOFmass)
        for i in range(0,massbins-1):  
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
           sigmaM= Ds*sigmaM
           if hist[i] >0:  print >> file, '%le %5d %7.4le %7.4lf %7.4lf %7.4lf %7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM), MF_fit(sigmaM, zseek))
	file.close()
print "\n##################### SO MF ########################################\n"
if case >= 1 and case<=4:
        tin= time.time()
        sodproperties= inputparams[16].rstrip()
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
        file= open(halofileroot+'.sod_mf','w')
        print>> file, "# SO Mass Msun/h # clusters dn/dln M (h/Mpc)^3 frac err 1/sigma(M) f(sigma)\n"
        meanmass, hist, binsize, bin_edges= calc_mf(SObinstart, SObinend, massbins, SOmass)
        for i in range(0,massbins-1):
           sigmaM= sigmam(Omegam, ns, meanmass[i], N, k, Tk, hubble)
           logsigmaM= logsigm(Omegam, ns, meanmass[i], N, k, Tk, sigmaM, hubble)
	   sigmaM= Ds*sigmaM
           if hist[i] >2:  print>> file, '%le %5d %7.4le %7.4lf %7.4lf %7.4lf' %(meanmass[i], hist[i], hist[i]/binsize/simvol,1/hist[i]**0.5, 1/sigmaM, hist[i]/binsize/simvol*meanmass[i]/(-rhoc*Omegam*logsigmaM))
        file.close()
	print "time takes= ",time.time()-tin
print "\n####################### c-M relation #################################\n"
if case >= 3:
        tin= time.time()      
        sodprofile= inputparams[17].rstrip()
	cts_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_count')
	radius_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_radius')
        overden_col= readheaders(halofileroot, sodprofile,'sod_halo_bin_rho_ratio')        
        print "\ndoing c-M calculation now.."
        conc, concerr, norm, massSO= conc_each_halo_lessmem(SOmass, SOradius, halofileroot, sodprofile, filenum, SObinstart, cts_col, radius_col, overden_col)
        cmean, cmeanerr, meanmass, count, variance= get_bin_cmean(massSO, conc, concerr, bin_edges, massbins)
	file= open(halofileroot+'.cM','w')
        print>>file, "# SO Mass Msun/h, # clusters, cmean, frac err, variance\n"
        for i in range(0,massbins-1):  
            if count[i] > 2:  print>>file, '%le %5d %3.4lf %3.4lf %3.4lf' %(meanmass[i], count[i],cmean[i], (cmeanerr[i]**2+cmean[i]**2/count[i])**0.5, (variance[i]/cmean[i]**2-1.0)**0.5)
        file.close()
	print "time takes= ",time.time()-tin
