# Britt Lundgren
# Feb. 28, 2017
# Cassie Crowe
# Feb 10, 2018
#
# DR12 QSOALS Pipeline
#
# Latest version re-evaluates the line EWs and FWHM using a Guassian profile
# Updated from v3.py by adding UNCA-specific directory structure
#
# Run with:
# python BOSSQALpipeline_UNCA.py [dirname]
#
# [dirname] = name of directory containing fits files
#

def between(value,low,high):
    if value>=low:
        if value<=high:
            return True
    if value<low:
        return False
    if value>high:
        return False
def smoothTriangle(data,degree,dropVals=False):
        """performs moving triangle smoothing with a variable degree."""
        """note that if dropVals is False, output length will be identical
        to input length, but with copies of data at the flanking regions"""
        triangle=np.array(range(degree)+[degree]+range(degree)[::-1])+1
        smoothed=[]
        for i in range(degree,len(data)-degree*2):
                point=data[i:i+len(triangle)]*triangle
                smoothed.append(sum(point)/sum(triangle))
        if dropVals: return smoothed
        #sp=sum(point)
        #print(sp,sum(triangle))
        smoothed=[smoothed[0]]*(degree+degree/2)+smoothed
        while len(smoothed)<len(data):smoothed.append(smoothed[-1])
        return smoothed
def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def twogauss(x,*p):
    A1, mu1, sigma1, A2, mu2, sigma2  = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))

def threegauss(x,*p):
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3 = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))+A3*np.exp(-(x-mu3)**2/(2.*sigma3**2))

def fourgauss(x,*p):
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3, A4, mu4, sigma4 = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))+A3*np.exp(-(x-mu3)**2/(2.*sigma3**2))+A4*np.exp(-(x-mu4)**2/(2.*sigma4**2))

def fivegauss(x,*p):
    A1, mu1, sigma1, A2, mu2, sigma2, A3, mu3, sigma3, A4, mu4, sigma4,  A5, mu5, sigma5 = p
    return A1*np.exp(-(x-mu1)**2/(2.*sigma1**2))+A2*np.exp(-(x-mu2)**2/(2.*sigma2**2))+A3*np.exp(-(x-mu3)**2/(2.*sigma3**2))+A4*np.exp(-(x-mu4)**2/(2.*sigma4**2))+A5*np.exp(-(x-mu5)**2/(2.*sigma5**2))


def gauss_fwhm(x, *p):
    A, mu, sigma = p
    return (2.355*sigma)

def readin_nspec(filename):
	lam_arr = []
	nflux_arr = []
	nfluxerr_arr = []
	cspec_arr = []
	for line in open(filename).readlines():
		cols = line.split()
		lam_arr.append(float(cols[0]))
		nflux_arr.append(float(cols[1]))
		nfluxerr_arr.append(float(cols[2]))
		cspec_arr.append(float(cols[3]))
        print "length of spectrum: %i" % (len(lam_arr))
	return(lam_arr,nflux_arr,nfluxerr_arr, cspec_arr)

def readzidin(filename):
    waves = []
    delwaves = []
    ews = []
    ewers = []
    fwhms = []
    sls = []
    for line in open(filename).readlines():
	if (line.startswith('#')>0):
	    continue
        if (line.startswith('-')>0):
            continue
	cols = line.split()
	if (len(cols)):
            if str(cols[0]).count('.')>1 or len(cols)<7:
                continue
	    elif (float(cols[0])!=9200.00):
                if 'nan' not in cols:
                    if (float(cols[0])<9200. and (float(cols[1]))<1000.) and float(cols[3])>0.:
                        waves.append(float(cols[0]))
                        delwaves.append(float(cols[1]))
                        ews.append(float(cols[2]))
                        ewers.append(float(cols[3]))
                        fwhms.append(float(cols[4]))
                        sls.append(float(cols[5]))
    return (waves,delwaves,ews,ewers,fwhms,sls)

def fix_ews(nspec, zid_input_file):

    zidin = zid_input_file+'.zid_in.dat'
    zidin_old = zid_input_file+'.zid_in_old.dat'
    os.system("mv %s %s" % (zidin, zidin_old))
    
    (waves,delwaves,ews,ewers,fwhms,sls) = readzidin(zidin_old)
    (lam_arr,nflux_arr,nfluxerr_arr, cspec_arr) = readin_nspec(nspec)
    
    fout = open(zidin, 'w')



    for i in range(0,len(waves)):
        #print "trying to fit line at lambda=%f with EW=%f" % (waves[i], ews[i])

        # bounding functions
        # lower bound:            lbound(boundary_value, parameter)
        # upper bound:            ubound(boundary_value, parameter)
        # lower and upper bounds: bound([low, high], parameter)
        # fixed parameter:        fixed(fixed_value, parameter)
        lbound = lambda p, x: 1e4*sqrt(p-x) + 1e-3*(p-x) if (x<p) else 0
        ubound = lambda p, x: 1e4*sqrt(x-p) + 1e-3*(x-p) if (x>p) else 0
        bound  = lambda p, x: lbound(p[0],x) + ubound(p[1],x)
        fixed  = lambda p, x: bound((p,p), x)
    

        # identify any of the nearest identified lines, separated by < dlam
        dlam = 15 # wavelength range across which lines are considered blended
        clumpwave = []
        clumpew = []
        clumpwave.append(float(waves[i]))
        clumpew.append(float(ews[i]))
        # number of lines to look for on each side = 1
        #if i-2>0 and abs(waves[i]-waves[i-2])<dlam:
        #    clumpwave.append(float(waves[i-2]))
        #    clumpew.append(float(ews[i-2]))
        if i-1>0 and abs(waves[i]-waves[i-1])<dlam:
            clumpwave.append(float(waves[i-1]))
            clumpew.append(float(ews[i-1]))
        if i+1<len(waves) and abs(waves[i+1]-waves[i])<dlam:
            clumpwave.append(float(waves[i+1]))
            clumpew.append(float(ews[i+1]))   
        #if i+2<len(waves) and abs(waves[i+2]-waves[i])<dlam:
        #    clumpwave.append(float(waves[i+2]))
        #    clumpew.append(float(ews[i+2]))   
 
        # for each clump, fit a Gaussian to each component 
        # first, for single cases 
        gaus = lambda p, x: p[0] * scipy.exp(-(x - p[1])**2 / (2. * p[2]**2))

        #print "clumpsize = %i" % (len(clumpwave))
        if len(clumpwave)==1:
            EWs = float(ews[i])
            wave = waves[i]
            p0 = [1., wave, EWs]
            goodfit=0
            # recalculate the EW of each line, fit with a single gaussian
            cutflux_zero = []
            cutfluxer_zero = []
            cutwave_zero = []
            counter2 = -1
            for m in range(0,len(lam_arr)):
                if lam_arr[m]>wave-40 and lam_arr[m]<wave+40:
                    counter2+=1
                    cutwave_zero.append(lam_arr[m])
                    cutflux_zero.append(-1.*float(nflux_arr[m])+1)
                    cutfluxer_zero.append(float(nfluxerr_arr[m]))
                    if lam_arr[m]<wave and lam_arr[m+1]>=wave:
                        centerwave_index = counter2            
            try:
                errfn = lambda p, x, y: gaus(p,x) - y + bound([wave-0.5,wave+0.5],p[1])
                coeff,var_matrix = scipy.optimize.leastsq(errfn, p0, args=(cutwave_zero, cutflux_zero))
                hist_fit = gaus(coeff,cutwave_zero)
                gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                gfwhm = gauss_fwhm(cutwave_zero, *coeff)

                #coeff, var_matrix = curve_fit(gauss,cutwave_zero, cutflux_zero, p0=p0, maxfev=10000)
                #hist_fit = gauss(cutwave_zero, *coeff)
                #gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                #gfwhm = gauss_fwhm(cutwave_zero, *coeff)
                if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-p0[1])<2. and gauss_ew/p0[2]<5.:
                    labelit='good'
                else:
                    labelit='reject'
                #pl.clf()
                yrange=np.linspace(-0.5,1., 10)
                xrange1 = []
                for m in range(0,len(yrange)):
                    xrange1.append(wave)
                #pl.plot(cutwave_zero, cutflux_zero, 'k-', label=labelit)
                #pl.plot(cutwave_zero, hist_fit, 'r-')
                #pl.plot(xrange1,yrange,'c-')
                #pl.legend(numpoints=1)
                #pl.savefig('gaussfit_test_'+str(zid_input_file)+'_'+str(int(waves[i]))+'_fit.png')
                noise=0
                if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-p0[1])<2. and gauss_ew/p0[2]<5.:
                    for m in range(0,len(cutflux_zero)):
                        if cutwave_zero[m]>wave-gfwhm and cutwave_zero[m]<wave+gfwhm:
                            # noise+=(cutfluxer_zero[m]**2.)
                            # fixing an error in the error calculation 07.08.14
                            noise+=((cutfluxer_zero[m]*(cutwave_zero[m]-cutwave_zero[m-1]))**2.)

                    if noise>0:
                        gauss_ewer = sqrt(noise)
                        ew_snr = gauss_ew/gauss_ewer
                        goodfit=1
                        #print >> fout, "%10.2f %10.4f %10.4f %10.4f %10.4f %10.2f 0.0" % (coeff[1], delwaves[i], gauss_ew, gauss_ewer, gfwhm, ew_snr)
                        #print "%10.2f %10.4f %10.4f %10.4f %10.4f %10.2f 0.0" % (coeff[1], delwaves[i], gauss_ew, gauss_ewer, gfwhm, ew_snr)

                else:
                    # pl.clf()
                    # pl.plot(cutwave_zero, cutflux_zero, 'k-', label='cent = %.2f' % (wave))
                    # pl.plot(cutwave_zero, hist_fit, 'r-')
                    # pl.legend(numpoints=1)
                    # pl.savefig('gaussfit_test_'+str(i)+'_unfit.png')
                    
                    # measure noise in 3pix radius around centroid
                    for m in range(centerwave_index-3,centerwave_index+3):
                        # sig+=cutflux_zero[m]
                        # noise+=(cutfluxer_zero[m]**2.) # changed this 07.08.14
                        noise+=((cutfluxer_zero[m]*(cutwave_zero[m]-cutwave_zero[m-1]))**2.)

            except:
                gauss_ew = -1.
                gauss_ewer = -1.
                ew_snr=-1.
                gfwhm = -1.
        else:
            EWs = []
            wave = []
            p0 = []
            params = []
            for m in range(0,len(clumpew)):
                EWs.append(clumpew[m])
                wave.append(clumpwave[m])
                params.append(1.)
                params.append(clumpwave[m])
                params.append(clumpew[m])
            goodfit=0
            # recalculate the EW of each line, fit with multiple gaussians
            cutflux_zero = []
            cutfluxer_zero = []
            cutwave_zero = []
            counter2 = -1
            for m in range(0,len(lam_arr)):
                if lam_arr[m]>waves[i]-40 and lam_arr[m]<waves[i]+40:
                    counter2+=1
                    cutwave_zero.append(lam_arr[m])
                    cutflux_zero.append(-1.*float(nflux_arr[m])+1)
                    cutfluxer_zero.append(float(nfluxerr_arr[m]))
                    if lam_arr[m]<waves[i] and lam_arr[m+1]>=waves[i]:
                        centerwave_index = counter2            
            try: 
                if len(EWs)==2:
                    gaus2 = lambda p, x: p[0] * exp(-(x - p[1])**2 / (2. * p[2]**2)) + p[3] * exp(-(x - p[4])**2 / (2. * p[5]**2))
                    errfn2 = lambda p, x, y: gaus2(p,x) - y + bound([wave[0]-0.5,wave[0]+0.5],p[1]) + bound([wave[1]-0.5,wave[1]+0.5],p[4])
                    coeff,var_matrix = scipy.optimize.leastsq(errfn2, params, args=(cutwave_zero, cutflux_zero))
                    coeff2 = coeff[0:3]
                    hist_fit = gauss(cutwave_zero, *coeff2)
                    gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                    gfwhm = gauss_fwhm(cutwave_zero, *coeff2)

                    #coeff, var_matrix = curve_fit(twogauss,cutwave_zero, cutflux_zero, p0=params, maxfev=10000)
                    #coeff2 = coeff[0:3]
                    #hist_fit = gauss(cutwave_zero, *coeff2)
                    #gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                    #gfwhm = gauss_fwhm(cutwave_zero, *coeff2)
                    
                    if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-params[1])<2. and gauss_ew/params[2]<5.:
                        labelit='good'
                    else:
                        labelit='reject'
                    #pl.clf()
                    yrange=np.linspace(-0.5,1.,10)
                    xrange1 = []
                    xrange2 = []
                    for m in range(0,len(yrange)):
                        xrange1.append(EWs[0])
                        xrange2.append(EWs[1])
                    #pl.plot(cutwave_zero, cutflux_zero, 'k-', label=labelit)
                    #pl.plot(cutwave_zero, twogauss(cutwave_zero, *coeff), 'b-')
                    #pl.plot(cutwave_zero, hist_fit, 'r-')
                    ##pl.plot(xrange1, yrange, 'c-')
                    ##pl.plot(xrange2, yrange, 'c-')
                    #pl.legend(numpoints=1)
                    #pl.savefig('twogaussfit_test_'+str(zid_input_file)+'_'+str(int(waves[i]))+'_fit.png')

                elif len(EWs)==3:
                    gaus3 = lambda p, x: p[0] * scipy.exp(-(x - p[1])**2 / (2. * p[2]**2)) + p[3] * scipy.exp(-(x - p[4])**2 / (2. * p[5]**2)) + p[6] * scipy.exp(-(x - p[7])**2 / (2. * p[8]**2))
                    errfn3 = lambda p, x, y: gaus3(p,x) - y + bound([wave[0]-0.5,wave[0]+0.5],p[1]) + bound([wave[1]-0.5,wave[1]+0.5],p[4]) + bound([wave[2]-0.5,wave[2]+0.5],p[7])
                    coeff,var_matrix = scipy.optimize.leastsq(errfn3, params, args=(cutwave_zero, cutflux_zero))
                    coeff2 = coeff[0:3]
                    hist_fit = gauss(cutwave_zero, *coeff2)
                    gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                    gfwhm = gauss_fwhm(cutwave_zero, *coeff2)

                #coeff, var_matrix = curve_fit(threegauss,cutwave_zero, cutflux_zero, p0=params, maxfev=10000)
                #coeff3 = coeff[0:3]
                #hist_fit = gauss(cutwave_zero, *coeff3)
                #gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                #gfwhm = gauss_fwhm(cutwave_zero, *coeff3)
                
                    if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-params[1])<2. and gauss_ew/params[2]<5.:
                        labelit='good'
                    else:
                        labelit='reject'
                    #pl.clf()
                    yrange=np.linspace(-0.5,1.,10)
                    xrange1 = []
                    xrange2 = []
                    xrange3 = []
                    for m in range(0,len(yrange)):
                        xrange1.append(EWs[0])
                        xrange2.append(EWs[1])
                        xrange3.append(EWs[2])
                    #pl.plot(cutwave_zero, cutflux_zero, 'k-', label=labelit)
                    #pl.plot(cutwave_zero, threegauss(cutwave_zero, *coeff), 'b-')
                    #pl.plot(cutwave_zero, hist_fit, 'r-')
                    ##pl.plot(xrange1, yrange, 'c-')
                    ##pl.plot(xrange2, yrange, 'c-')
                    ##pl.plot(xrange3, yrange, 'c-')
                    #pl.legend(numpoints=1)
                    #pl.savefig('threegaussfit_test_'+str(zid_input_file)+'_'+str(int(waves[i]))+'_fit.png')

                elif len(EWs)==4:
                    coeff, var_matrix = curve_fit(fourgauss,cutwave_zero, cutflux_zero, p0=params, maxfev=10000)
                    coeff4 = coeff[0:3]
                    hist_fit = gauss(cutwave_zero, *coeff4)
                    gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                    gfwhm = gauss_fwhm(cutwave_zero, *coeff4)
                    
                    if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-params[1])<2. and gauss_ew/params[2]<5.:
                        labelit='good'
                    else:
                        labelit='reject'
                    #pl.clf()
                    yrange=np.linspace(-0.5,1.,10)
                    xrange1 = []
                    xrange2 = []
                    xrange3 = []
                    for m in range(0,len(yrange)):
                        xrange1.append(EWs[0])
                        xrange2.append(EWs[1])
                        xrange3.append(EWs[2])
                    #pl.plot(cutwave_zero, cutflux_zero, 'k-', label=labelit)
                    #pl.plot(cutwave_zero, fourgauss(cutwave_zero, *coeff), 'b-')
                    #pl.plot(cutwave_zero, hist_fit, 'r-')
                    #pl.legend(numpoints=1)
                    #pl.savefig('fourgaussfit_test_'+str(zid_input_file)+'_'+str(int(waves[i]))+'_fit.png')
                
                elif len(EWs)==5:
                    coeff, var_matrix = curve_fit(fivegauss,cutwave_zero, cutflux_zero, p0=params, maxfev=10000)
                    coeff5 = coeff[0:3]
                    hist_fit = gauss(cutwave_zero, *coeff5)
                    gauss_ew = -1.*np.trapz(cutwave_zero, hist_fit)	
                    gfwhm = gauss_fwhm(cutwave_zero, *coeff5)
                    
                    if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-params[1])<2. and gauss_ew/params[2]<5.:
                        labelit='good'
                    else:
                        labelit='reject'
                    #pl.clf()
                    yrange=np.linspace(-0.5,1.,10)
                    xrange1 = []
                    xrange2 = []
                    xrange3 = []
                    for m in range(0,len(yrange)):
                        xrange1.append(EWs[0])
                        xrange2.append(EWs[1])
                        xrange3.append(EWs[2])
                    #pl.plot(cutwave_zero, cutflux_zero, 'k-', label=labelit)
                    #pl.plot(cutwave_zero, fivegauss(cutwave_zero, *coeff), 'b-')
                    #pl.plot(cutwave_zero, hist_fit, 'r-')
                    #pl.legend(numpoints=1)
                    #pl.savefig('fivegaussfit_test_'+str(zid_input_file)+'_'+str(int(waves[i]))+'_fit.png')
                
                    # sig=0
                noise=0
                if gfwhm>0 and gauss_ew>0 and abs(coeff[1]-params[1])<2. and gauss_ew/params[2]<5.:
                    for m in range(0,len(cutflux_zero)):
                        if cutwave_zero[m]>waves[i]-gfwhm and cutwave_zero[m]<waves[i]+gfwhm:
                            # sig+=cutflux_zero[m]
                            # noise+=(cutfluxer_zero[m]**2.) # changed this 07.08.14
                            noise+=((cutfluxer_zero[m]*(cutwave_zero[m]-cutwave_zero[m-1]))**2.)

                    # pl.clf()
                    # pl.plot(cutwave_zero, cutflux_zero, 'k-', label='cent = %.2f' % (wave))
                    # pl.plot(cutwave_zero, hist_fit, 'r-')
                    # pl.legend(numpoints=1)
                    # pl.savefig('gaussfit_test_'+str(i)+'_fit.png')
                    if noise>0:    
                        gauss_ewer = sqrt(noise)
                        ew_snr = gauss_ew/gauss_ewer
                        goodfit=1
                        #print >> fout, "%10.2f %10.4f %10.4f %10.4f %10.4f %10.2f 0.0" % (waves[i], delwaves[i], gauss_ew, gauss_ewer, gfwhm, ew_snr)
                        #print "%10.2f %10.4f %10.4f %10.4f %10.4f %10.2f 0.0" % (waves[i], delwaves[i], gauss_ew, gauss_ewer, gfwhm, ew_snr)

                else:
                    # pl.clf()
                    # pl.plot(cutwave_zero, cutflux_zero, 'k-', label='cent = %.2f' % (wave))
                    # pl.plot(cutwave_zero, hist_fit, 'r-')
                    # pl.legend(numpoints=1)
                    # pl.savefig('gaussfit_test_'+str(i)+'_unfit.png')
                    
                    # measure noise in 3pix radius around centroid
                    for m in range(centerwave_index-3,centerwave_index+3):
                        # sig+=cutflux_zero[m]
                        # noise+=(cutfluxer_zero[m]**2.) # changed this 07.08.14
                        noise+=((cutfluxer_zero[m]*(cutwave_zero[m]-cutwave_zero[m-1]))**2.)

            except:
                print "no good fit found!!! %f EW=%f" % (wave, EWs)
                gauss_ew = -1.
                gauss_ewer = -1.
                ew_snr=-1.
                gfwhm = -1.
                
    fout.close()
        
import os, sys
import math
from math import sqrt
import pylab as pl
from pylab import *
import matplotlib.pyplot as pl
import pyfits
from pyfits import writeto
import numpy as np
from numpy import median, std, sqrt
#import SDSScont
#from SDSScont import *
import BOSScont
from BOSScont import *
import linefinder_sv
from linefinder_sv import *
#import SDSSzidv071613
#import SDSSzidv030514
#import SDSSzidv070714
#import BOSSzidv030514
import BOSSzidv112614
import BOSSQALplotting_new_bulldogZOOM
#import BOSSQALplotting_new_bulldogM
import matplotlib.colors
from matplotlib.colors import colorConverter
import matplotlib.ticker
from matplotlib.ticker import *
import datetime
from datetime import datetime
import scipy
from scipy import optimize as opt
import scipy.stats as stats
from scipy.optimize import curve_fit, leastsq
import glob

if __name__ == '__main__':

    dirname = str(sys.argv[1])
    print("dirname is: "+dirname)
    catname = '/Users/sumrsch/Summer_Research/DR14Q_v4_4.fits'

    # make file containing the list of spectra
    files = glob.glob(dirname+'/*fits')

    spectra = []
    for file in files:
        fname = file.split('/')
        spectra.append(fname[len(fname)-1])

    # make a directory to contain pipeline output
    dirname2 = dirname+'.dir'
    if os.path.exists(dirname2):
        os.system("rm -r %s" % (dirname2))
    os.system("mkdir %s" % (dirname2))


    ## get list of plate,fiber,mjd for each spectrum
    print "Generating list of spectra in chosen directory..."
    pfmlistout =open(dirname+'/pfmlist.dat', 'w')
    platearr = []
    fiberarr = []
    mjdarr = []
    for m in range(0,len(spectra)):
        specname = spectra[m][:-5]
        pfm=specname.split('-')
        platearr.append(str(pfm[1]))
        mjdarr.append(str(pfm[2]))
        fiberarr.append(str(pfm[3]))
      
        #print int(pfm[1]),int(pfm[3]),int(pfm[2])
        #print >> pfmlistout, int(pfm[1]),int(pfm[3]),int(pfm[2])
    pfmlistout.close()

    # read in DR12Q Catalog
    qcat = '/Volumes/2T_EXT/DR14Q.fits'
    f = pyfits.open(qcat)
    tbdata2 = f[1].data
    #print tbdata2.size
    
    dr7plates = []
    dr7fibers = []
    dr7mjds = []
    dr12plates = []
    dr12fibers = []
    dr12mjds = []
    dr12ras = []
    dr12decs = []
    dr12mags = []
    dr12zPipe = []
    dr12zPCA = []
    dr12zVI = []
    dr12first_targ = []
    dr12first_flux = []
    dr12first_snr = []
    dr12first_sep = []
    dr12iabsmag = []
    dr12sdssj =[]
    dr12BALflag1= []
    dr12BALflag2= []
    dr12name = []
    dr12imags = []
    dr12counter=-1
    dr7names = []
    dr12z = []
    for m in range(0,len(platearr)):
        tbdata3 = tbdata2[tbdata2.field('PLATE')==int(platearr[m])]
        #print(platearr[m])
        #print len(tbdata3)
        #print len(tbdata3.field('FIBERID'))
        #print len(fiberarr)
        tbdata4 = tbdata3[tbdata3.field('FIBERID')==int(fiberarr[m])]
        print(fiberarr[m])
        #print len(tbdata4)
        tbdata = tbdata4[tbdata4.field('MJD')==int(mjdarr[m])]
        #print(len(tbdata))
        #print len(tbdata)
        #dr7plates.append(tbdata.field('PLATE_DUPLICATE')[0][0])
        dr7plates.append(-1)
        #dr7fibers.append(tbdata.field('FIBERID_DUPLICATE')[0][0])
        dr7fibers.append(-1)
        #dr7mjds.append(tbdata.field('MJD_DUPLICATE')[0][0])
        dr7mjds.append(-1)
        #dr7names.append(str(tbdata.field('PLATE_DUPLICATE')[0][0])+'-'+str(tbdata.field('FIBERID_DUPLICATE')[0][0])+'-'+str(tbdata.field('MJD_DUPLICATE')[0][0]))
        dr7names.append(str(-1)+'-'+str(-1)+'-'+str(-1))
        dr12plates.append(tbdata.field('PLATE')[0])
        dr12fibers.append(tbdata.field('FIBERID')[0])
        dr12mjds.append(tbdata.field('MJD')[0])
        dr12ras.append(tbdata.field('RA')[0])
        dr12decs.append(tbdata.field('DEC')[0])
        dr12zPipe.append(tbdata.field('Z_PIPE')[0])
        dr12zVI.append(tbdata.field('Z_VI')[0])##changed from Z_PCA for accuracy in DR14 catalog
        dr12first_targ.append(-1)
        dr12first_flux.append(-1)
        dr12first_snr.append(-1)
        dr12first_sep.append(-1)
        dr12iabsmag.append(tbdata.field('MI')[0])
        dr12sdssj.append(tbdata.field('SDSS_NAME')[0])
        dr12BALflag1.append(-1)
        dr12BALflag2.append(tbdata.field('BI_CIV')[0])
        #p=D12[D12['PLATE_DR7']==dr7plates[m]]
        #mj=p[p['MJD_DR7']==dr7mjds[m]]
        #f=mj[mj['FIBERID_DR7']==dr7fibers[m]]
        #dr12z.append(f['Z_PCA'])
        spec = 'spec-'+str(tbdata.field('PLATE')[0])+'-'+str(tbdata.field('MJD')[0])+'-'+str(tbdata.field('FIBERID')[0]).zfill(4)+'.fits'
        dr12imags.append(tbdata.field('PSFMAG')[0][3])
        dr12name.append(str(dr12plates[0])+'-'+str(dr12fibers[0])+'-'+str(dr12mjds[0]))
        
#        tbdata3 = tbdata2[tbdata2.field('PLATE')==int(platearr[m])]
#        tbdata4 = tbdata3[tbdata3.field('FIBERID')==int(fiberarr[m])]
#        tbdata = tbdata4[tbdata4.field('MJD')==int(mjdarr[m])]
#        dr7plates.append(tbdata.field('PLATE_DR7')[0])
#        dr7fibers.append(tbdata.field('FIBERID_DR7')[0])
#        dr7mjds.append(tbdata.field('MJD_DR7')[0])
#        dr7names.append(str(tbdata.field('PLATE_DR7')[0])+'-'+str(tbdata.field('FIBERID_DR7')[0])+'-'+str(tbdata.field('MJD_DR7')[0]))
#        dr12plates.append(tbdata.field('PLATE')[0])
#        dr12fibers.append(tbdata.field('FIBERID')[0])
#        dr12mjds.append(tbdata.field('MJD')[0])
#        dr12ras.append(tbdata.field('RA')[0])
#        dr12decs.append(tbdata.field('DEC')[0])
#        dr12zPipe.append(tbdata.field('Z_PIPE')[0])
#        dr12zPCA.append(tbdata.field('Z_PCA')[0])
#        dr12first_targ.append(tbdata.field('FIRST_MATCHED')[0])
#        dr12first_flux.append(tbdata.field('FIRST_FLUX')[0])
#        dr12first_snr.append(tbdata.field('FIRST_SNR')[0])
#        dr12first_sep.append(tbdata.field('SDSS2FIRST_SEP')[0])
#        dr12iabsmag.append(tbdata.field('MI')[0])
#        dr12sdssj.append(tbdata.field('SDSS_NAME')[0])
#        dr12BALflag1.append(tbdata.field('BAL_FLAG_VI')[0])
#        dr12BALflag2.append(tbdata.field('BI_CIV')[0])

        spec = 'spec-'+str(tbdata.field('PLATE')[0])+'-'+str(tbdata.field('MJD')[0])+'-'+str(tbdata.field('FIBERID')[0]).zfill(4)+'.fits'
        dr12imags.append(tbdata.field('PSFMAG')[0][3])
        dr12name.append(str(dr12plates[0])+'-'+str(dr12fibers[0])+'-'+str(dr12mjds[0]))
        
    newspeclist = []
    fibernamearr = []
    platenamearr = []
    newdatlist = []
    nspeclist = []
    plotlist = []
    nplotlist = []
    if len(dr12zVI)!= len(spectra):
      print "WARNING: FITS data and catalog of unequal length!"
      print "Only objects found in catalog will be fit..."
    # only use spectra in catalog
    counter=-1
    for i in range(0,len(dr12sdssj)):
        counter+=1
        fibername = str(dr12fibers[i]).zfill(4)
        platename = str(dr12plates[i])
        fibernamearr.append(fibername)
        platenamearr.append(platename)
        newspeclist.append(str('spec-'+str(platename)+'-'+str(dr12mjds[i])+'-'+fibername+'.fits'))
        newdatlist.append(str('spec_'+str(dr12plates[i])+'_'+str(fibername)+'_'+str(dr12mjds[i])+'.dat'))
        nspeclist.append(str('nspec_'+str(dr12plates[i])+'_'+str(fibername)+'_'+str(dr12mjds[i])+'.dat'))
        plotlist.append(str('spec_'+str(dr12plates[i])+'_'+str(fibername)+'_'+str(dr12mjds[i])))
        nplotlist.append(str('nspec_'+str(dr12plates[i])+'_'+str(fibername)+'_'+str(dr12mjds[i])))

        # convert each spectrum to ascii and fit a continuum
        specdirname = str(dirname+'/'+newspeclist[counter])
        print specdirname
        specname = newdatlist[counter]
        specfile = open(specname, 'w')
        
        # Read in spectrum from FITS file    
        print "Reading in FITS file: %s" % (specdirname)
        f = pyfits.open(specdirname)
        tbdata = f[1].data
        for m in range(0,len(tbdata.field('ivar'))):
            if tbdata.field('ivar')[m]>0.:
                fitserror=sqrt(1./tbdata.field('ivar')[m])
            else:
                fitserror=99.0
            fitsflux = tbdata.field('flux')[m]
            fitswave = 10**(tbdata.field('loglam')[m])
            print >> specfile, "%f %f %f" % (fitswave, fitsflux, fitserror)
        npix = len(tbdata.field('ivar'))
        specfile.close()

        #zbest = float(dr12zPCA[i])  #*# commented out for accuracy errors - ccrowe1 6/8/17 #*#

        #*# added due to accuracy errors in DR14 catalog - ccrowe1 6/8/17
        if dr12zVI[i]>0:
            zbest = float(dr12zVI[i])
            print('zbest is VI')
        elif dr12zVI[i]<0:
            zbest = float(dr12zPipe[i])
            print('zbest is PIPE')
        #if(abs(float(dr12VI[i])-float(dr12zPipe))):
            
        print("Z_Pipe is "+str(dr12zPipe[i]))
        print("Z_VI is "+str(dr12zVI[i]))
        #*#

         
        print "Sending %s to continuum-fitter...."  % (newdatlist[counter])
      
        BOSScont.contfit(int(float(dr12plates[i])),int(float(dr12fibers[i])),int(float(dr12mjds[i])),zbest,npix)
      
        # find lines in normalized spectrum
        print "Sending %s with redshift %f to line-finder...."  % (newdatlist[counter],zbest)
        linefinder_sv.linefind(int(dr12plates[i]),int(fibernamearr[i]), int(dr12mjds[i]), zbest)

        print "finished line-finding"

        # correct EWs using Gaussian profile
        zidin = 'nspec_'+str(dr12plates[i])+'_'+fibernamearr[i]+'_'+str(dr12mjds[i])
        fix_ews(nspeclist[i],zidin)
      

        # calculate average error/pix between major emission lines
        #
        # 12xx to 12yy (the main typical emission width), 12yy to NV, NV to SI IV, Si IV to C IV, 
        # C IV to C III, C III to Mg II and lam > Mg II. We could make the last Mg II to Ca II, then Ca II-8200A, 8200-9000.
        emlines = [1215.7, 1240.81, 1399.8, 1549.48, 1908.73, 2799.117, 3969., 8200.]
        mederr = []
        for k in range(0,len(emlines)-1):
            sumerr = []
            for line in open(nspeclist[i]).readlines():
                cols = line.split()
                if float(cols[0])>emlines[k]*(1.+zbest) and float(cols[0])<=emlines[k+1]*(1.+zbest):
                    sumerr.append(float(cols[2]))
            if len(sumerr)>0:
                mederr.append(2*median(sumerr))
            else:
                mederr.append(-1.0)
          
        # build systems from identified lines
        zid_input = 'nspec_'+str(platearr[i])+'_'+fibernamearr[i]+'_'+str(dr12mjds[i])

        # set BAL_ya = '-' as placeholder for Yusra's BAL info
        BAL_ya = '-'
        galz=-1
        #def brittzid(qsoname, zem, imag, ra, dec, mederr,first_targ,first_flux,first_snr,first_sep,iabMag,sdssj,zS10,BALflag, plate, fiber, mjd):

        BOSSzidv112614.brittzid(zid_input, zbest, dr12imags[i], dr12ras[i], dr12decs[i], mederr, dr12first_targ[i],dr12first_flux[i],dr12first_snr[i],dr12first_sep[i],dr12iabsmag[i],dr12sdssj[i], zbest, dr12BALflag1[i], int(dr12plates[i]), int(fiberarr[i]), int(dr12mjds[i]))

        # generate diagnostic plots
#        BOSSQALplotting_new_bulldogM.plotit(nspeclist[i],newdatlist[i],zbest,dr7names[i],dr12imags[i],dr12sdssj[i], galz)
        # if first time plotting, re-run!
 #       if i==0:
  #          BOSSQALplotting_new_bulldogM.plotit(nspeclist[i],newdatlist[i],zbest,dr7names[i],dr12imags[i],dr12sdssj[i], galz)

            # generate ZOOM diagnostic plots
        BOSSQALplotting_new_bulldogZOOM.plotit(nspeclist[i],newdatlist[i],zbest,dr7names[i],dr12imags[i],dr12sdssj[i], galz)
        # if first time plotting, re-run!
        if i==0:
            BOSSQALplotting_new_bulldogZOOM.plotit(nspeclist[i],newdatlist[i],zbest,dr7names[i],dr12imags[i],dr12sdssj[i], galz)

         ###plot zoom in on C4
        zscaled=dr12zPipe[i]
        zoomedinx=1350.0*(1+zbest)
        zoomedinxmax=1600.0*(1+zbest)
        redmin=(zoomedinx/1548.2)-1
        redmax=(zoomedinxmax/1548.2)-1

        nospecname = nspeclist[i].split('spec')[1]
        nplotname = 'norm_'+nospecname[:-4]
        plotname = 'DR10_'+nospecname[:-4]
        catname = nspeclist[i][:-4]+'CIV.brittzid_cattable.dat'
        zoomname = nspeclist[i][:-4]+'CIV.zoomin.png'
        # read in data arrays
        wave = []
        flux = []
        error = []
        # plot normalized spectrum, for visual inspection
        indexnum=-1
        probindices=[]
        x1=[]
        y1=[]
        y12=[]
        for line3 in open(newdatlist[i]).readlines():
            cols3 = line3.split()
            indexnum+=1
            #print line3.strip()
            #print(cols3[4],cols3[0])
            if float(cols3[4])<100:
                if float(cols3[0])>zoomedinx and float(cols3[0])<zoomedinxmax:
                    wave.append(float(cols3[0]))
                    x1.append(float(cols3[0]))
                    flux.append(float(cols3[1]))
                    y1.append(float(cols3[1]))
                    error.append(float(cols3[2]))
                    y12.append(float(cols3[2]))
                    
                else:
                    probindices.append(indexnum)
        
        print(zoomedinx,zoomedinxmax)
        plt.clf()
        #x1=[]
    	#y1=[]
    	#y12=[]
    	#y13=[]
    	#for m in range(0,waveliminds[1]):
            
            
            
            #y13.append(cntflux[m])
        smoothfact = 2
    	tempy1 = smoothTriangle(y1,smoothfact)#flux
        tempx1 = smoothTriangle(x1, smoothfact)#wave
        tempy12 = smoothTriangle(y12, smoothfact)#error
        #tempy13 = smoothTriangle(y13, smoothfact)
        print(len(tempx1),len(tempy1))
        civpeak=0
        sivpeak=0
        civp=[]
        sivp=[]
        for x in range(len(tempx1)):
            if between(tempx1[x],1520.087*(1+zbest),1593.387*(1+zbest)):
                civp.append(tempy1[x])
            if between(tempx1[x],1340.255*(1+zbest),1441.055*(1+zbest)):
                sivp.append(tempy1[x])
        civpeak=max(civp)
        peak=civpeak
        print("*-*-PEAK: -*-*")#copy line 926 here for missing SIV data
        if(len(sivp)>len(civp)/2):
            sivpeak=max(sivp)
            print(civpeak,sivpeak)
            if(between(sivpeak,civpeak,civpeak+(civpeak/2))):
                if(sivpeak>civpeak):
                    print("SiIV "+str(sivpeak))
                    peak=sivpeak
                else:
                    print("CIV"+str(civpeak))
                    peak=civpeak
        else:
            print("CIV "+str(civpeak))
            peak=civpeak
        if(len(sivp)<len(civp)/2):
            print("CIV "+str(civpeak))
            peak=civpeak
        medh=(median(tempy1))*8
        medtext=(median(tempy1))*7
        med=median(tempy1)
        print(medh,medtext,med)
        #sortflux=flux
        #sortflux.sort()
        #maxflux=sortflux[-7]
        #maxflux=8
        print("PeAk "+str(peak))
        plt.plot(tempx1,tempy1,'k-')#wave vs flux
        MIN=min(tempy1)-.35
        tempy12a=[x+MIN for x in tempy12]
        plt.plot(tempx1,tempy12a, 'b-.')#wave vs error
        ylab = 'F$_{\lambda}$ [10$^{-17}$ erg s$^{-1}$ cm$^{-2} \AA^{-1}$]'
        plt.ylabel(ylab)
        plt.xlabel("Wavelength(\AA)")
        ax=plt.subplot(111)
        majticx=np.arange(zoomedinx,zoomedinxmax,100)
        minticx=np.arange(zoomedinx,zoomedinxmax,20)
        majticy=np.arange(MIN,peak*1.5,2)
        minticy=np.arange(MIN,peak*1.5,1)
        if max(flux)>16:
            majticy=np.arange(MIN,peak*1.5,4)
            minticy=np.arange(MIN,peak*1.5,2)
        if max(flux)>32:
            majticy=np.arange(MIN,peak*1.5,8)
            minticy=np.arange(MIN,peak*1.5,4)
        ax.set_xticks(majticx)
        ax.set_xticks(minticx,minor=True)
        ax.set_yticks(majticy)
        ax.set_yticks(minticy,minor=True)
        ax.grid(which='minor',linestyle='dotted',alpha=.2)
        ax.grid(which='major',linestyle='dotted',alpha=.5)
        #plt.axis([zoomedinx,zoomedinxmax, 0, maxflux*1.1])
        #print("***MAX AND MIN FLUX***")
        #print(max(flux),min(flux),len(tempy12a))
        #print(min(tempy1),MIN)
        plt.axis([zoomedinx,zoomedinxmax, MIN, peak*1.5])
        if(peak>10):
            plt.axis([zoomedinx,zoomedinxmax, MIN, peak*1.3])
        B=0.006671
        NormSize=(2799.117*(1+zbest))-((1+zbest)*2799.117*sqrt((1-B)/(B+1)))
        plt.title("BI_CIV: "+str(dr12BALflag2[i])+" Normal Trough Size: "+str(round(NormSize,2))+" A",loc='right')
        #plt.title('Date observed: '+str(mjd_to_date(int(dr12mjds[i]))), loc='right')
        plt.title('DR14: '+newspeclist[i],loc='left')
        plt.title('Z_qso: '+str(round(zbest,2)),loc='center')
        #print("max flux: "+str(maxflux),"label flux: "+str(maxflux*1.025))
        #print("max flux: "+str(flux),"label flux: "+str(flux*1.025))
        #plt.text(1544.187*(1+zbest),maxflux*1.025,'CIV',fontsize=8)
        if(peak<=10):
            plt.text(1549.187*(1+zbest),peak*1.375,'CIV',fontsize=8)
            plt.axvline(1549.187*(1+zbest), color='r', linestyle='-.')
        if(peak>10):
            plt.text(1549.187*(1+zbest),peak*1.175,'CIV',fontsize=8)
            plt.axvline(1549.187*(1+zbest), color='r', linestyle='-.')
        #plt.text(1390.755*(1+zbest),maxflux*1.025,'SI IV',fontsize=8)
        plt.text(1398.755*(1+zbest),peak*1.375,'SI IV',fontsize=8)
        plt.axvline(1398.755*(1+zbest), color='r', linestyle='-.')
        if(peak>10):
            plt.text(1398.755*(1+zbest),peak*1.175,'SI IV',fontsize=8)
        plt.axvline(1398.755*(1+zbest), color='r', linestyle='-.')
        if(min(tempx1)<=5578.5 and max(tempx1)>=5578.5):
            plt.text(5578.5,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(5578.5, color='m', linestyle=':')
        if(min(tempx1)<=5894.6 and max(tempx1)>=5894.6):
            plt.text(5894.6,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(5894.6, color='m', linestyle=':')
        if(min(tempx1)<=6301.7 and max(tempx1)>=6301.8):
            plt.text(6301.7,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(6301.7, color='m', linestyle=':')
        if(min(tempx1)<=7246.0 and max(tempx1)>=7246.0):
            plt.text(7246.0,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(7246.0, color='m', linestyle=':')
        plt.savefig('/Users/sumrsch/Summer_Research/TestPull/VIz1_57-2_66/BI!=0CIV/'+zoomname, dpi=300)

        ###plot zoom in on MgII
        zscaled=dr12zPipe[i]
        zoomedinx=(redmin+1)*2796.0
        zoomedinxmax=(redmax+1)*2796.0

        nospecname = nspeclist[i].split('spec')[1]
        nplotname = 'norm_'+nospecname[:-4]
        plotname = 'DR10_'+nospecname[:-4]
        catname = nspeclist[i][:-4]+'MgII.brittzid_cattable.dat'
        zoomname = nspeclist[i][:-4]+'MgII.zoomin.png'
        # read in data arrays
        wave = []
        flux = []
        error = []
    # plot normalized spectrum, for visual inspection
        indexnum=-1
        mgFlux=0
        probindices=[]
        x1=[]
        y1=[]
        y12=[]
        for line3 in open(newdatlist[i]).readlines():
            cols3 = line3.split()
            indexnum+=1
            #print line3.strip()
            if float(cols3[4])<100:
                if float(cols3[0])>zoomedinx and float(cols3[0])<zoomedinxmax:
                    wave.append(float(cols3[0]))
                    x1.append(float(cols3[0]))
                    flux.append(float(cols3[1]))
                    y1.append(float(cols3[1]))
                    error.append(float(cols3[2]))
                    y12.append(float(cols3[2]))
                else:
                    probindices.append(indexnum)
        print("__**__zoomedinx__**__")
        print(zoomedinx,zoomedinxmax)
        plt.clf()
        #x1=[]
    	#y1=[]
    	#y12=[]
    	#y13=[]
    	#for m in range(0,waveliminds[1]):
            
            
            
            #y13.append(cntflux[m])
        smoothfact = 2
    	tempy1 = smoothTriangle(y1,smoothfact)#flux
        tempx1 = smoothTriangle(x1, smoothfact)#wave
        neivpeak=0
        neivp=[]
        mgiip=[]
        for x in range(len(tempx1)):
            if between(tempx1[x],2389.5*(1+zbest),2489.5*(1+zbest)):
                neivp.append(tempy1[x])
            if between(tempx1[x],2749.117*(1+zbest),2839.117*(1+zbest)):
                mgiip.append(tempy1[x])
        print(len(neivp),len(mgiip))
        print("*-*-PEAK: -*-*")
        neivpeak=max(neivp)
        if(len(mgiip)<len(neivp)/2):
            print("NeIV "+str(neivpeak))
            peak=neivpeak
        else:
            mgiipeak=max(mgiip)
            if(between(neivpeak,mgiipeak,mgiipeak+5)):
                print("NeIV "+str(neivpeak))
                peak=neivpeak
            else:
                print("MGII "+str(mgiipeak))
                peak=max(mgiip)
       
        #print("<<<MAX tempx1>>>")
        #print(max(tempx1))
        a=0
        shift=2800.8*(1+zbest)
        for a in range(len(tempx1)):
            if tempx1[a]<float(shift):
                mgFlux=tempy1[a]
        #print("<<<mgFLUX>>>")
        #print(mgFlux)
        tempy12 = smoothTriangle(y12, smoothfact)#error
        #tempy13 = smoothTriangle(y13, smoothfact)
        medh=(median(tempy1))*8
        medtext=(median(tempy1))*7
        med=median(tempy1)
        ##print(medh,medtext,med)
        #sortflux=flux
        #sortflux.sort()
        #maxflux=sortflux[-7]
        #maxflux=8
        plt.plot(tempx1,tempy1,'k-')#wave vs flux
        MIN=min(tempy1)-.35
        if(abs(min(tempy1)-med)>8):
            MIN=min(tempy1)-8
        tempy12a=[x+MIN for x in tempy12]
        plt.plot(tempx1,tempy12a, 'b-.')#wave vs error
        ylab = 'F$_{\lambda}$ [10$^{-17}$ erg s$^{-1}$ cm$^{-2} \AA^{-1}$]'
        plt.ylabel(ylab)
        plt.xlabel("Wavelength(\AA)")
        ax=plt.subplot(111)
        majticx=np.arange(zoomedinx,zoomedinxmax,100)
        minticx=np.arange(zoomedinx,zoomedinxmax,20)
        majticy=np.arange(MIN,peak*1.5,2)
        minticy=np.arange(MIN,peak*1.5,1)
        ax.set_xticks(majticx)
        ax.set_xticks(minticx,minor=True)
        ax.set_yticks(majticy)
        ax.set_yticks(minticy,minor=True)
        ax.grid(which='minor',linestyle='dotted',alpha=.2)
        ax.grid(which='major',linestyle='dotted',alpha=.5)
        B=0.006671
        NormSize=(2799.117*(1+zbest))-((1+zbest)*2799.117*sqrt((1-B)/(B+1)))
        plt.title("BI_CIV: "+str(dr12BALflag2[i])+" Normal Trough Size: "+str(round(NormSize,2))+" A",loc='right')
        #plt.title('Date observed: '+str(mjd_to_date(int(dr12mjds[i]))), loc='right')
        plt.title('DR14: '+newspeclist[i],loc='left')
        plt.title('Z_qso: '+str(round(zbest,2)),loc='center')
        #print("max flux: "+str(max(flux)),"label flux: "+str(maxflux*1.025))
        plt.axvline(2799.117*(1+zbest), color='r', linestyle='-.')
        plt.axis([zoomedinx,zoomedinxmax, MIN, peak*1.5])
        if((2440.187*(1+zbest))>zoomedinx):
            plt.text(2439.5*(1+zbest),peak*1.375,'[NeIV]',fontsize=8)
            plt.axvline(2439.5*(1+zbest), color='r', linestyle='-.')
        plt.text(2799.117*(1+zbest),peak*1.375,'MgII',fontsize=8)
        if(min(tempx1)<=5578.5 and max(tempx1)>=5578.5):
            plt.text(5578.5,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(5578.5, color='m', linestyle=':')
        if(min(tempx1)<=5894.6 and max(tempx1)>=5894.6):
            plt.text(5894.6,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(5894.6, color='m', linestyle=':')
        if(min(tempx1)<=6301.7 and max(tempx1)>=6301.8):
            plt.text(6301.7,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(6301.7, color='m', linestyle=':')
        if(min(tempx1)<=7246.0 and max(tempx1)>=7246.0):
            plt.text(7246.0,peak*1.275,'   Sky Line',fontsize=8)
            plt.axvline(7246.0, color='m', linestyle=':')
        #if max(tempy1)>10+tempy1[0]:
         #   MaxF=10+tempy1[0]
          #  plt.axis([zoomedinx,zoomedinxmax, MIN, MaxF])
           # majticy=np.arange(MIN,MaxF,2)
            #minticy=np.arange(MIN,MaxF,1)
#            plt.text(2420.187*(1+zbest),MaxF,'[NeIV]',fontsize=8)
 #           plt.text(2800.755*(1+zbest),MaxF,'MgII',fontsize=8)
  #      if min(tempy1)<tempy1[0]-6:
   #         MinF=tempy1[0]-6
    #        plt.axis([zoomedinx,zoomedinxmax, MinF, mgFlux*2])
     #       majticy=np.arange(MinF,mgFlux*1.7,2)
      #      minticy=np.arange(MinF,mgFlux*1.7,1)
       #     plt.text(2420.187*(1+zbest),mgFlux*1.725,'[NeIV]',fontsize=8)
        #    plt.text(2800.755*(1+zbest),mgFlux*1.725,'MgII',fontsize=8)
         #   if max(tempy1)>10+tempy1[0]:
          #      plt.axis([zoomedinx,zoomedinxmax, MinF, MaxF])
           #     majticy=np.arange(MIN,MaxF,2)
            #    minticy=np.arange(MIN,MaxF,1)
             #   plt.text(2420.187*(1+zbest),MaxF,'[NeIV]',fontsize=8)
              #  plt.text(2800.755*(1+zbest),MaxF,'MgII',fontsize=8)
#        else:
 #           majticy=np.arange(min(tempy1),mgFlux*1.7,2)
  #          minticy=np.arange(min(tempy1),mgFlux*1.7,1)
   #         plt.axis([zoomedinx,zoomedinxmax, MIN, mgFlux*2])
    #        plt.text(2420.187*(1+zbest),mgFlux*1.725,'[NeIV]',fontsize=8)
     #       plt.text(2800.755*(1+zbest),mgFlux*1.725,'MgII',fontsize=8)
        plt.grid(True)
        plt.savefig('/Users/sumrsch/Summer_Research/TestPull/VIz1_57-2_66/BI!=0MGII/'+zoomname, dpi=300)
        # clean up
        #
        # determine directory to contain pipeline output for each plate
        basename = str(dr12plates[i])+'_'+fibernamearr[i]+'_'+str(dr12mjds[i])
         
        # do not erase if plate directory already in place
        if os.path.exists(dirname2):
            print "plate directory %s exists" % (dirname2)
        else:
            os.system("mkdir %s" % (dirname2))

        # within the dir for each plate, set up a directory for each indiv spectrum
        specdir = str(dirname2)+'/spec_'+str(dr12plates[i])+'_'+fibernamearr[i]+'_'+str(dr12mjds[i])+'.dir'
        import shutil
        shutil.move(specdirname,'/Users/sumrsch/Summer_Research/TestPull/VIz1_57-2_66/PullHold')
        # do not erase if fiber directory already in place
        if os.path.exists(specdir):
            print "fiber directory %s exists" % (specdir)
        else:
            os.system("mkdir %s" % (specdir)) 

        os.system("mv *%s* %s" % (basename, specdir))
        
