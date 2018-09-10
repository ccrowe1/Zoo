# Cassie Crowe
# July 12, 2018
#
# PullDownFits.py
#
# Pulls down fits files based on what redshift range you're looking at and what type of redshift you want to use.
#
# To run: python PullDownFits.py, put in what type of redshift you want to query with, the range and the SDSS SAS user and pw
# SAS user: sdss
# SAS pw: 2.5-meter
#
# Note: FITS files should download size should be around 200 KB. If it's lower (like 4 KB) you might have to download from an earlier SAS site version.


def between(value,low,high):  ## function test if given value is between two other values
    if value>=low:
        if value<=high:
            return True
    if value<low:
        return False
    if value>high:
        return False

import numpy as np
from astropy.io import fits
import os, sys
import urllib

if __name__ == '__main__':

    Ztype=str(raw_input("Please input which redshift type you would like to use (VI, Pipe, PCA, or MG): "))
    minz=float(raw_input("Please input your minimum redshift value: "))
    maxz=float(raw_input("Please input your maximum redshift value: "))

    hdulist = fits.open('//Users/sumrsch/Summer_Research/DR14Q_v4_4.fits') ## location of catalog
    tbdata = hdulist[1].data
    tb1=tbdata[tbdata["MJD"]>56870] ## all new observations

    name=tb1["SDSS_NAME"]
    MJD14=tb1["MJD"]
    PLATE=tb1["PLATE"]
    FIBER=tb1["FIBERID"]
    BI=tb1["BI_CIV"]                     ## balnicity measure by pipeline
    Zvi=tb1["Z_VI"]                      ## visually inspected redshift
    Zpipe=tb1["Z_PIPE"]                  ## redshift measured by pipeline
    Zpca=tb1["Z_PCA"]                    ## PCA redshift measurement
    Zmg=tb1["Z_MGII"]                    ## redshift based on Magnesium II (often most accurate but not always taken)
    pdup=tbdata["PLATE_DUPLICATE"]
    mdup=tbdata["MJD_DUPLICATE"]
    fdup=tbdata["FIBERID_DUPLICATE"]
    p=[]               ## plate array for spec that are within range
    f=[]               ## fiber array for spec that are within range
    m=[]               ## MJD array for spec that are within range
    bi=[]              ## balnicity array ...
    zv=[]              ## visually inspected Z array ...
    zp=[]              ## Pipeline Z array ...
    zpca=[]            ## PCA Z array...
    zmg=[]             ## MGII Z array ...
    n=[]               ## SDSS Name array ...
    zt=[]              ## type of Z selected array ...
    pd=[]              ## plate dup array ...
    md=[]              ## MJD dup array...
    fd=[]              ## fiber dup array for spec that are within range

    if Ztype=="VI":
        zt=Zvi        ## if user chooses to use visually inspected Z values then puts all values in zt array
    if Ztype=="Pipe":
        zt=Zpipe      ## if user chooses pipeline Z...
    if Ztype=="PCA":
        zt=Zpca       ## if user chooses PCA...
    if Ztype=="MG":
        zt=Zmg        ## if user chooses MGII...
    
    for i in range(len(tb1)):
        #if between(Zvi[i],1.57,3.29): ##range used originally
        if between(zt[i],minz,maxz):
            if BI[i]>0:
                if mdup[i][1]>56870 & mdup[i][1]!=MJD14[i]:  ## if there is a new duplicate and it's not the same as the MJD at the same index...
                    pd.append(pdup[i][1])                       ## ...append the plate duplicate to pd...
                    md.append(mdup[i][1])                       ## ...append the MJD dup to md...
                    fd.append(fdup[i][1])                       ## ...append the fiber dup to fd
                if mdup[i][1]<56870:                         ## if the dup isn't new...
                    pd.append(0)                                ## ...we're not concerned so append 0
                    md.append(0)
                    fd.append(0)
                n.append(name[i])  
                p.append(PLATE[i])
                f.append(FIBER[i])
                m.append(MJD14[i])
                bi.append(BI[i])
                zp.append(Zpipe[i])
                zv.append(Zvi[i])
                zpca.append(Zpca[i])
                zmg.append(Zmg[i])
                
    print("Number of spectra in your query: "+str(len(n)))
    
    
    for i in range(len(p)):
        if md[i]>0 & pd[i] in p:                             ## if there is a duplicate spec and it's in the non-dup arrays...
            if fd[i] in f:                                      ##...
                print("DUPLICATE: ",pd[i], md[i], fd[i])        ##...print what dup PMF is...
                continue                                        ##...then skip
        pfm=str(p[i])+'-'+str(m[i])+'-'+str(f[i]).zfill(4)
        url='https://data.sdss.org/sas/ebosswork/eboss/spectro/redux/v5_10_4/spectra/lite/'+str(p[i])+'/spec-'+pfm+'.fits'  ## This is where the fits are being pulled from, change as needed
        print(pfm)
        place=os.path.join('/Users/sumrsch/Summer_Research/TestPull/BI!=0/','spec-'+pfm+'.fits')  ## This is where the fits are being placed on the computer, change as needed
        urllib.urlretrieve(url,place)
