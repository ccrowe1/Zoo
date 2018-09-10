# Cassie Crowe
# July 15, 2018
#
# ZooManifestMaker.py
#
# Makes a CSV manifest for Zooniverse image upload
#
# To run: python ZooManifestMaker.py, put in the location for each image set


import numpy as np
import glob
import csv
from astropy.io import fits

hdu=fits.open("/Users/sumrsch/Summer_Research/DR14Q_v4_4.fits")
tbdata = hdu[1].data
tbdata=tbdata[tbdata["MJD"]>56837]
MJD14=tbdata["MJD"]
name=tbdata["SDSS_NAME"]
Zvi=tbdata["Z_VI"]
Zpipe=tbdata["Z_PIPE"]
PLATE=tbdata["PLATE"]
FIBER=tbdata["FIBERID"]

B=[]
z=[]

hold1=[]
hold2=[]
CIVpmf1=[]
CIVfile1=[]
MGIIpmf1=[]
MGIIfile1=[]
CIVfile2=[]
CIVpmf2=[]
MGIIpmf2=[]
MGIIfile2=[]

SetLocal=str(raw_input("Please input the full location of the first (CIV) image set: "))  ## ask for the first image set location
SetLocal2e=str(raw_input("Please input the full location of the second (MGII) image set: "))  ## ask for the second image set location

hold1=glob.glob(SetLocal) ## makes list of CIV image file names
hold2=glob.glob(SetLocal2) ## makes list of MGII image file names

for i in range(len(hold1)):      ## for-loop no.1: appends the PMF from CIV image filepath to CIVpmf1 array and filename to CIVfile1 
    s=hold1[i]
    CIVfile1.append(s[45:])
    CIVpmf1.append(s[51:66])
for i in range(len(hold2)):      ## for-loop no.2: appends filename to MGIIfile2 array
    t=hold2[i]
    MGIIfile2.append(t[46:])
for i in range(len(CIVpmf1)):    ## for-loop no.3: loops through to line up CIV and MGII at the same indexes in their arrays
    a=CIV1[i]
    for i in range(len(MGII2)):             
        if MGIIpmf2[i]==a:
            MGIIpmf1.append(MGIIpmf2[i])
    for i in range(len(CIVfile1)):
        b=MGIIfile2[i]
        if b[6:21]==a:
            MGIIfile1.append(MGIIfile2[i])


for i in range(len(CIV1)):       ## for-loop no.4: assigns plate, fiber, mjd, redshift and balnicity to their own arrays
    p=CIVfile1[i][6:10]
    f=CIVfile1[i][11:15]
    m=CIVfile1[i][16:21]
    zp=tbdata[tbdata["PLATE"]==int(p)]
    zf=zp[zp["FIBERID"]==int(f)]
    z.append(zf["Z_VI"])
    B.append(zf["BI_CIV"])
    
f=open("QSOflagManifestCALIB.csv","wb") ## name of final manifest
f.write('id,image1,image2,Plate,Fiber_ID,MJD,BI_CIV,z_vi\n') ## Zooniverse requires id and image1 for show one image.
## Add image2, image3 etc for additional image to be shown together and whatever other info to be visable to volunteers
writer=csv.writer(f,dialect='excel',delimiter=',')
count=1
#k=0
#print(len(hold1))
print("CIV1")
for i in range(len(CIV1)):      ## writes the final csv manifest
    #print(i)
    man=[]
    count=int(count)
    man=np.append(man,str(count))
    #print(man)
    #print(CIV1[k])
    man=np.append(man,CIVfile1[i])
    man=np.append(man,MGIIfile1[i])
    man=np.append(man,CIVfile1[i][6:10])
    man=np.append(man,CIVfile1[i][11:15])
    man=np.append(man,CIVfile1[i][16:21])
    man=np.append(man,B[i])
    man=np.append(man,z[i])
    print(i,man)
    writer.writerow(man)
    count+=1
