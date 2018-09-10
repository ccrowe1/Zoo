# Cassie Crowe
# July 15, 2018
#
# ZooManifestMaker.py
#
# Makes a manifest for Zooniverse image upload
#
# To run: python ZooManifestMaker.py, put in the location for each image set


import numpy as np
import glob
import csv
from astropy.io import fits

hold1=[]
hold2=[]
CIV1=[]
CIVpmf1=[]
MGII1=[]
MGIIpmf1=[]
CIV2=[]
CIVpmf2=[]
MGII2=[]
MGIIpmf2=[]
MGII1=[]
MGIIpmf1=[]

SetLocal=str(raw_input("Please input the full location of the first (CIV) image set: "))  ## ask for the first image set location
SetLocal2e=str(raw_input("Please input the full location of the second (MGII) image set: "))  ## ask for the second image set location

hold1=glob.glob(SetLocal) ## makes list of CIV image file names
hold2=glob.glob(SetLocal2) ## makes list of MGII image file names

for i in range(len(hold1)):
    s=hold1[i]
    CIVpmf1.append(s[45:])
    CIV1.append(s[51:66])
for i in range(len(hold2)):
    t=hold2[i]
    MGIIpmf2.append(t[46:])
for i in range(len(CIV1)):
    a=CIV1[i]
    for i in range(len(MGII2)):
        if MGII2[i]==a:
            MGII1.append(MGII2[i])
    for i in range(len(CIVpmf1)):
        b=MGIIpmf2[i]
        if b[6:21]==a:
            MGIIpmf1.append(MGIIpmf2[i])


            

hdu=fits.open("/Users/sumrsch/Summer_Research/DR14Q_v4_4.fits")
tbdata = hdu[1].data
tbdata=tbdata[tbdata["MJD"]>56837]
MJD14=tbdata["MJD"]
name=tbdata["SDSS_NAME"]
Zvi=tbdata["Z_VI"]
Zpipe=tbdata["Z_PIPE"]
PLATE=tbdata["PLATE"]
FIBER=tbdata["FIBERID"]
BI=tbdata["BI_CIV"]

B=[]
z=[]

for i in range(len(CIV1)):
    p=CIVpmf1[i][6:10]
    f=CIVpmf1[i][11:15]
    m=CIVpmf1[i][16:21]
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
for i in range(len(CIV1)):
    #print(i)
    man=[]
    count=int(count)
    man=np.append(man,str(count))
    #print(man)
    #print(CIV1[k])
    man=np.append(man,CIVpmf1[i])
    man=np.append(man,MGIIpmf1[i])
    man=np.append(man,CIVpmf1[i][6:10])
    man=np.append(man,CIVpmf1[i][11:15])
    man=np.append(man,CIVpmf1[i][16:21])
    man=np.append(man,B[i])
    man=np.append(man,z[i])
    print(i,man)
    writer.writerow(man)
    count+=1
