import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import struct

#===============================================================================================
#Minimum
#===============================================================================================
#read position saddle1
npix1=377

f1=open("./MicroField_88/MicroField_0.28_0.28_0.11_1.04_1.00.bin","rb")
pot=struct.unpack("d"*npix1**2, f1.read(8*npix1**2))
pot=np.array(pot).reshape(-1,npix1)

plt.figure(figsize=(5,5))
plt.imshow(pot, vmin=-80, vmax=80)
plt.colorbar()

npix1=298

f1=open("./MicroField_88/MicroField_0.27_0.27_0.09_0.38_1.00.bin","rb")
pot=struct.unpack("d"*npix1**2, f1.read(8*npix1**2))
pot=np.array(pot).reshape(-1,npix1)

plt.figure(figsize=(5,5))
plt.imshow(pot)
plt.colorbar()
plt.savefig('test.png', dpi=450)

#===============================================================================================
#read stars XY
npix1=22872

f1=open("./MicroField/MicroLensCoorXY_0.40_0.40_0.24_0.74_1.00.bin","rb")
XY=struct.unpack("d"*npix1*2, f1.read(8*npix1*2))
XY=np.array(XY).reshape(npix1,2)


#===============================================================================================
#read stars Mass
npix1=246792

f1=open("./MicroField_100/Lens_Mass_2.14_2.14_2.14_1.04_1.00.bin","rb")
MASS=struct.unpack("d"*npix1, f1.read(8*npix1))
MASS=np.array(MASS)

plt.hist(np.log10(MASS_15), bins=100, histtype = 'step', label='15', density=True)
plt.hist(np.log10(MASS_88), bins=100, histtype = 'step', label='88', density=True)
plt.hist(np.log10(MASS), bins=100, histtype='step', label='100', density=True)
plt.semilogy()
plt.legend()

#===============================================================================================
#Saddle
#===============================================================================================
#read position saddle1
npix1=3048

f1=open("./MicroField/MicroField_0.55_0.55_0.46_1.00_1.00.bin","rb")
pot=struct.unpack("d"*npix1**2, f1.read(8*npix1**2))
pot=np.array(pot).reshape(-1,npix1)

plt.figure(figsize=(5,5))
plt.imshow(pot)
plt.colorbar()

#===============================================================================================
#read stars XY
npix1=92852

f1=open("./MicroField/MicroLensCoorXY_0.55_0.55_0.46_1.00_1.00.bin","rb")
XY=struct.unpack("d"*npix1*2, f1.read(8*npix1*2))
XY=np.array(XY).reshape(npix1,2)


#===============================================================================================
#read stars Mass
npix1=92852

f1=open("./MicroField/Lens_Mass_0.55_0.55_0.46_1.00_1.00.bin","rb")
MASS=struct.unpack("d"*npix1, f1.read(8*npix1))
MASS=np.array(MASS)



#===============================================================================================
#read ChariberIMF
npix1=1412
f1=open("./SampleMethod/Sample.bin","rb")
pot=struct.unpack("d"*npix1, f1.read(8*npix1))
plt.loglog(np.arange(0.08,1.5,0.01), pot)
plt.loglog(np.arange(0.08,1.5,0.01), np.arange(0.08,1.5,0.01) ** (-2.35))
np.sum(np.arange(0.01,1.5,0.01)[0:7] * pot[0:7]) / np.sum(np.arange(0.01,1.5,0.01) * pot)
np.sum(np.arange(0.01,1.5,0.01)[0:7] * np.arange(0.01,1.5,0.01)[0:7] ** (-2.35)) / np.sum(np.arange(0.01,1.5,0.01) * np.arange(0.01,1.5,0.01) ** (-2.35))



#===============================================================================================
#read ChariberIMF
RemnantMass = np.loadtxt('./SampleMethod/Remnant_MF.csv', delimiter=',')[0]
RemnantPDF = np.loadtxt('./SampleMethod/Remnant_MF.csv', delimiter=',')[1]
npix1=240763 + 6029
f1=open("./SampleMethod/SampleTest.bin","rb")
SampleIMF=struct.unpack("d"*npix1, f1.read(8*npix1))

IMF = np.array(SampleIMF[0:240763])
Remnant = np.array(SampleIMF[240763:])
print(str(np.sum(IMF>10))+'\n')
print(np.sum(IMF>10)/len(IMF))
plt.figure()
plt.hist(IMF[np.where(IMF>10)], bins=100)
plt.semilogy()

plt.figure()
plt.hist(IMF, bins=10000, density=True)
plt.plot(np.arange(np.min(IMF), np.max(IMF)), np.arange(np.min(IMF), np.max(IMF))**(-2.35) / (1 / (-1.35) * np.max(IMF)**(-1.35) - 1 / (-1.35) * np.min(IMF)**(-1.35)))
plt.loglog()

plt.loglog(RemnantMass, RemnantPDF)
plt.hist(Remnant, bins=100, density=True)
