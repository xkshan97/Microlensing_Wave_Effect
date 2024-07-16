import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
import struct
from tqdm import tqdm
import pandas as pd
import seaborn as sns
"""
新的数据放在/data/pubdata/Data3_cxc_to_sxk里了:
还是原来的参数和分辨率
lenspos_848.bin :  星的位置,848*2
map_266505625.bin : X->Y映射, 266505625*6, 像平面23231*23231个点的位置(0 1),对应的源的位置(2 3),放大率(4),引力势(5)
mag_for_each_src.bin : 1000000 int, 源平面上1000*1000个源对应的像的数目
Roots目录下data_i=i1_j=j1.bin,  记录编号i=i1,  j=j1的源成像的信息,  (N_image+1)*5,  最后一行记录该源的坐标和总放大率,  即data[N_image][0] data[N_image][1]  data[N_image][2] 
圆圈是5.5倍的像平面分辨率，会丢掉放大率0.01以下的像
"""


'''
#read map
M_L = 1
z_L = 0.5
M_sun = 2 * 10**30
G = 6.67 * 10**(-11)
c = 3 * 10**8
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
npix=20046**2
column = 1
f=open("./Phi.dat","rb")
data=struct.unpack("d"*npix*column, f.read(8*npix*column))
f.close()
# data = np.array(data).reshape(-1,6)
PhiMine = np.array(data)
data = 0

PhiMicroMine = PhiMine 
PhiMicroMine = PhiMicroMine.reshape(-1, 20046)
# Phi = np.flip(Phi,axis=0)
PhiMine = 0




plt.imshow(PhiMicroMine)
plt.colorbar()








#taylor
npix=20046**2
column = 1
f=open("./TaylorPhi.dat","rb")
data=struct.unpack("d"*npix*column, f.read(8*npix*column))
f.close()
# data = np.array(data).reshape(-1,6)
PhiMineTaylor = np.array(data)
data = 0

# print(len(PhiMine))
PhiMineTaylor = PhiMineTaylor.reshape(-1, 20046) 
PhiMicroMineTaylor = PhiMineTaylor - (1/2*0.4*(X1**2 + Y1**2) - 1/2 * 0.4 *(Y1**2 - X1**2))
PhiMineTaylor = 0
plt.imshow(PhiMicroMineTaylor)
plt.colorbar()


DeltaPhiMine = (PhiMicroMineTaylor - PhiMicroMine)/PhiMicroMine
plt.imshow(DeltaPhiMine)
plt.colorbar()
'''

#read Mnimum map
npix=20000**2
column = 1
f=open("/data/pubdata/Data_cxc_to_sxk/saddle_large/phi_stars_and_sheet_500.bin","rb")
data=struct.unpack("d"*npix*column, f.read(8*npix*column))
f.close()
# data = np.array(data).reshape(-1,6)
# ImagPosx = np.array(data[0::5])
# ImagPosy = np.array(data[1::5])
# SourPosx = np.array(data[2::5])
# SourPosy = np.array(data[3::5])
PhiMicro = np.array(data)
data = 0
PhiMicro = PhiMicro.reshape(-1, 20000) 
# ImagPosx = ImagPosx.reshape(-1, 20046)

plt.imshow(PhiMicro, vmin=-350, vmax=300)
plt.colorbar()

plt.savefig('PhiMicro_0.27.bin', dpi=450)


npix=19099*2
f=open("/data/pubdata/Data_cxc_to_sxk/saddle_large/lens_pos_19099.bin","rb")
data=struct.unpack("d"*npix, f.read(8*npix))
f.close()
LensPos = np.array(data)
data = 0
LensPos = LensPos.reshape(-1, 2) 
X = LensPos[:,0]
Y = LensPos[:,1]
plt.figure(figsize=(5,5))
plt.scatter(X, Y, s=0.1)



#read Mnimum map
npix=2000**2
f=open("./Creat_Micro_Minus.bin","rb")
data=struct.unpack("d"*npix, f.read(8*npix))
f.close()
# data = np.array(data).reshape(-1,6)
# ImagPosx = np.array(data[0::5])
# ImagPosy = np.array(data[1::5])
# SourPosx = np.array(data[2::5])
# SourPosy = np.array(data[3::5])
PhiMicro_little = np.array(data)
data = 0
PhiMicro_little = PhiMicro_little.reshape(-1, 2000) 
# ImagPosx = ImagPosx.reshape(-1, 20046)

plt.imshow(PhiMicro_little)
plt.colorbar()

X1Set = np.arange(-500 + 0.05/2, 500, 0.05)
X1Set_little = np.arange(-50 + 0.05/2, 50, 0.05)
np.where(X1Set<=X1Set_little[0])[0][-1]
index = 9000
Cut_map = np.empty([len(X1Set_little), len(X1Set_little)])
for i in range(len(X1Set_little)):
    for j in range(len(X1Set_little)):
        Cut_map[i][j] = PhiMicro[9000+i,9000+j]
        
  
plt.imshow(Cut_map)
plt.colorbar()

Error_map = np.empty([len(X1Set_little), len(X1Set_little)])
for i in range(len(X1Set_little)):
    for j in range(len(X1Set_little)):
        Error_map[i][j] = np.abs(Cut_map[i,j] - PhiMicro_little[i,j])

plt.imshow(Error_map)
plt.colorbar()    

 

#read Mnimum Layers
npix=3001 * 3001
f=open("./ResultMicro/LayersFile_min.bin","rb")
data=struct.unpack("l"*npix, f.read(8*npix))
f.close()
# data = np.array(data).reshape(-1,6)
# ImagPosx = np.array(data[0::5])
# ImagPosy = np.array(data[1::5])
# SourPosx = np.array(data[2::5])
# SourPosy = np.array(data[3::5])
LayersNum = np.array(data)
data = 0
LayersNum = LayersNum.reshape(-1, 3001) 
# ImagPosx = ImagPosx.reshape(-1, 20046)

plt.imshow(LayersNum)
plt.xlim(1300, 1700)
plt.ylim(1300, 1700)
plt.colorbar()
plt.savefig('./ResultMicro/Layer.png', dpi=450)



f=open("./Test/TFour.bin","rb")
data=struct.unpack("d"*6, f.read(8*6))
f.close()
TFour = np.array(data)


f=open("./Test/DSDT.bin","rb")
data=struct.unpack("d"*int(TFour[5]), f.read(8*int(TFour[5])))
f.close()
DSDT = np.array(data)


f=open("./Test/Time.bin","rb")
data=struct.unpack("d"*int(TFour[5]), f.read(8*int(TFour[5])))
f.close()
Time = np.array(data)




plt.plot(Time, DSDT)
plt.plot([TFour[0], TFour[0]], [0,TFour[4]], label='t1')
plt.plot([TFour[1], TFour[1]], [0,TFour[4]], label='t2')
plt.plot([TFour[2], TFour[2]], [0,TFour[4]], label='t3')
plt.plot([TFour[3], TFour[3]], [0,TFour[4]], label='t4')
plt.xlim(-186.991, -186.9875)
plt.xlim(-0.081, -0.079)

