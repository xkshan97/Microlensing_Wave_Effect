'''
Author: your name
Date: 2021-10-16 16:04:33
LastEditTime: 2021-10-22 17:34:46
LastEditors: Please set LastEditors
Description: In User Settings Edit
FilePath: /code/mock_/fromC2py.py
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
from tqdm import tqdm
import struct
from  multiprocessing import Process,Pool


M_L = 0.35491022739065014591
z_L = 0.5
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 4.250996065574553
gamma = 1.2136718743511312
mu = 1/((1-kappa)**2 - gamma**2)**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = 2*np.pi*mu/coeffi
cw = coeffi/2/np.pi
Cut = [0]


#adaptive
for yi in range(100):
    # yi = 0
    
    f1=open('./ResultMaximum/TimeLength_max_' + str(yi) + '.bin',"rb")
    TimeLengthFile_max=struct.unpack("l"*1, f1.read(8*1))
    f1.close()
    TimeLengthFile_max = np.array(TimeLengthFile_max)
    fileArea = "./ResultMaximum/adptive_Area_max_" + str(yi) + ".bin"
    fileTime = "./ResultMaximum/adptive_Time_max_" + str(yi) + ".bin"
    
    f1=open(fileArea,"rb")
    Area=struct.unpack("d"*TimeLengthFile_max[0], f1.read(8*TimeLengthFile_max[0]))
    f1.close()
    Area = np.array(Area)

    f2=open(fileTime,"rb")
    time=struct.unpack("d"*TimeLengthFile_max[0], f2.read(8*TimeLengthFile_max[0]))
    f2.close()
    time = np.array(time)
    
    
    index_none_zero = np.where(Area>0)[0][-1]
    time = time[0:index_none_zero+1]
    Area = Area[0:index_none_zero+1]
    
    
    
    Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
    time_raw = time[0:-1]

    time_zero_point = time_raw[-1]
    time_raw = time_raw - time_zero_point
    

    delta_time_raw = time_raw[1] - time_raw[0] 
    
    # plt.figure(1)
    # plt.semilogx(-time_raw,Ft_raw)
    # plt.semilogx([-time_raw[0],-time_raw[-1]],[constant,constant],label='smooth')
    # plt.xlabel('time')
    # plt.ylabel('dA/dt')
    # plt.xlim(10**(-4),10**(0))
    # plt.legend()
    
    
    """去掉尾巴"""
    length = len(time_raw)
    time_raw = time_raw[length*1//5::]
    Ft_raw = Ft_raw[length*1//5::]
    
    
    

    """without apodization"""
    Ft_raw[0:-1] = Ft_raw[0:-1] - constant
    Ft_raw[-1] = Ft_raw[-1] - constant/2
    
    
    plt.figure(1)
    plt.semilogx(-time_raw, Ft_raw,label='$F(t)-F_{smooth}(t)$')
    plt.xlabel('time')
    plt.ylabel('F(t)')
    plt.legend()
    plt.grid()
    """以上"""

    time_new = time_raw
    Ft_new = Ft_raw
    # """加窗归零"""
    # lengthnew = len(time_new)
    # window = np.hanning(lengthnew*2//5)
    # try:
    #     Ft_new[lengthnew*4//5+1::] = Ft_new[lengthnew*4//5+1::] * window[lengthnew*1//5::] 
    # except ValueError:
    #     Ft_new[lengthnew*4//5::] = Ft_new[lengthnew*4//5::] * window[lengthnew*1//5::] 
        
    
    np.savetxt('./ResultMaximum/time_new_adaptive_' + str(yi) + '.csv',time_new,delimiter=',')
    np.savetxt('./ResultMaximum/Ft_new_adaptive_' + str(yi) + '.csv',Ft_new,delimiter=',')
     
    
    
    """自己写傅里叶变换"""
    omegafreq = np.arange(0.1 * 2 * np.pi,2000*2*np.pi,2*np.pi)
    freq = omegafreq/2/np.pi
    Freal = []
    Fimag = []
    for i in tqdm(range(len(omegafreq))):
        Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
        Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
        
    Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
    Ff = Ff* omegafreq / complex(0,1)*cw
    Ffsgn = -constant * cw 
    Ff = Ff + Ffsgn
    
   
    Ffabs = np.abs(Ff)
    ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
    # Ffabs = signal.savgol_filter(Ffabs,3,1)
    
    
    plt.figure(2)
    
    plt.semilogx(freq, Ffabs / mu,label = 'Numerical Result') #//FIXME
    # plt.plot(freq, FfabsOld)
    # plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
    # plt.plot(freq, Ffabs1, label='Cut')
    
    plt.grid()
    # plt.legend()
    # plt.ylim(1,3)
    plt.xlabel('f[Hz]')
    plt.ylabel('|F(f)|')
    plt.xlim(0.1, 3000)
    # plt.savefig("Ffembed.png",dpi=450)
    """相位"""
    plt.figure(3)
    
    plt.semilogx(freq, np.unwrap(ThetaF), '--' , label = 'Numerical Result') #//FIXME
    # plt.plot(freq, ThetaF, label='Cut')
    # plt.loglog(test1,test2)
    # plt.ylim(-0.2,0.2)
    plt.grid()
    # plt.legend()
    plt.xlabel('f[Hz]')
    plt.ylabel(r'$\theta_F$')
    plt.xlim(0.1, 3000)


    np.savetxt('./ResultMaximum/Ffabs_adaptive_' + str(yi) + '.csv',Ffabs,delimiter=',')
    np.savetxt('./ResultMaximum/ThetaF_adaptive_' + str(yi) + '.csv',ThetaF,delimiter=',')
