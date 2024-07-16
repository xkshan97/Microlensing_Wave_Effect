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
from scipy.optimize import curve_fit
import struct
from  multiprocessing import Process,Pool


f0=open('./MicroField/AveMassAndNum_0.55_0.55_0.46_1.00_1.00.bin',"rb")
AveMassAndNum=struct.unpack("d"*3, f0.read(8*3))
f0.close()


f0=open('./ResultSaddle/X1020New_0.55_0.55_0.46_1.00_1.00.bin',"rb")
X10X20New=struct.unpack("d"*2, f0.read(8*2))
f0.close()

M_L = AveMassAndNum[0]
z_L = 1
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 0.55
gamma = 0.55
mu = abs(1/((1-kappa)**2 - gamma**2))**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = mu/coeffi
cw = coeffi/2/np.pi
x_step = 0.001 #像平面分辨率


"""adaptive"""
for yi in range(60):
    
    f1=open('./ResultSaddle/TimeLength_sad_0.55_0.55_0.46_1.00_1.00.bin',"rb")
    TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
    f1.close()
    TimeLengthFile_min = np.array(TimeLengthFile_min)
    fileArea = "./ResultSaddle/adptive_Area_sad_0.55_0.55_0.46_1.00_1.00.bin"
    fileTime = "./ResultSaddle/adptive_Time_sad_0.55_0.55_0.46_1.00_1.00.bin"
    
    f1=open(fileArea,"rb")
    Area=struct.unpack("d"*TimeLengthFile_min[0], f1.read(8*TimeLengthFile_min[0]))
    f1.close()
    Area = np.array(Area)

    f2=open(fileTime,"rb")
    time=struct.unpack("d"*TimeLengthFile_min[0], f2.read(8*TimeLengthFile_min[0]))
    f2.close()
    time = np.array(time)
    
    
    
    
    Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
    time_raw = time[0:-1]

    

    delta_time_raw = time_raw[1] - time_raw[0] 
    
    # plt.figure(1)
    # plt.semilogx(time_raw,Ft_raw)
    # # plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
    # plt.xlabel('time')
    # plt.ylabel('dA/dt')
    # # plt.xlim(10**(-2),10**(-1))
    # plt.legend()
    
    
    index0 = np.where(time_raw ==0)[0]
    if len(index0) != 0:
        print(index0)
        time_raw = (time_raw[1::] + time_raw[0:-1])/2
        Ft_raw = Ft_raw[0:-1]
        plt.plot(time_raw, Ft_raw)
    else:
        pass
    
    
    length = len(time_raw)
    time_raw = time_raw[length//5:length//5*4] #去掉尾巴处的下降部分
    Ft_raw = Ft_raw[length//5:length//5*4]
   
    # np.savetxt('./ResultSaddle/time_raw'+file+'.csv',time_raw,delimiter=',')
    # np.savetxt('./ResultSaddle/Ft_raw'+file+'.csv',Ft_raw,delimiter=',')
    
    x10 = X10X20New[0]
    x20 = X10X20New[1]
    mur = 1 - kappa + gamma
    mut = kappa + gamma - 1
    Axisa = np.sqrt(1/coeffi/mur)
    Axisb = np.sqrt(1/coeffi/mut)
    index1 = np.where(time_raw < 0)
    index2 = np.where(time_raw >= 0)
    F1 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisa) + np.log(np.abs(time_raw[index1])) - 2 * np.log(x10 + np.sqrt(2 * Axisa**2 * np.abs(time_raw[index1]) + x10**2)))
    F2 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisb) + np.log(np.abs(time_raw[index2])) - 2 * np.log(x20 + np.sqrt(2 * Axisb**2 * np.abs(time_raw[index2]) + x20**2))) 
    Ft_theory = np.append(F1, F2)
    
    
    # step = 0.723605
    # X1Start = -1103.76 + step / 2
    # X1Tmp = X1Start
    # X1Set = []
    # i = 0
    # while X1Tmp < 1103.76:
    #     X1Set.append(X1Tmp)
    #     i += 1
    #     X1Tmp += step
    # X1Set = np.array(X1Set)
    # X10New = 317.922
    # X20New = 953.765
    # X1IndexSet = np.where((X1Set > - X10New)&(X1Set < X10New))[0]
    # X2IndexSet = np.where((X1Set > - X20New)&(X1Set < X20New))[0]
    # x11 = - X1Set[X1IndexSet[0]] + step / 2
    # x12 = X1Set[X1IndexSet[-1]] + step / 2
    # x21 = - X1Set[X2IndexSet[0]] + step / 2
    # x22 = X1Set[X2IndexSet[-1]] + step / 2
    
    # mur = 1 - kappa + gamma
    # mut = kappa + gamma - 1
    # Axisa = np.sqrt(1/coeffi/mur)
    # Axisb = np.sqrt(1/coeffi/mut)
    # index1 = np.where(time_raw < 0)
    # index2 = np.where(time_raw >= 0)
    # # F1 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisa) + np.log(np.abs(time_raw[index1])) - 2 * np.log(x10 + np.sqrt(2 * Axisa**2 * np.abs(time_raw[index1]) + x10**2)))
    # # F2 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisb) + np.log(np.abs(time_raw[index2])) - 2 * np.log(x20 + np.sqrt(2 * Axisb**2 * np.abs(time_raw[index2]) + x20**2))) 
    # F1 = mu / coeffi * (-2 + x11/np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x11 ** 2) + (2 * Axisa ** 2 * np.abs(time_raw[index1]))/(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x11 * (x11 + np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x11 ** 2))) + x12/np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x12 ** 2) + (2 * Axisa ** 2 * np.abs(time_raw[index1]))/(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x12 * (x12 + np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x12 ** 2))) - np.log(4) - 4 * np.log(Axisa) - 2 * np.log(np.abs(time_raw[index1])) + 2 * np.log(x11 + np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x11 ** 2)) + 2 * np.log(x12 + np.sqrt(2 * Axisa ** 2 * np.abs(time_raw[index1]) + x12 ** 2)))
    # F2 = mu / coeffi * (-2 + x21/np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x21 ** 2) + (2 * Axisb ** 2 * np.abs(time_raw[index2]))/(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x21 * (x21 + np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x21 ** 2))) + x22/np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x22 ** 2) + (2 * Axisb ** 2 * np.abs(time_raw[index2]))/(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x22 * (x22 + np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x22 ** 2))) - np.log(4) - 4 * np.log(Axisb) - 2 * np.log(np.abs(time_raw[index2])) + 2 * np.log(x21 + np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x21 ** 2)) + 2 * np.log(x22 + np.sqrt(2 * Axisb ** 2 * np.abs(time_raw[index2]) + x22 ** 2)))
    # Ft_theory = np.append(F1, F2)

    Ft_subtract = Ft_raw - Ft_theory
    
    time_new = time_raw
    # """加窗归零"""
    # lengthnew = len(time_new)
    # window = np.hanning(lengthnew*4//10)
    # try:
    #     Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
    #     Ft_subtract[lengthnew*8//10+1::] = Ft_subtract[lengthnew*8//10+1::] * window[lengthnew*2//10::] 
    # except ValueError:
    #     Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
    #     Ft_subtract[lengthnew*8//10::] = Ft_subtract[lengthnew*8//10::] * window[lengthnew*2//10::] 
        
        
    np.savetxt('./To_Meena/Saddle_time.csv',time_new,delimiter=',')
    np.savetxt('./To_Meena/Saddle_Ft_new.csv',Ft_subtract,delimiter=',')
    
    plt.figure(1)
    plt.plot(time_new,Ft_raw)
    # plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
    plt.xlabel('time')
    plt.ylabel('dA/dt')
    # plt.xlim(10**(-2),10**(-1))
    plt.legend()
    plt.grid()
    plt.savefig('./To_Meena/Saddle_Ft.png',dpi=450)

    


    """自己写傅里叶变换"""
    omegafreq = np.arange(0.1*2*np.pi,2000*2*np.pi,2*np.pi)
    freq = omegafreq/2/np.pi
    Freal = []
    Fimag = []
    for i in tqdm(range(len(omegafreq))):
        Freal.append(np.sum(Ft_subtract*np.cos(omegafreq[i]*time_new))*delta_time_raw)
        Fimag.append(np.sum(Ft_subtract*np.sin(omegafreq[i]*time_new))*delta_time_raw)
        
    Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
    Ff = Ff* omegafreq / complex(0,1)*cw
    Ffsgn = - np.sign(omegafreq)* 2*np.pi*constant * cw * complex(0,1)
    Ff = Ff + Ffsgn

    

    Ffabs = np.abs(Ff) #取模   
    ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs)) #计算相位
    # Ffabs = signal.savgol_filter(Ffabs, 9, 1) #这个是把振幅做平滑平均
    # ThetaF = signal.savgol_filter(ThetaF, 119, 1)
    Ff = Ffabs * np.exp(complex(0,1)*ThetaF) #用平滑后的振幅重组复数Ff
    # freq = freq[1::] #去除0频率
    # Ff = Ff[1::]
    Ffabs = np.abs(Ff)
    
    
    plt.figure(2)
    plt.semilogx(freq, Ffabs, label = 'Sign') #//FIXME
    # plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
    # plt.plot(freq, Ffabs1, label='Cut')
    # plt.plot(freq, FfTheoryAbs, '--', color = 'red',label='Geometric')
    
    # plt.xlim(1,2000)
    # plt.loglog(test1,test2)
    # plt.ylim(0,10)
    plt.grid()
    # plt.legend()
    plt.xlabel('f[Hz]')
    plt.ylabel('|F(f)|')
    plt.savefig("./To_Meena/Saddle_Ffabs.png",dpi=450)
    
    plt.figure(3)
    plt.semilogx(freq, ThetaF, label = 'Sign') #//FIXME
    # plt.semilogx(freq, FfTheoryTheta, label='Geometric')
    # plt.xlim(1,2000)
    # plt.loglog(test1,test2)
    # plt.ylim(-2.3,-1)
    plt.grid()
    # plt.legend()
    plt.xlabel('f[Hz]')
    plt.ylabel(r'$\theta_F$')
    plt.savefig("./To_Meena/Saddle_ThetaF.png",dpi=450)


    
    np.savetxt('./To_Meena/Saddle_Ffabs.csv',Ffabs,delimiter=',')
    np.savetxt('./To_Meena/Saddle_ThetaF.csv',ThetaF,delimiter=',')


plt.figure(2)
plt.grid()

plt.figure(3)
plt.grid()
