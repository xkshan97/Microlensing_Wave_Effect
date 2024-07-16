import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.interpolate import interp1d
import scipy.fft
from numpy import polyfit, poly1d
from tqdm import tqdm
import struct
from  multiprocessing import Process,Pool

class Diffrac(object):
    def __init__(self):
        
        self.M_sun = 1.9884099 * 10**30

        self.G = 6.6743*10**(-11)

        self.c = 2.9979246*10**8
        

    def SIS(self, M_L, z_L, y1, y2, fileArea, fileTime, binlength, freq):

        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        sqrtmu = 1
        constant = 2*np.pi/coeffi * sqrtmu
        cw = coeffi/2/np.pi
        print(constant)
        
        f1=open(fileArea,"rb")
        Area=struct.unpack("d"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)

        f2=open(fileTime,"rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)
        
    
        Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
        time_raw = time[0:-1]
    
        time_raw = time_raw - time_raw[0]
        

        delta_time_raw = time_raw[1] - time_raw[0]
        
        plt.semilogx(time_raw,Ft_raw)
        plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
        plt.xlabel('time')
        plt.ylabel('dA/dt')
        # plt.xlim(10**(-2),10**(-1))
        plt.legend()
        
        
        """去掉尾巴"""
        length = len(time_raw)
        time_raw = time_raw[0:length*4//5]
        Ft_raw = Ft_raw[0:length*4//5]
        


        """without apodization"""
        Ft_raw[1::] = Ft_raw[1::] - constant
        Ft_raw[0] = Ft_raw[0] - constant/2
        
        
        
        plt.semilogx(time_raw, Ft_raw,label='$F(t)-F_{smooth}(t)$')
        plt.xlabel('time')
        plt.ylabel('F(t)')
        plt.legend()
        plt.grid()
        """以上"""

        time_new = time_raw
        Ft_new = Ft_raw
        #保存时域
        savetmp1 = fileArea.split('Read')[0]
        savetmp2 = savetmp1.split('/')[-1]
        np.savetxt('./Ftnew/FtnewSIS_' + y1 + '_' + y2 + '.csv', Ft_new, delimiter=',')
        np.savetxt('./Ftnew/timenewSIS_' + y1 + '_' + y2 + '.csv', time_new, delimiter=',')
    
        """加窗归零"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*2//5)
        try:
            Ft_new[lengthnew*4//5+1::] = Ft_new[lengthnew*4//5+1::] * window[lengthnew*1//5::] 
        except ValueError:
            Ft_new[lengthnew*4//5::] = Ft_new[lengthnew*4//5::] * window[lengthnew*1//5::] 
            
        
            
        """自己写傅里叶变换"""
        omegafreq = freq * 2 * np.pi
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = constant * cw 
        Ff = Ff + Ffsgn
        
    
        # """以上"""
        Ffabs = np.abs(Ff)
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
        # plt.figure(1)
        
        # plt.semilogx(freq, Ffabs,label = 'Numerical Result') #//FIXME
    
        # plt.plot([freq[0],freq[-1]],[sqrtmu,sqrtmu],label='Macro magnification')
    
        # plt.grid()
        # plt.legend()
    
        # plt.xlabel('f[Hz]')
        # plt.ylabel('|F(f)|')
        # # plt.savefig("Ffembed.png",dpi=450)
        # """相位"""
        # plt.figure(2)
        
        # plt.plot(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
    
        # plt.grid()
        # plt.legend()
        # plt.xlabel('f[Hz]')
        # plt.ylabel(r'$\theta_F$')
        return Ffabs, ThetaF
    
    
    def Minimum(self, M_L, z_L, kappa, gamma, fileArea, fileTime, binlength, freq, kappa_star):

        
        mu = 1/((1-kappa)**2 - gamma**2)**0.5 #sqrt(mu)
        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        constant = 2*np.pi*mu/coeffi
        cw = coeffi/2/np.pi
        
        print(constant)
        
        f1=open(fileArea,"rb")
        Area=struct.unpack("d"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)
        
        f2=open(fileTime,"rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)
        
        index_none_zero = np.where(Area>0)[0][0]
        time = time[index_none_zero::]
        Area = Area[index_none_zero::]
        
        
        
        Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
        time_raw = time[0:-1]
    
        time_raw = time_raw - time_raw[0]
        

        delta_time_raw = time_raw[1] - time_raw[0] 
    
        
        plt.semilogx(time_raw,Ft_raw)
        plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
        plt.xlabel('time')
        plt.ylabel('dA/dt')
        # plt.xlim(10**(-2),10**(-1))
        plt.legend()
        
        
        """去掉尾巴"""
        length = len(time_raw)
        time_raw = time_raw[0:length*4//5]
        Ft_raw = Ft_raw[0:length*4//5]
        


        """without apodization"""
        Ft_raw[1::] = Ft_raw[1::] - constant
        Ft_raw[0] = Ft_raw[0] - constant/2
        
        
        
        # plt.semilogx(time_raw, Ft_raw,label='$F(t)-F_{smooth}(t)$')
        # plt.xlabel('time')
        # plt.ylabel('F(t)')
        # plt.legend()
        # plt.grid()
        """以上"""

        time_new = time_raw
        Ft_new = Ft_raw
        #保存时域
        savetmp1 = fileArea.split('Read')[0]
        savetmp2 = savetmp1.split('/')[-1]
        np.savetxt('./Ftnew/FtnewMinimum_'+kappa_star+'.csv', Ft_new, delimiter=',')
        np.savetxt('./Ftnew/timenewMinimum_'+kappa_star+'.csv', time_new, delimiter=',')
        """加窗归零"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*2//5)
        try:
            Ft_new[lengthnew*4//5+1::] = Ft_new[lengthnew*4//5+1::] * window[lengthnew*1//5::] 
        except ValueError:
            Ft_new[lengthnew*4//5::] = Ft_new[lengthnew*4//5::] * window[lengthnew*1//5::] 
            
        
            
        """自己写傅里叶变换"""
        omegafreq = freq * 2 * np.pi
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = constant * cw 
        Ff = Ff + Ffsgn
        
    
        # """以上"""
        Ffabs = np.abs(Ff)
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
        # plt.figure(1)
        
        # plt.semilogx(freq, Ffabs,label = 'Numerical Result') #//FIXME
    
        # plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
    
        # plt.grid()
        # plt.legend()
    
        # plt.xlabel('f[Hz]')
        # plt.ylabel('|F(f)|')
        # # plt.savefig("Ffembed.png",dpi=450)
        # """相位"""
        # plt.figure(2)
        
        # plt.plot(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
    
        # plt.grid()
        # plt.legend()
        # plt.xlabel('f[Hz]')
        # plt.ylabel(r'$\theta_F$')
        return Ffabs, ThetaF
    
    def Saddle(self, M_L, z_L, kappa, gamma, fileArea, fileTime, binlength, x10, x20, freq,kappa_star):
        
        mu = abs(1/((1-kappa)**2 - gamma**2))**0.5
        coeffi = 4 * self.G * M_L * self.M_sun * (1 + z_L) / self.c**3
        constant = mu/coeffi
        cw = coeffi/2/np.pi
        
        print(constant)
        
        f1=open(fileArea,"rb")
        Area=struct.unpack("d"*binlength, f1.read(8*binlength))
        f1.close()
        Area = np.array(Area)
        
        f2=open(fileTime,"rb")
        time=struct.unpack("d"*binlength, f2.read(8*binlength))
        f2.close()
        time = np.array(time)
        
        
        Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
        time_raw = time[0:-1]
    
        time_raw = time_raw #- time_raw[0]
        

        delta_time_raw = time_raw[1] - time_raw[0] 
        

        
        # plt.savefig('./test.png')
       
        #防止出现t=0，因为ln|t|。
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

        plt.plot(time_raw,Ft_raw)
        # plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant])
        plt.xlabel('time')
        plt.ylabel('dA/dt')
        plt.grid()
        
        mur = 1 - kappa + gamma
        mut = kappa + gamma - 1
        Axisa = np.sqrt(1/coeffi/mur)
        Axisb = np.sqrt(1/coeffi/mut)
        index1 = np.where(time_raw < 0)
        index2 = np.where(time_raw >= 0)
        F1 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisa) + np.log(np.abs(time_raw[index1])) - 2 * np.log(x10 + np.sqrt(2 * Axisa**2 * np.abs(time_raw[index1]) + x10**2)))
        F2 = - 2 * mu / coeffi * (np.log(2) + 2 * np.log(Axisb) + np.log(np.abs(time_raw[index2])) - 2 * np.log(x20 + np.sqrt(2 * Axisb**2 * np.abs(time_raw[index2]) + x20**2))) 
        Ft_theory = np.append(F1, F2)
        plt.plot(time_raw, Ft_theory, '--')
        Ft_subtract = Ft_raw - Ft_theory
        
        time_new = time_raw
        #保存时域
        savetmp1 = fileArea.split('Read')[0]
        savetmp2 = savetmp1.split('/')[-1]
        np.savetxt('./Ftnew/FtnewSaddle_'+kappa_star+'.csv', Ft_subtract, delimiter=',')
        np.savetxt('./Ftnew/timenewSaddle_'+kappa_star+'.csv', time_new, delimiter=',')
        np.savetxt('./Ftnew/FtTheorySaddle_'+kappa_star+'.csv', Ft_theory, delimiter=',')
        """加窗归零"""
        lengthnew = len(time_new)
        window = np.hanning(lengthnew*4//10)
        try:
            Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
            Ft_subtract[lengthnew*8//10+1::] = Ft_subtract[lengthnew*8//10+1::] * window[lengthnew*2//10::] 
        except ValueError:
            Ft_subtract[0:lengthnew*2//10] = Ft_subtract[0:lengthnew*2//10] * window[0:lengthnew*2//10]
            Ft_subtract[lengthnew*8//10::] = Ft_subtract[lengthnew*8//10::] * window[lengthnew*2//10::] 
            


        """自己写傅里叶变换"""
        omegafreq = 2 * np.pi * freq
        Freal = []
        Fimag = []
        for i in range(len(omegafreq)):
            Freal.append(np.sum(Ft_subtract*np.cos(omegafreq[i]*time_new))*delta_time_raw)
            Fimag.append(np.sum(Ft_subtract*np.sin(omegafreq[i]*time_new))*delta_time_raw)
            
        Ff = np.array(Freal) + complex(0,1)*np.array(Fimag)
        Ff = Ff* omegafreq / complex(0,1)*cw
        Ffsgn = - np.sign(omegafreq)* 2*np.pi*constant * cw * complex(0,1)
        Ff = Ff + Ffsgn


        Ffabs = np.abs(Ff) #取模   
        ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs)) #计算相位
        
        return Ffabs, ThetaF 
    