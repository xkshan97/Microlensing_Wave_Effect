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
import copy

M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
            
AcceptLensIndex = np.loadtxt('../Paper4_CE_Modify/SampleResult/AcceptLensIndex.csv', delimiter=',')
imagenum = np.loadtxt('../Paper4_CE_Modify/SampleResult/imagenumber.csv', delimiter=',')
SNR_network = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult/SNR_network_only_macro.csv',delimiter=',')
kappa_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/kappa.csv', delimiter=',')
gamma_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/gamma.csv', delimiter=',')
kappa_s_set = np.loadtxt("../Paper4_CE_Modify/SampleResult/kappa_s.csv", delimiter=',')
timedelay = np.loadtxt('../Paper4_CE_Modify/SampleResult/timedelayDays.csv', delimiter=',')
lens_z_set = np.loadtxt('../Paper4_CE_Modify/SampleResult/lens_z.csv', delimiter=',')
magnification = np.loadtxt('../Paper4_CE_Modify/SampleResult/magnification.csv', delimiter=',')
Total_inject_time = np.loadtxt('../Paper4_CE_Modify/Lensed_SampleResult_Micro_quadruple/Total_inject_time.csv', delimiter=',')

IndexSumimage = 0
for ImageIndex in tqdm(range(100)):
    IndexInAccept = AcceptLensIndex[ImageIndex]
    LensRedshift = lens_z_set[ImageIndex]
    for subImageIndex in range(int(imagenum[ImageIndex])):
        if SNR_network[IndexSumimage] >= 12:
            kappa = kappa_set[IndexSumimage]
            gamma = gamma_set[IndexSumimage]
            kappaStar_stellar = kappa_s_set[IndexSumimage]
            
            # MainDiffraction(kappa, gamma, 1.1 * kappaStar_stellar, LensRedshift, 1, 200);
            
            if 1.1 * kappaStar_stellar > kappa:
                kappaStar_stellar = kappa / 1.1



            # f0=open('./MicroField_88/AveMassAndNum_' + str(round(kappa, 2)) + '_' + str(round(gamma, 2)) + '_' + str(round(1.1 * kappaStar_stellar, 2)) + '_' + str(round(LensRedshift, 2)) + '_1.00.bin',"rb")
            f0=open('./MicroField/AveMassAndNum_%.2f_%.2f_%.2f_%.2f_1.00.bin'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),"rb")
            AveMassAndNum=struct.unpack("d"*3, f0.read(8*3))
            f0.close()

            M_L = AveMassAndNum[0]
            z_L = LensRedshift
            mu = 1/((1-kappa)**2 - gamma**2)**0.5
            coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
            constant = 2*np.pi*mu/coeffi
            cw = coeffi/2/np.pi


            
            f1=open('./ResultMinimum/TimeLength_min_%.2f_%.2f_%.2f_%.2f_1.00.bin'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),"rb")
            TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
            f1.close()
            TimeLengthFile_min = np.array(TimeLengthFile_min)
            fileArea = "./ResultMinimum/adptive_Area_min_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
            fileTime = "./ResultMinimum/adptive_Time_min_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
            
            f1=open(fileArea,"rb")
            Area=struct.unpack("d"*TimeLengthFile_min[0], f1.read(8*TimeLengthFile_min[0]))
            f1.close()
            Area = np.array(Area)

            f2=open(fileTime,"rb")
            time=struct.unpack("d"*TimeLengthFile_min[0], f2.read(8*TimeLengthFile_min[0]))
            f2.close()
            time = np.array(time)
            
            
            index_none_zero = np.where(Area>0)[0][0]
            time = time[index_none_zero::]
            Area = Area[index_none_zero::]
            
            
            
            Ft_raw = Area[0:-1] / (time[1::] - time[0:-1])
            time_raw = time[0:-1]
            time_zero = time_raw[0]

            time_raw = time_raw - time_zero
            #下面两个变量是为了后面与Paper4作比较的
            time_raw_new = copy.deepcopy(time_raw)
            Ft_raw_new = copy.deepcopy(Ft_raw)
            

            delta_time_raw = time_raw[1] - time_raw[0] 
            
            # plt.figure(1)
            # plt.semilogx(time_raw,Ft_raw)
            # plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant],label='smooth')
            # plt.xlabel('time')
            # plt.ylabel('dA/dt')
            # # plt.xlim(10**(-2),10**(-1))
            # plt.legend()
            
            
            """去掉尾巴"""
            length = len(time_raw)
            time_raw = time_raw[0:length*4//5]
            Ft_raw = Ft_raw[0:length*4//5]
            
            

            """without apodization"""
            Ft_raw[1::] = Ft_raw[1::] - constant
            Ft_raw[0] = Ft_raw[0] - constant/2
            
            
            plt.figure(1)
            plt.semilogx(time_raw, Ft_raw,label='$F(t)-F_{smooth}(t)$')
            plt.plot([time_raw[0], time_raw[-1]], [0, 0])
            plt.xlabel('time')
            plt.ylabel('F(t)')
            # plt.legend()
            plt.grid()
            plt.savefig('./To_Meena/Sta_Minimum_Ft.png',dpi=450)
            """以上"""

            time_new = time_raw
            Ft_new = Ft_raw
            # # """加窗归零"""
            # lengthnew = len(time_new)
            # window = np.hanning(lengthnew*2//5)
            # try:
            #     Ft_new[lengthnew*4//5+1::] = Ft_new[lengthnew*4//5+1::] * window[lengthnew*1//5::] 
            # except ValueError:
            #     Ft_new[lengthnew*4//5::] = Ft_new[lengthnew*4//5::] * window[lengthnew*1//5::] 
                
            
            np.savetxt('./To_Meena/Sta_Minimum_time_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),time_new,delimiter=',')
            np.savetxt('./To_Meena/Sta_Minimum_Ft_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ft_new,delimiter=',')
            
            
            
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
            Ffsgn = constant * cw 
            Ff = Ff + Ffsgn
            
        
            Ffabs = np.abs(Ff)
            ThetaF = np.real(-complex(0,1)*np.log(Ff/Ffabs))
            # Ffabs = signal.savgol_filter(Ffabs,3,1)
            
            
            """以上"""
            
            plt.figure(2)
            
            plt.semilogx(freq, Ffabs,label = 'Numerical Result') #//FIXME
            # plt.plot(freq, FfabsOld)
            # plt.plot([freq[0],freq[-1]],[mu,mu],label='Macro magnification')
            # plt.plot(freq, Ffabs1, label='Cut')
        
            plt.grid()
            # plt.legend()
            # plt.ylim(1,3)
            plt.xlabel('f[Hz]')
            plt.ylabel('|F(f)|')
            plt.xlim(0.1, 3000)
            plt.savefig("./To_Meena/Sta_Minimum_Ffabs.png",dpi=450)
            """相位"""
            plt.figure(3)
            
            plt.semilogx(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME
            # plt.semilogx(freq, FfTheoryTheta, label='Geometric')
            # plt.plot(freq, ThetaF, label='Cut')
            # plt.loglog(test1,test2)
            # plt.ylim(-0.2,0.2)
            plt.grid()
            # plt.legend()
            plt.xlabel('f[Hz]')
            plt.ylabel(r'$\theta_F$')
            plt.xlim(0.1, 3000)
            plt.savefig("./To_Meena/Sta_Minimum_ThetaF.png",dpi=450)


            np.savetxt('./To_Meena/Sta_Minimum_Ffabs_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ffabs,delimiter=',')
            np.savetxt('./To_Meena/Sta_Minimum_ThetaF_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),ThetaF,delimiter=',')

        IndexSumimage += 1
