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
from tqdm import tqdm
import struct

M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8


# f0=open('./MicroField_88/AveMassAndNum_' + str(round(kappa, 2)) + '_' + str(round(gamma, 2)) + '_' + str(round(1.1 * kappaStar_stellar, 2)) + '_' + str(round(LensRedshift, 2)) + '_1.00.bin',"rb")
f0=open('./test/Coordinate_Microlens_at_Minimum.bin',"rb")
AveMassAndNum=struct.unpack("d"*16738*2, f0.read(8*16738*2))
f0.close()
M_L = AveMassAndNum[0]

coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3

            if 1 - kappa - gamma > 0:                
                constant = 2*np.pi*mu/coeffi
                cw = coeffi/2/np.pi


                
                f1=open('./ResultMinimum_88/TimeLength_min_%.2f_%.2f_%.2f_%.2f_1.00.bin'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),"rb")
                TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
                f1.close()
                TimeLengthFile_min = np.array(TimeLengthFile_min)
                fileArea = "./ResultMinimum_88/adptive_Area_min_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
                fileTime = "./ResultMinimum_88/adptive_Time_min_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
                
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
                plt.plot([time_raw[0], time_raw[-1]], [0, 0], 'k--')
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
                    
                
                np.savetxt('./To_Meena_88/Sta_Minimum_time_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),time_new,delimiter=',')
                np.savetxt('./To_Meena_88/Sta_Minimum_Ft_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ft_new,delimiter=',')
                
                
                
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
                
                """相位"""
                plt.figure(3)
                plt.semilogx(freq, ThetaF, '--' , label = 'Numerical Result') #//FIXME


                np.savetxt('./To_Meena_88/Sta_Minimum_Ffabs_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ffabs,delimiter=',')
                np.savetxt('./To_Meena_88/Sta_Minimum_ThetaF_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),ThetaF,delimiter=',')
            else:
                
                f0=open('./ResultSaddle_88/X1020New_%.2f_%.2f_%.2f_%.2f_1.00.bin'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),"rb")
                X10X20New=struct.unpack("d"*2, f0.read(8*2))
                f0.close()

                
                constant = mu/coeffi
                cw = coeffi/2/np.pi
                
                
                f1=open('./ResultSaddle_88/TimeLength_sad_%.2f_%.2f_%.2f_%.2f_1.00.bin'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),"rb")
                TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
                f1.close()
                TimeLengthFile_min = np.array(TimeLengthFile_min)
                fileArea = "./ResultSaddle_88/adptive_Area_sad_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
                fileTime = "./ResultSaddle_88/adptive_Time_sad_%.2f_%.2f_%.2f_%.2f_1.00.bin"%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift)
                
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
                    
                    
                np.savetxt('./To_Meena_88/Sta_Saddle_time_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),time_new,delimiter=',')
                np.savetxt('./To_Meena_88/Sta_Saddle_Ft_new_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ft_subtract,delimiter=',')
                
                plt.figure(4)
                plt.plot(time_new,Ft_subtract)
                plt.plot([time_new[0], time_new[-1]], [0, 0], 'k--')

                


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
                
                
                plt.figure(5)
                plt.semilogx(freq, Ffabs, label = 'Sign') #//FIXME
                
                plt.figure(6)
                plt.semilogx(freq, ThetaF, label = 'Sign') #//FIXME


                
                np.savetxt('./To_Meena_88/Sta_Saddle_Ffabs_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),Ffabs,delimiter=',')
                np.savetxt('./To_Meena_88/Sta_Saddle_ThetaF_%.2f_%.2f_%.2f_%.2f.csv'%(kappa, gamma, 1.1*kappaStar_stellar, LensRedshift),ThetaF,delimiter=',')

                
        IndexSumimage += 1

plt.figure(1)
plt.xlabel('time')
plt.ylabel('F(t)')
plt.grid()
plt.savefig('./To_Meena_88/Sta_Minimum_Ft.png',dpi=450)

plt.figure(2)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('|F(f)|')
plt.savefig("./To_Meena_88/Sta_Minimum_Ffabs.png",dpi=450)

plt.figure(3)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel(r'$\theta_F$')
plt.savefig("./To_Meena_88/Sta_Minimum_ThetaF.png",dpi=450)


plt.figure(4)
plt.xlabel('time')
plt.ylabel('dA/dt')
plt.grid()
plt.savefig('./To_Meena_88/Sta_Saddle_Ft.png',dpi=450)

plt.figure(5)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('|F(f)|')
plt.savefig("./To_Meena_88/Sta_Saddle_Ffabs.png",dpi=450)

plt.figure(6)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel(r'$\theta_F$')
plt.savefig("./To_Meena_88/Sta_Saddle_ThetaF.png",dpi=450)
