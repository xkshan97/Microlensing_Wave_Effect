'''
Author: your name
Date: 2021-10-16 16:04:33
LastEditTime: 2023-08-03 00:27:04
LastEditors: xk_shan xk_shan@mail.bnu.edu.cn
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

import numpy as np
import matplotlib.pyplot as plt
from GetMicroDiffrac import Diffrac
import struct
from scipy.interpolate import interp1d
from astropy.io import fits
import mpmath
from scipy.optimize import curve_fit
from tqdm import tqdm



def Theory(kappa,w, x1, x2, T0, k):
    mu = 1/(1 - kappa)**2
    return -(mu/w)*2**(-2 - (complex(0,1)*w)/2)*np.exp(1/2*complex(0,1)*w*(2*T0/k+ (x1**2 + x2**2)/np.sqrt(mu)))*np.abs(mu)**(1/2 - (complex(0,1)*w)/4)*np.abs(w)**(-1 + (complex(0,1)*w)/2)*mpmath.gamma(1 - (complex(0,1)*w)/2)*(-4*w**2/mu*mpmath.hyper([1/2 - (complex(0,1)*w)/4, 1 - (complex(0,1)*w)/4], [1/2, 1/2, 1], -1/16/mu*w**2*(x1**2 + x2**2)**2)*(mpmath.cosh(np.pi*w/4)*mpmath.sign(w/mu) + mpmath.sinh(np.pi*w/4)) + (2*complex(0,1) + w)*(x1**2 + x2**2)*np.abs(mu)**(-3/2)*np.abs(w)**3*mpmath.hyper([1 - (complex(0,1)*w)/4, 3/2 - (complex(0,1)*w)/4], [1, 3/2, 3/2], -1/16/mu*w**2*(x1**2 + x2**2)**2)*(mpmath.cosh(np.pi*w/4) + mpmath.sign(w/np.sqrt(mu))*mpmath.sinh(np.pi*w/4)))

#Minimum
kappa_star = "0.06"
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 0.7
gamma = -0.25
y1 = ""
y2 = ""
M_L = 1
z_L = 0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
sqrt_mu = np.sqrt(1 / (1 - kappa - gamma) / (1 - kappa + gamma))

#adptive
f1=open('./ResultMicro/MicroTimeLength_min.bin',"rb")
TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
f1.close()
TimeLengthFile_min = np.array(TimeLengthFile_min)
fileArea = "./ResultMicro/Microadptive_Area_min.bin"
fileTime = "./ResultMicro/Microadptive_Time_min.bin"
lens_freq = np.arange(0.1, 1000, 0.5)
Ff_abs_adp, ThetaF_adp = Diffrac().Minimum(M_L, z_L, kappa, gamma, fileArea, fileTime, TimeLengthFile_min[0], lens_freq, kappa_star)

Ftnew_adp = np.loadtxt('./Ftnew/FtnewMinimum_' + kappa_star + '.csv', delimiter=',')
timenew_adp = np.loadtxt('./Ftnew/timenewMinimum_' + kappa_star + '.csv', delimiter=',')
timenew_adp = timenew_adp
Func_Ft_adp = interp1d(timenew_adp, Ftnew_adp)

#No-shear 理论对比

f2=open(fileTime,"rb")
time=struct.unpack("d"*TimeLengthFile_min[0], f2.read(8*TimeLengthFile_min[0]))
f2.close()
time = np.array(time)

w = 4*G*M_L*M_sun*(1 + z_L) * lens_freq * 2 * np.pi/c**3
ResTheory = []
FfTheoryAbs = []
FfTheoryTheta = [] 
for i in range(len(w)):
    ResTheory.append(Theory(kappa, w[i], 0.1, 0, -time[0], coeffi))
    FfTheoryAbs.append(float(np.abs(ResTheory[i])))
    FfTheoryTheta.append(np.real(-complex(0,1)*np.log(complex(float(ResTheory[i].real), float(ResTheory[i].imag))/FfTheoryAbs[i])))


#网格数据
f2=open('./ResultMinimum_100/LayersFile_min_0.28_0.28_0.11_1.04_1.00.bin',"rb")
LayerNum_Res01=struct.unpack("l"*97969, f2.read(8*97969))
f2.close()
LayerNum_Res01 = np.array(LayerNum_Res01)
# x1_set01 = np.arange(-500 + 0.723594 / 2, 500, 0.723594)
# x2_set01 = x1_set01
LayerNum_Res01 = LayerNum_Res01.reshape(int(np.sqrt(len(LayerNum_Res01))), -1)
plt.imshow(LayerNum_Res01, cmap='gnuplot2', vmin=0, vmax=7)

plt.xticks(np.arange(0,len(x1_set01), 200), x1_set01[::200].astype(int))
plt.yticks(np.arange(0,len(x1_set01), 200), x1_set01[::200].astype(int))
plt.colorbar()
# plt.xlim(590, 780)
# plt.ylim(600, 830)

plt.savefig('Fig/Micro_min_Layers.png', dpi=450)
plt.close()

#zoom in
f2=open('./ResultMicro/LayersFile_min.bin',"rb")
LayerNum_Res01=struct.unpack("l"*1909924, f2.read(8*1909924))
f2.close()
LayerNum_Res01 = np.array(LayerNum_Res01)
x1_set01 = np.arange(-500 + 0.723594 / 2, 500, 0.723594)
x2_set01 = x1_set01
LayerNum_Res01 = LayerNum_Res01.reshape(len(x1_set01), -1)
plt.imshow(LayerNum_Res01, cmap='gnuplot2')

plt.xticks(np.arange(0,len(x1_set01), 50), x1_set01[::50].astype(int))
plt.yticks(np.arange(0,len(x1_set01), 50), x1_set01[::50].astype(int))
plt.colorbar()
plt.xlim(550, 830)
plt.ylim(550, 830)

plt.savefig('Fig/Micro_min_Layers_zoom_in.png', dpi=450)
plt.close()


# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogx(timenew_adp, Ftnew_adp , label='adaptive')
plt.plot(timenew_adp, [0]*len(timenew_adp))
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-5), 10**(-1))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-100000,40000)
plt.grid()
plt.xlabel('t[s]')
plt.ylabel('$\\frac{|F(t) - F_s(t)|}{F_s(t)}$')
plt.legend()
plt.savefig('./Fig/Min_F_t_error.png', dpi=450)
plt.close()


# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogx(lens_freq, Ff_abs_adp, label='adaptive')
# plt.semilogx(lens_freq, FfTheoryAbs, '--',label='Analytic')
plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$\\frac{||F(f)| - |F(f)_\mathrm{theory}||}{|F(f)_\mathrm{theory}|}$')
plt.legend()
plt.savefig('./Fig/Min_F_f.png', dpi=450)
plt.close()



# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogx(lens_freq, ThetaF_adp, label='adaptive')
# plt.semilogx(lens_freq, , '--',label='Analytic')
# plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$|\\theta_f - \\theta_\mathrm{theory}|$')
plt.legend()
plt.savefig('./Fig/Min_F_theta.png', dpi=450)
plt.close()


#Saddle
kappa_star = "0"
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 0.8
gamma = 0.25
y1 = ""
y2 = ""
M_L = 1
z_L = 0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
sqrt_mu = np.sqrt(1 / np.abs(1 - kappa - gamma) / np.abs(1 - kappa + gamma))

#adptive
f1=open('./ResultMicro/TimeLength_sad.bin',"rb")
TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
f1.close()
TimeLengthFile_min = np.array(TimeLengthFile_min)
fileArea = "./ResultMicro/adptive_Area_sad.bin"
fileTime = "./ResultMicro/adptive_Time_sad.bin"
lens_freq = np.arange(0.1, 1000, 0.5)
Ff_abs_adp, ThetaF_adp = Diffrac().Saddle(M_L, z_L, kappa, gamma, fileArea, fileTime, TimeLengthFile_min[0], 216.081, 648.244, lens_freq, kappa_star)

Ftnew_adp = np.loadtxt('./Ftnew/FtnewSaddle_' + kappa_star + '.csv', delimiter=',')
timenew_adp = np.loadtxt('./Ftnew/timenewSaddle_' + kappa_star + '.csv', delimiter=',')
timenew_adp = timenew_adp
Func_Ft_adp = interp1d(timenew_adp, Ftnew_adp)

FtTheory_adp = np.loadtxt('./Ftnew/FtTheorySaddle_' + kappa_star + '.csv', delimiter=',')



X10New = 137.662# + 20
X20New = 412.985 + 20
X10New = X20New
#网格数据
f2=open('./ResultMicro/LayersFile_sad.bin',"rb")
LayerNum_Res01=struct.unpack("l"*7639696, f2.read(8*7639696))
f2.close()
LayerNum_Res01 = np.array(LayerNum_Res01)
x1_set01 = np.arange(-1000. + 0.723603 / 2, 1000, 0.723603)
x2_set01 = x1_set01
x1range = np.where((x1_set01<X10New)&(x1_set01>-X10New))[0]
x2range = np.where((x2_set01<X20New)&(x2_set01>-X20New))[0]

LayerNum_Res01 = LayerNum_Res01.reshape(len(x1_set01), -1).T
pixel_num = 0
for i in range(len(LayerNum_Res01)):
    for j in range(len(LayerNum_Res01[i])):
        if i > x1range[0] and i < x1range[-1] and j > x2range[0] and j < x2range[-1]:
            
            pixel_num += 4 ** LayerNum_Res01[i][j]
uniform_pixel_num = (X20New * 2 / (0.723603/2**7)) * (X10New * 2 / (0.723603/2**7))
print(uniform_pixel_num / pixel_num)

plt.figure(figsize=(10,10))
plt.xticks(np.arange(0,len(x1_set01), 400), x1_set01[::400].astype(int))
plt.yticks(np.arange(0,len(x1_set01), 400), x1_set01[::400].astype(int))
plt.imshow(LayerNum_Res01, cmap='gnuplot2')
plt.colorbar()
plt.xlim(x1range[0], x1range[-1])
plt.ylim(x2range[0], x2range[-1])
plt.savefig('./Fig/Micro_sad_Layers.png', dpi=450)
plt.close()

#zoom in
f2=open('./ResultMicro/LayersFile_sad.bin',"rb")
LayerNum_Res01=struct.unpack("l"*7639696, f2.read(8*7639696))
f2.close()
LayerNum_Res01 = np.array(LayerNum_Res01)
x1_set01 = np.arange(-1000. + 0.723603 / 2, 1000, 0.723603)
x2_set01 = x1_set01
LayerNum_Res01 = LayerNum_Res01.reshape(len(x1_set01), -1)

plt.xticks(np.arange(0,len(x1_set01), 50), x1_set01[::50].astype(int))
plt.yticks(np.arange(0,len(x1_set01), 50), x1_set01[::50].astype(int))
plt.imshow(LayerNum_Res01, cmap='gnuplot2')
plt.colorbar()
plt.xlim(1260, 1500)
plt.ylim(1260, 1500)
plt.savefig('./Fig/Micro_sad_Layers_zoom_in.png', dpi=450)
plt.close()




#Saddle均匀
kappa_star = "0"
M_sun = 1.9884099 * 10**30
G = 6.6743*10**(-11)
c = 2.9979246*10**8
kappa = 0.875
gamma = 0.325
y1 = ""
y2 = ""
M_L = 100
z_L = 0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
sqrt_mu = np.sqrt(1 / np.abs(1 - kappa - gamma) / np.abs(1 - kappa + gamma))


f1=open('./ResultMicro_Uniform/TimeLength_sad.bin',"rb")
TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
f1.close()
TimeLengthFile_min = np.array(TimeLengthFile_min)
fileArea = "./ResultMicro_Uniform/adptive_Area_sad.bin"
fileTime = "./ResultMicro_Uniform/adptive_Time_sad.bin"
lens_freq = np.arange(0.1, 2000, 0.5)
Ff_abs_uniform, ThetaF_uniform = Diffrac().Saddle(M_L, z_L, kappa, gamma, fileArea, fileTime, TimeLengthFile_min[0], 187.487, 281.23, lens_freq, kappa_star)

Ftnew_uniform = np.loadtxt('./Ftnew/FtnewSaddle_' + kappa_star + '.csv', delimiter=',')
timenew_uniform = np.loadtxt('./Ftnew/timenewSaddle_' + kappa_star + '.csv', delimiter=',')
Func_Ft_uniform = interp1d(timenew_uniform, Ftnew_uniform)

FtTheory_uniform = np.loadtxt('./Ftnew/FtTheorySaddle_' + kappa_star + '.csv', delimiter=',')


# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogy(timenew_uniform, np.abs(Ftnew_uniform) / FtTheory_uniform, label='Uniform $10^{-5}$')
plt.semilogy(timenew_adp, np.abs(Ftnew_adp) / FtTheory_adp, label='adaptive 0.1')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('t[s]')
plt.ylabel('$\\frac{|F(t) - F_s(t)|}{F_s(t)}$')
plt.legend()
plt.savefig('./Fig/Sad_F_t_error.png', dpi=450)
plt.close()

# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogy(lens_freq, np.abs(ThetaF_uniform - ThetaF_adp))
# plt.semilogx(lens_freq, Ff_abs_adp, '--', label='adaptive 0.1')
# plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$|\\theta_\mathrm{adp} - \\theta_\mathrm{Uni}|$')
plt.legend()
plt.savefig('./Fig/Sad_theta_f.png', dpi=450)
plt.close()


# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogy(lens_freq, np.abs(Ff_abs_uniform - Ff_abs_adp) / Ff_abs_uniform)
# plt.semilogx(lens_freq, Ff_abs_adp, '--', label='adaptive 0.1')
# plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$\\frac{||F(f)| - \sqrt{\mu}|}{\sqrt{\mu}}$')
plt.legend()
plt.savefig('./Fig/Sad_F_f.png', dpi=450)
plt.close()




# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogx(lens_freq, Ff_abs_adp, label='adaptive')
plt.semilogx(lens_freq, Ff_abs_uniform, label='Uniform')
# plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$|F(f)|$')
plt.legend()
plt.savefig('./Fig/abs_Sad_F_f.png', dpi=450)
plt.close()


# plt.semilogx(timenew, Ftnew, label='tradition')
plt.semilogx(lens_freq, ThetaF_adp, label='adaptive')
plt.semilogx(lens_freq, ThetaF_uniform, label='Uniform')
# plt.semilogx(lens_freq, , '--',label='Analytic')
# plt.plot(lens_freq, [sqrt_mu]*len(lens_freq), '--', label='Theory')
# plt.scatter(0.0004734142, Func_Ft(0.0004734142))
# plt.xlim(4*10**(-4), 6*10**(-4))
# plt.xlim(10**(-4), 0.09)
# plt.ylim(-10,10)
plt.grid()
plt.xlabel('f[Hz]')
plt.ylabel('$|\\theta_f|$')
plt.legend()
plt.savefig('./Fig/abs_sad_F_theta.png', dpi=450)
plt.close()




#窗函数方法 没改完因为程序没运行完
cw = coeffi/2/np.pi
f1=open('./ResultMicro_Uniform/TimeLength_sad.bin',"rb")
TimeLengthFile_min=struct.unpack("l"*1, f1.read(8*1))
f1.close()
TimeLengthFile_min = np.array(TimeLengthFile_min)
fileArea = "./ResultMicro_Uniform/adptive_Area_sad.bin"
fileTime = "./ResultMicro_Uniform/adptive_Time_sad.bin"

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

time_raw = time_raw #- time_raw[0]


delta_time_raw = time_raw[1] - time_raw[0] 
mu = abs(1/((1-kappa)**2 - gamma**2))**0.5
coeffi = 4 * G * M_L * M_sun * (1 + z_L) / c**3
constant = mu/coeffi

"""延拓对数图像"""

length = len(time_raw)

def func2(t, b): #理论表达式
    
    a = -2*constant
    return a*np.log(abs(t)) + b
    
        

"""后面延拓"""
indexend = np.arange(length*7//10,length*9//10) #取这个范围的时间和dS/dt来拟合公式func2，得到参数b的值。
popt, pcov = curve_fit(func2, time_raw[indexend], Ft_raw[indexend])

print(popt)
# a = popt[0] 
b = popt[0]
delta_t = time_raw[1] - time_raw[0]
time_new2 = [time_raw[indexend[-1]]] #新的时间延迟从上面的0.13秒开始

Ft_new2 = [func2(time_new2[0],b)]
i_new = 0
while Ft_new2[i_new] >= 0: #直到dS/dt=0，否则一直往外延拓
    time_new2.append(time_new2[i_new]+delta_t)
    Ft_new2.append(func2(time_new2[-1],b))
    i_new += 1
    
time_new2 = np.array(time_new2)
Ft_new2 = np.array(Ft_new2)
time_raw = np.append(time_raw[0:indexend[-1]],time_new2) #把延拓的部分和前面部分粘合在一起。
Ft_raw = np.append(Ft_raw[0:indexend[-1]],Ft_new2)
plt.figure(1)
plt.plot(time_raw, Ft_raw)
"""以上"""
"""理论值"""
# time_theory1 = [delta_t]
# Ft_theory1 = [func2(time_theory1[0],b)]
# i_theory1 = 0
# while Ft_theory1[i_theory1] >= 0:
#     time_theory1.append(time_theory1[i_theory1] + delta_t)
#     Ft_theory1.append(func2(time_theory1[-1],b))
#     i_theory1 += 1

# plt.loglog(time_theory1, Ft_theory1)
# plt.loglog(time_raw,Ft_raw)

"""前面延拓"""
def func1(t,b):
    
    a = -2*constant
    return a*np.log(abs(t)) + b

indexstart = indexend = np.arange(length*1//10,length*3//10)

popt, pcov = curve_fit(func1, time_raw[indexstart], Ft_raw[indexstart])

print(popt)
# a = popt[0] 
b = popt[0]
delta_t = time_raw[1] - time_raw[0]
time_new1 = [time_raw[indexstart[0]]-delta_t]

Ft_new1 = [func1(time_new1[0],b)]
i_new = 0
while Ft_new1[i_new] >= 0:
    time_new1.append(time_new1[i_new]-delta_t)
    Ft_new1.append(func1(time_new1[-1],b))
    i_new += 1

time_new1.reverse()
Ft_new1.reverse()

time_new1 = np.array(time_new1)
Ft_new1 = np.array(Ft_new1)
time_raw = np.append(time_new1,time_raw[indexstart[0]::])
Ft_raw = np.append(Ft_new1,Ft_raw[indexstart[0]::])
plt.figure(2)
plt.plot(time_raw, Ft_raw)

# plt.semilogx(time_raw,Ft_raw)
# plt.semilogx([time_raw[0],time_raw[-1]],[constant,constant])
# plt.xlabel('time')
# plt.ylabel('dA/dt')



"""前后等长"""
index0 = np.where((delta_time_raw/2>= time_raw)&(time_raw>=-delta_time_raw/2))[0][0]
lengthraw = len(Ft_raw)
lengtha = index0
lengthb = lengthraw - lengtha
if(lengtha > lengthb):
    lengthtmp = lengtha - lengthb
    time_raw = time_raw[lengthtmp::]
    Ft_raw = Ft_raw[lengthtmp::]
else:
    lengthtmp = lengthb - lengtha
    time_raw = time_raw[0:lengthraw - lengthtmp]
    Ft_raw = Ft_raw[0:lengthraw - lengthtmp]



"""with apodization"""
windows = np.hanning(len(Ft_raw)) #* 2
Ft_raw = windows * Ft_raw 
# plt.plot(time_raw,windows)
# plt.semilogx(time_raw, Ft_raw)


time_new = time_raw
Ft_new = Ft_raw

# plt.semilogx(time_new, Ft_new)

FrealWindow = []
FimagWindow = []
lens_freq = np.arange(0.1, 2000, 0.5)
omegafreq = lens_freq[0:200] * 2 * np.pi
for i in tqdm(range(len(omegafreq))):
    FrealWindow.append(np.sum(Ft_new*np.cos(omegafreq[i]*time_new))*delta_time_raw)
    FimagWindow.append(np.sum(Ft_new*np.sin(omegafreq[i]*time_new))*delta_time_raw)
    
FfWindow = np.array(FrealWindow) + complex(0,1)*np.array(FimagWindow)
FfWindow = FfWindow* omegafreq[0:len(FfWindow)] / complex(0,1)*cw


FfabsWindow = np.abs(FfWindow)
ThetaFWindow = np.real(-complex(0,1)*np.log(FfWindow/FfabsWindow)) 
# Ffabs = signal.savgol_filter(Ffabs[1::], 9, 1)

# FfabsWindow = FfabsWindow/np.average(FfabsWindow[0])*mu

plt.figure(3)
# plt.semilogx(freq, FfTheoryAbs,label = 'Geometric Result')

plt.semilogx(lens_freq[0:len(FfWindow)], FfabsWindow,label = 'Window Numerical Result') #//FIXME
plt.semilogx(lens_freq, Ff_abs_adp, '--', label='adaptive')
plt.xlim(1,2000)
# plt.loglog(test1,test2)
plt.ylim(0.01,3.5)
plt.grid()
plt.legend()
plt.xlabel('f[Hz]')
plt.ylabel('|F(f)|')
plt.savefig('Fig/windows_sad_Ff.png', dpi=450)
plt.close()
# plt.savefig("Ffembed.png",dpi=450)
"""相位"""
plt.semilogx(lens_freq[0:len(FfWindow)], ThetaFWindow , label = 'Numerical Result') #//FIXME
plt.semilogx(lens_freq, ThetaF_adp, '--',label='adaptive')
plt.xlim(1,2000)
# plt.loglog(test1,test2)
# plt.ylim(-0.2,0.2)
plt.grid()
plt.legend()
plt.xlabel('f[Hz]')
plt.ylabel('$\\theta_F$')
plt.savefig('./Fig/windows_sad_ThetaF.png', dpi=450)
plt.close()

np.savetxt('./Fig/Windows_sad_Ff.csv', FfabsWindow, delimiter=',')
np.savetxt('./Fig/Windows_sad_ThetaF.csv', ThetaFWindow, delimiter=',')






#验证单网格近似的
# f1=open('./Test/TFour.bin',"rb")
# TFour=struct.unpack("d"*6, f1.read(8*6))
# f1.close()
# TFour = np.array(TFour)

# f2=open('./Test/DSDT.bin',"rb")
# DSDT=struct.unpack("d"*int(TFour[5]), f2.read(8*int(TFour[5])))
# f2.close()
# DSDT = np.array(DSDT)

# f3=open('./Test/Time.bin',"rb")
# Time=struct.unpack("d"*int(TFour[5]), f3.read(8*int(TFour[5])))
# f3.close()
# Time = np.array(Time)

# plt.plot(Time, DSDT)
# for i in range(4):
#     plt.plot([TFour[i], TFour[i]], [0,TFour[4]], label='t'+str(i))
# plt.plot([TFour[2], TFour[3]], [TFour[4], TFour[4]], '--', linewidth = 0.5, label = 'Platform')
# plt.legend()
# plt.xlabel('t[s]')
# plt.ylabel('dS/dt')
# plt.grid()
# plt.ylim(0, 1500)