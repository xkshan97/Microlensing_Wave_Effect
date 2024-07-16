import matplotlib.pyplot as plt
import numpy as np



"""Minimum"""
Ffabs_adaptive_min = np.loadtxt('./ResultMinimum/Ffabs_adaptive.csv', delimiter=',')
Theta_F_adaptive_min = np.loadtxt('./ResultMinimum/ThetaF_adaptive.csv', delimiter=',')

Ffabs_direct_min = np.loadtxt('./ResultMinimum/Ffabs_uniform_Direct.csv', delimiter=',')
Theta_F_direct_min = np.loadtxt('./ResultMinimum/ThetaF_uniform_Direct.csv', delimiter=',')

Ffabs_uniform_min = np.loadtxt('./ResultMinimum/Ffabs_uniform.csv', delimiter=',')
Theta_F_uniform_min = np.loadtxt('./ResultMinimum/ThetaF_uniform.csv', delimiter=',')

omegafreq = np.arange(0.1*2*np.pi,2000*2*np.pi,2*np.pi)
freq = omegafreq/2/np.pi

"""minimum |Ff|"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, Ffabs_adaptive_min, label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('|F(f)|')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(Ffabs_adaptive_min - Ffabs_direct_min) / Ffabs_adaptive_min)
# sub_ax.semilogy(freq, np.abs(Ffabs_direct_min - Ffabs_uniform_min) / Ffabs_uniform_min, '--')
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\\frac{\Delta |F|}{|F|}$')
plt.savefig('./Fig/Ff_min.png', dpi=450)

"""Theta F"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, Theta_F_adaptive_min, label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('$\\theta_F$')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(Theta_F_adaptive_min - Theta_F_direct_min))
# sub_ax.semilogy(freq, np.abs(Theta_F_direct_min - Theta_F_uniform_min))
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\Delta |\\theta_F|$')
plt.savefig('./Fig/Theta_F_min.png', dpi=450)



"""Saddle"""
Ffabs_adaptive_sad = np.loadtxt('./ResultSaddle/Ffabs_adaptive.csv', delimiter=',')
Theta_F_adaptive_sad = np.loadtxt('./ResultSaddle/ThetaF_adaptive.csv', delimiter=',')

Ffabs_direct_sad = np.loadtxt('./ResultSaddle/Ffabs_uniform_Direct.csv', delimiter=',')
Theta_F_direct_sad = np.loadtxt('./ResultSaddle/ThetaF_uniform_Direct.csv', delimiter=',')

Ffabs_uniform_sad = np.loadtxt('./ResultSaddle/Ffabs_uniform.csv', delimiter=',')
Theta_F_uniform_sad = np.loadtxt('./ResultSaddle/ThetaF_uniform.csv', delimiter=',')

omegafreq = np.arange(0.1*2*np.pi,2000*2*np.pi,2*np.pi)
freq = omegafreq/2/np.pi




"""saddle |Ff|"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, Ffabs_adaptive_sad, label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('|F(f)|')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(Ffabs_adaptive_sad - Ffabs_direct_sad) / Ffabs_direct_sad)
# sub_ax.semilogy(freq, np.abs(Ffabs_direct_sad - Ffabs_uniform_sad) / Ffabs_uniform_sad, '--')
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\\frac{\Delta |F|}{|F|}$')
plt.savefig('./Fig/Ff_sad.png', dpi=450)

"""Theta F"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, Theta_F_adaptive_sad, label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('$\\theta_F$')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(np.unwrap(Theta_F_adaptive_sad) - np.unwrap(Theta_F_direct_sad)))
# sub_ax.semilogy(freq, np.abs(np.unwrap(Theta_F_direct_sad) - np.unwrap(Theta_F_uniform_sad)))
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\Delta |\\theta_F|$')
plt.savefig('./Fig/Theta_F_sad.png', dpi=450)



"""Maximum"""
Ffabs_adaptive_max = np.loadtxt('./ResultMaximum/Ffabs_adaptive.csv', delimiter=',')
Theta_F_adaptive_max = np.loadtxt('./ResultMaximum/ThetaF_adaptive.csv', delimiter=',')

Ffabs_direct_max = np.loadtxt('./ResultMaximum/Ffabs_uniform_Direct.csv', delimiter=',')
Theta_F_direct_max = np.loadtxt('./ResultMaximum/ThetaF_uniform_Direct.csv', delimiter=',')

Ffabs_uniform_max = np.loadtxt('./ResultMaximum/Ffabs_uniform.csv', delimiter=',')
Theta_F_uniform_max = np.loadtxt('./ResultMaximum/ThetaF_uniform.csv', delimiter=',')

omegafreq = np.arange(0.1*2*np.pi,2000*2*np.pi,2*np.pi)
freq = omegafreq/2/np.pi




"""Maximum |Ff|"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, Ffabs_adaptive_max, label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('|F(f)|')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(Ffabs_adaptive_max - Ffabs_direct_max) / Ffabs_direct_max)
# sub_ax.semilogy(freq, np.abs(Ffabs_direct_max - Ffabs_uniform_max) / Ffabs_uniform_max, '--')
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\\frac{\Delta |F|}{|F|}$')
plt.savefig('./Fig/Ff_max.png', dpi=450)

"""Theta F"""
plt.figure(figsize=(6.5,5))
plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
grid = plt.GridSpec(4, 1, hspace = 0)
main_ax = plt.subplot(grid[0:3,0])

main_ax.semilogx(freq, np.unwrap(Theta_F_adaptive_max), label = 'Adaptive')
# main_ax.semilogx(freq, Ffabs_uniform_min, '--', label = 'Uniform')

main_ax.grid()
main_ax.legend()
main_ax.set_ylabel('$\\theta_F$')

sub_ax = plt.subplot(grid[3,0], sharex=main_ax)
sub_ax.semilogy(freq, np.abs(np.unwrap(Theta_F_adaptive_max) - np.unwrap(Theta_F_direct_max)))
# sub_ax.semilogy(freq, np.abs(np.unwrap(Theta_F_direct_max) - np.unwrap(Theta_F_uniform_max)))
sub_ax.grid()
sub_ax.set_ylim(10**(-10), 1)
# sub_ax.legend()
sub_ax.set_xlabel('$f$[Hz]')
sub_ax.set_ylabel('$\Delta |\\theta_F|$')
plt.savefig('./Fig/Theta_F_max.png', dpi=450)