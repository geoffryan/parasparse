from os import system
import numpy as np
import matplotlib.pyplot as plt

num = 21
minsize = 5
maxsize = 15

n_sine = np.logspace(minsize,maxsize,num=num,base=2).astype(int)
n_grav2d = np.logspace(minsize,maxsize,num=num,base=2).astype(int)
n_grav3d = np.logspace(minsize,maxsize,num=num,base=2).astype(int)

Xsine = np.logspace(np.log2(n_sine[0]), np.log2(n_sine[-1]), num=200, base=2)
Xgrav2d = np.logspace(np.log2(n_grav2d[0]), np.log2(n_grav2d[-1]), num=200, base=2)
Xgrav3d = np.logspace(np.log2(n_grav2d[0]), np.log2(n_grav2d[-1]), num=200, base=2)

system("rm -f time_sine_solve.out")
system("rm -f time_grav2d_solve.out")
system("rm -f time_grav3d_solve.out")

for i in range(num):
	system("./ParaSparse sine solve "+str(n_sine[i]))
	system("./ParaSparse grav2d solve "+str(n_grav2d[i]))
	system("./ParaSparse grav3d solve "+str(n_grav3d[i]))
	
print "Simulation done, calculating."
	
N_sine, size_sine, iters_sine, dot_sine, ticks_sine, time_sine = np.loadtxt("time_sine_solve.out", unpack=True)
N_grav2d, size_grav2d, iters_grav2d, dot_grav2d, ticks_grav2d, time_grav2d = np.loadtxt("time_grav2d_solve.out", unpack=True)
N_grav3d, size_grav3d, iters_grav3d, dot_grav3d, ticks_grav3d, time_grav3d = np.loadtxt("time_grav3d_solve.out", unpack=True)

fit_sine_iters = np.polyfit(np.log(N_sine[5:]), np.log(iters_sine[5:]), 1) 
fit_sine_time = np.polyfit(np.log(N_sine[5:]), np.log(time_sine[5:]), 1)
fit_sine_ticks = np.polyfit(np.log(N_sine[5:]), np.log(ticks_sine[5:]), 1)
fit_grav2d_iters = np.polyfit(np.log(N_grav2d[5:]), np.log(iters_grav2d[5:]), 1) 
fit_grav2d_time = np.polyfit(np.log(N_grav2d[5:]), np.log(time_grav2d[5:]), 1)
fit_grav2d_ticks = np.polyfit(np.log(N_grav2d[5:]), np.log(ticks_grav2d[5:]), 1)
fit_grav3d_iters = np.polyfit(np.log(N_grav3d[5:]), np.log(iters_grav3d[5:]), 1) 
fit_grav3d_time = np.polyfit(np.log(N_grav3d[5:]), np.log(time_grav3d[5:]), 1)
fit_grav3d_ticks = np.polyfit(np.log(N_grav3d[5:]), np.log(ticks_grav3d[5:]), 1)

p_sine_iters = np.poly1d(fit_sine_iters)
p_sine_time = np.poly1d(fit_sine_time)
p_sine_ticks = np.poly1d(fit_sine_ticks)
p_grav2d_iters = np.poly1d(fit_grav2d_iters)
p_grav2d_time = np.poly1d(fit_grav2d_time)
p_grav2d_ticks = np.poly1d(fit_grav2d_ticks)
p_grav3d_iters = np.poly1d(fit_grav3d_iters)
p_grav3d_time = np.poly1d(fit_grav3d_time)
p_grav3d_ticks = np.poly1d(fit_grav3d_ticks)

Ysine_iters = np.exp(p_sine_iters(np.log(Xsine)))
Ysine_time = np.exp(p_sine_time(np.log(Xsine)))
Ysine_ticks = np.exp(p_sine_ticks(np.log(Xsine)))
Ygrav2d_iters = np.exp(p_grav2d_iters(np.log(Xgrav2d)))
Ygrav2d_time = np.exp(p_grav2d_time(np.log(Xgrav2d)))
Ygrav2d_ticks = np.exp(p_grav2d_ticks(np.log(Xgrav2d)))
Ygrav3d_iters = np.exp(p_grav3d_iters(np.log(Xgrav3d)))
Ygrav3d_time = np.exp(p_grav3d_time(np.log(Xgrav3d)))
Ygrav3d_ticks = np.exp(p_grav3d_ticks(np.log(Xgrav3d)))

print "Plotting!"

ax = plt.subplot(331)
plt.loglog(N_sine, iters_sine, '.')
plt.loglog(Xsine, Ysine_iters)
plt.ylabel("Iterations")
plt.text(0.3, 0.9, "O = {0:g}".format(fit_sine_iters[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(334)
plt.loglog(N_sine, time_sine, '.')
plt.loglog(Xsine, Ysine_time)
plt.ylabel("Time (s)")
plt.text(0.3, 0.9, "O = {0:g}".format(fit_sine_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(337)
plt.loglog(N_sine, ticks_sine, '.')
plt.loglog(Xsine, Ysine_ticks)
plt.ylabel("Ticks")
plt.text(0.3, 0.9, "O = {0:g}".format(fit_sine_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Sine")

ax = plt.subplot(332)
plt.loglog(N_grav2d, iters_grav2d, '.')
plt.loglog(Xgrav2d, Ygrav2d_iters)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav2d_iters[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(335)
plt.loglog(N_grav2d, time_grav2d, '.')
plt.loglog(Xgrav2d, Ygrav2d_time)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav2d_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(338)
plt.loglog(N_grav2d, ticks_grav2d, '.')
plt.loglog(Xgrav2d, Ygrav2d_ticks)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav2d_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Grav 2D")

ax = plt.subplot(333)
plt.loglog(N_grav3d, iters_grav3d, '.')
plt.loglog(Xgrav3d, Ygrav3d_iters)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav3d_iters[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(336)
plt.loglog(N_grav3d, time_grav3d, '.')
plt.loglog(Xgrav3d, Ygrav3d_time)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav3d_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(339)
plt.loglog(N_grav3d, ticks_grav3d, '.')
plt.loglog(Xgrav3d, Ygrav3d_ticks)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav3d_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Grav 3D")

plt.show()