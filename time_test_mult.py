from os import system
import numpy as np
import matplotlib.pyplot as plt

num = 21
minsize = 5
maxsize = 15

n_sine = np.logspace(minsize,maxsize,num=num,base=2).astype(int)
n_grav1d = np.logspace(minsize,maxsize,num=num,base=2).astype(int)
n_grav2d = np.logspace(minsize,maxsize,num=num,base=2).astype(int)

Xsine = np.logspace(np.log2(n_sine[0]), np.log2(n_sine[-1]), num=200, base=2)
Xgrav1d = np.logspace(np.log2(n_grav1d[0]), np.log2(n_grav1d[-1]), num=200, base=2)
Xgrav2d = np.logspace(np.log2(n_grav2d[0]), np.log2(n_grav2d[-1]), num=200, base=2)

system("rm time_sine_mult.out")
system("rm time_grav1d_mult.out")
system("rm time_grav2d_mult.out")

for i in range(num):
	system("./ParaSparse sine multiply "+str(n_sine[i]))
	system("./ParaSparse grav1d multiply "+str(n_grav1d[i]))
	system("./ParaSparse grav2d multiply "+str(n_grav2d[i]))
	
print "Simulation done, calculating."
	
N_sine, ticks_sine, time_sine = np.loadtxt("time_sine_mult.out", unpack=True)
N_grav1d, ticks_grav1d, time_grav1d = np.loadtxt("time_grav1d_mult.out", unpack=True)
N_grav2d, ticks_grav2d, time_grav2d = np.loadtxt("time_grav2d_mult.out", unpack=True)

fit_sine_time = np.polyfit(np.log(N_sine[5:]), np.log(time_sine[5:]), 1)
fit_sine_ticks = np.polyfit(np.log(N_sine[5:]), np.log(ticks_sine[5:]), 1)
fit_grav1d_time = np.polyfit(np.log(N_grav1d[5:]), np.log(time_grav1d[5:]), 1)
fit_grav1d_ticks = np.polyfit(np.log(N_grav1d[5:]), np.log(ticks_grav1d[5:]), 1)
fit_grav2d_time = np.polyfit(np.log(N_grav2d[5:]), np.log(time_grav2d[5:]), 1)
fit_grav2d_ticks = np.polyfit(np.log(N_grav2d[5:]), np.log(ticks_grav2d[5:]), 1)

p_sine_time = np.poly1d(fit_sine_time)
p_sine_ticks = np.poly1d(fit_sine_ticks)
p_grav1d_time = np.poly1d(fit_grav1d_time)
p_grav1d_ticks = np.poly1d(fit_grav1d_ticks)
p_grav2d_time = np.poly1d(fit_grav2d_time)
p_grav2d_ticks = np.poly1d(fit_grav2d_ticks)

Ysine_time = np.exp(p_sine_time(np.log(Xsine)))
Ysine_ticks = np.exp(p_sine_ticks(np.log(Xsine)))
Ygrav1d_time = np.exp(p_grav1d_time(np.log(Xgrav1d)))
Ygrav1d_ticks = np.exp(p_grav1d_ticks(np.log(Xgrav1d)))
Ygrav2d_time = np.exp(p_grav2d_time(np.log(Xgrav2d)))
Ygrav2d_ticks = np.exp(p_grav2d_ticks(np.log(Xgrav2d)))

print "Plotting!"

ax = plt.subplot(231)
plt.loglog(N_sine, time_sine, '.')
plt.loglog(Xsine, Ysine_time)
plt.ylabel("Time (s)")
plt.text(0.3, 0.9, "O = {0:g}".format(fit_sine_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(234)
plt.loglog(N_sine, ticks_sine, '.')
plt.loglog(Xsine, Ysine_ticks)
plt.ylabel("Ticks")
plt.text(0.3, 0.9, "O = {0:g}".format(fit_sine_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Sine")

ax = plt.subplot(232)
plt.loglog(N_grav1d, time_grav1d, '.')
plt.loglog(Xgrav1d, Ygrav1d_time)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav1d_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(235)
plt.loglog(N_grav1d, ticks_grav1d, '.')
plt.loglog(Xgrav1d, Ygrav1d_ticks)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav1d_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Grav 1D")

ax = plt.subplot(233)
plt.loglog(N_grav2d, time_grav2d, '.')
plt.loglog(Xgrav2d, Ygrav2d_time)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav2d_time[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
ax = plt.subplot(236)
plt.loglog(N_grav2d, ticks_grav2d, '.')
plt.loglog(Xgrav2d, Ygrav2d_ticks)
plt.text(0.3, 0.9, "O = {0:g}".format(fit_grav2d_ticks[0]), horizontalalignment='center', 
		verticalalignment='center', transform = ax.transAxes)
plt.xlabel("Grav 2D")

plt.show()