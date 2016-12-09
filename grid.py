import matplotlib.pyplot as plt

#execfile("radoöan5.py")

plt.pcolormesh(x, y, rwdata, color='blue',vmin=0.1,vmax=10, alpha=0.3)
plt.scatter(x,y, color='blue', lw=1, alpha=0.9, label='RADOLAN Grid')
plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),color='black',vmin=0.1,vmax=10, alpha=0.3)
plt.scatter(gpm_x,gpm_y, alpha=0.9, lw=5, color='black', label='GPROF Grid')
plt.xlim(-85,-35); plt.ylim(-4140, -4080)
plt.legend(loc='upper right',scatterpoints = 1)
plt.xlabel('x [km]'); plt.ylabel('y [km]')
plt.show()


plt.pcolormesh(x, y, rwdata, color='blue',vmin=0.1,vmax=10, alpha=0.3)
plt.scatter(x,y, color='blue', lw=1, alpha=0.9, label='RADOLAN Grid')
plt.pcolormesh(gpm_x, gpm_y,np.ma.masked_invalid(gprof_pp[latstart:latend]),color='black',vmin=0.1,vmax=10, alpha=0.3)
plt.scatter(gpm_x,gpm_y, alpha=0.9, lw=5, color='black', label='GPROF Grid')
plt.xlim(-360,-300); plt.ylim(-4460, -4380)
plt.legend(loc='upper right',scatterpoints = 1)
plt.xlabel('x [km]'); plt.ylabel('y [km]')
plt.show()