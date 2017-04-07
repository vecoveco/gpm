import satlib as sl
import glob
from netCDF4 import Dataset
import numpy as np

pfad_3d = sorted(glob.glob(('/automount/ags/velibor/gpmdata/nicht3dComposit/*.nc')))

comp3d = Dataset(pfad_3d[0], mode='r')

x3 = comp3d.variables['x'][:]
y3 = comp3d.variables['y'][:]
z3 = comp3d.variables['z'][:]

zh3 = comp3d['cband_oase_zh']

pfad = sorted(glob.glob(('/automount/ags/velibor/gpmdata/dpr/*.HDF5')))[69]

print pfad

glat, glon, gpp, gtime = sl.read_dpr(pfad, 'NS')
gtime =  sl.get_time_of_gpm(glon, glat, gtime)


print glat.shape, glon.shape, gpp.shape

blon, blat, bpp = sl.cut_the_swath(glon, glat, gpp, eu=0)

print blon.shape, blat.shape, bpp.shape


print gtime


x,y,rwdata, rn = sl.read_rado(gtime, r_pro='rx')

print x.shape, y.shape, rwdata.shape

gx, gy = sl.proj_gpm2radolan(blon, blat)

rrr = sl.ipol_rad2gpm(x,y,gx,gy,rwdata)
print rrr.shape

res_bin = sl.ipol_rad2gpm(x,y,gx,gy,rn)
res_bin[res_bin!=0]= 1

bpp = bpp*res_bin


res_bin[res_bin!=0]= 1 # Randkorrektur

rand_y_unten = -4658.6447242655722
rand_y_oben = -3759.6447242655722
rand_x_rechts = 375.5378330781441


rrr[np.where(gy < rand_y_unten)] = np.nan
rrr[np.where(gy > rand_y_oben)] = np.nan
rrr[np.where(gx > rand_x_rechts)] = np.nan

res_bin[np.where(gy < rand_y_unten)] = np.nan
res_bin[np.where(gy > rand_y_oben)] = np.nan
res_bin[np.where(gx > rand_x_rechts)] = np.nan
res_bin[res_bin == 0] = np.nan #check nur 1 un NaN

ggg = bpp * res_bin

## Dynamischer Threshold
#THref = np.nanmax([np.nanmin(rrr),np.nanmin(ggg)])
THref = 0.1
## Nur Niederschlagsrelevante
rrr[rrr < THref]=np.nan
ggg[ggg < THref]=np.nan



import matplotlib.pyplot as plt
from pcc import get_miub_cmap as mappo

plt.subplot(1,3,1)
plt.pcolormesh(x,y,rwdata,vmin=0, cmap=mappo())
plt.xlim(-200,200)
plt.ylim(-4200,-4000)
plt.subplot(1,3,2)
plt.pcolormesh(gx,gy,np.ma.masked_invalid(ggg), vmin=0, cmap=mappo())
plt.xlim(-200,200)
plt.ylim(-4200,-4000)
plt.subplot(1,3,3)
plt.pcolormesh(x3,y3,zh3[3,:,:], vmin=0, cmap=mappo())
plt.xlim(-200,200)
plt.ylim(-4200,-4000)
plt.tight_layout()
plt.show()



