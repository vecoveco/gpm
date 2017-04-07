import satlib as sl
import glob

pfad = sorted(glob.glob(('/automount/ags/velibor/gpmdata/dpr/*.HDF5')))[34]

print pfad

glat, glon, gpp, gtime = sl.read_dpr(pfad, 'NS')

print glat.shape, glon.shape, gpp.shape

blon, blat, bpp = sl.cut_the_swath(glon, glat, gpp, eu=0)

print blon.shape, blat.shape, bpp.shape

gtime =  sl.get_time_of_gpm(blon, blat, gtime)

print gtime

ht, mt = gtime[14:16], str(int(round(float(gtime[17:19])/5.0)*5.0))
print ht, mt
if mt == '0':
    mt = '00'
elif mt == '5':
    mt = '05'
elif mt =='60':
    mt = '55'

x,y,rwdata, rn = sl.read_rado(gtime, r_pro='rx')

