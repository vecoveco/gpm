'''

By Kai in wradlib from Schwaller und Morris 2011

@ http://wradlib.org/wradlib-docs/latest/notebooks/match3d/wradlib_match_workflow.html#Conventions

"Schwaller, MR, and Morris, KR. 2011. A ground validation network for the Global Precipitation Measurement mission. J. Atmos. Oceanic Technol., 28, 301-319."

'''

import wradlib as wradlib
import matplotlib.pyplot as pl
import matplotlib as mpl
from matplotlib.collections import PatchCollection
from matplotlib.colors import from_levels_and_colors
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.cm as cm
import warnings

#warnings.filterwarnings('ignore')
#try:
#    get_ipython().magic("matplotlib inline")
#except:
#    pl.ion()
import numpy as np
import datetime as dt
from osgeo import osr
import wradlib as wrl
import datetime as dt
import numpy as np
from osgeo import gdal
# flake8: noqa


# Import from external libraries
#    these functions and objects will finally have to be moved to wradlib!
# use relative import
#try:
#    from io_func import *
#except:
#    from . io_func import *
def read_gpm(filename):

    pr_data = wrl.io.read_generic_hdf5(filename)

    lon = pr_data['NS/Longitude']['data']
    lat = pr_data['NS/Latitude']['data']

    year = pr_data['NS/ScanTime/Year']['data']
    month = pr_data['NS/ScanTime/Month']['data']
    dayofmonth = pr_data['NS/ScanTime/DayOfMonth']['data']
    dayofyear = pr_data['NS/ScanTime/DayOfYear']['data']
    hour = pr_data['NS/ScanTime/Hour']['data']
    minute = pr_data['NS/ScanTime/Minute']['data']
    second = pr_data['NS/ScanTime/Second']['data']
    secondofday = pr_data['NS/ScanTime/SecondOfDay']['data']
    millisecond = pr_data['NS/ScanTime/MilliSecond']['data']
    date_array = zip(year, month, dayofmonth,
                     hour, minute, second,
                     millisecond.astype(np.int32) * 1000)
    pr_time = np.array(
        [dt.datetime(d[0], d[1], d[2], d[3], d[4], d[5], d[6]) for d in
         date_array])

    sfc = pr_data['NS/PRE/landSurfaceType']['data']
    pflag = pr_data['NS/PRE/flagPrecip']['data']

    bbflag = pr_data['NS/CSF/flagBB']['data']
    zbb = pr_data['NS/CSF/heightBB']['data']
    print(zbb.dtype)
    bbwidth = pr_data['NS/CSF/widthBB']['data']
    qbb = pr_data['NS/CSF/qualityBB']['data']
    qtype = pr_data['NS/CSF/qualityTypePrecip']['data']
    ptype = pr_data['NS/CSF/typePrecip']['data']

    quality = pr_data['NS/scanStatus/dataQuality']['data']
    refl = pr_data['NS/SLV/zFactorCorrected']['data']

    # Check for bad data
    if max(quality) != 0:
        raise ValueError('GPM contains Bad Data')

    pflag = pflag.astype(np.int8)

    # Determine the dimensions
    ndim = refl.ndim
    if ndim != 3:
        raise ValueError('GPM Dimensions do not match! '
                         'Needed 3, given {0}'.format(ndim))

    tmp = refl.shape
    nscan = tmp[0]
    nray = tmp[1]
    nbin = tmp[2]

    # Reverse direction along the beam
    # TODO: Why is this reversed?
    refl = refl[::-1]

    # Change pflag=1 to pflag=2 to be consistent with 'Rain certain' in TRMM
    pflag[pflag == 1] = 2

    # Simplify the precipitation types
    ptype = (ptype/1e7).astype(np.int16)

    # Simplify the surface types
    imiss = (sfc == -9999)
    sfc = (sfc/1e2).astype(np.int16) + 1
    sfc[imiss] = 0

    # Set a quality indicator for the BB and precip type data
    # TODO: Why is the `quality` variable overwritten?

    quality = np.zeros((nscan, nray), dtype=np.uint8)

    i1 = ((qbb == 0) | (qbb == 1)) & (qtype == 1)
    quality[i1] = 1

    i2 = ((qbb > 1) | (qtype > 2))
    quality[i2] = 2

    gpm_data = {}
    gpm_data.update({'nscan': nscan, 'nray': nray, 'nbin': nbin,
                     'date': pr_time, 'lon': lon, 'lat': lat,
                     'pflag': pflag, 'ptype': ptype, 'zbb': zbb,
                     'bbwidth': bbwidth, 'sfc': sfc, 'quality': quality,
                     'refl': refl})

    return gpm_data
def _get_tilts2(dic):
    i = 0
    for k in dic.keys():
        if 'dataset' in k:
            i += 1
    return i/5

def read_gr2(filename, loaddata=True):

    gr_data = wrl.io.read_generic_hdf5(filename)
    dat = gr_data['what']['attrs']['date'].decode()
    tim = gr_data['what']['attrs']['time'].decode()
    date = dt.datetime.strptime(dat + tim, "%Y%d%m%H%M%S")
    source = gr_data['what']['attrs']['source']

    lon = gr_data['where']['attrs']['lon']
    lat = gr_data['where']['attrs']['lat']
    alt = gr_data['where']['attrs']['height']

    if gr_data['what']['attrs']['object'].decode() == 'PVOL':
        ntilt = int(_get_tilts2(gr_data))
        print("ntilt:", ntilt)
    else:
        raise ValueError('GR file is no PPI/Volume File')

    ngate = np.zeros(ntilt, dtype=np.int16)
    nbeam = np.zeros(ntilt)
    elang = np.zeros(ntilt)
    r0 = np.zeros(ntilt)
    dr = np.zeros(ntilt)
    a0 = np.zeros(ntilt)

    for i in np.arange(0, ntilt, dtype=np.uint8):
        a0[i] = gr_data['dataset{0}/how'.format(i+1)]['attrs']['astart']
        elang[i] = gr_data['dataset{0}/where'.format(i+1)]['attrs']['elangle']
        ngate[i] = gr_data['dataset{0}/where'.format(i+1)]['attrs']['nbins']
        r0[i] = gr_data['dataset{0}/where'.format(i+1)]['attrs']['rstart']
        dr[i] = gr_data['dataset{0}/where'.format(i+1)]['attrs']['rscale']
        nbeam[i] = gr_data['dataset{0}/where'.format(i+1)]['attrs']['nrays']

    if ((len(np.unique(r0)) != 1) |
            (len(np.unique(dr)) != 1) |
            (len(np.unique(a0)) != 1) |
            (len(np.unique(nbeam)) != 1) |
            (nbeam[0] != 360)):
        raise ValueError('GroundRadar Data layout dos not match')

    gr_dict = {}
    gr_dict.update({'source': source, 'date': date, 'lon': lon, 'lat': lat,
                    'alt': alt, 'ngate': ngate, 'nbeam': nbeam, 'ntilt': ntilt,
                    'r0': r0, 'dr': dr, 'a0': a0, 'elang': elang})
    if not loaddata:
        return gr_dict

    sdate = []
    refl = []
    for i in np.arange(0, ntilt, dtype=np.uint8):
        dat = gr_data['dataset{0}/what'.format(i+1)]['attrs']['startdate'].decode()
        tim = gr_data['dataset{0}/what'.format(i+1)]['attrs']['starttime'].decode()
        date = dt.datetime.strptime(dat + tim, "%Y%d%m%H%M%S")
        sdate.append(date)
        quantity = gr_data['dataset{0}/data1/what'.format(i+1)]['attrs']['quantity']
        factor = gr_data['dataset{0}/data1/what'.format(i+1)]['attrs']['gain']
        offset = gr_data['dataset{0}/data1/what'.format(i+1)]['attrs']['offset']
        if quantity.decode() == 'DBZH':
            dat = gr_data['dataset{0}/data1/data'.format(i+1)]['data'] * factor + offset
            refl.append(dat)

    sdate = np.array(sdate)
    refl = np.array(refl)

    gr_dict.update({'sdate': sdate, 'refl': refl})

    return gr_dict

# Space-born precipitation radar parameters
pr_pars = {"trmm": {
   "zt": 402500.,  # orbital height of TRMM (post boost)   APPROXIMATION!
   "dr": 250.,     # gate spacing of TRMM
    }, "gpm": {
   "zt": 407000.,  # orbital height of GPM                 APPROXIMATION!
   "dr": 125.      # gate spacing of GPM
}}

# Set parameters for this procedure
bw_pr = 0.71                  # PR beam width
platf = "gpm"                 # PR platform/product: one out of ["gpm", "trmm_2a23", "trmm_2a25"]
zt = pr_pars[platf]["zt"]     # PR orbit height (meters)
dr_pr = pr_pars[platf]["dr"]  # PR gate length (meters)
ee = 2                        # Index that points to the GR elevation angle to be used

# define GPM data set
gpm_file = wradlib.util.get_wradlib_data_file('/automount/ags/velibor/gpmdata/dpr/2A.GPM.DPR.V6-20160118.20141007-S015721-E032951.003445.V04A.HDF5')

# define matching ground radar file
#gr2gpm_file = wradlib.util.get_wradlib_data_file('hdf5/IDR66_20141206_094829.vol.h5')
gr2gpm_file = wradlib.util.get_wradlib_data_file('/automount/radar-archiv/scans/2014/2014-10/2014-10-07/ppi_1p5deg/2014-10-07--02:30:00,00.mvol')

# define TRMM data sets
#trmm_2a23_file = wradlib.util.get_wradlib_data_file('trmm/2A-RW-BRS.TRMM.PR.2A23.20100206-S111422-E111519.069662.7.HDF')
#trmm_2a25_file = wradlib.util.get_wradlib_data_file('trmm/2A-RW-BRS.TRMM.PR.2A25.20100206-S111422-E111519.069662.7.HDF')

# define matching ground radar file
#gr2trmm_file = wradlib.util.get_wradlib_data_file('hdf5/IDR66_20100206_111233.vol.h5')


##Read and organize the data
# read spaceborn PR data
if platf == "gpm":
    pr_data = read_gpm(gpm_file)
#elif platf == "trmm":
#    try:
#        pr_data = read_trmm(trmm_2a23_file, trmm_2a25_file)
#        print("test")
#    except IOError:
#        # Mostly needed on Windows Anaconda (netcdf4 without hdf4 support)
#        pr_data = read_trmm_gdal(trmm_2a23_file, trmm_2a25_file)
else:
    raise("Invalid platform")
# read matching GR data
if platf == "gpm":
    gr_data = read_gr2(gr2gpm_file)
#elif platf=="trmm":
#    gr_data = read_gr2(gr2trmm_file)
else:
    raise("Invalid platform")

##Extract relevant GR data and meta-data
# number of rays in gr sweep
nray_gr = gr_data['nbeam'].astype("i4")[ee]
# number of gates in gr beam
ngate_gr = gr_data['ngate'].astype("i4")[ee]
# number of sweeps
nelev = gr_data['ntilt']
# elevation of sweep (degree)
elev = gr_data['elang'][ee]
# gate length (meters)
dr_gr = gr_data['dr'][ee]
# reflectivity array of sweep
ref_gr = gr_data['refl'][ee]
# sweep datetime stamp
date_gr = gr_data['sdate'][ee]
# range of first gate
r0_gr = gr_data['r0'][ee]
# azimuth angle of first beam
a0_gr = gr_data['a0'][ee]
# Longitude of GR
lon0_gr = gr_data['lon']
# Latitude of GR
lat0_gr = gr_data['lat']
# Altitude of GR (meters)
alt0_gr = gr_data['alt']
# Beam width of GR (degree)
bw_gr = 1.

##Extract relevant PR data and meta-data

# Longitudes of PR scans
pr_lon = pr_data['lon']
# Latitudes of PR scans
pr_lat = pr_data['lat']
# Precip flag
pflag = pr_data['pflag']
# Number of scans on PR data
nscan_pr= pr_data['nscan']
# Number of rays in one PR scan
nray_pr = pr_data['nray']
# Number of gates in one PR ray
ngate_pr = pr_data['nbin']

##Set fundamental georeferencing parameters

# Calculate equivalent earth radius
wgs84 = wradlib.georef.get_default_projection()
re1 = wradlib.georef.get_earth_radius(lat0_gr, wgs84) * 4./3.
print("eff. Earth radius 1:", re1)
a = wgs84.GetSemiMajor()
b = wgs84.GetSemiMinor()
print("SemiMajor, SemiMinor:", a, b)

# Set up aeqd-projection gr-centered
rad = wradlib.georef.proj4_to_osr(('+proj=aeqd +lon_0={lon:f} ' +
                                   '+lat_0={lat:f} +a={a:f} ' +
                                   '+b={b:f}').format(lon=lon0_gr,
                                                      lat=lat0_gr,
                                                      a=a, b=b))
re2 = wradlib.georef.get_earth_radius(lat0_gr, rad) * 4./3.
print("eff. Earth radius 2:", re2)


##Georeference GR data

# TODO: Seperate the insides of wradlib.georef.polar2lonlatalt_n

# create gr range and azimuth arrays
rmax_gr = r0_gr + ngate_gr * dr_gr
r_gr = np.arange(0, ngate_gr) * dr_gr + dr_gr/2.
az_gr = np.arange(0, nray_gr) - a0_gr
print("Range/Azi-Shape:", r_gr.shape, az_gr.shape)

# create gr lonlat grid
gr_polargrid = np.meshgrid(r_gr, az_gr)
gr_lon, gr_lat, gr_alt = wradlib.georef.polar2lonlatalt_n(gr_polargrid[0], gr_polargrid[1], elev, (lon0_gr, lat0_gr, alt0_gr ))
gr_ll = np.dstack((gr_lon, gr_lat, gr_alt))
print("LonLatAlt-Grid-Shape", gr_ll.shape)

# reproject to xyz
gr_xyz = wradlib.georef.reproject(gr_ll, projection_source=wgs84, projection_target=rad)
print("XYZ-Grid-Shape:", gr_xyz.shape)

# get radar domain (outer ring)
gr_domain = gr_xyz[:,-1,0:2]
gr_domain = np.vstack((gr_domain, gr_domain[0]))
print("Domain-Shape:", gr_domain.shape)

##Georeference PR data


pr_x, pr_y = wradlib.georef.reproject(pr_lon, pr_lat,
                                      projection_source=wgs84,
                                      projection_target=rad)
pr_xy = np.dstack((pr_x, pr_y))
print("PR-GRID-Shapes:", pr_x.shape, pr_y.shape, pr_xy.shape)


##Subset relevant PR data

# Create ZonalData for spatial subsetting (inside GR range domain)

# get precip indexes
precip_mask = (pflag == 2) & wradlib.zonalstats.get_clip_mask(pr_xy, gr_domain, rad)

# get iscan/iray boolean arrays
#iscan = precip_idx.nonzero()[0]
#iray = precip_idx.nonzero()[1]
#print(iscan, iray)
print(precip_mask.shape, pr_xy.shape)
print(pflag.shape, pr_xy.reshape(-1, pr_xy.shape[-1]).shape)
pl.imshow(precip_mask)

##PR Parallax Correction
print("NRAY", nray_pr)
print("NBIN", ngate_pr)

# Approximation!
alpha = abs(-17.04 + np.arange(nray_pr) * bw_pr)

print(alpha.shape)

# Correct for parallax, get 3D-XYZ-Array
#   xyzp_pr: Parallax corrected xyz coordinates
#   r_pr_inv: range array from ground to PR platform
#   zp: PR bin altitudes
xyp_pr, r_pr_inv, z_pr = wradlib.georef.correct_parallax(pr_xy, ngate_pr, dr_pr, alpha)
xyzp_pr = np.concatenate((xyp_pr, np.repeat(z_pr[np.newaxis, ..., np.newaxis], xyp_pr.shape[0], axis=0)),
                   axis=-1)
print(pr_xy.shape)
print("PR_XYP:", xyp_pr.shape, xyzp_pr.shape, r_pr_inv.shape, z_pr.shape)

#Compute spherical coordinates of PR bins with regard to GR
# TODO: Do we have to consider refraction in sat2pol?

r_pr, elev_pr, az_pr = wradlib.georef.sat2pol(xyzp_pr, alt0_gr, re1)
mask = (elev_pr > (1.0 - bw_gr/2.)) & (elev_pr < (1.0 + bw_gr/2.))
##pl.figure()
##pl.pcolormesh(mask[90,:,:].T)
print(r_pr.shape, elev_pr.shape, az_pr.shape)
print(r_pr, az_pr)


##Compute PR and GR pulse volumes

# PR pulse volumes

# Range of PR bins
dists = wradlib.georef.dist_from_orbit(zt, alpha, r_pr_inv)
print(dists.shape, dists)
## Original IDL code...
##    rt=zt/COS(!dtor*alpha)-range
##    volp=(1.e-9)*!pi*(rt*!dtor*bwt/2.)^2*drt
## Translated to Python
vol_pr2  = np.pi * dr_pr * (dists * np.radians(bw_pr / 2.))**2
##fig = pl.figure(figsize=(12,4))
##pm = pl.pcolor(vol_pr.T)
##pl.colorbar(pm)

# Or using wradlib's native function
vol_pr = wradlib.qual.pulse_volume(dists, dr_pr, bw_pr)
#vol_pr = np.pi * dr_pr * (dists ** 2) * (np.tan(np.radians(bw_pr/2.))) ** 2

# Evaluate difference between both approaches
print("Min. difference (m3):", (vol_pr - vol_pr2).min())
print("Max. difference (m3): ", (vol_pr - vol_pr2).max())
print("Average rel. difference (%):", round(np.mean(vol_pr-vol_pr2)*100./np.mean(np.mean(vol_pr2)), 4))

# Verdict: differences are negligble - use wradlibs's native function!

# GR pulse volumes
#   along one beam
vol_gr = wradlib.qual.pulse_volume(r_gr, dr_gr, bw_gr)
#   with shape (nray_gr, ngate_gr)
vol_gr = np.repeat(vol_gr, nray_gr).reshape((nray_gr,ngate_gr), order="F")

##Median Brightband Width/Height


print(z_pr.shape)
print(pr_data['zbb'].shape, pr_data['bbwidth'].shape, pr_data['quality'].shape)
#qual =  pr_data['quality'].copy()
#qual[...] = 0.
ratio, ibb = wradlib.qual.get_bb_ratio(pr_data['zbb'], pr_data['bbwidth'], pr_data['quality'], z_pr)
zbb = pr_data['zbb'].copy()
zbb[~ibb] = np.nan
##pl.pcolormesh(ratio[60,:,:].T, vmin=-1, vmax=2)
##pl.colorbar()


##Convert PR Ku reflectivities to S-band

# REVERSE!!!
refp = pr_data['refl'][:,:,::-1]
print(refp.max(), refp.min(), np.nanmedian(refp))
print("REFP:", refp.shape)

refp_ss = np.zeros_like(refp) * np.nan
refp_sh = np.zeros_like(refp) * np.nan

a_s, a_h = (wradlib.trafo.ku2s.snow, wradlib.trafo.ku2s.hail)

ia = (ratio >= 1)
print("SHAPE:", refp[ia].shape, a_s[:,10].shape, wradlib.util.calculate_polynomial(refp[ia], a_s[:,10]).shape)
refp_ss[ia] = refp[ia] + wradlib.util.calculate_polynomial(refp[ia], a_s[:,10])
refp_sh[ia] = refp[ia] + wradlib.util.calculate_polynomial(refp[ia], a_h[:,10])
ib = (ratio <= 0)
refp_ss[ib] = refp[ib] + wradlib.util.calculate_polynomial(refp[ib], a_s[:,0])
refp_sh[ib] = refp[ib] + wradlib.util.calculate_polynomial(refp[ib], a_h[:,0])
im = (ratio > 0) & (ratio < 1)
ind = np.round(ratio[im] * 10).astype(np.int)
#print("W:", a_s[:,ind].shape)
refp_ss[im] = refp[im] + wradlib.util.calculate_polynomial(refp[im], a_s[:,ind])
refp_sh[im] = refp[im] + wradlib.util.calculate_polynomial(refp[im], a_h[:,ind])

# Jackson Tan's fix for C-band
is_cband = False
if (is_cband):
    deltas = (refp_ss - refp) * 5.3 / 10.0
    refp_ss = refp + deltas
    deltah = (refp_sh - refp) * 5.3 / 10.0
    refp_sh = refp + deltah

refp_ss[refp < 0] = np.nan
out = np.ma.masked_invalid(refp_ss)
##pl.figure()
##pl.pcolormesh(out[60,:,:].T, vmin=0, vmax=60)
##pl.colorbar()
##pl.figure()
##pl.pcolormesh(refp[60,:,:].T, vmin=0, vmax=60)
##pl.colorbar()
##pl.figure()
##pl.pcolormesh(ratio[60,:,:].T, vmin=-1, vmax=2)
##pl.colorbar()



##Convert S-band GR reflectivities to Ku-band

# Convert S-band GR reflectivities to Ku-band using method of Liao and Meneghini (2009)
ref_gr_ku = np.zeros_like(ref_gr) * np.nan

# Which zbb value should we take here???
#    Q'n'Dirty: just take the mean of all PR profiles
#    TODO: Consider zbb for each profile during the matching process

# Snow
ia = ( gr_xyz[...,2] >= np.nanmean(zbb) )
#ref_gr_ku[ia] = wradlib.trafo.ku2s.snow[0] + wradlib.trafo.ku2s.snow[1]*ref_gr[ia] + wradlib.trafo.ku2s.snow[2]*ref_gr[ia]**2
ref_gr_ku[ia] = wradlib.util.calculate_polynomial(ref_gr[ia], wradlib.trafo.s2ku.snow)
# Rain
ib = ( gr_xyz[...,2] < np.nanmean(zbb) )
#ref_gr_ku[ib] = wradlib.trafo.ku2s.rain[0] + wradlib.trafo.ku2s.rain[1]*ref_gr[ia] + wradlib.trafo.ku2s.rain[2]*ref_gr[ia]**2
ref_gr_ku[ib] = wradlib.util.calculate_polynomial(ref_gr[ib], wradlib.trafo.s2ku.rain)

# Jackson Tan's fix for C-band
is_cband = False
if (is_cband):
    delta = (ref_gr_ku - ref_gr) * 5.3/10.0
    ref_gr_ku = ref_gr + delta

## Identify which PR rays actually intersect with the GR sweep


# First assumption: no valid PR bins (all False)
valid = np.asarray(elev_pr, dtype=np.bool)==False
print(valid.shape, precip_mask.shape)
# PR is inside GR range and is precipitating
iscan = precip_mask.nonzero()[0]
iray = precip_mask.nonzero()[1]
valid[iscan,iray] = True
# PR bins intersect with GR sweep
valid = valid & (elev_pr >= elev-bw_gr/2.) & (elev_pr <= elev+bw_gr/2.)
#valid = precip_mask[..., np.newaxis] & (elev_pr >= elev-bw_gr/2.) & (elev_pr <= elev+bw_gr/2.)
# Number of matching PR bins per profile
nvalids = np.sum(valid, axis=2)
# scan and ray indices for profiles with at least one valid bin
vscan, vray = np.where(nvalids>0)
# number of profiles with at least one valid bin
nprof = len(vscan)

# Lots of containers to store samples (only for one GR sweep angle!)
x = np.zeros(nprof)*np.nan        # x coordinate of sample
y = np.zeros(nprof)*np.nan        # y coordinate of sample
z = np.zeros(nprof)*np.nan        # z coordinate of sample
dz = np.zeros(nprof)*np.nan       # depth of sample
ds = np.zeros(nprof)*np.nan       # width of sample
rs = np.zeros(nprof)*np.nan       # range of sample from GR
refpr1 = np.zeros(nprof)*np.nan     # PR reflectivity
refpr2 = np.zeros(nprof)*np.nan     # PR reflectivity (S-band, snow)
refpr3 = np.zeros(nprof)*np.nan     # PR reflectivity (S-band, hail)
refgr1 = np.zeros(nprof)*np.nan     # GR reflectivity
refgr2 = np.zeros(nprof)*np.nan     # GR reflectivity (Ku-band)
ntotpr = np.zeros(nprof,dtype="i4")# total number of PR bins in sample
nrej1 = np.zeros(nprof,dtype="i4")# number of rejected PR bins in sample
ntotgr = np.zeros(nprof,dtype="i4")# total number of GR bins in sample
nrej2 = np.zeros(nprof,dtype="i4")# number of rejected GR bins in sample
iref1 = np.zeros(nprof)*np.nan    # path-integrated PR reflectivity
iref2 = np.zeros(nprof)*np.nan    # path-integrated GR reflectivity
stdv1 = np.zeros(nprof)*np.nan    # std. dev. of PR reflectivity in sample
stdv2 = np.zeros(nprof)*np.nan    # std. dev. of GR reflectivity in sample
volpr = np.zeros(nprof)*np.nan     # total volume of PR bins in sample
volgr = np.zeros(nprof)*np.nan     # total volume of GR bins in sample

# Loop over relevant PR profiles
for ii, (ss, rr)  in enumerate(zip(vscan,vray)):
    # Index and count valid bins in each profile
    ip = np.where(valid[ss,rr])[0]
    numbins = len(ip)
    ntotpr[ii]=numbins
    if numbins == 0:
        continue
    # Compute the mean position of these bins
    x[ii]=np.mean(xyzp_pr[ss,rr,ip,0])
    y[ii]=np.mean(xyzp_pr[ss,rr,ip,1])
    z[ii]=np.mean(xyzp_pr[ss,rr,ip,2])

    #np.stack((pr_xp, pr_yp,
    #          np.repeat(zp[np.newaxis, ...], pr_xp.shape[0], axis=0)),
    #          axis=3)

    # Thickness of the layer
    dz[ii]=(numbins * dr_pr) * np.cos( np.radians(alpha[rr]) )

    # PR averaging volume
    volpr[ii]=np.sum(vol_pr2[rr,ip])

    # Note mean TRMM beam diameter
    ds[ii]=np.radians(bw_pr) * np.mean( ( (zt-z[ii]) / np.cos( np.radians(alpha[rr]) ) ) )

    # Note distance from radar
    s=np.sqrt(x[ii]**2+y[ii]**2)
    rs[ii]=(re2+z[ii]) * np.sin(s/re2) / np.cos(np.radians(elev))

    # This should not be required because we applied ZonalData
    ### Check that sample is within radar range
    ##if r[ii,jj]+ds[ii,jj]/2. gt rmax then continue

    ## THIS IS THE ORIGINAL IDL CODE - IS THIS A BUG???
    ##ref1[ii,jj]=MEAN(refp1,/nan)
    ##ref3[ii,jj]=MEAN(refp2,/nan)
    ##ref4[ii,jj]=MEAN(refp3,/nan)

    # Simple linear average of reflectivity
    #   - we can become fancier in the next step
    # ATTENTION: NEED TO FLIP ARRAY
    refpr1[ii]=np.nanmean(np.flipud(refp)   [ss,rr,ip])
    refpr2[ii]=np.nanmean(np.flipud(refp_ss)[ss,rr,ip])
    refpr3[ii]=np.nanmean(np.flipud(refp_sh)[ss,rr,ip])

    ## Not sure why we need this...
    ### Note the number of rejected bins
    ##nrej1[ii,jj]=ROUND(TOTAL(FINITE(refp1,/nan)))
    ##if FINITE(stdv1[ii,jj]) eq 0 and np-nrej1[ii,jj] gt 1 then STOP

    # SHOULD WE USE ZONALDATA INSTEAD? COULD BE MORE ACCURATE, BUT ALSO SLOWER
    # WE COULD BASICALLY START A NEW LOOP HERE AND RUN ZONALDATA BEFORE

    # Compute the horizontal distance to all the GR bins
    d = np.sqrt((gr_xyz[...,0]-x[ii])**2 + (gr_xyz[...,1]-y[ii])**2)

    # Find all GR bins within the SR beam
    aa, bb = np.where(d <= ds[ii]/2.)

    # Store the number of bins
    ntotgr[ii] = len(aa)

    if len(aa) == 0:
        continue

    # Extract the relevant GR bins

    # Compute the GR averaging volume
    volgr[ii]=np.sum(vol_gr[aa,bb])

    # Average over those bins that exceed the reflectivity threshold
    #   IDL code does exponential distance and volume weighting
    #   Let's try simple mean first,
    #   THEN ZonalStats!
    refgr1[ii]=np.nanmean(ref_gr[aa,bb])
    refgr2[ii]=np.nanmean(ref_gr_ku[aa,bb])


fig = pl.figure(figsize=(12,5))
ax = fig.add_subplot(121, aspect="equal")
pl.scatter(refgr1, refpr1, marker="+", c="black")
pl.plot([0,60],[0,60], linestyle="solid", color="black")
pl.xlim(10,50)
pl.ylim(10,50)
pl.xlabel("GR reflectivity (dBZ)")
pl.ylabel("PR reflectivity (dBZ)")
ax = fig.add_subplot(122)
pl.hist(refgr1[refpr1>-10], bins=np.arange(-10,50,5), edgecolor="None", label="GR")
pl.hist(refpr1[refpr1>-10], bins=np.arange(-10,50,5), edgecolor="red", facecolor="None", label="PR")
pl.xlabel("Reflectivity (dBZ)")
pl.legend()
pl.show()

fig = pl.figure(figsize=(12,12))
ax = fig.add_subplot(121, aspect="equal")
pl.scatter(x, y, c=refpr1, cmap=pl.cm.jet, vmin=0, vmax=50, edgecolor="None")
pl.title("PR reflectivity")
pl.grid()
ax = fig.add_subplot(122, aspect="equal")
pl.scatter(x, y, c=refgr1, cmap=pl.cm.jet, vmin=0, vmax=50, edgecolor="None")
pl.title("GR reflectivity")
pl.grid()
pl.show()