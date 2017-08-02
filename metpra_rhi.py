### ---- Read and plot rhi

import numpy as np
import matplotlib.pyplot as plt
import wradlib
import h5py
import glob




ZP = '20170615'
print 'Zeitpunkt: ',ZP

year, m, d = ZP[0:4], ZP[4:6], ZP[6:8]
ye = ZP[2:4]



if year < '2015':
    print 'archive'
    sc = 'radar-archiv'
else:
    sc = 'radar'
################################################### Read RADOLAN GK Koordinaten

iii = 0
pfad = glob.glob('/automount/'+sc+'/scans/'+str(year)+'/'+str(year)+'-'+str(m)+'/'+
        str(year)+'-'+str(m)+'-'+str(d)+'/n_rhi_1/*')

print len(pfad)

#len(pfad)
for xi in range(len(pfad)):
    pfad_rhi = pfad[xi]
    print pfad_rhi

    boxpol_filename = wradlib.util.get_wradlib_data_file(pfad_rhi)

    ppi=h5py.File(boxpol_filename,'r')
    data, attrs = wradlib.io.read_GAMIC_hdf5(boxpol_filename)

    #print data[u'SCAN0'].keys()
    #print attrs['VOL'].keys()
    #print attrs['SCAN0'].keys()
    azim =  attrs['SCAN0'][u'azimuth']
    #print azim


    zh = data['SCAN0'][u'ZH']['data']
    phidp = data['SCAN0'][u'PHIDP']['data']
    rhohv = data['SCAN0'][u'RHOHV']['data']
    zv = data['SCAN0'][u'ZV']['data']
    zdr = data['SCAN0'][u'ZDR']['data']
    kdp = data['SCAN0'][u'KDP']['data']

    r = attrs['SCAN0']['r']
    az = attrs['SCAN0']['az']
    #print az
    th = attrs['SCAN0']['el']

    lon_ppi = attrs['VOL']['Longitude']
    lat_ppi = attrs['VOL']['Latitude']
    alt_ppi = attrs['VOL']['Height']

    from mpl_toolkits.axisartist.grid_finder import FixedLocator, DictFormatter
    # reading in GAMIC hdf5 file

    mask_ind = np.where(data <= np.nanmin(zh))
    zh[mask_ind] = np.nan
    ma = np.ma.array(zh, mask=np.isnan(zh))


    fft = 15
    import pcc as pcc
    cmap2 = pcc.get_miub_cmap()

    radpar = [zh, zv, zdr,phidp, rhohv, kdp]
    vmini = [-30, 0, -1, -100, .7,-0.75]
    vmaxi = [60, 50, 3, 50, .99, 2]
    name = ['Zh (dBZ)','Zv (dBZ)','Zdr (db)','PhiDP','RhoHV (/)', r"kdp ( $^\circ$ $km^{-1}$)"]
    #cc = [cmap2, cmap2, cmap2, 'jet', 'jet', 'jet']
    iii=0
    #fig = plt.figure(figsize=(16,12))
    import matplotlib.gridspec as gridspec
    gs = gridspec.GridSpec(5, 5)

    cgax, caax, paax, pm, xxx, yyy = pcc.pcc_plot_cg_rhi(radpar[iii], r=r,
                                                         th=th,
                                                         rf=1e3,
                                                         subplot=gs[1:4,:],

                                                         autoext=True,
                                                         vmin=vmini[iii],
                                                         vmax= vmaxi[iii],
                                                         cmap=cmap2)
    cgax.set_ylim(0,12)
    cbar = plt.gcf().colorbar(pm, shrink=0.8, extend='both')#, pad=0.05)
    plt.gca().invert_xaxis()
    cbar.set_label(name[iii],fontsize=fft)
    caax.set_xlabel('x (km)',fontsize=fft)
    caax.set_ylabel('z (km)',fontsize=fft)

    cbar.ax.tick_params(labelsize=fft)

    caax.tick_params(labelsize=fft)

    tit = pfad_rhi[69:73]+'-'+ pfad_rhi[73:75]+'-'+ pfad_rhi[75:77]+' T'+pfad_rhi[77:79]+':'+pfad_rhi[79:81]+':'+pfad_rhi[81:83] + ' UTC - ZH - '+str(round(azim,1))+ r'$^\circ$ az - Bonn'

    plt.title(tit, fontsize=fft)


    gh = cgax.get_grid_helper()
    locs = [0.]
    gh.grid_finder.grid_locator1 = FixedLocator(locs)
    gh.grid_finder.tick_formatter1 = DictFormatter(dict([(i, r"${0:.0f}^\circ$".format(i)) for i in locs]))
    plt.yticks(fontsize=fft)
    plt.xticks(fontsize=fft)



    plt.tight_layout()
    plt.savefig('/automount/ftp/velibor/metprarhi/rhi_boxpol_'+ str(pfad_rhi[69:83]) + '.png' )
    plt.close()
    #plt.show()

