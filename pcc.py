def dscat(xdat, ydat):
    import matplotlib.pyplot as plt, numpy as np, numpy.random, scipy
    #data definition
    N = 1000
    #xdat, ydat = np.random.normal(size=N), np.random.normal(size=N)
    #xdat, ydat = np.sin(np.arange(1, N)),np.cos(np.arange(1, N)+10)
    #xdat, ydat = np.array([1,2,3,4,4,4,4,4,4,4]), np.array([1,2,3,4,4,4,4,4,4,4])


    #histogram definition
    xyrange = [[np.nanmin(xdat),np.nanmax(xdat)],[np.nanmin(ydat),np.nanmax(ydat)]] # data range
    bins = [100,100] # number of bins
    thresh = 0.1  #density threshold

    #linear regression
    from scipy import stats
    mask = ~np.isnan(xdat) & ~np.isnan(xdat)
    slope, intercept, r_value, p_value, std_err = stats.linregress(xdat[mask], ydat[mask])
    line = slope * xdat + intercept

    from scipy import stats
    #xdat = xdat.reshape(xdat.shape[0]*xdat.shape[1],)
    #ydat = ydat.reshape(ydat.shape[0]*ydat.shape[1],)

    # histogram the data
    hh, locx, locy = scipy.histogram2d(xdat, ydat, range=xyrange, bins=bins)
    posx = np.digitize(xdat, locx)
    posy = np.digitize(ydat, locy)

    #select points within the histogram
    ind = (posx > 0) & (posx <= bins[0]) & (posy > 0) & (posy <= bins[1])
    hhsub = hh[posx[ind] - 1, posy[ind] - 1] # values of the histogram where the points are
    xdat1 = xdat[ind][hhsub < thresh] # low density points
    ydat1 = ydat[ind][hhsub < thresh]
    hh[hh < thresh] = np.nan # fill the areas with low density by NaNs

    plt.imshow(np.flipud(hh.T),cmap='jet',extent=np.array(xyrange).flatten(), interpolation='none', origin='upper')
    plt.colorbar()
    plt.plot(xdat1, ydat1, '.',color='darkblue')
    plt.plot(xdat, line, 'r', label = 'Korrelation \n' + str(round(r_value, 3))
                                           + r'$\pm$' + str(round(std_err, 3)))
    print xdat, xdat1
    plt.legend(loc='upper left', ncol=1, fancybox=True, shadow=True,
               fontsize='20')
    plt.xlim(xyrange[0])
    plt.ylim(xyrange[1])
    plt.grid()




def rmse(predictions, targets):

    import numpy as np
    # RMS sqrt(1/n SUM/d_i - pi)^2
    return np.sqrt(((predictions - targets) ** 2).mean())


def histo(data1, data2, bino):

    import matplotlib.pyplot as plt
    #Todo: in arbeit
    plt.hist(data1, bins=bino)
    plt.hist(data2, bins=bino)
    plt.show()


#Heidke Skill Score
# Lit: J. Tan A Novel Approach to Indetify Source of Errors in IMERG for GPM Ground Validation
#                   |  Estimate  |  Reference
#-----------------------------------------------
#Hit H :            |  Yes      |    Yes
#Miss M :           |  No       |    Yes
#F False:           |  Yes      |    No
#C Correctnegative :|  No       |    No


def skill_score(estimate, reference, th=None):
    #contigency tables Lit: Tang et al. 2015
    import numpy as np
    #reshapen von arrays
    try:
        estimate = estimate.reshape(estimate.shape[0]*estimate.shape[1])
        reference = reference.reshape(reference.shape[0]*reference.shape[1])
    except IndexError:
        estimate = estimate
        reference = reference

    if th == None:
        th = 0.1 # GMI 0.1, ka 0.2 ku 0.5 Hou et al 2014
    reference, estimate = np.nan_to_num(reference), np.nan_to_num(estimate)

    #Hit
    H_pos = np.array(np.flatnonzero((estimate > th) & (reference > th)))
    #Miss
    M_pos = np.array(np.flatnonzero((estimate <= th) & (reference > th)))
    #FalseAlarm
    F_pos = np.array(np.flatnonzero((estimate > th) & (reference <= th)))
    #Correctnegative
    C_pos = np.array(np.flatnonzero((estimate <= th) & (reference <= th)))
    #Samplesize
    N = len(estimate)

    H, M, F, C = len(H_pos), len(M_pos), len(F_pos), len(C_pos)

    H, M, F, C, N = float(H), float(M), float(F), float(C), float(N)

    E = 1/N * (((H + M)*(H + F)) +
                       ((C + M) * (C + F)))

    POD = H/(H + M)
    FAR = F/(H + F)
    BID = (H +F)/(H + M)
    HSS = (H + C - E) / (N - E)

    result = {'H': H,
              'M': M,
              'F': F,
              'C': C,
              'H_pos': H_pos,
              'M_pos': M_pos,
              'F_pos': F_pos,
              'C_pos': C_pos,
              'N': N,
              'POD': POD,
              'FAR': FAR,
              'BID': BID,
              'HSS': HSS}

    #for key, value in result.items():
    #    print(key + ':', value)

    return result


def scores(H_tab, M_tab, F_tab, C_tab, N_tab):
    #score Lit: Tan et al 2016

    E_tab = 1/N_tab * (((H_tab + M_tab)(H_tab + F_tab)) +
                       ((C_tab + M_tab) * (C_tab + F_tab)))

    POD = H_tab/(H_tab + M_tab)
    FAR = F_tab/(H_tab + F_tab)
    BID = (H_tab +F_tab)/(H_tab + M_tab)
    HSS = (H_tab + C_tab - E_tab) / (N_tab - E_tab)

    result = {'POD': POD,
              'FAR': FAR,
              'BID': BID,
              'HSS': HSS}

    return result


def plot_score(estimate, reference, scoreval):
    import numpy as np
    #reshapen von arrays
    try:
        estimate = estimate.reshape(estimate.shape[0]*estimate.shape[1])
        reference = reference.reshape(reference.shape[0]*reference.shape[1])
    except IndexError:
        estimate = estimate
        reference = reference

    import matplotlib.pyplot as plt

    plt.scatter(estimate[scoreval['H_pos']],reference[scoreval['H_pos']],
                color='black', label='Hit')
    plt.scatter(estimate[scoreval['M_pos']],reference[scoreval['M_pos']],
                color='blue', label='Miss')
    plt.scatter(estimate[scoreval['F_pos']],reference[scoreval['F_pos']],
                color='red', label='False')
    plt.scatter(estimate[scoreval['C_pos']],reference[scoreval['C_pos']],
                color='green', label='Correct Negative')
    plt.grid()
    plt.legend(loc='upper right')
    plt.figtext(0.005,0.005,' POD: '+ str(scoreval['POD'])+ ' , ' + ' FAR: ' +
              str(scoreval['FAR']) + ' BID: '+ str(scoreval['BID'])+
              ' , ' + ' HSS: ' + str(scoreval['HSS'])+
                '-- H: '+ str(scoreval['H'])+ ' , ' + ' M: ' +
              str(scoreval['M']) + ' F: '+ str(scoreval['F'])+
              ' , ' + ' C: ' + str(scoreval['C']) +' N: ' + str(scoreval['N'])
                ,fontsize=10)

    #plt.show()




#Todo: Zeitschleife so bauen das ein Zeitvector ztime raus kommt!
def zeitschleife(sY,sm,sd,sH,sM,sS, eY,em,ed,eH,eM,eS):
    import datetime as dt
    import numpy as np
    ztime = []
    def datespan(startDate, endDate, delta=dt.timedelta(days=1)):
        currentDate = startDate
        while currentDate < endDate:
            yield currentDate
            currentDate += delta

    start = dt.datetime(sY, sm, sd, sH, sM,sS)
    ende = dt.datetime(eY, em, ed, eH, eM,eS)
    ztime=[]

    for timestamp in datespan(start,ende,delta=dt.timedelta(minutes=5)):
        ztime.append(timestamp.strftime('%Y%m%d%H%M%S'))


        #print timestamp.strftime('%Y%m%d%H%M%S')
        #return timestamp.strftime('%Y%m%d%H%M%S')
    return ztime

#### Idee
#Todo: Korrelation von Radar und Satellit Daten in bestimmten Bereichen
# Z.B. Nur bestimmte Niederschlagsintensiteten!
#RR1_3 = A[mask]
#GR1_3 = B[mask]
#
#Rmin, Rmax = 2, 6#
#
#RR1_3[RR1_3 < Rmin] = np.nan
#GR1_3[GR1_3< Rmax] = np.nan
#
#RR1_3[RR1_3 > Rmin] = np.nan
#GR1_3[GR1_3> Rmax] = np.nan
#
#plt.scatter(RR1_3,GR1_3)
#from pcc import rmse
#plt.title(str(rmse(RR1_3,GR1_3)))
#plt.show()

#ag = gpm_gpm.reshape(gpm_gpm.shape[0]*gpm_gpm.shape[1])
#ar = rrr.reshape(rrr.shape[0]*rrr.shape[1])
def plot_ocean(ax):
    import wradlib
    from osgeo import osr
    import os
    import wradlib as wrl
    import numpy as np

    filename = os.path.join('/automount/db01/python/data/NED/10m/physical/10m_physical/ne_10m_ocean.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    inLayer.SetSpatialFilterRect(1, 45, 19, 56.5)
    borders, keys = wradlib.georef.get_shape_coordinates(inLayer)
    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)
    print borders.shape

    for j in range(borders.shape[0]):
        bu = np.array(borders[j].shape)

        a = np.array(bu.shape)

        if a==1:
            for i in range(0,borders[j].shape[0],1):
                bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
        if a==2:    #richtig
            bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()
            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

        bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

        wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
    ax.autoscale()

def plot_borders(ax):

    from osgeo import osr
    import wradlib as wrl
    import wradlib
    import numpy as np
    wgs84 = osr.SpatialReference()
    wgs84.ImportFromEPSG(4326)
    india = osr.SpatialReference()
    # asia south albers equal area conic
    india.ImportFromEPSG(102028)

    proj_gk = osr.SpatialReference()
    proj_gk.ImportFromEPSG(31466)
    proj_ll = osr.SpatialReference()
    proj_ll.ImportFromEPSG(4326)
    gk3 = wradlib.georef.epsg_to_osr(31467)
    proj_stereo = wrl.georef.create_osr("dwd-radolan")
    proj_wgs = osr.SpatialReference()
    proj_wgs.ImportFromEPSG(4326)

    # country list
    countries = ['Germany']#,'France','Denmark', 'Netherlands', 'Poland']
    # open the input data source and get the layer
    filename = wradlib.util.get_wradlib_data_file('/automount/db01/python/data/NED/10m/cultural/10m_cultural/10m_cultural/ne_10m_admin_0_countries.shp')
    dataset, inLayer = wradlib.io.open_shape(filename)
    # iterate over countries, filter accordingly, get coordinates and plot
    for item in countries:
        #print item
        # SQL-like selection syntax
        fattr = "(name='"+item+"')"
        inLayer.SetAttributeFilter(fattr)
        # get borders and names
        borders, keys = wradlib.georef.get_shape_coordinates(inLayer, key='name')

        for j in range(borders.shape[0]):
            bu = np.array(borders[j].shape)
            a = np.array(bu.shape)

            if a==1:
                for i in range(0,borders[j].shape[0],1):
                    bordx, bordy = wrl.georef.reproject(borders[j][i][:,0], borders[j][i][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                    bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                    wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
            if a==2:    #richtig
                bordx, bordy = wrl.georef.reproject(borders[j][:,0], borders[j][:,1], projection_source=proj_wgs, projection_target=proj_stereo)
                bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

                wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)

            bord_xy = np.vstack((bordx.ravel(), bordy.ravel())).transpose()

            wradlib.vis.add_lines(ax, bord_xy, color='black', lw=2, zorder=3)
    ax.autoscale()


def boxpol_pos():
    pos_boxpol = {'lat_ppi' : 50.730519999999999, 'lon_ppi' : 7.071663
    ,'gky_ppi' : -4235.233235191105, 'gkx_ppi' : -216.64772430049572}
    return pos_boxpol