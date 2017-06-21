"""

Reading MRR Data

"""
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob


yy = '2015'
mm = '11'
dd = '08'
hh = '02'
min = '40'

pfad_mrr = glob.glob('/automount/mrr/mrr2/'+yy+'/'+yy+'-'+mm+'/'+yy+'-'+mm+'-'+dd+'/AveData_mrr2_'+yy+mm+dd+hh+min+'*.ave.gz')[0]



df = pd.read_csv(pfad_mrr, compression='gzip', header=1, delim_whitespace=True,index_col=False)
df = df.set_index(u'H')


h = np.arange(150,4800,150)


plt.plot(df.loc['Z'].values,h, label='Ref. in dBZ', color='blue', linestyle='-', lw=2)
plt.plot(df.loc['PIA'].values,h,label='PIA in dB', color='blue', linestyle='-.', lw=2)
plt.plot(df.loc['z'].values,h,label='att. Ref in dBZ', color='blue', linestyle='--', lw=2)
plt.plot(df.loc['TF'].values,h,label='TF', color='grey')
plt.plot(df.loc['RR'].values,h,label='RR', linestyle='-', color='black',lw=2)
plt.plot(df.loc['LWC'].values,h,label='Liquid Water Content')
plt.plot(df.loc['W'].values,h,label='Fallgeschwindigkeit')
plt.grid()
plt.legend(loc='lower right')
plt.xlabel('Reflectivity in dBZ')
plt.ylabel('Hight in m')
plt.title('MRR - ' + pfad_mrr[44:44+28])
plt.ylim(0,6000)
plt.xlim(0,50)
plt.show()


#ToDo : Zeitlicher VErlauf vom MRR einfugen

"""
def plot_mrr1(t_start='2016-6-12 20:55', time_delta=20):

    import numpy as np
    import pymysql as mysql
    import datetime as dt
    from datetime import datetime
    import matplotlib
    import matplotlib.pyplot as pl

    #datetime_object = datetime.strptime('Jun 1 2005  1:33PM', '%b %d %Y %I:%M%p')
    t_start = datetime.strptime(t_start,'%Y-%m-%d %H:%M')
    t_stop  = t_start + dt.timedelta(minutes=time_delta)#time_delta

    cnx = mysql.connect(user='wetterdaten', passwd='miub321', host='mysql', db='MIUB_messdaten')
    cursor = cnx.cursor()
    # Humicap 1:
    query1 = ('SELECT Base.time, '
              'h.h_1, h.h_2, h.h_3, h.h_4, h.h_5, h.h_6, h.h_7, h.h_8, h.h_9, h.h_10, '
              'h.h_11, h.h_12, h.h_13, h.h_14, h.h_15, h.h_16, h.h_17, h.h_18, h.h_19, h.h_20, '
              'h.h_21, h.h_22, h.h_23, h.h_24, h.h_25, h.h_26, h.h_27, h.h_28, h.h_29, h.h_30, h.h_31, '
              'z.h_1, z.h_2, z.h_3, z.h_4, z.h_5, z.h_6, z.h_7, z.h_8, z.h_9, z.h_10, '
              'z.h_11, z.h_12, z.h_13, z.h_14, z.h_15, z.h_16, z.h_17, z.h_18, z.h_19, z.h_20, '
              'z.h_21, z.h_22, z.h_23, z.h_24, z.h_25, z.h_26, z.h_27, z.h_28, z.h_29, z.h_30, z.h_31 '
              'FROM MIUB_messdaten.MRR  AS Base '
              'LEFT OUTER JOIN MIUB_messdaten.MRR_z AS z ON Base.mes_id = z.mes_id '
              'LEFT OUTER JOIN MIUB_messdaten.MRR_heights AS h ON Base.mes_id = h.mes_id '
              'WHERE Base.Geraete_idGeraete = 40 '
              'AND Base.time BETWEEN ' + t_start.strftime('"%Y-%m-%d %H:%M:%S" ') +
              'AND' + t_stop.strftime(' "%Y-%m-%d %H:%M:%S" ')+
              'ORDER BY Base.time ASC'
             )
    #print(query1)
    cursor.execute(query1)
    MRR1 = np.asarray(cursor.fetchall())

    x = MRR1[:,0]
    y = MRR1[0,1:32]
    z = np.array(MRR1[:,32:],dtype=float)

    fig=pl.figure()
    pl.pcolormesh(x,y,z.T,vmin=-10,vmax=35)
    #pl.pcolormesh(x,y,np.ma.masked_where(np.isnan(z.T),z.T))
    pl.colorbar(label='z in dBz')
    pl.ylim(min(y),max(y))
    pl.title('MRR1 from ' + t_start.isoformat() + ' to ' + t_stop.isoformat())
    pl.ylabel('Height in m')
    fig.autofmt_xdate()
    pl.show()


plot_mrr1(t_start='2016-6-12 20:55', time_delta=20)