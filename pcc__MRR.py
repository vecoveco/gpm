"""

Reading MRR Data

"""

def plot_mrr1(t_start='2016-10-07 02:35', time_delta=10):
    """
    Creat by Jo. Beer

    """

    import numpy as np
    import pymysql as mysql
    import datetime as dt
    from datetime import datetime
    import matplotlib
    import matplotlib.pyplot as plt
    from pcc import get_miub_cmap
    my_cmap2 = get_miub_cmap()


    # Bonn 40
    # Roleber 42
    # Bergheim 41
    station = [40, 41, 42]
    name_mrr = ['Bonn', 'Bergheim', 'Roleber']
    fig=plt.figure(figsize=(18,10))
    #datetime_object = datetime.strptime('Jun 1 2005  1:33PM', '%b %d %Y %I:%M%p')
    t_start = datetime.strptime(t_start,'%Y-%m-%d %H:%M')
    t_stop  = t_start + dt.timedelta(minutes=time_delta)#time_delta

    for iii in range(3):



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
                  'WHERE Base.Geraete_idGeraete = '+str(station[iii]) +' '
                  'AND Base.time BETWEEN ' + t_start.strftime('"%Y-%m-%d %H:%M:%S" ') +
                  'AND' + t_stop.strftime(' "%Y-%m-%d %H:%M:%S" ')+
                  'ORDER BY Base.time ASC'
                 )
        #print(query1)
        cursor.execute(query1)
        MRR1 = np.asarray(cursor.fetchall())

        if MRR1.shape==(0,):
            pass
        else:
            x = MRR1[:,0]
            y = MRR1[0,1:32]
            z = np.array(MRR1[:,32:],dtype=float)


            plt.subplot(3,1,iii+1)
            plt.pcolormesh(x,y,z.T,vmin=-10,vmax=35, cmap=my_cmap2)
            # Idee zeige BB
            #plt.plot(x,y[np.argmax(z.T,axis=0)],lw=1.5,color='grey')

            #pl.pcolormesh(x,y,np.ma.masked_where(np.isnan(z.T),z.T))
            plt.colorbar(label='Z in DBz')
            plt.ylim(min(y),max(y))
            plt.title(name_mrr[iii]+' MRR1 from ' + t_start.isoformat() + ' to ' + t_stop.isoformat())
            plt.ylabel('Height in m')
            plt.grid()
            fig.autofmt_xdate()
            import matplotlib.dates as mdates
            #fig.autofmt_xdate = mdates.DateFormatter('%Y-%m-%d')
    plt.tight_layout()
    plt.show()


#plot_mrr1(t_start='2014-07-29 22:00', time_delta=60)
plot_mrr1(t_start='2016-02-09 10:00', time_delta=60)