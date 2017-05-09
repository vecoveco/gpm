import pandas as pd
import matplotlib.pyplot as plt
import numpy as np



a_bd = pd.read_csv('/automount/ags/velibor/text/overpass_dpr_boxpolrr.csv', sep=',', header=0)
a_bd = a_bd.set_index('time')
a_rd = pd.read_csv('/automount/ags/velibor/text/dpr_stat.csv', sep=',', header=0)
a_rd = a_rd.set_index('time')

'''
a = pd.read_csv('/home/velibor/shkgpm/texte/dpr_stat_rr.csv', sep=',')
b= a.set_index('time')

y = a.values

for i in range(18):
    plt.subplot(6,3,i+1)
    plt.plot(y[:,i+1])
    plt.title(a.columns[i+1])
    plt.grid()

plt.tight_layout()
plt.show()


plt.subplot(2,2,1)
plt.plot(y[:,[1]], label='r_value')
plt.plot(y[:,[2]], label='std_err')
plt.plot(y[:,[8]], label='HR')
plt.plot(y[:,[9]], label='POD')
plt.plot(y[:,[10]], label='FAR')
plt.plot(y[:,[11]], label='BID')
plt.plot(y[:,[12]], label='HSS')
plt.xlabel('overpass')
plt.grid()
plt.legend()

plt.subplot(2,2,2)
plt.plot(y[:,4]/y[:,3]*100, label='H in %')
plt.plot(y[:,5]/y[:,3]*100, label='M in %')
plt.plot(y[:,6]/y[:,3]*100, label='F in %')
plt.plot(y[:,7]/y[:,3]*100, label='C in %')
plt.legend()
plt.xlabel('overpass')
plt.grid()

plt.subplot(2,2,3)
plt.plot(y[:,[13]], label='bias')
plt.plot(y[:,[14]], label='RMSE')
plt.legend()
plt.grid()
plt.xlabel('overpass')


plt.subplot(2,2,4)
plt.plot(y[:,[15]], label='Mean GPM in dBZ')
plt.plot(y[:,[16]], label='Mean RADOLAN in dBZ')
plt.plot(y[:,[17]], label='Median GPM in dBZ')
plt.plot(y[:,[18]], label='Median RADOLAN in dBZ')
plt.legend()
plt.xlabel('overpass')
plt.grid()
plt.show()


plt.scatter(np.arange(0,len(y),1), y[:,4], c=y[:,[1]],marker='o',s=50)
plt.colorbar()
plt.xlim(0,len(y))
plt.title("DPR Overpasses with Number of Hits and Correlation")
#plt.hlines(0.5)
plt.show()

def plot_stat(b):

    c = b.copy()
    c.index=c.index.to_datetime()

    ff=15
    fig = plt.figure(figsize=(10,8))
    fig.autofmt_xdate()
    c['r_value'].plot(style='.', color='black')
    pd.rolling_mean(c['r_value'],40, center=False,
                    min_periods=2).plot(linewidth=3, color='red',
                                        label='Rollin Median \n Window 40')

    plt.xticks(rotation=45)
    plt.ylabel('Correlation', fontsize=ff)
    plt.legend(loc='lower right')
    plt.grid()

plt.show()


plot_stat(a_rd)
plot_stat(a_rd)
plt.show()





c = a_rd.copy()
c.index=c.index.to_datetime()

ff=15
fig = plt.figure(figsize=(10,8))
fig.autofmt_xdate()
c['r_value'].plot(style='.', color='black')
pd.rolling_mean(c['r_value'],40, center=False,
                min_periods=2).plot(linewidth=3, color='red',
                                    label='Rollin Median \n Window 40')

plt.xticks(rotation=45)
plt.ylabel('Correlation', fontsize=ff)
plt.legend(loc='lower right')
plt.grid()
plt.show()
'''

df = a_rd
df2 = a_bd

time = np.array(df.index)#.strftime("%Y%m%d-%H%M%S")
time2 = np.array(df2.index)#.strftime("%Y%m%d-%H%M%S")


cor = df['r_value'].values
err = df['std_err'].values
Z = df['H'].values
#plt.plot(cor, 'o')
ff = 15

def statplot(cor, err, Z):
    from pcc import get_my_cmap
    plt.errorbar(range(0,len(cor)),cor, yerr=err, fmt='o', zorder=1, color='grey')
    plt.scatter(range(0,len(cor)),cor,c=Z, vmin=100, vmax=1000, cmap=get_my_cmap(), zorder=2)
    cb = plt.colorbar(shrink=0.5)
    cb.set_label("Hits in #",fontsize=ff)
    cb.ax.tick_params(labelsize=ff)
    #plt.show()
    plt.ylim(-1,1)
    plt.xlim(0,len(cor))
    plt.grid()
    plt.ylabel('Correlation',fontsize=ff)
    plt.xlabel('overpasses with time',fontsize=ff)
    plt.title('Overpass statistics betwen DPR and RADOLAN',fontsize=ff)

statplot(cor, err, Z)
plt.show()

def plot_ostat(V, tit):
    fig = plt.figure(figsize=(16,10))
    ff = 20
    fig.suptitle('Number of '+ tit +' for each Overpass', fontsize = ff)
    ax1 = plt.subplot(121)
    plt.plot(V,'ok')
    plt.ylabel('#')
    plt.xlabel('Overpass')
    plt.grid(color='blue')

    ax2 = plt.subplot(122)
    plt.hist(V, bins=len(V), orientation='horizontal')
    plt.grid(color='blue')
    plt.xlabel('freq')
    plt.show()
    #plt.savefig('/automount/ftp/velibor/GPM/Stat/Overpassstat_dpr-rado_'+str(tit)+'.png')
    #plt.close()


'''
for i in range(17):
    print i
    df.columns[i]
    try:
        plot_ostat(df[df.columns[i]].values, df.columns[i])
    except:
        print str(df.columns[i]) + 'NICHT BERECHNET!'

plt.close()
'''
#df.values[:,1][np.where(df['N']<300)]= np.nan
#plot_ostat(df['HSS'].values, 'HSST')
#plot_ostat(df['BID'].values, 'BID')
#plot_ostat(df['FAR'].values, 'FAR')
#plot_ostat(df['r_value'].values, 'Corr')
#plot_ostat(df['POD'].values, 'POD')


#df[np.where(df['N']<300)] = np.nan#

A = df['r_value'].values.copy()
B = df['std_err'].values.copy()
C = df['H'].values.copy()

T = 100
A[np.where(C<T)] = np.nan
B[np.where(C<T)] = np.nan
C[np.where(C<T)] = np.nan



from pcc import get_my_cmap
plt.errorbar(range(0,len(A)),A, B, fmt='o', zorder=1, color='grey')
plt.scatter(range(0,len(A)),A,c=C,s=300,linewidths=0.001, vmin=T, vmax=1000, cmap=get_my_cmap(), zorder=2)
cb = plt.colorbar(shrink=0.5)
cb.set_label("Hits in #",fontsize=ff)
cb.ax.tick_params(labelsize=ff)
plt.ylim(-1,1)
plt.xlim(0,len(A))
plt.grid()
plt.ylabel('Correlation',fontsize=ff)
plt.xlabel('overpasses with time',fontsize=ff)
plt.title('Overpass statistics betwen DPR and RADOLAN, \n Threshold: N = '+str(T),fontsize=ff)
plt.show()

dataframes = [df, df2]
names = ['RADOLAN','BoXPol']

for j in range(2):
    DF = dataframes[j]
    hits = DF['H'].values.copy()

    th = np.arange(np.nanmax(hits))
    anzahl_overpass = np.zeros(len(th))

    for i in range(len(th)):
        anzahl_overpass[i] = len(hits[hits >th[i]])

    plt.plot(th, anzahl_overpass, label=names[j])
    plt.legend(loc='upper right')
    plt.grid(color='blue')
    plt.xlabel('Threshold in hits')
    plt.ylabel('Overpasses')

plt.show()