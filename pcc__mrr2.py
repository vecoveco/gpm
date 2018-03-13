
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob


yy = '2016'
mm = '10'
dd = '07'
hh = '02'
min = '35'

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

