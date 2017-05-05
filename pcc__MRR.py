"""

Reading MRR Data

"""


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from io import StringIO


pfad_mrr = '/automount/mrr/mrr2/2014/2014-10/2014-10-07/AveData_mrr2_20141007023628.ave.gz'


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
plt.ylabel('Hight in m')
plt.title('MRR - ' + pfad_mrr[44:44+28])
plt.show()





