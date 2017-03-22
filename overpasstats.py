"""

Program fuer Statistik er Overpasses Daten

"""

import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

a = pd.read_csv('/automount/user/velibor/SHKGPM/prog/output3radolandpr.csv', sep=',')


fr = [1,2,3]
tit = ['Maximal detected Rainrate','Minimal detected Rainrate','Mean Rainrate']

for ii in fr:

    plt.subplot(3,1,ii)
    plt.plot(a.index, a[[ii]])
    plt.title(tit[ii-1])
    plt.ylabel('Rainrate in mm/h \n')
    plt.xlim(0, a.shape[0])
    plt.grid()


plt.xlabel('timesteps of overpasses')
plt.show()


