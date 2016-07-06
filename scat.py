
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
from scipy import stats
import pandas as pd

cor_jo = pd.read_csv('/user/velibor/SHKGPM/data/plot/cor_jo.csv')
std_jo = pd.read_csv('/user/velibor/SHKGPM/data/plot/std1_jo.csv')
cor_v = pd.read_csv('/user/velibor/SHKGPM/data/plot/cor_v.csv')
std_v = pd.read_csv('/user/velibor/SHKGPM/data/plot/std1_v.csv')

print (cor_jo.shape)
print (cor_v[[1]].values)

j = cor_jo[[1]].values
v = cor_v[[1]].values

v = v[0:6]

A = j
B = v
mask = ~np.isnan(B) & ~np.isnan(A)
slope, intercept, r_value, p_value, std_err = stats.linregress(B[mask], A[mask])
line = slope*B+intercept
plt.plot(B,line,'r-',B,A,'o')
maxAB = np.nanmax([np.nanmax(A),np.nanmax(B)])
plt.xlim(0,maxAB)
plt.ylim(0,maxAB)
plt.xlabel("Correlation Raw")
plt.ylabel("Correlation Attenuation Corrected")
plt.grid(True)
plt.title("Scatterplot, cor: " + str(r_value))
plt.show()

