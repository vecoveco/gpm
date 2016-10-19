import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mmm
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats
import matplotlib.cm as cm


'''
x = np.array([-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0])  # Werte
L=len(x)  # bins

S=(L+1)/2  # Summer der Betrege

W = x/ S
W = W / L

y = 2*x/ (L*(L+1))

A = x
print (A)
print (abs(A))
print (np.sum(abs(A)))



plt.plot(A)
plt.plot(W)
plt.plot(y, 'ok')

plt.grid()
plt.show()
'''

# test for 1 dimension in space and two value dimensions
trg = np.arange(10)[:,None]
src = np.linspace(0,20,40)[:,None]
vals = np.hstack((np.sin(src), 10.+np.sin(src)))
# here we introduce missing values only in the second dimension
vals[3:5,1] = np.nan
ipol_result = wradlib.ipol.interpolate(src, trg, vals, wradlib.ipol.Idw, nnearest=2)
import matplotlib.pyplot as plt
plt.plot(src, vals, 'ro')
plt.plot(trg, ipol_result, 'b+')

plt.show()