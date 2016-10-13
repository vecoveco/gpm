import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import wradlib
import glob
import math
import pandas as pd
from scipy import stats


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