import numpy as np
import matplotlib.pyplot as plt
a1 = np.arange(0,10,1)
a2 = np.arange(0,10,1)
aa1, aa2 = np.meshgrid(a1,a2)

A = np.sin(aa1) + np.cos(aa2)

plt.pcolormesh(aa1,aa2,A)
plt.scatter(aa1,aa2)

plt.show()



b1 = np.arange(0,5,1)
b2 = np.arange(0,5,1)
bb1, bb2 = np.meshgrid(b1,b2)

B = np.sin(bb1)+ np.cos(bb2)

plt.pcolormesh(bb1,bb2,B)
plt.scatter(bb1,bb2)
plt.show()


plt.plot