print('start')

import math
print('imported math')

import matplotlib 
matplotlib.use('TkAgg')
print('imported matplotlib')

import matplotlib.pyplot as plt
print('matplotlib.pyplot imported')

import numpy as np
print('imported numpy')

import scipy.interpolate
print('imported scypy')


idummy,zsl, rsl,d2,d3,d4,d5,d6 = np.genfromtxt(r'sls.DAT', unpack=True)
print('data2 read')

fig = plt.figure()

plt.xlabel("z")
plt.ylabel("r")
plt.legend(loc='best')
plt.title('streamlines cT = 1')
plt.plot(zsl,rsl)
plt.show()
#plt.savefig('streamline-opt-vK-fit.png')
#print('ready')
#plt.close()

