import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

data=np.genfromtxt('sobol.dat',dtype=None,delimiter=" ")

dataX=data[0:,0]
dataY=data[0:,1]
dataZ=data[0:,2]
dataW=data[0:,3]

figXY = plt.figure()
ax = figXY.add_subplot(111)
t = ax.scatter(dataX, dataY)
figXY.show()


figXZ = plt.figure()
ax = figXZ.add_subplot(111)
t = ax.scatter(dataX, dataZ)
figXZ.show()


figXW = plt.figure()
ax = figXW.add_subplot(111)
t = ax.scatter(dataX, dataW)
figXW.show()


figYZ = plt.figure()
ax = figYZ.add_subplot(111)
t = ax.scatter(dataY, dataZ)
figYZ.show()


figYW = plt.figure()
ax = figYW.add_subplot(111)
t = ax.scatter(dataY, dataW)
figXW.show()

figZW = plt.figure()
ax = figZW.add_subplot(111)
t = ax.scatter(dataZ, dataW)
figZW.show()