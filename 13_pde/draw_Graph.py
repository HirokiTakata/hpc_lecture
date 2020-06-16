import numpy as np
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

nx = 41
ny = 41
X = np.zeros((ny,nx))
Y = np.zeros((ny,nx))
u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
read_file = open("output.txt","r")
lines = read_file.readlines()
l = 0
for i in range(nx):
	for j in range(ny):	
		X[i][j] = lines[l*4]
		Y[i][j] = lines[l*4+1]
		u[i][j] = lines[l*4+2]
		v[i][j] = lines[l*4+3]
		l+=1
		
	

fig = pyplot.figure(figsize = (11,7), dpi=100)
pyplot.quiver(X, Y, u, v)
pyplot.savefig("output.png")
