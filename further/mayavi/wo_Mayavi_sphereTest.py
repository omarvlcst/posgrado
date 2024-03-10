import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d
from matplotlib import cm

fig = plt.figure()
ax = fig.gca(projection='3d')

theta_coord = np.linspace(np.pi/180, np.pi, 180)
phi_coord = np.linspace(2*np.pi/360, 2*np.pi, 360)
T,P = np.meshgrid(theta_coord,phi_coord)
X, Y, Z = np.cos(P)*np.sin(T), np.sin(P)*np.sin(T), np.cos(T)


ax.scatter3D(X, Y, Z,s=0.01)
plt.show()
