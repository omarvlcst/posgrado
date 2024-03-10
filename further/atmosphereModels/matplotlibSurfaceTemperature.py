from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np  

gamma = 0.5*np.pi # magnetization angle, in radians
Tequat5 = 3.0e5   # temperature in the magnetic equator * 10^5 kelvin
Tpolar5 = 7.0e5   # temperature in the hotspots (magnetic poles) * 10^5 kelvin
chi04 = (Tequat5/Tpolar5)**4

def tsup_norm(x,y,z):
    rsq = np.square(x) + np.square(y) + np.square(z)
    f = x*np.sin(gamma) + z*np.cos(gamma)
    f1 = 4*f**2
    f2 = 3*f + rsq
    ts4 = f1/f2 * (1.0-chi04) + chi04
    return Tpolar5*np.abs(ts4)**(0.25)

def sphere_surface(r):
    u = np.linspace(2 * np.pi/720, 2 * np.pi, 720)
    v = np.linspace(np.pi/360, np.pi, 360)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))
    return x,y,z

x,y,z = sphere_surface(Tpolar5)
tsup = tsup_norm(x,y,z)

norm=matplotlib.colors.SymLogNorm(1,vmin=tsup.min(),vmax=tsup.max())
colors=plt.cm.coolwarm(norm(tsup))

fig = plt.figure()
ax = fig.gca(projection='3d')
# plot the surface
surf = ax.plot_surface(x,y,z, facecolors=colors,
                       linewidth=0, antialiased=False)
# colorbar
sm = plt.cm.ScalarMappable(cmap=plt.cm.coolwarm, norm=norm)
sm.set_array(tsup)
fig.colorbar(sm, shrink=0.5, aspect=5)

ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

plt.show()
