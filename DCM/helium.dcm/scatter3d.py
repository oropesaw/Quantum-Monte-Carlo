# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


data = np.loadtxt("psi.data", dtype=float)
xs = data[:,0]
ys = data[:,1]
zs = data[:,2]
ax.scatter(xs, ys, zs, marker=".")

plt.title("Función de onda del átomo en el estado fundamental")
ax.set_xlabel(r"$r_{1}$", fontsize=18)
ax.set_ylabel(r"$r_{2}$", fontsize=18)
ax.set_zlabel(r"$\psi(r_{1}, r_{2})$", fontsize=18)

plt.show()


