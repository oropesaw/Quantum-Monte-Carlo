# -*- coding: utf-8 -*-

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


data = np.loadtxt("data.txt", dtype=float)
xs = data[:,4]
ys = data[:,3]
zs = data[:,0]
ax.scatter(xs, ys, zs, marker=".")


ax.set_xlabel(r"$\alpha$", fontsize=18)
ax.set_ylabel(r"$\beta$", fontsize=18)
ax.set_zlabel(r"$<E>$", fontsize=18)

plt.show()


