import numpy as np
import matplotlib.pyplot as plt

def graphy_energy(string_name):
    data = np.loadtxt(string_name, dtype=float)
    x = data[:,0]
    y = data[:,1]
    
    
    plt.title(r"$<N(t)>\approx${:.0f}".format(y.mean()),fontsize=18)
    plt.plot(x, y, marker=".", color="k")
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlabel(r"$t$", fontsize=18)
    plt.ylabel(r"$N(t)$", fontsize=18)

    plt.show()
    
def main():
    graphy_energy("time.txt")

if __name__ == '__main__':
    main()

