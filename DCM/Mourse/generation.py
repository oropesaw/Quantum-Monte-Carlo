import numpy as np
import matplotlib.pyplot as plt

def graphy_energy(string_name):
    data = np.loadtxt(string_name, dtype=float)
    y = data[:,0]
    x = data[:,3]
    yerr = data[:,1]
    
    lines = {'linestyle': 'None'}
    plt.rc('lines', **lines)
    plt.title("""Valor medio de la energía del estado fundamental en función del tamaño del paso de tiempo""")
    plt.errorbar(x, y, yerr=yerr, marker="o", color="k")
    plt.grid(True, linestyle="--", color="0.5")
    plt.xlabel(r"$\Delta t$", fontsize=18)
    plt.ylabel(r"$<E>$", fontsize=18)

    plt.show()
    
def main():
    graphy_energy("data.txt")

if __name__ == '__main__':
    main()
