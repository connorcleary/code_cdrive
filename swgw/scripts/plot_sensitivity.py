from pkgutil import extend_path
import numpy as np
import matplotlib.pyplot as plt

mixing_zone = np.genfromtxt(r'.\results\sensitivity\mixing_zone_size.csv', delimiter=',')
toe = np.genfromtxt(r'.\results\sensitivity\extent.csv', delimiter=',')
inflow = -1*np.genfromtxt(r'.\results\sensitivity\inflow.csv', delimiter=',')
outflow = np.genfromtxt(r'.\results\sensitivity\outflow.csv', delimiter=',')

qinflow = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 5.702)  # m3/day
dmcoef = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 0.57024)  # m2/day  Could also try 1.62925 as another case of the Henry problem
hk = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 864.0)  # m/day

for result, name in zip([mixing_zone, toe, inflow, outflow], ['mixing_zone', 'toe', 'inflow', 'outflow']):
    f, axs = plt.subplots(1, 3, sharey=True)
    f.suptitle(f'Sensitivity of {name}')

    for i, (param, value) in enumerate(zip(['qinflow', 'dmcoef', 'hk'], [qinflow, dmcoef, hk])):
        axs[i].scatter(value, result[i,:])
        axs[i].set_xlabel(param)
        axs[i].set_ylabel(name)

    plt.savefig(f'.\\results\\sensitivity\\{name}.png')

