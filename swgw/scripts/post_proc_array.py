from importlib_metadata import files
import numpy as np
import flopy
import glob 
import matplotlib.pyplot as plt

def load_files(prefix):
    # load all files with that prefix
    matrix = np.zeros((20, 50, 100))
    hk = np.linspace(864, 9000, 20)
    f_list = glob.glob(f"**/{prefix}*.*", recursive=True)

    for idx, file in enumerate(f_list):
        matrix[idx, :, :] = np.loadtxt(file)

    return matrix, hk

def plot_mixing_zone(matrix, hk):
    f, ax = plt.subplots()
    mixing_zone_size = np.zeros(20)
    for idx, real in enumerate(range(matrix.shape[0])):
        mix = 0
        for i in range(matrix.shape[1]):
            for j in range(matrix.shape[2]):
                if 3.5 <= matrix[real, i, j] <= 31.5:
                    mix += 1
        mixing_zone_size[idx] = 2*mix/5000
    ax.scatter(hk, mixing_zone_size, marker='x')
    ax.set_title("Area of mixing zone")
    ax.set_xlabel('K high')
    ax.set_ylabel('Area (m^2)')
    ax.set_aspect(0.5/ax.get_data_ratio())
    plt.show()

def plot_subset(matrix, hk, name):

    f,ax = plt.subplots(1, 3, sharey=True)
    for idx, plot in enumerate([0, 10, 19]):
        pcm = ax[idx].pcolormesh(np.flipud(matrix[plot]))
        ax[idx].set_ylabel("z (m)")
        ax[idx].set_xlabel("x (m)")
        ax[idx].set_aspect('equal')
        ax[idx].set_title(f"K_high = {hk[plot]:.2f}")
    f.suptitle(name)
    plt.colorbar(pcm,  ax=ax.ravel().tolist(), location="bottom")

    plt.show()
        

def mixing_zone_volume(concentration, delV):
    '''
    Calculate the volume of the mixing zone
    '''
    mix = 0
    for i in range(concentration.shape[0]):
        for j in range(concentration.shape[1]):
            for k in range(concentration.shape[2]):
                if 3.5 <= concentration[i][j][k] <= 31.5:
                    mix += 1
    return mix*delV

def fresh_water_volume(concentration, delV, colq):
    '''
    Calculate the volume of fresh water seaward of the pumping locationY
    '''
    fresh = 0
    for i in range(concentration.shape[0]):
        for j in range(concentration.shape[1]):
            for k in range(colq, concentration.shape[2]):
                if 3.5 > concentration[i][j][k]:
                    fresh += 1
    return fresh*delV

def well_salinity(concentration, lay, row, col):
    '''
    Find the concentration at the well over time
    '''
    pass


def main():
    # load concentrations 
    concentrations, hk = load_files("concentration")
    heads, _ = load_files("head")

    #plot_subset(concentrations, hk, "Saline concentration")
    plot_mixing_zone(concentrations, hk)

    pass


if __name__ == "__main__":
    main()