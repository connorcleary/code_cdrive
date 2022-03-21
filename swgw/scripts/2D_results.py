import numpy as np
import matplotlib.pyplot as plt
import os
import post_proc_utils as proc

model_mults = ["1", "2", "3", "4"]
repos = [f".\\results\\pirot2D_{mult}" for mult in model_mults]

def load_data(model_mults, repos):
    conc_het = {}
    conc_hom = {} 
    qx_het = {}
    qx_hom = {}
    qy_het = {}
    qy_hom = {}
    qz_het = {} 
    qz_hom = {}
    head_het = {}
    head_hom = {}

    for mult, repo in zip(model_mults, repos):
        conc_het[mult] = np.genfromtxt(os.path.join(repo, 'concentration_heterogenous'))
        conc_hom[mult] = np.genfromtxt(os.path.join(repo, 'concentration_homogenous'))
        qx_het[mult] = np.genfromtxt(os.path.join(repo, 'qx_heterogenous'))
        qx_hom[mult] = np.genfromtxt(os.path.join(repo, 'qx_homogenous'))
        qy_het[mult] = np.genfromtxt(os.path.join(repo, 'qy_heterogenous'))
        qy_hom[mult] = np.genfromtxt(os.path.join(repo, 'qy_homogenous'))
        qz_het[mult] = np.genfromtxt(os.path.join(repo, 'qz_heterogenous'))
        qz_hom[mult] = np.genfromtxt(os.path.join(repo, 'qz_homogenous'))
        head_het[mult] = np.genfromtxt(os.path.join(repo, 'head_heterogenous'))
        head_hom[mult] = np.genfromtxt(os.path.join(repo, 'head_homogenous'))
        
    return conc_het, conc_hom, qx_het, qx_hom, qy_het, qy_hom, qz_het, qz_hom, head_het, head_hom

def plot_for_het_hom(y_het, y_hom, x, title , xlabel, ylabel):
    fig, ax = plt.subplots()
    het = ax.scatter(x, y_het, c='red', label='heterogenous')
    hom = ax.scatter(x, y_hom, c='blue', label = 'homogenous')
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    #ax.set_ylim(bottom=0)
    ax.legend()
    plt.show()

    fig, ax = plt.subplots()
    difference = [100*het/hom for het, hom in zip(y_het, y_hom)]
    diff = ax.scatter(x, difference)
    ax.set_title("Difference in" + title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Percentage difference (%)")
    plt.show()

    fig, ax = plt.subplots()
    difference = [het-hom for het, hom in zip(y_het, y_hom)]
    diff = ax.scatter(x, difference)
    ax.set_title("Difference in" + title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Absolute difference (m)")
    plt.show()

# import files and compile arrays of data
def results(model_mults, repos):
    conc_het, conc_hom, qx_het, qx_hom, qy_het, qy_hom, qz_het, qz_hom, head_het, head_hom = load_data(model_mults, repos)
    aquifer_lengths = [800*mult for mult in range(1, 5)]
   
    # mixing_zones_het = []
    # mixing_zones_hom = []
    # for chet, chom in zip(conc_het.values(), conc_hom.values()):
    #     mixing_zones_het.append(proc.find_mixing_zone(chet))
    #     mixing_zones_hom.append(proc.find_mixing_zone(chom))
    # plot_for_het_hom(mixing_zones_het, mixing_zones_hom, aquifer_lengths, title = "Mixing Zone Area", xlabel = "Aquifer Length (m)", ylabel = "Area (m^2)")

    toe_pen_het = []
    toe_pen_hom = []
    for chet, chom in zip(conc_het.values(), conc_hom.values()):
        toe_pen_het.append(proc.find_toe_penetration(chet))
        toe_pen_hom.append(proc.find_toe_penetration(chom))
    plot_for_het_hom(toe_pen_het, toe_pen_hom, aquifer_lengths, title = "Toe Penetration", xlabel = "Aquifer Length (m)", ylabel = "Penetration (m)")


if __name__ == "__main__":
    results(model_mults, repos)
    pass
