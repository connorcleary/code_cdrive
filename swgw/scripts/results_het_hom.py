import matplotlib.pyplot as plt
import numpy as np 
import flopy
import post_proc_utils as proc

def single_results(modelname):
    pass

def results(modelname):

    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'
    heterogenous = {}
    homogenous = {}
    for case_dict, case_name in zip([heterogenous, homogenous],['heterogenous', 'homogenous']):
        case_dict['qx'] = np.genfromtxt(f"{results_location}\qx_{case_name}")
        case_dict['qy'] = np.genfromtxt(f"{results_location}\qy_{case_name}")
        case_dict['qz'] = np.genfromtxt(f"{results_location}\qz_{case_name}")
        case_dict['head'] = np.genfromtxt(f"{results_location}\head_{case_name}")
        case_dict['concentration'] = np.genfromtxt(f"{results_location}\concentration_{case_name}")

    fig, axs = plt.subplots(2, 2, figsize=(18, 9))
    axs[0][0] = proc.plot_conc(axs[0][0], swt, heterogenous)
    axs[1][0] = proc.plot_head(axs[1][0], swt, heterogenous)
    axs[0][1] = proc.plot_conc(axs[0][1], swt, homogenous)
    axs[1][1] = proc.plot_head(axs[1][1], swt, homogenous)
    plt.savefig(f"{results_location}\\comparison_{modelname}.jpg", dpi=1200)

def results_3D(modelname):
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'

    qx, qy, qz, head, concentration = proc.load_results_3D(modelname, "heterogenous")
    qx_eq, qy_eq, qz_eq, head_eq, concentration_eq = proc.load_results_3D(modelname, "homogenous")

    nrow = swt.dis.nrow
    for iplot in range(int(np.floor(int(nrow)/5))):
        rows = [irow for irow in range(iplot*5, np.min([nrow, iplot*5+5]))]
        nhetsub = len(rows)
        f, axs = plt.subplots(nhetsub+1,2, sharex=True, figsize=(8, 12))
        axs[0][0] = proc.plot_head(axs[0][0], swt, head_eq[:, :, :], 0)
        axs[0][1] = proc.plot_conc(axs[0][1], swt, {'qx':qx_eq[:, :, :], 'qy':qy_eq[:, :, :], 'qz':qz_eq[:, :, :], 'concentration':concentration_eq[:, :, :]}, 0)
        for i, row in zip(range(1,nhetsub+1), rows): 
            axs[i][0] = proc.plot_head(axs[i][0], swt, head[:, :, :], row)
            axs[i][1] = proc.plot_conc(axs[i][1], swt, {'qx':qx[:, :, :], 'qy':qy[:, :, :], 'qz':qz[:, :, :], 'concentration':concentration[:, :, :]}, row)

        plt.savefig(f"{results_location}\\comparison_{modelname}_rows_{rows[0]}-{rows[-1]}.jpg", dpi=1200)

def results_3D_single(modelname, lay, row, col):
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'

    qx, qy, qz, head, concentration = proc.load_results_3D(modelname, "heterogenous")
    f, axs = plt.subplots(1, 2, figsize=(12,9))
    row = 0
    axs[0] = proc.plot_head(axs[0], swt, head[:, :, :], vmin=-10, vmax=0.6)
    axs[1] = proc.plot_conc(axs[1], swt, {'qx':qx[:, :, :], 'qy':qy[:, :, :], 'qz':qz[:, :, :], 'concentration':concentration[:, :, :]}, row, vmax=35, vmin=0)
    plt.savefig(f"{results_location}\\pumping_{modelname}_law{lay}_col{col}_row{row}.jpg", dpi=1200)

def results_pumping_recovery(modelname, lay, row, col):
    
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'
    qx, qy, qz, head, concentration = proc.load_results_3D(modelname, "heterogenous")
    f, axs = plt.subplots(3, 1, figsize = (9, 12))
    axs[0] = proc.plot_conc(axs[0], swt, {'qx':qx.item().get('steady')[:, :, :], 'qy':qy.item().get('steady')[:, :, :], 'qz':qz.item().get('steady')[:, :, :], 'concentration':concentration.item().get('steady')[:, :, :]}, row, vmax=35, vmin=0)
    axs[0].set_title("Steady state concentration")
    axs[1] = proc.plot_conc(axs[1], swt, {'qx':qx.item().get('pumping')[:, :, :], 'qy':qy.item().get('pumping')[:, :, :], 'qz':qz.item().get('pumping')[:, :, :], 'concentration':concentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1].set_title("Concentration after pumping")
    axs[2] = proc.plot_conc(axs[2], swt, {'qx':qx.item().get('recovery')[:, :, :], 'qy':qy.item().get('recovery')[:, :, :], 'qz':qz.item().get('recovery')[:, :, :], 'concentration':concentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2].set_title("Concentration after recovery")
    axs[1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    plt.savefig(f"{results_location}\\pumping_and_recovery{modelname}.jpg", dpi=1200)
    plt.show()



def probability_of_saline(modelname):
    _, _, _, _, concentration = proc.load_results_3D(modelname, "heterogenous")
    _, _, _, _, concentration_eq = proc.load_results_3D(modelname, "homogenous")
    concentration_counts = np.zeros((concentration.shape[0],concentration.shape[2]))
    for row in range(concentration.shape[1]):
        for lay in range(concentration.shape[0]):
            for col in range(concentration.shape[2]):
                if concentration[lay, row, col] > 0.35:
                    concentration_counts[lay, col] += 1
    homogenous_y = np.linspace(50, 0, 50)
    homogenous_interface = []

    for  lay in range(concentration.shape[0]):
        for col in range(concentration.shape[2]):
            if concentration_eq[lay, 0, col] > 0.35:
                homogenous_interface.append(col)
                break

    # x = np.linspace(0, 800, 80)
    # y = np.linspace(0, 25, 50)
    # X, Y = np.meshgrid(x, y)
    conc_prob = concentration_counts/concentration.shape[1]
    f, ax = plt.subplots()
    prob = ax.pcolor(np.flipud(conc_prob), cmap="coolwarm")
    # ax.plot(X, Y, homogenous_interface, homogenous_y, color='black')
    ax.plot(homogenous_interface, homogenous_y, color='black')
    cb = plt.colorbar(prob)
    ax.set_title("Probability of salinization")
    ax, 
    plt.show()
    
    pass

if __name__=="__main__":
    results("pirot2D_check_k_eff")