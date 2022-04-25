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
    plt.show()

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
    f, axs = plt.subplots(3, 1, figsize = (18, 12))
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

def results_pumping_recovery_comparison(modelname, lay, row, col):
    
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'
    qx, qy, qz, head, concentration = proc.load_results_3D(modelname, "heterogenous")
    eqqx, eqqy, eqqz, eqhead, eqconcentration = proc.load_results_3D(modelname, "homogenous")
    f, axs = plt.subplots(3, 2, figsize = (18, 12))
    axs[0][0] = proc.plot_conc(axs[0][0], swt, {'qx':qx.item().get('steady')[:, :, :], 'qy':qy.item().get('steady')[:, :, :], 'qz':qz.item().get('steady')[:, :, :], 'concentration':concentration.item().get('steady')[:, :, :]}, row, vmax=35, vmin=0)
    axs[0][0].set_title("Steady state concentration")
    axs[1][0] = proc.plot_conc(axs[1][0], swt, {'qx':qx.item().get('pumping')[:, :, :], 'qy':qy.item().get('pumping')[:, :, :], 'qz':qz.item().get('pumping')[:, :, :], 'concentration':concentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1][0].set_title("Concentration after pumping")
    axs[2][0] = proc.plot_conc(axs[2][0], swt, {'qx':qx.item().get('recovery')[:, :, :], 'qy':qy.item().get('recovery')[:, :, :], 'qz':qz.item().get('recovery')[:, :, :], 'concentration':concentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2][0].set_title("Concentration after recovery")
    axs[1][0].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][0].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[0][1] = proc.plot_conc(axs[0][1], swt, {'qx':eqqx.item().get('steady')[:, :, :], 'qy':eqqy.item().get('steady')[:, :, :], 'qz':eqqz.item().get('steady')[:, :, :], 'concentration':eqconcentration.item().get('steady')[:, :, :]}, row, vmax=35, vmin=0)
    axs[0][1].set_title("Steady state concentration")
    axs[1][1] = proc.plot_conc(axs[1][1], swt, {'qx':eqqx.item().get('pumping')[:, :, :], 'qy':eqqy.item().get('pumping')[:, :, :], 'qz':eqqz.item().get('pumping')[:, :, :], 'concentration':eqconcentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1][1].set_title("Concentration after pumping")
    axs[2][1] = proc.plot_conc(axs[2][1], swt, {'qx':eqqx.item().get('recovery')[:, :, :], 'qy':eqqy.item().get('recovery')[:, :, :], 'qz':eqqz.item().get('recovery')[:, :, :], 'concentration':eqconcentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2][1].set_title("Concentration after recovery")
    axs[1][1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    plt.savefig(f"{results_location}\\pumping_and_recovery_comparison_{modelname}.jpg", dpi=1200)
    plt.show()

def results_pumping_recovery_comparison_with_heterogeneity(modelname, lay, row, col, hk, eqhk):
    
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'
    qx, qy, qz, head, concentration = proc.load_results_3D(modelname, "heterogenous")
    eqqx, eqqy, eqqz, eqhead, eqconcentration = proc.load_results_3D(modelname, "homogenous")

    f, axs = plt.subplots(4, 2, figsize = (18, 18))

    axs[0][0] = proc.plot_conc(axs[0][0], swt, {'qx':qx.item().get('steady')[:, :, :], 'qy':qy.item().get('steady')[:, :, :], 'qz':qz.item().get('steady')[:, :, :], 'concentration':concentration.item().get('steady')[:, :, :]}, row, vmax=35, vmin=0)
    axs[0][0].set_title("Steady state concentration")
    axs[1][0] = proc.plot_conc(axs[1][0], swt, {'qx':qx.item().get('pumping')[:, :, :], 'qy':qy.item().get('pumping')[:, :, :], 'qz':qz.item().get('pumping')[:, :, :], 'concentration':concentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1][0].set_title("Concentration after pumping")
    axs[2][0] = proc.plot_conc(axs[2][0], swt, {'qx':qx.item().get('recovery')[:, :, :], 'qy':qy.item().get('recovery')[:, :, :], 'qz':qz.item().get('recovery')[:, :, :], 'concentration':concentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2][0].set_title("Concentration after recovery")
    axs[1][0].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][0].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[0][1] = proc.plot_conc(axs[0][1], swt, {'qx':eqqx.item().get('steady')[:, :, :], 'qy':eqqy.item().get('steady')[:, :, :], 'qz':eqqz.item().get('steady')[:, :, :], 'concentration':eqconcentration.item().get('steady')[:, :, :]}, row, vmax=35, vmin=0)
    axs[0][1].set_title("Steady state concentration")
    axs[1][1] = proc.plot_conc(axs[1][1], swt, {'qx':eqqx.item().get('pumping')[:, :, :], 'qy':eqqy.item().get('pumping')[:, :, :], 'qz':eqqz.item().get('pumping')[:, :, :], 'concentration':eqconcentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1][1].set_title("Concentration after pumping")
    axs[2][1] = proc.plot_conc(axs[2][1], swt, {'qx':eqqx.item().get('recovery')[:, :, :], 'qy':eqqy.item().get('recovery')[:, :, :], 'qz':eqqz.item().get('recovery')[:, :, :], 'concentration':eqconcentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2][1].set_title("Concentration after recovery")
    axs[1][1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)

    x = np.linspace(0, 800, 80)
    y = np.linspace(-25, 0, 50)

    cmhk = axs[3][0].pcolormesh(x, y, np.flipud(np.log(hk)), cmap="coolwarm")
    cmeqhk = axs[3][1].pcolormesh(x, y, np.ones_like(hk)*np.log(eqhk), cmap="coolwarm", vmax=np.max(np.log(hk)), vmin= np.min(np.log(hk)))
    axs[3][1].set_aspect(10)
    axs[3][0].set_aspect(10)
    plt.colorbar(cmhk,ax=axs[3][0])
    plt.colorbar(cmeqhk,ax=axs[3][1])

    plt.savefig(f"{results_location}\\pumping_and_recovery_comparison_{modelname}.jpg", dpi=1200)
    plt.show()    

def results_single_staged(modelname, realization, hk, lay, row, col):

    qx, qy, qz, head, concentration =  proc.load_results_3D(modelname, realization, stress_period="transient")
    qx_s, qy_s, qz_s, head_s, concentration_s =  proc.load_results_3D(modelname, realization, stress_period="steady")

    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}_steady.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{modelname}'

    f, axs = plt.subplots(4, 1, figsize = (8, 12))
    axs[0], _ = proc.plot_conc(axs[0], swt, {'qx':qx_s[:, :, :], 'qy':qy_s[:, :, :], 'qz':qz_s[:, :, :], 'concentration':concentration_s[:, :, :]}, row, vmax=35, vmin=0)
    axs[0].set_title("Steady state concentration")
    axs[1], _ = proc.plot_conc(axs[1], swt, {'qx':qx.item().get('pumping')[:, :, :], 'qy':qy.item().get('pumping')[:, :, :], 'qz':qz.item().get('pumping')[:, :, :], 'concentration':concentration.item().get('pumping')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[1].set_title("Concentration after pumping")
    axs[2], _ = proc.plot_conc(axs[2], swt, {'qx':qx.item().get('recovery')[:, :, :], 'qy':qy.item().get('recovery')[:, :, :], 'qz':qz.item().get('recovery')[:, :, :], 'concentration':concentration.item().get('recovery')[:, :, :]}, row, vmax=35, vmin=0, zorder=1)
    axs[2].set_title("Concentration after recovery")
    axs[1].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2].scatter(col*10+5, lay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)

    x = np.linspace(0, 800, 80)
    y = np.linspace(-25, 0, 50)
    cmhk = axs[3].pcolormesh(x, y, np.flipud(np.log(hk)), cmap="coolwarm")
    axs[3].set_aspect(10)
    plt.colorbar(cmhk,ax=axs[3])

    plt.savefig(f"{results_location}\\pumping_and_recovery{modelname}{realization}.jpg", dpi=300)
    # plt.show()

def compare_steady_states(name, realizations, hks):
    nreal = len(realizations)
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{name}\{name}_steady.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{name}'

    f, axs = plt.subplots(3, nreal, sharex=True, sharey=True, figsize=(15, 5), constrained_layout=True)
    for i, (realization, hk) in enumerate(zip(realizations, hks)):
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
        x = np.linspace(0, 800, 80)
        y = np.linspace(-25, 0, 50)
        cmhk = axs[0][i].pcolormesh(x, y, np.flipud(np.log(hk)), cmap="coolwarm", vmax=-1, vmin=-13)
        axs[0][i].set_aspect(10)
        axs[1][i], conc = proc.plot_conc(axs[1][i], swt, {'qx':qx[:, :, :], 'qy':qy[:, :, :], 'qz':qz[:, :, :], 'concentration':concentration[:, :, :]}, row=0, vmax=35, vmin=0)
        axs[2][i], head = proc.plot_head(axs[2][i], swt, head_steady[:, :, :], vmin=-0.1, vmax=0.6);
        

    axs[0][0].set_title("Heterogenenous")
    axs[0][1].set_title("Homogeneous (equiv. head)")
    axs[0][2].set_title("Homogeneous (equiv. flux)")
    f.suptitle("Comparison of conductivity, concentration, and head", fontsize=16)

    cb1 = plt.colorbar(cmhk,ax=axs[0][:], location="right", shrink=0.75, pad=0.02)
    cb2 = plt.colorbar(conc,ax=axs[1][:], location="right", shrink=0.75, pad=0.02)
    cb3 = plt.colorbar(head,ax=axs[2][:], location="right", shrink=0.75, pad=0.02)
    cb1.ax.set_title('log(K)', fontsize = 'small')
    cb2.ax.set_title('C (kg/m^3)', fontsize = 'small')
    cb3.ax.set_title('h (m)', fontsize = 'small')
    cb1.ax.locator_params(nbins=6)
    cb2.ax.locator_params(nbins=4)
    cb3.ax.locator_params(nbins=4)
    axs[1][0].set_ylabel("Depth (m)")
    axs[2][1].set_xlabel("Distance (m)")
    f.set_constrained_layout_pads(w_pad=0.01, wspace=0.0)
    plt.show()

def compare_transient_response(name, realizations, hks, row, qlay, qcol):
    nreal = len(realizations)
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{name}\{name}_steady.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    results_location = f'.\\results\\{name}'

    f, axs = plt.subplots(3, nreal, sharex=True, sharey=True, figsize=(15, 5), constrained_layout=True)
    for i, (realization, hk) in enumerate(zip(realizations, hks)):
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
        qx_t, qy_t, qz_t, head_t, concentration_t =  proc.load_results_3D(name, realization, stress_period="transient")
        x = np.linspace(0, 800, 80)
        y = np.linspace(-25, 0, 50)
        cmhk = axs[0][i].pcolormesh(x, y, np.flipud(np.log(hk)), cmap="coolwarm", vmax=-1, vmin=-13)
        axs[0][i].set_aspect(10)
        axs[1][i], conc = proc.plot_conc(axs[1][i], swt, {'qx':qx[:, :, :], 'qy':qy[:, :, :], 'qz':qz[:, :, :], 'concentration':concentration[:, :, :]}, row=0, vmax=35, vmin=0)
        axs[2][i], conct = proc.plot_conc(axs[2][i], swt, {'qx':qx_t.item().get('pumping'), 'qy':qy_t.item().get('pumping'), 'qz':qz_t.item().get('pumping'), 'concentration':concentration_t.item().get('pumping')}, row=0, vmax=35, vmin=0)
        

    axs[0][0].set_title("Heterogenenous")
    axs[0][1].set_title("Homogeneous (equiv. head)")
    axs[0][2].set_title("Homogeneous (equiv. flux)")
    f.suptitle("Comparison of conductivity, steady concentration and transient concentration", fontsize=16)

    cb1 = plt.colorbar(cmhk,ax=axs[0][:], location="right", shrink=0.75, pad=0.02)
    cb2 = plt.colorbar(conc,ax=axs[1][:], location="right", shrink=0.75, pad=0.02)
    cb3 = plt.colorbar(conct,ax=axs[2][:], location="right", shrink=0.75, pad=0.02)
    cb1.ax.set_title('log(K)', fontsize = 'small')
    cb2.ax.set_title('C (kg/m^3)', fontsize = 'small')
    cb3.ax.set_title('C (kg/m^3)', fontsize = 'small')
    cb1.ax.locator_params(nbins=6)
    cb2.ax.locator_params(nbins=4)
    cb3.ax.locator_params(nbins=4)
    axs[1][0].set_ylabel("Depth (m)")
    axs[2][1].set_xlabel("Distance (m)")
    axs[1][0].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][0].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[1][1].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[1][2].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][1].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    axs[2][2].scatter(qcol*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
    f.set_constrained_layout_pads(w_pad=0.01, wspace=0.0)
    
    plt.savefig(f"{results_location}\\transient_comparison_row{row}.png", dpi=300)

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

def plot_metric_over_layers(modelname, metric, title):
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}_transient.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    delL = -swt.dis.botm[0,0,0]
    nlay = swt.dis.nlay
    (qlay, _, _, _) = swt.wel.stress_period_data['1'][-1]


    stoch = True
    real_color = "grey"
    try:
        rows = metric.shape[1]
    except:
        rows = 1

    layers=np.linspace(0, -nlay*delL, nlay)

    if rows == 1:
        real_color = "green"
        stoch = False

    f, ax = plt.subplots()
    ax.set_title(title)
    ax.set_ylabel("Depth (m)")
    for row in range(rows):
        ax.plot(metric[:,row], layers, color=real_color, lw=0.5)

    if stoch:
        average = np.average(metric, 1)
        ax.plot(average, layers, color="green", zorder=2)

    ax.axhline(y=-qlay*delL, color = "red", linestyle=":")

    plt.show()

def plot_heatmap(X, Y, heatmap, qlay):

    f, ax = plt.subplots()
    cm = ax.pcolormesh(X, Y, heatmap, cmap="hot_r")
    cb = plt.colorbar(cm)

    ax.axhline(y=qlay*(Y[1]-Y[0]), color = "blue", linestyle=":")
    plt.show()

def compute_steadystate_statistics(name, realizations, rows):

    n = len(rows)
    m = len(realizations)
    toe_position = np.zeros((m, n))
    mixing_area = np.zeros((m, n))
    centre_of_mass = np.zeros((m, n))
    fresh_sgd_flux = np.zeros((m, n))
    

    for i, realization in enumerate(realizations):
        for row in rows:
            realization_name = f"row{row}{realization}"
            qx, qy, qz, head, concentration = proc.load_results_3D(name, realization_name, stress_period="steady")
            toe_position[i, row] = proc.find_toe_penetration(concentration[:,0,:])
            mixing_area[i, row] = proc.find_mixing_zone(concentration[:,0,:])
            centre_of_mass[i, row] = proc.find_mixing_com(concentration[:,0,:], 800)
            fresh_sgd_flux[i, row] = proc.find_sgd(concentration[:,0,:], qz[:,0,:], 0.5)
 

    # for i in range(len(realizations)):
    #     for metric, name in zip([toe_position, mixing_area, centre_of_mass, fresh_sgd_flux], ['toe_position', 'mixing_area', 'centre_of_mass', 'fresh_sgd_flux']):
    #         print(f"mean of {name} for realization {realizations[i]}: {np.mean(metric[i, :])}")
    #         print(f"std of {name} for realization {realizations[i]}: {np.std(metric[i, :])}")

    for metric, name, unit in zip([toe_position, mixing_area, centre_of_mass, fresh_sgd_flux], \
                                ['toe_position', 'mixing_area', 'centre_of_mass', 'fresh_sgd_flux'],\
                                ['toe position (m)', 'mixing_area (m^2)', 'centre of mass position (m)', 'sgd flux (m^3/day)']):
        f, ax = plt.subplots()
        ax.boxplot([metric[0, :], metric[1, :], metric[2, :]], labels=["heterogenous", "equivalent flux", "equivalent head"])
        if name == 'centre_of_mass': ax.axhline(400, c='r', zorder = -1, linestyle=':')
        ax.set_title(name)
        plt.show()

    # f, ax = plt.subplots()
    # colors = ['r', 'b', 'g']
    # labels=["heterogenous", "equivalent flux", "equivalent head"]
    # for i in range(len(realizations)):
    #     ax.scatter(centre_of_mass[i,:], fresh_sgd_flux[i,:], c=colors[i], label=labels[i])

    # ax.set_title("Fresh SGD vs COM of mixing zone")
    # ax.set_ylabel("sgd (m^3/day)")
    # ax.set_xlabel("com position (m)")
    # ax.legend()
    # plt.show()

def transient_box_plots(name, realizations, rows, qcol, qlay):
    n = len(rows)
    m = len(realizations)
    toe_position = np.zeros((m, n))
    mixing_area = np.zeros((m, n))
    centre_of_mass = np.zeros((m, n))
    fresh_sgd_flux = np.zeros((m, n))
    fresh = np.zeros((m, n))
    wel = np.zeros((m, n))

    for i, realization in enumerate(realizations):
        for row in rows:
            realization_name = f"row{row}{realization}"
            qx, qy, qz, head, concentration = proc.load_results_3D(name, realization_name, stress_period="transient")
            toe_position[i, row] = proc.find_toe_penetration(concentration.item().get('pumping')[:,0,:])
            mixing_area[i, row] = proc.find_mixing_zone(concentration.item().get('pumping')[:,0,:])
            centre_of_mass[i, row] = proc.find_mixing_com(concentration.item().get('pumping')[:,0,:], 800)
            fresh_sgd_flux[i, row] = proc.find_sgd(concentration.item().get('pumping')[:,0,:], qz.item().get('pumping')[:,0,:], 0.5)
            fresh[i, row] = proc.find_fresh_water_volume(concentration.item().get('pumping')[:,0,:], 5, qcol)
            wel[i, row] = concentration.item().get('pumping')[qlay,0,qcol]

    # for i in range(len(realizations)):
    #     for metric, name in zip([toe_position, mixing_area, centre_of_mass, fresh_sgd_flux], ['toe_position', 'mixing_area', 'centre_of_mass', 'fresh_sgd_flux']):
    #         print(f"mean of {name} for realization {realizations[i]}: {np.mean(metric[i, :])}")
    #         print(f"std of {name} for realization {realizations[i]}: {np.std(metric[i, :])}")

    for metric, name, unit in zip([toe_position, mixing_area, centre_of_mass, fresh_sgd_flux, fresh, wel], \
                                ['toe position', 'mixing area', 'mixing centre of mass position', 'fresh sgd flux', "freshwater volume", "well salinity"],\
                                ['distance inland (m)', 'area (m^2)', 'distance offshore (m)', 'flux (m^3/day)', 'area (m^2)', "salinity (kg/m^3)"]):
        f, ax = plt.subplots()
        ax.boxplot([metric[0, :], metric[1, :], metric[2, :]], labels=["heterogenous", "equivalent flux", "equivalent head"])
        ax.set_title(name)
        ax.set_ylabel(unit)
        plt.show()

def plot_steady_time_evolutions(name, realizations):

    com, toe, mix = proc.load_steady_metric_evolutions(name)
    com_mean = np.mean(com, axis=1)
    com_lower = np.percentile(com, 25, axis = 1)
    com_upper = np.percentile(com, 75, axis = 1)
    toe_mean = np.mean(toe, axis=1)
    toe_lower = np.percentile(toe, 25, axis = 1)
    toe_upper = np.percentile(toe, 75, axis = 1)
    mix_mean = np.mean(mix, axis=1)
    mix_lower = np.percentile(mix, 25, axis = 1)
    mix_upper = np.percentile(mix, 75, axis = 1)
    
    c = ["r", "b", "g"]
    f, axs = plt.subplots(3, sharex=True, constrained_layout=True)
    t = 1/365*1e6*np.logspace(0, 3.56, base=10, num=30).astype(int)
    axs[0].set_xscale('log')
    axs[1].set_xscale('log')
    axs[2].set_xscale('log')

    
    for i in range(len(realizations)):
        axs[0].plot(t, com_mean[i], c=c[i])
        axs[1].plot(t, toe_mean[i], c=c[i], label= realizations[i])
        axs[2].plot(t, mix_mean[i], c=c[i])
        axs[0].fill_between(t, com_lower[i], com_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[1].fill_between(t, toe_lower[i], toe_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[2].fill_between(t, mix_lower[i], mix_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")

    axs[2].set_xlabel("log(years)")
    axs[0].set_ylabel("distance offshore (m)")
    axs[1].set_ylabel("distance onshore (m)")
    axs[2].set_ylabel("area (m^2)")

    axs[0].set_title("Centre of mass of mixing zone")
    axs[1].set_title("Toe position")
    axs[2].set_title("Area of mixing zone")

    axs[1].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.show()
    # axs[1].
    # axs[2].
    pass

def plot_transient_time_evolutions(name, realizations):

    com, toe, mix, fresh, wel = proc.load_transient_metric_evolutions(name)
    com_mean = np.mean(com, axis=1)
    com_lower = np.percentile(com, 25, axis = 1)
    com_upper = np.percentile(com, 75, axis = 1)
    toe_mean = np.mean(toe, axis=1)
    toe_lower = np.percentile(toe, 25, axis = 1)
    toe_upper = np.percentile(toe, 75, axis = 1)
    mix_mean = np.mean(mix, axis=1)
    mix_lower = np.percentile(mix, 25, axis = 1)
    mix_upper = np.percentile(mix, 75, axis = 1)
    fresh_mean = np.mean(fresh, axis=1)
    fresh_lower = np.percentile(fresh, 25, axis = 1)
    fresh_upper = np.percentile(fresh, 75, axis = 1)
    wel_mean = np.mean(wel, axis=1)
    wel_lower = np.percentile(wel, 25, axis = 1)
    wel_upper = np.percentile(wel, 75, axis = 1)
    
    c = ["r", "b", "g"]
    f, axs = plt.subplots(5, sharex=True, constrained_layout=True)
    t = np.linspace(0, 1000, 50).astype(int)*3650/2/1000

    
    for i in range(len(realizations)):
        axs[0].plot(t, com_mean[i], c=c[i])
        axs[1].plot(t, toe_mean[i], c=c[i])
        axs[2].plot(t, mix_mean[i], c=c[i], label= realizations[i])
        axs[3].plot(t, fresh_mean[i], c=c[i])
        axs[4].plot(t, wel_mean[i], c=c[i])
        
        axs[0].fill_between(t, com_lower[i], com_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[1].fill_between(t, toe_lower[i], toe_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[2].fill_between(t, mix_lower[i], mix_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[3].fill_between(t, fresh_lower[i], fresh_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")
        axs[4].fill_between(t, wel_lower[i], wel_upper[i], alpha=0.1, color=c[i], label="Q25-Q75")

    axs[4].set_xlabel("time (days)")
    axs[0].set_ylabel("distance offshore (m)")
    axs[1].set_ylabel("distance onshore (m)")
    axs[2].set_ylabel("area (m^2)")
    axs[3].set_ylabel("area (m^2)")
    axs[4].set_ylabel("salinity (kg/m^3)")

    axs[0].set_title("Centre of mass of mixing zone")
    axs[1].set_title("Toe position")
    axs[2].set_title("Area of mixing zone")
    axs[3].set_title("Area of freshwater")
    axs[4].set_title("Well salinity")

    axs[2].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.show()

def plot_well_ensemble(name, realizations):

    _, _, _, _, wel = proc.load_transient_metric_evolutions(name)
    wel_mean = np.mean(wel, axis=1)

    c = ["r", "b", "g"]
    f, axs = plt.subplots(3, sharex=True, constrained_layout=True)
    t = np.linspace(0, 1000, 50).astype(int)*3650/2/1000

    for i in range(len(realizations)):
        axs[i].plot(t, wel_mean[i], c=c[i], lw=2, label=f'{realizations[i]} mean')
        axs[i].plot(t, wel[i].T, c=c[i], alpha = 0.2, lw=0.5)

    axs[1].set_ylabel("salinity (kg/m^3)")
    axs[2].set_xlabel("time (days)")
    axs[0].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    axs[1].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    axs[2].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.suptitle("Well salinity under constant pumping")
    plt.show()
        

if __name__=="__main__":
    results("pirot2D_check_k_eff")