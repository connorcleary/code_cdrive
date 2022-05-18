from ctypes import c_ssize_t
from lib2to3.pgen2.pgen import DFAState
import matplotlib.pyplot as plt
import numpy as np 
import flopy
from pyrsistent import v
from zmq import TCP_MAXRT
import post_proc_utils as proc
import matplotlib.cm as cm
import plot_helpers as plth
from matplotlib.gridspec import GridSpec

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

def compare_steady_states_2D_3D(name2D, name3D, realizations, hk, row):

    f, axs = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(6, 7), constrained_layout=True)

    x = np.linspace(-400, 400, 80)
    y = np.linspace(-25, 0, 50)

    qx2D, _, qz2D, _, concentration2D =  proc.load_results_3D(name2D, f"row{row}{realizations[0]}", stress_period="steady")
    qx3D, _, qz3D, _, concentration3D =  proc.load_results_3D(name3D, realizations[1], stress_period="steady")

    cmhk = axs[0].pcolormesh(x, y, np.flipud(np.log(hk[:,row,:])), cmap="coolwarm", vmax=-1, vmin=-13, )
    cm2D = axs[1].pcolormesh(x, y, np.flipud(concentration2D[:,0,:]), cmap="viridis", vmax=35, vmin=0)
    cm3D = axs[2].pcolormesh(x, y, np.flipud(concentration3D[:,row,:]), cmap="viridis", vmax=35, vmin=0)
    axs[1].quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                    plth.sample_grid_data(np.flipud(qx2D[:,0,:]), 3,3), -plth.sample_grid_data(np.flipud(qz2D[:,0,:]), 3,3), color="white")
    axs[2].quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                    plth.sample_grid_data(np.flipud(qx3D[:,row,:]), 3,3), -plth.sample_grid_data(np.flipud(qz3D[:,row,:]), 3,3), color="white")

    axs[0].set_title("Horizontal hydraulic conductivity")
    axs[1].set_title("2D aquifer model")
    axs[2].set_title("3D aquifer model")
    f.suptitle(f"Concentration in row {row}", fontsize=16)
    axs[1].set_ylabel("Depth (m)")
    axs[2].set_xlabel("Distance offshore (m)")

    cb1 = plt.colorbar(cmhk,ax=axs[0], location="right", shrink=0.75, pad=0.02)
    cb2 = plt.colorbar(cm2D,ax=axs[1:], location="right", shrink=0.75/2, pad=0.02)

    cb1.ax.set_title('log10[hk]', fontsize = 'small')
    cb2.ax.set_title('C (kg/m^3)', fontsize = 'small')

    for ax in axs[:]: ax.set_aspect(10)

    #plt.show()

    results_location = f'.\\results\\{name3D}'
    plt.savefig(f"{results_location}\\steady2Dvs3Drow{row}", dpi=300)

def steady_states_gridspec_evolution3D(name3D, realization, hk, row):
    
    f = plt.figure(constrained_layout=True, figsize=(7,7))
    gs = GridSpec(4, 2, figure=f, width_ratios=[1,2])

    axleg = f.add_subplot(gs[0, 0])
    axtime = f.add_subplot(gs[1:, 0])
    axhk = f.add_subplot(gs[0, 1])
    axmax = f.add_subplot(gs[1, 1])
    axmin = f.add_subplot(gs[2, 1])
    axss = f.add_subplot(gs[3, 1])

    axhk.set_aspect(10)
    axmax.set_aspect(10)
    axmin.set_aspect(10)
    axss.set_aspect(10)

    t = 1/365*1e6*np.logspace(0, 3.56, base=10, num=30).astype(int)
    # load time evolutions
    com, toe, mix = proc.load_steady_metric_evolutions(name3D)
    concentration, qx, qz = proc.load_time_evolution_3D(name3D, realization,"steady")

    imax = np.atleast_1d(np.argmax(com[0,row,:]))[0]
    imin = np.atleast_1d(np.argmin(com[0,row,imax:]))[0] + imax
    iss = com.shape[-1] -1
    for i in range(com[0, row, imin:].shape[0]):
        if np.average(np.absolute(concentration[i+imin,:,:,:]-concentration[i+imin-1,:,:,:])) < 0.05:
            iss = i + imin
            break
    imax = 13


    nmax = int(t[imax]//(1e6/365))
    nmin = int(t[imin]//(1e6/365))
    nss = int(t[iss]//(1e6/365))

    cmax = concentration[nmax,:,:,:]
    cmin = concentration[nmin,:,:,:] 
    css = concentration[nss,:,:,:]
    qxmax = qx[nmax,:,:,:]
    qxmin = qx[nmin,:,:,:] 
    qxss  = qx[nss,:,:,:]
    qzmax  = qz[nmax,:,:,:]
    qzmin  = qz[nmin,:,:,:]
    qzss = qz[nss,:,:,:]

    x = np.linspace(-400, 400, 80)
    y = np.linspace(-25, 0, 50)

    cmhk = axhk.pcolormesh(x, y, np.flipud(np.log(hk[:,row,:])), cmap="coolwarm", vmax=-1, vmin=-13)
    cmmax = axmax.pcolormesh(x, y, np.flipud(cmax[:,row,:]), cmap="viridis", vmax=35, vmin=0)
    cmmin = axmin.pcolormesh(x, y, np.flipud(cmin[:,row,:]), cmap="viridis", vmax=35, vmin=0)
    cmss = axss.pcolormesh(x, y, np.flipud(css[:,row,:]), cmap="viridis", vmax=35, vmin=0)

    axmax.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                    plth.sample_grid_data(np.flipud(qxmax[:,0,:]), 3,3), -plth.sample_grid_data(np.flipud(qzmax[:,0,:]), 3,3), color="white")
    axmin.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                    plth.sample_grid_data(np.flipud(qxmin[:,row,:]), 3,3), -plth.sample_grid_data(np.flipud(qzmin[:,row,:]), 3,3), color="white")
    axss.quiver(plth.sample_grid_data(x, 3, 3), plth.sample_grid_data(y, 3, 3), 
                    plth.sample_grid_data(np.flipud(qxss[:,row,:]), 3,3), -plth.sample_grid_data(np.flipud(qzss[:,row,:]), 3,3), color="white")

    
    axtime.invert_yaxis()
    axtime.plot(com[0,row,:], t)
    axtime.axhline(y=t[imax], color = "red", linestyle=(0,(1,2)), label="maximum")
    axtime.axhline(y=t[imin], color = "blue", linestyle=(1,(1,2)), label="minimum")
    axtime.axhline(y=t[iss], color = "green", linestyle=(2,(1,2)), label="steady")

    axtime.set_yscale('log')

    axhk.set_title("Horizontal hydraulic conductivity")
    axmax.set_title("Maximum distance offshore")
    axmin.set_title("Minimum distance offshore")
    axss.set_title("Steady state")
    axtime.set_title("Time evolution")

    axmin.set_ylabel("Depth (m)")
    axss.set_xlabel("Distance offshore (m)")
    axtime.set_xlabel("Distance offshore (m)")
    axtime.set_ylabel("time (log[years])", labelpad=10, rotation=270)

    h, l = axtime.get_legend_handles_labels()
    axleg.axis("off")
    axleg.legend(h, l, loc="lower center", borderaxespad=0)

    f.suptitle("Saline concentrations at different states")

    cb1 = plt.colorbar(cmhk,ax=axhk, location="right", shrink=0.75, pad=0.05)
    cb2 = plt.colorbar(cmmax,ax=[axmax, axmin, axss], location="right", shrink=0.75/3, pad=0.05)

    cb1.ax.set_title('log10[hk]', fontsize = 'small')
    cb2.ax.set_title('C (kg/m^3)', fontsize = 'small')

    results_location = f'.\\results\\{name3D}'
    plt.savefig(f"{results_location}\\steady3D_evolution{row}", dpi=300)
    

def probability_of_saline(modelname):
    """
        For 3D
    """
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
    plt.show()
    
    pass

def plot_metric_over_layers(modelname, metric, realizations, title):
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

    f, axs = plt.subplots(1, len(realizations), constrained_layout=True, sharey=True)
    f.suptitle(title)
    axs[0].set_ylabel("Depth (m)")
    for i, realization in enumerate(realizations):
        metric_real = metric[:,:,i]
        axs[i].set_xlabel("movement inland (m)")
        for row in range(rows):
            axs[i].plot(metric_real[:,row], layers, color=real_color, lw=0.5, label="_nolegend_")
        axs[i].plot(metric_real[:,0], layers, color=real_color, lw=0.5, label="single row")

        if stoch:
            average = np.average(metric_real, 1)
            axs[i].plot(average, layers, color="green", zorder=2, label="average")

        axs[i].axhline(y=-qlay*delL, color = "red", linestyle=":", label="well layer")
        axs[i].set_title(realization)

    
    axs[i].legend(loc='center left')
    plt.show()

def plot_metric_as_heatmap_layered(modelname, metric, realizations, title):
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

    f, axs = plt.subplots(1, len(realizations), constrained_layout=True, sharey=True)
    f.suptitle("Probability of lateral movement of mixing zone")
    axs[0].set_ylabel("Depth (m)")
    metric_real0 = metric[0,:,:]
    range_dist = range(int(np.min(metric_real0)), int(np.max(metric_real0)+10), 10)
    for i, realization in enumerate(realizations):
        metric_real = metric[i,:,:]
        probability = np.zeros((nlay, len(range_dist)))
        axs[i].set_xlabel("movement inland (m)")
        for row in range(rows):
            for lay in range(nlay):
                for j, distance in enumerate(range_dist):
                    if (metric_real[row, lay] >= distance and distance >= 0) or (metric_real[row, lay] <= distance and distance < 0):
                        probability[lay, j] += 1

        prob = axs[i].pcolor(range_dist, layers, probability/rows, cmap="coolwarm")
    # ax.plot(X, Y, homogenous_interface, homogenous_y, color='black')
        axs[i].set_title(realization)

        axs[i].axhline(y=-qlay*delL, color = "red", linestyle=":", label="well layer")

    cb = plt.colorbar(prob, ax=axs[:], location="bottom", shrink=0.2)
    axs[i].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.show()

def plot_heatmaps(X, Y, modelname, heatmaps, qlay, realizations):

    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{modelname}\{modelname}_transient.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    delL = -swt.dis.botm[0,0,0]
    nlay = swt.dis.nlay
    (qlay, _, _, _) = swt.wel.stress_period_data['1'][-1]

    f, axs = plt.subplots(1, len(realizations), constrained_layout=True, sharey=True)
    f.suptitle("Probable area salinized")
    axs[0].set_ylabel("Depth (m)")

    for i, realization in enumerate(realizations):
        axs[i].set_xlabel("x (m)")
        axs[i].set_aspect(10)
        cm = axs[i].pcolormesh(X, Y, heatmaps[i], cmap="coolwarm")
        axs[i].set_title(realization)
        axs[i].scatter(15*10+5, qlay*(-0.5) - 0.25, c='red', edgecolors='black', zorder=2)
        # axs[i].axhline(y=qlay*(Y[1]-Y[0]), color = "red", linestyle=":", label="well layer")

    cb = plt.colorbar(cm, ax=axs[:], location="bottom", shrink=0.2)
    axs[i].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
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
        if name == 'centre_of_mass': ax.axhline(0, c='r', zorder = -1, linestyle=':')
        ax.set_title(name)
        plt.show()

    f, ax = plt.subplots()
    colors = ['r', 'b', 'g']
    labels=["heterogenous", "equivalent flux", "equivalent head"]
    for i in range(len(realizations)):
        ax.scatter(centre_of_mass[i,:], fresh_sgd_flux[i,:], c=colors[i], label=labels[i])

    ax.set_title("Fresh SGD vs COM of mixing zone")
    ax.set_ylabel("sgd (m^3/day)")
    ax.set_xlabel("com position (m)")
    ax.legend()
    plt.show()

def steady_boxplots_2D_vs_3D(name2D, name3D, realizations, rows):

    if type(rows) is int:
        n = rows
    else:
        n = len(rows)

    m = len(realizations)
    toe_position = np.zeros((2, m, n))
    mixing_area = np.zeros((2, m, n))
    centre_of_mass = np.zeros((2, m, n))
    fresh_sgd_flux = np.zeros((2, m, n))

    for i, model in enumerate([name2D, name3D]):
        for j, realization in enumerate(realizations):
            if i == 0:
                realization_name = f"pirot_basic{realization}"
                _, _, qz, _, concentration = proc.load_results_3D(model, realization_name, stress_period="steady")
            for k in range(rows):
                if i == 1:
                    realization_name = f"row{k}{realization}"
                    row = 0
                    
                    _, _, qz, _, concentration = proc.load_results_3D(model, realization_name, stress_period="steady")
                else:
                    row = k 
                
                toe_position[i, j, k] = proc.find_toe_penetration(concentration[:,row,:])
                mixing_area[i, j, k] = proc.find_mixing_zone(concentration[:,row,:])
                centre_of_mass[i, j, k] = proc.find_mixing_com(concentration[:,row,:], 800)
                fresh_sgd_flux[i, j, k] = proc.find_sgd(concentration[:,row,:], qz[:,row,:], 0.5)
            
    f, axs = plt.subplots(4, 1, sharex=True, constrained_layout=True)

    for i, (metric, name, unit) in enumerate(zip([toe_position, mixing_area, centre_of_mass, fresh_sgd_flux], \
                                ['toe_position', 'mixing_area', 'centre_of_mass_position', 'fresh_sgd_flux'],\
                                ['distance onshore (m)', 'area (m^2)', 'distance offshore (m)', 'flux (m^3/day)'])):

        bp3D = axs[i].boxplot([metric[0, 0, :], metric[0, 1, :], metric[0, 2, :]], positions=np.array(range(len(realizations)))*2.0-0.4, sym='', widths=0.6)
        bp2D = axs[i].boxplot([metric[1, 0, :], metric[1, 1, :], metric[1, 2, :]], positions=np.array(range(len(realizations)))*2.0+0.4, sym='', widths=0.6)
        plth.set_box_color(bp3D, "b")
        plth.set_box_color(bp2D, "r")

        axs[i].set_title(name)
        axs[i].set_ylabel(unit)

    realizations[0] = "heterogenous"

    axs[i].set_xticks(range(0, len(realizations) * 2, 2), realizations)
    axs[i].set_xlim(-2, len(realizations)*2)

    axs[1].plot([], c='b', label='3D model')
    axs[1].plot([], c='r', label='2D ensemble')
    axs[1].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')

    plt.show()

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

def plot_steady_time_evolution_ensemble(name, realizations="heterogenous", hk=None):

    # colors = plt.cm.jet(np.linspace(0,1,40))
    # colors = np.concatenate((colors, colors, colors, colors))
    norm = plt.Normalize()
    colors = plt.cm.coolwarm(norm(np.log10(hk)))

    com, toe, mix = proc.load_steady_metric_evolutions(name)
    com_mean = np.mean(com, axis=1)
    com_median = np.median(com, axis=1)
    com_lower = np.percentile(com, 25, axis = 1)
    com_upper = np.percentile(com, 75, axis = 1)
    toe_mean = np.mean(toe, axis=1)
    toe_median = np.median(toe, axis=1)
    toe_lower = np.percentile(toe, 25, axis = 1)
    toe_upper = np.percentile(toe, 75, axis = 1)
    mix_mean = np.mean(mix, axis=1)
    mix_median = np.median(mix, axis=1)
    mix_lower = np.percentile(mix, 25, axis = 1)
    mix_upper = np.percentile(mix, 75, axis = 1)
    
    c = ["r", "b", "g"]
    f, axs = plt.subplots(3, len(realizations), sharex=True, sharey="row", constrained_layout=True)
    t = 1/365*1e6*np.logspace(0, 3.56, base=10, num=30).astype(int)

    commin = np.min(com.flatten()) -50
    commax = np.max(com.flatten()) +50
    toemin = np.min(toe.flatten()) -50
    toemax = np.max(toe.flatten()) +50
    mixmin = np.min(mix.flatten()) -50
    mixmax = np.max(mix.flatten()) +50

    for real in range(len(realizations)):

        axs[0][real].set_xscale('log')
        axs[1][real].set_xscale('log')
        axs[2][real].set_xscale('log')

        axs[0][real].set_ylim(commin, commax)
        axs[1][real].set_ylim(toemin, toemax)
        axs[2][real].set_ylim(mixmin, mixmax)

        i = 0
        for j in range(40):
            axs[0][real].plot(t, com[real,j,:].T, c=colors[j], alpha=1, linewidth=0.35)
            axs[1][real].plot(t, toe[real,j,:], c=colors[j], label="_nolegend_", alpha=1, linewidth=0.35)
            axs[2][real].plot(t, mix[real,j,:], c=colors[j],label="_nolegend_", alpha=1, linewidth=0.35)

        axs[0][real].plot(t, com_mean[real], c=c[2])
        axs[1][real].plot(t, toe_mean[real], c=c[2], label= f"{realizations[real]} mean")
        axs[2][real].plot(t, mix_mean[real], c=c[2], label= f"ensemble mean")
        axs[0][real].plot(t, com_median[real], c=c[2], linestyle=":")
        axs[1][real].plot(t, toe_median[real], c=c[2], linestyle=":", label= f"{realizations[real]} median")
        axs[2][real].plot(t, mix_median[real], c=c[2], linestyle=":", label= f"ensemble median")

        axs[2][real].set_xlabel("log(years)")

    axs[0][0].set_ylabel("distance offshore (m)")
    axs[1][0].set_ylabel("distance onshore (m)")
    axs[2][0].set_ylabel("area (m^2)")

    # axs[0][1].set_title("Centre of mass of mixing zone")
    # axs[1][1].set_title("Toe position")
    # axs[2][1].set_title("Area of mixing zone")

    axs[0][0].set_title("Heterogenous")
    axs[0][1].set_title("Equivalent flux")
    axs[0][2].set_title("Equivalent head")

    for ax, row in zip(axs[:,0], ["Centre of mass of mixing zone", "Toe position", "Area of mixing zone"]):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center',rotation=90)

    axs[2][1].legend(bbox_to_anchor=(0.5, -0.3), loc='upper center', ncol=3)

    cb = f.colorbar(cm.ScalarMappable(norm, cmap="coolwarm"), ax=axs[:,2], shrink=0.3, extend="both")
    cb.ax.set_title('log10[Kh_eff]', fontsize = 'small', pad=10)
    plt.show()
    # axs[1].
    # axs[2].
    pass

def plot_steady_time_evolution_ensemble_2D_vs_3D(name2D, name3D, realizations="heterogenous", hk=None):

    # colors = plt.cm.jet(np.linspace(0,1,40))
    # colors = np.concatenate((colors, colors, colors, colors))
    norm = plt.Normalize()
    colors = plt.cm.coolwarm(norm(np.log10(hk)))

    com, toe, mix = proc.load_steady_metric_evolutions(name2D)
    com_3D, toe_3D, mix_3D = proc.load_steady_metric_evolutions(name3D)
    com = np.concatenate([[com[0,:,:]], com_3D])
    toe = np.concatenate([[toe[0,:,:]], toe_3D])
    mix = np.concatenate([[mix[0,:,:]], mix_3D])

    com_mean = np.mean(com, axis=1)
    com_median = np.median(com, axis=1)
    com_lower = np.percentile(com, 25, axis = 1)
    com_upper = np.percentile(com, 75, axis = 1)
    toe_mean = np.mean(toe, axis=1)
    toe_median = np.median(toe, axis=1)
    toe_lower = np.percentile(toe, 25, axis = 1)
    toe_upper = np.percentile(toe, 75, axis = 1)
    mix_mean = np.mean(mix, axis=1)
    mix_median = np.median(mix, axis=1)
    mix_lower = np.percentile(mix, 25, axis = 1)
    mix_upper = np.percentile(mix, 75, axis = 1)
    
    c = ["r", "b", "g"]
    f, axs = plt.subplots(3, len(realizations), sharex=True, sharey="row", constrained_layout=True)
    t = 1/365*1e6*np.logspace(0, 3.56, base=10, num=30).astype(int)

    commin = np.min(com_3D.flatten()) -10
    commax = np.max(com_3D.flatten()) +10
    toemin = np.min(toe.flatten()) -50
    toemax = np.max(toe.flatten()) +50
    mixmin = np.min(mix.flatten()) -50
    mixmax = np.max(mix.flatten()) +50

    for real in range(len(realizations)):

        axs[0][real].set_xscale('log')
        axs[1][real].set_xscale('log')
        axs[2][real].set_xscale('log')

        axs[0][real].set_ylim(commin, commax)
        axs[1][real].set_ylim(toemin, toemax)
        axs[2][real].set_ylim(mixmin, mixmax)

        i = 0
        for j in range(40):
            axs[0][real].plot(t, com[real,j,:].T, c=colors[j], alpha=1, linewidth=0.35)
            axs[1][real].plot(t, toe[real,j,:].T, c=colors[j], label="_nolegend_", alpha=1, linewidth=0.35)
            axs[2][real].plot(t, mix[real,j,:].T, c=colors[j],label="_nolegend_", alpha=1, linewidth=0.35)

        axs[0][real].plot(t, com_mean[real], c=c[2])
        axs[1][real].plot(t, toe_mean[real], c=c[2], label= f"{realizations[real]} mean")
        axs[2][real].plot(t, mix_mean[real], c=c[2], label= f"ensemble mean")
        axs[0][real].plot(t, com_median[real], c=c[2], linestyle=":")
        axs[1][real].plot(t, toe_median[real], c=c[2], linestyle=":", label= f"{realizations[real]} median")
        axs[2][real].plot(t, mix_median[real], c=c[2], linestyle=":", label= f"ensemble median")

        axs[2][real].set_xlabel("log(years)")

    axs[0][0].set_ylabel("distance offshore (m)")
    axs[1][0].set_ylabel("distance onshore (m)")
    axs[2][0].set_ylabel("area (m^2)")

    # axs[0][1].set_title("Centre of mass of mixing zone")
    # axs[1][1].set_title("Toe position")
    # axs[2][1].set_title("Area of mixing zone")

    axs[0][0].set_title("2D")
    axs[0][1].set_title("3D")

    for ax, row in zip(axs[:,0], ["Centre of mass of mixing zone", "Toe position", "Area of mixing zone"]):
        ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - 5, 0),
                xycoords=ax.yaxis.label, textcoords='offset points',
                size='large', ha='right', va='center',rotation=90)

    axs[2][1].legend(bbox_to_anchor=(0.5, -0.3), loc='upper center', ncol=3)

    cb = f.colorbar(cm.ScalarMappable(norm, cmap="coolwarm"), ax=axs[:,1], shrink=0.3, extend="both")
    cb.ax.set_title('log10[Kh_eff]', fontsize = 'small', pad=10)
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
        axs[i].axhline(0.35, c='k', alpha=0.5, zorder = -1, linestyle=':', label="Potable maximum")

    axs[1].set_ylabel("salinity (kg/m^3)")
    axs[2].set_xlabel("time (days)")
    axs[0].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    axs[1].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    axs[2].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    plt.suptitle("Well salinity under constant pumping")
    plt.show()
        
def plot_effective_conductivities(Kh, Kv):

    f, axs = plt.subplots(1, 3, constrained_layout=True, sharey=True)
    axs[0].boxplot([Kh, Kv], labels=["Kh_eff", "Kv_eff"], vert=True)
    axs[0].set_yscale("log")
    axs[1].plot(range(len(Kh)), np.array([Kh, Kv]).T)
    axs[1].set_yscale("log")
    axs[1].legend(labels=["Kh_eff", "Kv_eff"])
    f.suptitle("Effective conductivities across all rows")
    axs[0].set_ylabel("Conductivity (log[m/day])")
    axs[1].set_xlabel("Row number")
    axs[2].set_yscale("log")
    scat = axs[2].scatter(Kv, Kh, c=range(len(Kh)))
    scatcb = plt.colorbar(scat, ax=axs[2])
    scatcb.ax.set_title("row")
    axs[2].set_xlabel("Kv_eff")

    plt.show()

def compare_wel_salinity_with_effective_conductivity(Kh, Kv, wel, realizations=["heterogeneous"]):
    
    f, axs = plt.subplots(1,2, constrained_layout=True, sharey=True)
    colors = ["r", "b", "g"]

    for i in range(len(realizations)):
        axs[0].scatter(Kh, wel[i, :], color=colors[i])
        axs[1].scatter(Kv, wel[i, :], color=colors[i], label=realizations[i])

    axs[0].set_xscale("log")
    axs[1].set_xscale("log")
    axs[0].set_yscale("log")
    axs[1].set_yscale("log")
    axs[0].set_ylabel("salinity (log[kg/m^3])")
    axs[0].set_xlabel("Kh_eff (log[m/day])")
    axs[1].set_xlabel("Kv_eff (log[m/day])")
    axs[1].legend(bbox_to_anchor=(1.04, 0.5), loc='center left')
    f.suptitle("Effective conducitivities vs well salinity")

    plt.show()

if __name__=="__main__":
    results("pirot2D_check_k_eff")