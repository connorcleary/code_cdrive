import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
import os

def find_mixing_zone(conc):
    
    mix = 0
    for i in range(conc.shape[0]):
        for j in range(conc.shape[1]):
            if 3.5 <= conc[i, j] <= 31.5:
                mix += 1
    # need to input array shape as well
    mixing_zone = 5*mix

    return mixing_zone


def find_toe_penetration(conc):
    for i in range(conc.shape[1]):
        if conc[-1][i] > 0.5:
            return (conc.shape[1]/4 - i)*10

def find_sgd(names):
    return sgd_fresh, sgd_saline

def plot_conc(ax, swt, results_dict):

    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(results_dict['concentration'])
    pmv.plot_vector(results_dict['qx'], results_dict['qz'], -results_dict['qz'], color="white", kstep=3, hstep=3)
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Concentrations")


def plot_head(ax, swt, results_dict):

    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(results_dict['head'])
    contours = pmv.contour_array(results_dict['head'], colors="white")
    ax.clabel(contours, fmt="%2.2f")
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Heads")


def extract_results(swt):

    shape = swt.lpf.hk.shape
    if len(shape) == 2:
        nrow = 1
        nlay, ncol = swt.lpf.hk.shape
    else:
        nlay, nrow, ncol = swt.lpf.hk.shape

    ws = swt._model_ws
    ucnobj = bf.UcnFile(os.path.join(ws, "MT3D001.UCN"), model=swt)
    times = ucnobj.get_times()
    concentration = ucnobj.get_data(totim=times[-1])

    cbbobj = bf.CellBudgetFile(os.path.join(ws, f'{swt._BaseModel__name}.cbc'))
    qx = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
    try:
        qy = cbbobj.get_data(text="flow front face", totim=times[-1])[0]
    except:
        qy = np.zeros((nlay, nrow, ncol), dtype=float)
    qz = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]

    headobj = bf.HeadFile(os.path.join(ws, f'{swt._BaseModel__name}.hds'))
    head = headobj.get_data(totim=times[-1])

    return concentration, qx, qy, qz, head

def save_results(swt, realization, concentration, qx, qy, qz, head):

    ws = os.path.join(f'.\\results\\{swt._BaseModel__name}')
    if not os.path.exists(ws):
        os.mkdir(ws)

    np.savetxt(os.path.join(ws, f"qx_{realization}"), qx[:,0,:])
    np.savetxt(os.path.join(ws, f"qy_{realization}"), qy[:,0,:])
    np.savetxt(os.path.join(ws, f"qz_{realization}"), qz[:,0,:])
    np.savetxt(os.path.join(ws, f"head_{realization}"), head[:,0,:])
    np.savetxt(os.path.join(ws, f"concentration_{realization}"), concentration[:,0,:])


