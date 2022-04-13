import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
import os
import post_proc_array as proc_functions

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

def plot_conc(ax, swt, results_dict, row=0, **kwargs):

    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": row})
    arr = pmv.plot_array(results_dict['concentration'], **kwargs)
    pmv.plot_vector(results_dict['qx'], -results_dict['qz'], -results_dict['qz'], kstep=3, hstep=3, normalize=False, color="white")
    plt.colorbar(arr, shrink=0.5, ax=ax)

    ax.set_title("Simulated Concentrations")

    return ax


def plot_head(ax, swt, results_dict, row=0, **kwargs):

    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": row})

    if type(results_dict) is dict:
        arr = pmv.plot_array(results_dict['head'], **kwargs)
        contours = pmv.contour_array(results_dict['head'], colors="white")
    else: 
        contours = pmv.contour_array(results_dict, colors="white")
        arr = pmv.plot_array(results_dict, **kwargs)
        

    ax.clabel(contours, fmt="%2.2f")
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Heads")

    return ax

def extract_results(swt, nstp, pump=False, rec=False):

    shape = swt.lpf.hk.shape
    if len(shape) == 2:
        nrow = 1
        nlay, ncol = swt.lpf.hk.shape
    else:
        nlay, nrow, ncol = swt.lpf.hk.shape

    ws = swt._model_ws
    ucnobj = bf.UcnFile(os.path.join(ws, "MT3D001.UCN"), model=swt)
    times = ucnobj.get_times()
    cbbobj = bf.CellBudgetFile(os.path.join(ws, f'{swt._BaseModel__name}.cbc'))
    headobj = bf.HeadFile(os.path.join(ws, f'{swt._BaseModel__name}.hds'))
    
    if not pump:
        concentration = ucnobj.get_data(totim=times[-1])
        qx = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
        try:
            qy = cbbobj.get_data(text="flow front face", totim=times[-1])[0]
        except:
            qy = np.zeros((nlay, nrow, ncol), dtype=float)
        qz = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]
        head = headobj.get_data(totim=times[-1])

    else:
        concentration = {"pumping": ucnobj.get_data(kstpkper=(nstp-1,0))
        }
        head = {"pumping": headobj.get_data(kstpkper=(nstp-1,0))
        }
        qx = {"pumping": cbbobj.get_data(text="flow right face", kstpkper=(nstp-1,0))[0]
        }
        qy = {}
        try:
            qy["pumping"] = cbbobj.get_data(text="flow front face", kstpkper=(nstp-1,0))[0]
        except:
            qy["pumping"] = np.zeros((nlay, nrow, ncol), dtype=float) 
        qz = {"pumping": cbbobj.get_data(text="flow lower face", kstpkper=(nstp-1,0))[0]
        }

        if rec:
            concentration["recovery"] = ucnobj.get_data(kstpkper=(nstp-1,1))
            head["recovery"] = headobj.get_data(kstpkper=(nstp-1,1))
            qx["recovery"] = cbbobj.get_data(text="flow right face", kstpkper=(nstp-1,1))[0]
            qz["recovery"] = cbbobj.get_data(text="flow lower face", kstpkper=(nstp-1,1))[0]
            try:
                qy["recovery"] = cbbobj.get_data(text="flow front face", kstpkper=(nstp-1,1))[0]
            except:
                qy["recovery"] = np.zeros((nlay, nrow, ncol), dtype=float) 


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

def save_results_3D(swt, realization, concentration, qx, qy, qz, head, stress_period=None, name=None):
    
    ws = os.path.join(f'.\\results\\{name}')
    if not os.path.exists(ws):
        os.mkdir(ws)

    with open(os.path.join(ws, f"qx_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(qx))
    with open(os.path.join(ws, f"qy_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(qy))
    with open(os.path.join(ws, f"qz_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(qz))
    with open(os.path.join(ws, f"head_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(head))
    with open(os.path.join(ws, f"concentration_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(concentration))

def save_equivalent_parameters(qx, qy, qz, Kh, Kv, name, realization):
    ws = os.path.join(f'.\\results\\{name}')
    if not os.path.exists(ws):
        os.mkdir(ws)

    with open(os.path.join(ws, f"parameters_{realization}.npy"), 'wb') as f: np.save(f, np.array([qx, qy, qz, Kh, Kv]))

def load_equivalent_parameters(name, realization):
    
    ws = os.path.join(f'.\\results\\{name}')
    with open(os.path.join(ws, f"parameters_{realization}.npy"), 'rb') as f: pars = np.load(f, allow_pickle=True)

    return pars

def load_results_3D(modelname, realization, stress_period=None):

    ws = os.path.join(f'.\\results\\{modelname}')

    with open(os.path.join(ws, f"qx_{realization}{stress_period}.npy"), 'rb') as f: qx = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qy_{realization}{stress_period}.npy"), 'rb') as f: qy = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"qz_{realization}{stress_period}.npy"), 'rb') as f: qz = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"head_{realization}{stress_period}.npy"), 'rb') as f: head = np.load(f, allow_pickle=True)
    with open(os.path.join(ws, f"concentration_{realization}{stress_period}.npy"), 'rb') as f: concentration = np.load(f, allow_pickle=True)

    return qx, qy, qz, head, concentration

def get_time_evolution(swt, nstp):

    ws = swt._model_ws
    ucnobj = bf.UcnFile(os.path.join(ws, "MT3D001.UCN"), model=swt)
    times = ucnobj.get_times()
    cbbobj = bf.CellBudgetFile(os.path.join(ws, f'{swt._BaseModel__name}.cbc'))
    headobj = bf.HeadFile(os.path.join(ws, f'{swt._BaseModel__name}.hds'))
    delV = -1*swt.dis.delc[0]*swt.dis.delr[0]*swt.dis.botm[0][0]
    
    concentration_data = ucnobj.get_alldata()[nstp:]
    # budget_data = cbbobj.get_alldata()[15:]
    # head_data = headobj.get_alldata()[15:]
    (qlay, qrow, qcol, _) = swt.wel.stress_period_data['1'][-1]

    mixing_zone_evolution =  list(map(proc_functions.mixing_zone_volume, concentration_data, 2*nstp*[delV]))
    fresh_volume_evolution = list(map(proc_functions.fresh_water_volume, concentration_data, 2*nstp*[delV], 2000*[15]))
    well_salinity = concentration_data[:,qlay,qrow,qcol]
    
    return mixing_zone_evolution, fresh_volume_evolution, well_salinity

def compare_time_evolutions(times, arrays, realizations, metrics, colors, modelname):
    '''
        Compare time evolutions for a number of metrics between a number of
        realizations
    '''
    f, axs = plt.subplots(len(metrics), sharex=True, figsize=(15, 12))
    for mdx, (array, metric) in enumerate(zip(arrays, metrics)):
        for series, realization in zip(array, realizations):
            if realization == "homogenous":
                linestyle = "--"
            else:
                linestyle = "-"
            axs[mdx].plot(times, series, label=realization, linestyle=linestyle, color=colors[mdx])

        axs[mdx].set_title(f"{metric} over time") 
        axs[mdx].set_ylabel(metric)
        axs[mdx].set_xlabel("Time (days)")
        axs[mdx].axvline(x=times[-1]/2, color="grey", linestyle=":")
        axs[mdx].legend()

    results_location = f'.\\results\\{modelname}'
    plt.savefig(f"{results_location}\\time_evolution_{modelname}.jpg", dpi=1200)
    plt.show()

