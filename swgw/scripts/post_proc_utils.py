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

def find_fresh_water_volume(concentration, delV, colq):
    '''
    Calculate the volume of fresh water seaward of the pumping locationY
    '''
    fresh = 0
    for i in range(concentration.shape[0]):
            for k in range(colq, concentration.shape[1]):
                if 0.35 > concentration[i][k]:
                    fresh += 1
    return fresh*5


def find_toe_penetration(conc):
    for i in range(conc.shape[1]):
        if conc[-1][i] > 0.5:
            return (conc.shape[1]/2 - i)*10

def find_sgd(conc, qz, onshore_proportion):
    sgd = 0
    for j in range(int(conc.shape[1]*onshore_proportion), conc.shape[1]):
            if 3.5 >= conc[0, j] and qz[0, j] < 0:
                sgd += np.max([float(0), -qz[0,j]])
    return sgd

def find_mixing_com(conc, Lx):
    x_tot = 0 
    mix_n = 0

    for i in range(conc.shape[0]):
        for j in range(conc.shape[1]):
            if 3.5 <= conc[i, j] <= 31.5:
                mix_n += 1 
                x_tot += j/conc.shape[1]*Lx-Lx/2

    return x_tot/mix_n

def plot_conc(ax, swt, results_dict, row=0, **kwargs):

    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": row})
    arr = pmv.plot_array(results_dict['concentration'], cmap="viridis", **kwargs)
    pmv.plot_vector(results_dict['qx'], -results_dict['qz'], -results_dict['qz'], kstep=3, hstep=3, normalize=False, color="white")
    # plt.colorbar(arr, shrink=0.5, ax=ax)

    # ax.set_title("Simulated Concentrations")

    return ax, arr


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
    # plt.colorbar(arr, shrink=0.5, ax=ax)
    # ax.set_title("Simulated Heads")

    return ax, arr

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

def load_time_evolution_3D(modelname, realization, stress_period=None):

    ws = os.path.join(f'.\\results\\{modelname}')

    with open(os.path.join(ws, f"concentration_time_evolution_{realization}{stress_period}.npy"), 'rb') as f: evol = np.load(f, allow_pickle=True)

    return evol

def get_time_evolution(swt, nstp, steady=False):

    ws = swt._model_ws
    ucnobj = bf.UcnFile(os.path.join(ws, "MT3D001.UCN"), model=swt)
    times = ucnobj.get_times()
    cbbobj = bf.CellBudgetFile(os.path.join(ws, f'{swt._BaseModel__name}.cbc'))
    headobj = bf.HeadFile(os.path.join(ws, f'{swt._BaseModel__name}.hds'))
    delV = -1*swt.dis.delc[0]*swt.dis.delr[0]*swt.dis.botm.array[0][0][0]
    
    concentration_data = ucnobj.get_alldata()[:]
    # budget_data = cbbobj.get_alldata()[15:]
    # head_data = headobj.get_alldata()[15:]
    # if not steady:
    #     (qlay, qrow, qcol, _) = swt.wel.stress_period_data['1'][-1]
    #     well_salinity = concentration_data[:,qlay,qrow,qcol]
    # else:
    #     well_salinity = None

    # mixing_zone_evolution =  list(map(proc_functions.mixing_zone_volume, concentration_data, 2*nstp*[delV]))
    # fresh_volume_evolution = list(map(proc_functions.fresh_water_volume, concentration_data, 2*nstp*[delV], nstp*2*[15]))
    
    return concentration_data

def save_concentration_time_evolution(modelname, realization, concentration_time_evolution, stress_period=None):
    ws = os.path.join(f'.\\results\\{modelname}')
    with open(os.path.join(ws, f"concentration_time_evolution_{realization}{stress_period}.npy"), 'wb') as f: np.save(f, np.array(concentration_time_evolution))

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

def horizontal_movement_of_groundwater(conc1, conc2):
    '''
    Find horizontal distance moved of saline front.
    '''
    original_extent = np.zeros((conc1.shape[0], conc1.shape[1]))
    new_extent = np.zeros((conc1.shape[0], conc1.shape[1]))
    for row in range(conc1.shape[1]):
        for lay in range(conc1.shape[0]):
            new_found = False
            old_found = False
            for col in range(conc1.shape[2]):
                if conc1[lay, row, col] >= 0.35 and not old_found:
                    original_extent[lay, row] = col
                    old_found = True     

                if conc2[lay, row, col] >= 0.35 and conc2[lay, row, col] <= 35 and not new_found:
                    new_extent[lay, row] = col
                    new_found = True
            
    return original_extent - new_extent

def probability_of_movement(movement, ncol_onshore, nlay):
    heatmap = np.zeros((nlay, ncol_onshore))
    for row in range(movement.shape[1]):
        for lay in range(movement.shape[0]):
            if movement[lay, row] > 0:
                heatmap[lay,-int(movement[lay, row])] += 1
            
    heatmap = heatmap/movement.shape[1]

    return heatmap

def probability_of_salinization(conc1, conc2):

    heatmap = np.zeros_like(conc1[:, 0, :])
    for row in range(conc1.shape[1]):
        for lay in range(conc1.shape[0]):
            for col in range(conc1.shape[2]):
                if conc1[lay, row, col] < 3.5 and conc2[lay, row, col] > 3.5:
                    heatmap[lay, col] += 1
                
    heatmap = heatmap/conc1.shape[1]
    return heatmap

def get_steady_state_time_evolutions(name, realizations, rows):

    concentration = [[] for i in range(len(realizations))]
    
    for i, realization in enumerate(realizations):
        for j in range(rows):
            concentration[i].append(load_time_evolution_3D(name, f"row{j}{realization}","steady"))

    # com = np.zeros((i+1, j+1, int(0.05*len(concentration[0][0])))) 
    # toe = np.zeros((i+1, j+1, int(0.05*len(concentration[0][0]))))
    # mix = np.zeros((i+1, j+1, int(0.05*len(concentration[0][0]))))
    com = np.zeros((i+1, j+1, 30))
    toe = np.zeros((i+1, j+1, 30))
    mix = np.zeros((i+1, j+1, 30))

    for i in range(len(realizations)):
        for j in range(rows):
            for kdx, k in enumerate(np.logspace(0, 3.56, base=10, num=30).astype(int)):
                com[i, j, kdx] = find_mixing_com(concentration[i][j][k][:,0,:], 800)
                toe[i, j, kdx] = find_toe_penetration(concentration[i][j][k][:,0,:])
                mix[i, j, kdx] = find_mixing_zone(concentration[i][j][k][:,0,:])

    results_location = f'.\\results\\{name}'
    np.save(f"{results_location}\\com_evolution_log", com)
    np.save(f"{results_location}\\toe_evolution_log", toe)
    np.save(f"{results_location}\\mix_evolution_log", mix)
    pass

def get_transient_time_evolutions(name, realizations, rows, qcol, qlay):

    concentration = [[] for i in range(len(realizations))]
    
    for i, realization in enumerate(realizations):
        for j in range(rows):
            concentration[i].append(load_time_evolution_3D(name, f"row{j}{realization}","transient"))

    com = np.zeros((i+1, j+1, 50))
    toe = np.zeros((i+1, j+1, 50))
    mix = np.zeros((i+1, j+1, 50))
    fresh = np.zeros((i+1, j+1, 50))
    wel = np.zeros((i+1, j+1, 50))

    for i in range(len(realizations)):
        for j in range(rows):
            for kdx, k in enumerate(np.linspace(0, 1000, 50).astype(int)):
                com[i, j, kdx] = find_mixing_com(concentration[i][j][k][:,0,:], 800)
                toe[i, j, kdx] = find_toe_penetration(concentration[i][j][k][:,0,:])
                mix[i, j, kdx] = find_mixing_zone(concentration[i][j][k][:,0,:])
                fresh[i, j, kdx] = find_fresh_water_volume(concentration[i][j][k][:,0,:], 5, qcol)
                wel[i, j, kdx] = concentration[i][j][k][qlay,0,qcol]

    results_location = f'.\\results\\{name}'
    np.save(f"{results_location}\\com_evolution_transient", com)
    np.save(f"{results_location}\\toe_evolution_transient", toe)
    np.save(f"{results_location}\\mix_evolution_transient", mix)
    np.save(f"{results_location}\\fresh_evolution_transient", fresh)
    np.save(f"{results_location}\\wel_evolution_transient", wel)
    pass

def load_steady_metric_evolutions(name):

    results_location = f'.\\results\\{name}'
    com = np.load(f"{results_location}\\com_evolution_log.npy", allow_pickle=True)
    toe = np.load(f"{results_location}\\toe_evolution_log.npy", allow_pickle=True)
    mix = np.load(f"{results_location}\\mix_evolution_log.npy", allow_pickle=True)

    return com, toe, mix

def load_transient_metric_evolutions(name):
    results_location = f'.\\results\\{name}'
    com = np.load(f"{results_location}\\com_evolution_transient.npy", allow_pickle=True)
    toe = np.load(f"{results_location}\\toe_evolution_transient.npy", allow_pickle=True)
    mix = np.load(f"{results_location}\\mix_evolution_transient.npy", allow_pickle=True)
    fresh = np.load(f"{results_location}\\fresh_evolution_transient.npy", allow_pickle=True)
    wel = np.load(f"{results_location}\\wel_evolution_transient.npy", allow_pickle=True)

    return com, toe, mix, fresh, wel