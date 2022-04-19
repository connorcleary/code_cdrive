import numpy as np
import flopy
from scipy.io import loadmat
import coastal_aquifer_model as coastal
import calculate_K_eff as k_eff
import post_proc_utils as proc
import results_het_hom as results
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt

flipped=False
Lx = 800.0
Ly = 10.0
Lz = 25.0
nlay = 50
nrow = 1
ncol = 80
head = 0.5
perlen = 1e9
dt = 1e6
nstp = 1000
q = 3e1
pumpT = 3650/2
recT = 3650/2
recQmult = 0
onshore_proportion=0.75
qrow = 0
qlay = 15
qcol = 15


def load_field():
    fields = loadmat(r'./fields/80x40x50aTest')
    Kh = fields['model_Kh']
    Kv = fields['model_Kv']

    return Kh, Kv

def run_model(swt):
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def calc_qinflow(qx, qy):
    qinflow_approx = np.sum(qx[:,:,0])
    qinflow = 0
    for xrow, zrow in zip(qx[:,:,0], qy[:,:,0]):
        for x, z in zip(xrow, zrow):
            qinflow += np.sqrt(x**2+z**2) 
    return qinflow, qinflow_approx

def run_entire_model(name):
    
    row = 0
    # swt = coastal.cam_with_pumping(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, 25, 0, 15, q, pumpT=36500)
    swt = coastal.cam_with_pumping_and_recovery(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, 14, 0, 20, q, pumpT=pumpT, recT=recT, recQmult=recQmult, onshore_proportion=onshore_proportion)
    # change hk
    Kh, Kv = load_field()
    swt.lpf.hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    swt.lpf.vka = np.transpose(Kv, (2, 0, 1))[:,0,:]

    run_model(swt)
    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=True, rec=True)
    proc.save_results_3D(swt, "heterogenous", concentration, qx, qy, qz, head_array)

    het_mix, het_fresh, het_wel = proc.get_time_evolution(swt, nstp)

    qinflow = np.zeros(3)
    _, qinflow[0] = calc_qinflow(qx["steady"], qz["steady"])
    _, qinflow[1] = calc_qinflow(qx["pumping"], qz["pumping"])
    _, qinflow[2] = calc_qinflow(qx["recovery"], qz["recovery"])
    Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, swt.lpf.hk, swt.lpf.vka, Lx, Ly, Lz, head, 0)

    swt = coastal.cam_with_pumping_and_recovery_homogenous(name, Lz, Lx, Ly, nlay, nrow, ncol, qinflow, perlen, dt, nstp, 14, 0, 20, q, pumpT=pumpT, recT=recT, recQmult=recQmult, onshore_proportion=onshore_proportion)
    swt.lpf.hk = Kh_eff
    swt.lpf.vka = Kv_eff

    run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=True, rec=True)
    proc.save_results_3D(swt, "homogenous", concentration, qx, qy, qz, head_array)

    # results.results_pumping_recovery_comparison(name, 25, 0, 15)
    results.results_pumping_recovery_comparison_with_heterogeneity(name, 14, 0, 20,  np.transpose(Kh, (2, 0, 1))[:,0,:], Kh_eff)

    hom_mix, hom_fresh, hom_wel = proc.get_time_evolution(swt, nstp)
    # proc.plot_time_evolution(swt)

    proc.compare_time_evolutions(np.linspace(0, pumpT+recT, nstp*2), [[het_mix, hom_mix], [het_fresh, hom_fresh], [het_wel, hom_wel]], ["heterogenous", "homogenous"], ["Mixing zone volume", "Freshwater volume", "Well salinity"], ["green", "blue", "red"], name)

def run_steady(name, k_type, realization, hk, vka, qinflow=None):  
    swt = coastal.cam_steady(name, Lz, Lx, Ly, nlay, nrow, ncol, perlen, dt, nstp,onshore_proportion, k_type, head=head, qinflow=qinflow)
    swt.lpf.hk = hk
    swt.lpf.vka = vka
    run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=False, rec=False)
    proc.save_results_3D(swt, realization, concentration, qx, qy, qz, head_array, "steady", name)

def run_transient(name, k_type, realization, hk, vka , start_conc, start_head, qinflow=None):
    swt = coastal.cam_transient(name, Lz, Lx, Ly, nlay, nrow, ncol, pumpT, recT, pumpT/1000, onshore_proportion, k_type, q, qlay, qrow, qcol, recQmult, start_head, start_conc, head=head, qinflow=qinflow)
    swt.lpf.hk = hk
    swt.lpf.vka = vka
    run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=True, rec=True)
    proc.save_results_3D(swt, realization, concentration, qx, qy, qz, head_array, "transient", name)

    mix, fresh, wel = proc.get_time_evolution(swt, nstp)

    return mix, fresh, wel

def run_single_2D():
    Kh, Kv = load_field()
    hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    vka = np.transpose(Kv, (2, 0, 1))[:,0,:]
    realization = "heterogenous_test"
    name = "steady_test"

    # run_steady(name, "heterogenous", "heterogenous_test", hk, vka)

    qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
    
    # # Kv_eff, Kh_eff, _ = k_eff.calculate_K_eff(name, hk, vka, Lx, Ly, Lz, head, 0)
    # # proc.save_equivalent_parameters(qx, qy, qz, Kh_eff, Kv_eff, name, realization)

    #  #pars = proc.load_equivalent_parameters(name, realization)

    mix, fresh, wel = run_transient(name, "heterogenous", realization, hk, vka , concentration, head_steady, qinflow=None)
    qx_t, qy_t, qz_t, head_t, concentration_t =  proc.load_results_3D(name, realization, stress_period="transient")

    # results.results_single_staged(name, realization, qlay, qrow, qcol)

    # proc.compare_time_evolutions(np.linspace(0, pumpT+recT, nstp*2), [[mix], [fresh], [wel]], ["heterogenous"], ["Mixing zone volume", "Freshwater volume", "Well salinity"], ["green", "blue", "red"], name)
    movement = proc.horizontal_movement_of_groundwater(concentration, concentration_t.item().get('pumping'))
    movement = movement*Lx/ncol

    results.plot_metric_over_layers("steady_test", movement, "Lateral movement of Saline Area")


def run_all_2D_realizations():
    name = "2D_ergodic"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))

    rows = hk_all.shape[1]

    concentration_steady = np.zeros((nlay, rows, ncol))
    concentration_transient = np.zeros((nlay, rows, ncol))
    for row in range(rows):
        realization = f"row{row}"
        # run_steady(name, "heterogenous", realization, hk_all[:,row,:], vka_all[:,row,:])
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
        # mix, fresh, wel = run_transient(name, "heterogenous", realization, hk_all[:,row,:], vka_all[:,row,:], concentration, head_steady, qinflow=None)
        qx_t, qy_t, qz_t, head_t, concentration_t =  proc.load_results_3D(name, realization, stress_period="transient")
        # results.results_single_staged(name, realization, hk_all[:,row,:], qlay, qrow, qcol)

        concentration_steady[:, row, :] = concentration[:, 0,:]
        concentration_transient[:, row, :] = concentration_t.item().get("pumping")[:, 0, :]

    movement = proc.horizontal_movement_of_groundwater(concentration_steady, concentration_transient)
    # movement_dist = movement*Lx/ncol
    # results.plot_metric_over_layers("steady_test", movement_dist, "Lateral movement of Saline Area")

    # heatmap = proc.probability_of_movement(movement, int(ncol*onshore_proportion), nlay)
    # X = np.linspace(Lx*onshore_proportion, 0, int(ncol*onshore_proportion))
    # Y = np.linspace(0, -Lz, nlay)
    # results.plot_heatmap(X, Y, heatmap, qlay)
    
    heatmap = proc.probability_of_salinization(concentration_steady, concentration_transient)
    X = np.linspace(0,Lx, ncol)
    Y = np.linspace(0, -Lz, nlay)

    results.plot_heatmap(X, Y, heatmap, qlay)



def main():
    run_all_2D_realizations()
    pass 

    


if __name__ == "__main__":
    main()