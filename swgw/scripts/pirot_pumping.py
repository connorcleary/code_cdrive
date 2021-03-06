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
Ly = 400.0
Lz = 25.0
nlay = 50
nrow = 40
ncol = 80
head = 0.6
perlen = 3.65e9
dt = 1e6
nstp = 1000
q = 1.5e1
pumpT = 3650/2
recT = 3650/2
recQmult = 0
onshore_proportion=0.5
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
    swt = coastal.cam_steady(name, Lz, Lx, Ly, nlay, nrow, ncol, perlen, dt, int(perlen/dt) ,onshore_proportion, k_type, head=head, qinflow=qinflow)
    swt.lpf.hk = hk
    swt.lpf.vka = vka
    # run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=False, rec=False)
    proc.save_results_3D(swt, realization, concentration, qx, qy, qz, head_array, "steady", name)

    concentration_data, qx_data, qz_data = proc.get_time_evolution(swt, nstp, steady=True)
    proc.save_concentration_time_evolution(name, realization, concentration_data, qx_data, qz_data, stress_period="steady")

def run_transient(name, k_type, realization, hk, vka , start_conc, start_head, qinflow=None):
    swt = coastal.cam_transient(name, Lz, Lx, Ly, nlay, nrow, ncol, pumpT, recT, pumpT/1000, onshore_proportion, k_type, q, qlay, qrow, qcol, recQmult, start_head, start_conc, head=head, qinflow=qinflow)
    swt.lpf.hk = hk
    swt.lpf.vka = vka
    run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, nstp, pump=True, rec=True)
    proc.save_results_3D(swt, realization, concentration, qx, qy, qz, head_array, "transient", name)

    concentration_data, mix, fresh, wel = proc.get_time_evolution(swt, nstp)
    proc.save_concentration_time_evolution(name, realization, concentration_data, stress_period="transient")

def compare_steady_states_2D():
    name = "2D_steady_state_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    concentration_steady = np.zeros((nlay, rows, ncol))
    concentration_transient = np.zeros((nlay, rows, ncol))
    for row in range(rows):
        realization = f"row{row}"
        #run_steady(name, "heterogenous", realization, hk_all[:,row,:], vka_all[:,row,:])
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
        Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
        _, qinflow_steady = calc_qinflow(qx, qz)
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="transient")
        _, qinflow_transient = calc_qinflow(qx, qz)

        run_steady(name, "heterogenous", realization+"equivalent_head", Kh_eff, Kv_eff)
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_head", stress_period="steady")
        run_transient(name, "heterogenous", realization, Kh_eff, Kv_eff, concentration, head_steady, qinflow=None)
    	

        run_steady(name, "homogenous", realization+"equivalent_flux", Kh_eff, Kv_eff, qinflow=qinflow_steady)
        qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_flux", stress_period="steady")
        run_transient(name, "heterogenous", realization, Kh_eff, Kv_eff, concentration, head_steady, qinflow=qinflow_transient)



    Kh_array = np.ones((nlay, ncol))*Kh_eff
    results.compare_steady_states(name, ['row0', 'row0equivalent_head', 'row0equivalent_flux'], [hk_all[:,row,:], Kh_array, Kh_array])

def compare_transient_2D():
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    for row in range(rows):
        Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
        Kh_array = np.ones((nlay, ncol))*Kh_eff
        results.compare_transient_response(name, [f'row{row}', f'row{row}equivalent_head', f'row{row}equivalent_flux'], [hk_all[:,row,:], Kh_array, Kh_array], row, qlay, qcol)

def plot_effective_conductivities():
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    Kv_all = np.zeros(rows)
    Kh_all = np.zeros(rows)

    for row in range(rows):
        Kv_all[row], Kh_all[row], qinflow_bad = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
        
    results.plot_effective_conductivities(Kh_all, Kv_all)

def compare_wel_conductivity():
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    Kv_all = np.zeros(rows)
    Kh_all = np.zeros(rows)

    for row in range(rows):
        Kv_all[row], Kh_all[row], qinflow_bad = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
    
    com, toe, mix, fresh, wel = proc.load_transient_metric_evolutions(name)
    realizations = ["heterogenous", "equivalent_flux", "equivalent_head"]
    wel = wel[:,:,-1]
    results.compare_wel_salinity_with_effective_conductivity(Kh_all, Kv_all, wel, realizations)
    pass

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
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))

    rows = hk_all.shape[1]

    concentration_steady = np.zeros((nlay, rows, ncol))
    concentration_transient = np.zeros((nlay, rows, ncol))
    for row in range(rows):
        
            realization = f"row{row}"
            run_steady(name, "heterogenous", realization, hk_all[:,row,:], vka_all[:,row,:])
            qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
            Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
            _, qinflow_steady = calc_qinflow(qx, qz)
            run_transient(name, "heterogenous", realization, hk_all[:,row,:], vka_all[:,row,:], concentration, head_steady, qinflow=None)
            qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="transient")
            _, qinflow_pumping = calc_qinflow(qx.item().get('pumping'), qz.item().get('pumping'))
            _, qinflow_recovery = calc_qinflow(qx.item().get('recovery'), qz.item().get('recovery'))

            run_steady(name, "heterogenous", realization+"equivalent_head", Kh_eff, Kv_eff)
            qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_head", stress_period="steady")
            run_transient(name, "heterogenous", realization+"equivalent_head", Kh_eff, Kv_eff, concentration, head_steady, qinflow=None)
            
            run_steady(name, "homogenous", realization+"equivalent_flux", Kh_eff, Kv_eff, qinflow=qinflow_steady)
            qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_flux", stress_period="steady")
            run_transient(name, "homogenous", realization+"equivalent_flux", Kh_eff, Kv_eff, concentration, head_steady, qinflow=[qinflow_pumping, qinflow_recovery])
        
            print("Something went wrong")

        # if row%5 == 0:
        #     results.results_single_staged(name, realization+"equivalent_head", hk_all[:,row,:], qlay, qrow, qcol)
        #     results.results_single_staged(name, realization+"equivalent_flux", hk_all[:,row,:], qlay, qrow, qcol)

        # concentration_steady[:, row, :] = concentration[:, 0,:]
        # concentration_transient[:, row, :] = concentration_t.item().get("pumping")[:, 0, :]

    # movement = proc.horizontal_movement_of_groundwater(concentration_steady, concentration_transient)
    # movement_dist = movement*Lx/ncol
    # results.plot_metric_over_layers("steady_test", movement_dist, "Lateral movement of Saline Area")

    # heatmap = proc.probability_of_movement(movement, int(ncol*onshore_proportion), nlay)
    # X = np.linspace(Lx*onshore_proportion, 0, int(ncol*onshore_proportion))
    # Y = np.linspace(0, -Lz, nlay)
    # results.plot_heatmap(X, Y, heatmap, qlay)
    
    # heatmap = proc.probability_of_salinization(concentration_steady, concentration_transient)
    # X = np.linspace(0,Lx, ncol)
    # Y = np.linspace(0, -Lz, nlay)

    # results.plot_heatmap(X, Y, heatmap, qlay)

def run_all_3D_realizations(name):

    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))

    rows = hk_all.shape[1]

    concentration_steady = np.zeros((nlay, rows, ncol))
    concentration_transient = np.zeros((nlay, rows, ncol))

    realization = f"pirot_basic_all_fresh"
    run_steady(name, "heterogenous", realization, hk_all, vka_all)
    qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="steady")
    Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, hk_all, vka_all, Lx, Ly, Lz, head, 0)
    _, qinflow_steady = calc_qinflow(qx, qz)
    # run_transient(name, "heterogenous", realization, hk_all, vka_all, concentration, head_steady, qinflow=None)
    # qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization, stress_period="transient")
    # _, qinflow_pumping = calc_qinflow(qx.item().get('pumping'), qz.item().get('pumping'))
    # _, qinflow_recovery = calc_qinflow(qx.item().get('recovery'), qz.item().get('recovery'))

    # run_steady(name, "heterogenous", realization+"equivalent_head", Kh_eff, Kv_eff)
    # qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_head", stress_period="steady")
    # run_transient(name, "heterogenous", realization+"equivalent_head", Kh_eff, Kv_eff, concentration, head_steady, qinflow=None)
    
    # run_steady(name, "homogenous", realization+"equivalent_flux", Kh_eff, Kv_eff, qinflow=qinflow_steady)
    # qx, qy, qz, head_steady, concentration =  proc.load_results_3D(name, realization+"equivalent_flux", stress_period="steady")
    # run_transient(name, "homogenous", realization+"equivalent_flux", Kh_eff, Kv_eff, concentration, head_steady, qinflow=[qinflow_pumping, qinflow_recovery])

def analyse_2D_steady_state():
    realization_suffices = ["", "equivalent_flux", "equivalent_head"]
    rows = range(40)
    name = "2D_ergodic_0.5"

    results.compute_steadystate_statistics(name, realization_suffices, rows)

def compare_distance_moved(name, realizations, rows):

    movement = np.zeros((len(realizations), rows, nlay))
    for i, realization in enumerate(realizations):
        for row in range(rows):
            realization_name = f"row{row}{realization}"
            _, _, _, _, concentration =  proc.load_results_3D(name, realization_name, stress_period="steady")
            _, _, _, _, concentration_t =  proc.load_results_3D(name, realization_name, stress_period="transient")
            movement_temp = proc.horizontal_movement_of_groundwater(concentration, concentration_t.item().get('pumping'))

            movement[i,row,:] = movement_temp[:,0]*Lx/ncol

    realizations[0] = "heterogenous"
    results.plot_metric_over_layers(name, movement.T, realizations, "Lateral movement of Saline Area")

def compare_probability_distance_moved(name, realizations, rows):

    movement = np.zeros((len(realizations), rows, nlay))
    for i, realization in enumerate(realizations):
        for row in range(rows):
            realization_name = f"row{row}{realization}"
            _, _, _, _, concentration =  proc.load_results_3D(name, realization_name, stress_period="steady")
            _, _, _, _, concentration_t =  proc.load_results_3D(name, realization_name, stress_period="transient")
            movement_temp = proc.horizontal_movement_of_groundwater(concentration, concentration_t.item().get('pumping'))

            movement[i,row,:] = movement_temp[:,0]*Lx/ncol

    realizations[0] = "heterogenous"

    results.plot_metric_as_heatmap_layered(name, movement, realizations, "title")

def compare_heat_maps(name, realizations, rows):
    heatmaps = [[]]*len(realizations)
    concentration_steady = np.zeros((nlay, rows, ncol))
    concentration_transient = np.zeros((nlay, rows, ncol))
    X = np.linspace(0, Lx, int(ncol))
    Y = np.linspace(0, -Lz, nlay)

    for i, realization in enumerate(realizations):
        concentration_steady = np.zeros((nlay, rows, ncol))
        concentration_transient = np.zeros((nlay, rows, ncol))
        for row in range(rows):
            realization_name = f"row{row}{realization}"
            _, _, _, _, concentration =  proc.load_results_3D(name, realization_name, stress_period="steady")
            _, _, _, _, concentration_t =  proc.load_results_3D(name, realization_name, stress_period="transient")
            concentration_steady[:, row, :] = concentration[:, 0,:]
            concentration_transient[:, row, :] = concentration_t.item().get("pumping")[:, 0, :]
            
        heatmaps[i] = proc.probability_of_salinization(concentration_steady, concentration_transient)
    
    realizations[0] = "heterogenous"
    results.plot_heatmaps(X, Y, name, heatmaps, qlay, realizations)

def compare_steady_evolutions_with_hk():
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    Kh_all = np.zeros(rows)

    for row in range(rows):
        _, Kh_all[row], _ = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
        
    results.plot_steady_time_evolution_ensemble(name, realizations=["heterogenous", "equivalent_flux", "equivalent_head"], hk=Kh_all)

def compare_3D_2D_steady_evolutions_with_hk():
    name = "2D_ergodic_0.5"
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))
    vka_all = np.transpose(Kv, (2, 0, 1))
    rows = hk_all.shape[1]

    Kh_all = np.zeros(rows)

    for row in range(rows):
        _, Kh_all[row], _ = k_eff.calculate_K_eff(name, hk_all[:,row,:], vka_all[:,row,:], Lx, Ly, Lz, head, 0)
        
    results.plot_steady_time_evolution_ensemble_2D_vs_3D("2D_ergodic_0.5" ,"3D_0.5", realizations=["2D", "3D",], hk=Kh_all)

def compare_3D_2D_steady_crosssections():
    name2D = "2D_ergodic_0.5"
    name3D = "3D_0.5"
    realizations = ["", "pirot_basic"]
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))

    for row in range(40):
        results.compare_steady_states_2D_3D(name2D, name3D, realizations, hk_all, row)

def compare3D_time_evolution():
    name = "3D_0.5_fresh"
    swt = flopy.seawat.swt.Seawat.load(f'.\\model_files\\{name}\{name}_steady.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    # concentration_time_evolution, qx, qz = proc.get_time_evolution(swt, nstp, steady=True)
    # proc.save_concentration_time_evolution(name, "pirot_basic_all_fresh", concentration_time_evolution, qx, qz, stress_period="steady")
    proc.get_steady_state_time_evolutions(name, ["pirot_basic_all_fresh"], 40, 3)
    Kh, Kv = load_field()
    hk_all = np.transpose(Kh, (2, 0, 1))

    for row in range(40):
        results.steady_states_gridspec_evolution3D(name, "pirot_basic_all_fresh", hk_all, row)

def main():
    # compare_steady_states_2D()
    # run_all_2D_realizations()
    # analyse_2D_steady_state()
    name = "3D_0.5_FRESH"
    # run_all_3D_realizations(name)
    Kh, Kv = load_field()
    realizations = ["", "equivalent_flux", "equivalent_head"]
    realizations= ["pirot_basic"]
    
    compare3D_time_evolution()
    #compare_3D_2D_steady_evolutions_with_hk()
    #compare_3D_2D_steady_crosssections()
    # results.steady_boxplots_2D_vs_3D("3D_0.5", "2D_ergodic_0.5", realizations, 40)
    ## proc.get_steady_state_time_evolutions(name, realizations, 40, 3)

    # results.plot_steady_time_evolutions(name, realizations)
    # results.plot_steady_time_evolution_ensemble(name)

    # proc.get_transient_time_evolutions(name, realizations, 40, qcol, qlay)
    # realizations = ["", "equivalent_flux", "equivalent_head"]
    # results.compute_steadystate_statistics(name, realizations, range(40))
    # results.plot_transient_time_evolutions(name, realizations)
    # results.transient_box_plots(name, realizations, range(40), qcol, qlay)
    # results.plot_well_ensemble(name, realizations)
    # compare_transient_2D()
    # plot_effective_conductivities()
    #compare_wel_conductivity()
    # compare_probability_distance_moved(name, realizations, 40)
    # compare_heat_maps(name, realizations, 40)
    # compare_steady_evolutions_with_hk()

    pass


if __name__ == "__main__":
    main()