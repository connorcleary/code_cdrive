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
head = 0.6
perlen = 1.5e9
dt = 1e7
nstp = 15
q = 1e0

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

def main():
    name = "pirot2D_pumping_recovery_2_MAR"
    row = 0
    # swt = coastal.cam_with_pumping(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, 25, 0, 15, q, pumpT=36500)
    swt = coastal.cam_with_pumping_and_recovery(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, 25, 0, 15, q, pumpT=36500, recT=36500, recQmult=2)
    # change hk
    Kh, Kv = load_field()
    swt.lpf.hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    swt.lpf.vka = np.transpose(Kv, (2, 0, 1))[:,0,:]

    run_model(swt)

    concentration, qx, qy, qz, head_array = proc.extract_results(swt, pump=True, rec=True)
    proc.save_results_3D(swt, "heterogenous", concentration, qx, qy, qz, head_array)

    results.results_pumping_recovery(name, 25, 0, 20)




if __name__ == "__main__":
    main()