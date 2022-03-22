from winreg import QueryValueEx
import numpy as np
import flopy
from scipy.io import loadmat
import coastal_aquifer_model as coastal
import calculate_K_eff as k_eff
import post_proc_utils as proc

'''
Test pristine conditions in 3D normal size pirot model
'''

flipped=False
Lx = 800.0
Ly = 10.0
Lz = 25.0
nlay = 50
nrow = 40
ncol = 80
head = 0.6
perlen = 1.5e9
dt = 1e7
nstp = 15


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
    name = "pirot3D"

    swt = coastal.build_coastal_aquifer_model(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp)
    # change hk
    Kh, Kv = load_field()
    swt.lpf.hk = np.transpose(Kh, (2, 0, 1))
    swt.lpf.vka = np.transpose(Kv, (2, 0, 1))

    # run_model(swt)
    concentration, qx, qy, qz, head_array = proc.extract_results(swt)
    qinflow, qinflow_approx = calc_qinflow(qx, qy)
    Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, swt.lpf.hk, swt.lpf.vka, Lx, Ly, Lz, head, 0)
    
    swt = coastal.change_to_homogenous(swt, nlay, nrow, ncol, qinflow=qinflow)
    swt.lpf.hk = Kh_eff
    swt.lpf.vka = Kv_eff

    pass

if __name__=="__main__":
    main()