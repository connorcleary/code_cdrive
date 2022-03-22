import numpy as np
from scipy.io import loadmat
import os
import coastal_aquifer_model as coastal
import calculate_K_eff as k_eff
import post_proc_utils as proc
import results_het_hom as results
 
mult = 5
flipped=False
Lx = 800.0*mult
Ly = 10.0
Lz = 25.0*mult
nlay = 50*mult
nrow = 1
ncol = 80*mult
head = 0.6*mult
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
    return np.sum(qx[:,:,0])
    # more accurate but negligable
    # return np.sum([np.sqrt(qx**2+qy**2) for qx, qy in zip(qx[:,:,0], qy[:,:,0])])

def main():

    name = "pirot2D_" + str(mult)

    swt = coastal.build_coastal_aquifer_model(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt/2, nstp)
    # change hk
    Kh, Kv = load_field()
    # swt.lpf.hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    # swt.lpf.vka = np.transpose(Kv, (2, 0, 1))[:,0,:]

    Kh_big = np.zeros((nlay, ncol))
    Kv_big = np.zeros((nlay, ncol))
    for i in range(mult):
        for j in range(mult):
            if not flipped:
                Kh_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kh, (2, 0, 1))[:,(10*i+j)%40,:]
                Kv_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kv, (2, 0, 1))[:,(10*i+j)%40,:]
            else:
                Kh_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kh, (2, 0, 1))[:,(10*i+j)%40,::(-1)**j]
                Kv_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kv, (2, 0, 1))[:,(10*i+j)%40,::(-1)**j]

    swt.lpf.hk = Kh_big
    swt.lpf.vka = Kv_big

    # swt.write_input()

    run_model(swt)
    concentration, qx, qy, qz, head_array = proc.extract_results(swt)# extract and save results
    proc.save_results(swt, "heterogenous", concentration, qx, qy, qz, head_array)
    # calculate hk and vka and qinflow
    
    Kv_eff, Kh_eff, qinflow_bad = k_eff.calculate_K_eff(name, swt.lpf.hk, swt.lpf.vka, Lx, Ly*10, Lz, head, 0)
    qinflow = calc_qinflow(qx, qy)
    # qinflow = 0.008917525


    swt = coastal.change_to_homogenous(swt, nlay, nrow, ncol, qinflow=qinflow)
    swt.dis.firstdt = 0.5e7
    swt.lpf.hk = Kh_eff
    swt.lpf.vka = Kv_eff
    run_model(swt)
    concentration, qx, qy, qz, head_array = proc.extract_results(swt)# extract and save results
    proc.save_results(swt, "homogenous", concentration, qx, qy, qz, head_array)
    # extract and save results
    print(qinflow)
    results.results(name)

if __name__=="__main__":
    main()