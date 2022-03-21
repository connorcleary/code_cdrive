import coastal_aquifer_model as coastal
import numpy as np
import flopy
from scipy.io import loadmat

'''
Test pristine conditions in 3D normal size pirot model
'''

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

def main():
    name = "pirot3D"

    swt = coastal.build_coastal_aquifer_model(name, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp)
    # change hk
    Kh, Kv = load_field()
    swt.lpf.hk = np.transpose(Kh, (2, 0, 1))
    swt.lpf.vka = np.transpose(Kv, (2, 0, 1))
