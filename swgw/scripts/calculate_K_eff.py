from warnings import filterwarnings
import numpy as np
import os
import flopy
import flopy.utils.binaryfile as bf
from scipy.io import loadmat

def build_model(Lx, Ly, Lz, modelname, hk, vka, direction, head):
    
    shape = hk.shape
    if len(shape) == 2:
        nrow = 1
        nlay, ncol = hk.shape
    else:    
        nlay, nrow, ncol = hk.shape

    gwf = flopy.modflow.Modflow(modelname, exe_name=r"C:\Users\ccl124\bin\MF2005.exe")

    delr = Lx/ncol
    delc = Ly/nrow
    delv = Lz/nlay

    top = 0.0 
    botm = np.linspace(top - delv, -Lz, nlay)

    dis = flopy.modflow.ModflowDis(
        gwf, 
        nlay, 
        nrow, 
        ncol, 
        delr=delr, 
        delc=delc, 
        top=top,
        botm=botm,
        itmuni=1,
        perlen= 0.001,
        tsmult = 1.1,
        nstp=100,
    )

    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    strt = np.ones((nlay, nrow, ncol), dtype=np.float32)

    if direction.lower()[0] == "h":
        ibound[:, :, 0] = -1
        ibound[:, :, -1] = -1
        strt[:, :, 0] = head
        strt[:, :, -1] = 0.0
    else:
        ibound[0, :, :] = -1
        ibound[-1, :, :] = -1
        strt[0, :, :] = head
        strt[-1, :, :] = 0.0

    bas = flopy.modflow.ModflowBas(
        gwf, 
        ibound=ibound, 
        strt=strt
    )

    lpf = flopy.modflow.ModflowLpf(
        gwf,
        hk=hk,
        vka=vka, 
        ipakcb=53,
        laytyp=1
        )

    spd = {(0, 0): ["print head", "print budget", "save head", "save budget"]}
    oc = flopy.modflow.ModflowOc(
        gwf, 
        stress_period_data=spd, 
        compact=True
    )

    pcg = flopy.modflow.ModflowPcg(gwf)

    return gwf

def run_model(gwf):
    gwf.write_input()
    success, buff = gwf.run_model()
    if not success:
        raise Exception("MODFLOW did not terminate normally.")




def extract_model(modelname):

    cbb = bf.CellBudgetFile(modelname + ".cbc")
    times = cbb.get_times()
    kstpkper_list = cbb.get_kstpkper()
    frf = cbb.get_data(text="FLOW RIGHT FACE", totim=times[-1])[0]
    flf = cbb.get_data(text="FLOW LOWER FACE", totim=times[-1])[0]
    return frf, flf

def calculate_K_eff(modelname, hk, vka, Lx, Ly, Lz, head, row=-1):

    # field = np.genfromtxt(os.path.join(r'.\fields', fieldname))

    for direction in ["horizontal", "vertical"]:
        gwf = build_model(Lx, Ly, Lz, modelname, hk, vka, direction, head)
        run_model(gwf)
        frf, flf = extract_model(modelname)
        if direction[0] == "v":
            Kv_eff = np.sum(flf[0,:,:])/(head/(Lz))/(Lx*Ly)
        else:
            qinflow = np.sum(frf[:,:,0])
            Kh_eff = np.sum(frf[:, :, 0])/(head/Lx)/(Lz*Ly)

    return Kv_eff, Kh_eff, qinflow

if __name__ == "__main__":
    fields = loadmat(r'./fields/80x40x50aTest')
    Kh = fields['model_Kh']
    Kv = fields['model_Kv']
    hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    vka = np.transpose(Kv, (2, 0, 1))[:,0,:]
    # Kv_eff, Kh_eff = calculate_K_eff("80x40x50aTest.mat", hk, vka, 800, 20, 25, 5.0, 0)
    # print(Kh_eff)
    # print(Kv_eff)