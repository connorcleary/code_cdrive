import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os

def build_coastal_aquifer_model(modelname, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp=15):

    model_ws = f".\\model_files\\{modelname}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = modelname
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    print(swt.namefile)

    delr = Lx/ncol
    delc = Ly/nrow
    delv = Lz/nlay

    top = 0.0
    botm = np.linspace(top - delv, -Lz, nlay)

    ipakcb = 53

    dis = flopy.modflow.ModflowDis(
        swt,
        nlay,
        nrow,
        ncol,
        nper=1,
        delr=delr,
        delc=delc,
        laycbd=0,
        top=top,
        botm=botm,
        perlen=perlen, # COULD BE A POTENTIAL PROBLEM
        nstp=nstp
    )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol/4):ncol] = -1

    strt = np.zeros((nlay, ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=1)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-8)

    oc = flopy.modflow.ModflowOc(
        swt,
        stress_period_data={(0, nstp-1): ["save head", "save budget"]},
        compact=True
    )

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol/4), ncol):
        for j in range(nrow):
            offshore_boundary_cells.append([0, j, i])

    # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # wel data
    wel_data = {}
    wel_sp1 = []
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []
    for cell in onshore_boundary_cells:
        #ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["WEL"]])
        ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["BAS6"]])
        #wel_sp1.append([cell[0], cell[1], cell[2], qinflow / nlay])
        chd_sp1.append([cell[0], cell[1], cell[2], head, head])
    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])

    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1
    wel_data[0] = wel_sp1

    #wel = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb=ipakcb)

    chd = flopy.modflow.ModflowChd(
        swt,
        stress_period_data=chd_data,
        ipakcb = ipakcb
    )

    sconc = np.zeros((nlay, nrow, ncol))
    for col in range(int(ncol/4), ncol):
        for row in range(nrow): 
            for lay in range(nlay):
                sconc[lay, row, col] = 35.0

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-5,
        prsity=0.35, # can change this
        sconc=sconc,#sconc, # can change this: to each area having different starti
        ifmtcn=0,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=dt
    )

    adv = flopy.mt3d.Mt3dAdv(swt, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=1, cclose=1e-7)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=dt,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()

    return swt 

def change_to_homogenous(swt, nlay, nrow, ncol, qinflow):

    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, -1] = -1
    swt.bas6.ibound = ibound

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol/4), ncol):
        for j in range(nrow):
            offshore_boundary_cells.append([0, j, i])

        # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # wel data
    wel_data = {}
    wel_sp1 = []
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []
    for cell in onshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["WEL"]])
        wel_sp1.append([cell[0], cell[1], cell[2], qinflow / (nlay*nrow)])
    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])

    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1
    wel_data[0] = wel_sp1

    wel = flopy.modflow.ModflowWel(
        swt, 
        stress_period_data=wel_data, 
        ipakcb=53
    )

    swt.chd.stress_period_data = chd_data
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data)
    swt.write_input()
    return swt