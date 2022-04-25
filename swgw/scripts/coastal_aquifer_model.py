from msilib import change_sequence
import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os

def build_coastal_aquifer_model(modelname, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt,  onshore_proportion, nstp=15):

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
    ibound[0, :, int(ncol*onshore_proportion):ncol] = -1

    strt = np.zeros((nlay, nrow, ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=1)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)

    oc = flopy.modflow.ModflowOc(
        swt,
        stress_period_data={(0, 0): ["save head", "save budget"], 
                            (0, nstp-1): ["save head", "save budget"]},
        compact=True
    )

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
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
        sconc= sconc, # can change this: to each area having different starti
        ifmtcn=0,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=dt
    )

    adv = flopy.mt3d.Mt3dAdv(swt, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
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

    # swt.write_input()

    return swt 

def change_to_homogenous(swt, nlay, nrow, ncol, qinflow, onshore_proportion):

    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
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
        ibound[cell[0], cell[1], cell[2]] = -1

    swt.bas6.ibound = ibound

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

def add_pumping_period(swt, lay, row, col, q, pumpT, recT=None, recMult=None):
    dis = swt.dis
    oc = swt.oc
    ssm = swt.ssm
    chd = swt.chd
    itype = flopy.mt3d.Mt3dSsm.itype_dict()

    if recMult:
        dis.nper = 3
        dis.perlen = [dis.perlen, pumpT, recT]
        dis.steady = [True, False, False]
        oc.stress_period_data[(1, dis.nstp-1)] =  ["save head", "save budget"]
        oc.stress_period_data[(2, dis.nstp-1)] =  ["save head", "save budget"]
        ssm.ssm_data[1] = ssm.ssm_data[0] # change
        chd.chd_data[1] = chd.chd_data[0]
        wel.wel_data[1] = wel.wel_data[0] # change
        ssm.ssm_data[2] = ssm.ssm_data[0] # change
        chd.chd_data[2] = chd.chd_data[0]
        wel.wel_data[2] = wel.wel_data[0] # change
    else: 
        dis.nper = 2
        dis.perlen.shape = (2,) 
        dis.nstp.shape = (2,)
        dis.tsmult.shape = (2,)
        dis.perlen[1] = pumpT
        dis.nstp[1] = 15
        dis.tsmult[1] = 1
        dis.steady.shape = (2,) 
        dis.steady[1] = [False]
        oc.stress_period_data[(1, dis.nstp[0]-1)] =  ["save head", "save budget"]
        new = np.array([(lay, row, col, 0.0, itype["WEL"])], dtype=ssm.stress_period_data[0].dtype.descr).view(np.recarray)
        ssm.stress_period_data.data[1] = np.append(ssm.stress_period_data[0], new)
        chd.stress_period_data.data[1] = chd.stress_period_data[0]
        wel_data = {}
        wel_sp1 = []
        wel_sp1 = np.append(wel_sp1, (lay, row, col, -q))
        wel_data[1] = wel_sp1
        
        wel = flopy.modflow.ModflowWel(
            swt, 
            stress_period_data=wel_data, 
            ipakcb=53
        )

    return swt

def cam_with_pumping(modelname, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, qlay, qrow, qcol, q, pumpT):
    
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
        nper=2,
        delr=delr,
        delc=delc,
        laycbd=0,
        top=top,
        botm=botm,
        perlen=[perlen, pumpT], # COULD BE A POTENTIAL PROBLEM
        nstp=nstp,
        steady=[True, False]
    )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol/4):ncol] = -1

    strt = np.zeros((nlay, nrow, ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laytyp[0] = 1
    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=0)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)

    oc = flopy.modflow.ModflowOc(
        swt,
        stress_period_data={(0, 0): ["save head", "save budget"], 
                            (0, nstp-1): ["save head", "save budget"],
                            (1, nstp-1): ["save head", "save budget"]},
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
    # ssm_data[1] = ssm_sp1
    ssm_data[1] = np.append(ssm_sp1, [[qlay, qrow, qcol, 35.0, itype['WEL']]], axis=0)
    chd_data[0] = chd_sp1
    chd_data[1] = chd_sp1

    wel_data = {}
    wel_sp1 = []
    wel_sp1 = np.append(wel_sp1, (qlay, qrow, qcol, -q))
    wel_data[1] = wel_sp1
        
    wel = flopy.modflow.ModflowWel(
        swt, 
        stress_period_data=wel_data, 
        ipakcb=53
    )  

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
        nprs=0,
        prsity=0.35, # can change this
        sconc= sconc, # can change this: to each area having different starti
        ifmtcn=0,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=[dt, dt/100000]
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=500)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=pumpT/15,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()

    return swt 

def cam_with_pumping_and_recovery(modelname, Lz, Lx, Ly, nlay, nrow, ncol, head, perlen, dt, nstp, qlay, qrow, qcol, q, pumpT, recT, recQmult, onshore_proportion=0.25):
    
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
        nper=3,
        itmuni=4,
        delr=delr,
        delc=delc,
        laycbd=0,
        top=top,
        botm=botm,
        perlen=[perlen, pumpT, recT], # COULD BE A POTENTIAL PROBLEM
        nstp=nstp,
        steady=[True, False, False]
    )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, 0] = -1
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol*onshore_proportion):ncol] = -1

    strt = np.zeros((nlay, nrow, ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laytyp[0] = 1
    laywet=np.zeros(50)
    laywet[0] = 1
    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)

    oc_spd = {}
    for kper in range(3):
        for kstp in range(nstp):
            oc_spd[(kper, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)
    

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
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
    # ssm_data[1] = ssm_sp1
    ssm_data[1] = np.append(ssm_sp1, [[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    ssm_data[2] = np.append(ssm_sp1, [[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    chd_data[0] = chd_sp1
    chd_data[1] = chd_sp1
    chd_data[2] = chd_sp1

    wel_data = {}
    wel_sp1 = []
    wel_sp1 = np.append(wel_sp1, (qlay, qrow, qcol, -q))
    wel_data[1] = wel_sp1
    wel_sp2 = []
    wel_sp2 = np.append(wel_sp2, (qlay, qrow, qcol, q*recQmult))
    wel_data[2] = wel_sp2
         
    wel = flopy.modflow.ModflowWel(
        swt, 
        stress_period_data=wel_data, 
        ipakcb=53
    )  

    chd = flopy.modflow.ModflowChd(
        swt,
        stress_period_data=chd_data,
        ipakcb = ipakcb
    )

    sconc = np.zeros((nlay, nrow, ncol))
    for col in range(int(ncol*onshore_proportion), ncol):
        for row in range(nrow): 
            for lay in range(nlay):
                sconc[lay, row, col] = 35.0

    times=[t for t in np.linspace(perlen/nstp, perlen, nstp)]\
            + [t for t in np.linspace(perlen+pumpT/nstp, perlen+pumpT, nstp)]\
            + [t for t in np.linspace(perlen+pumpT+recT/nstp, perlen+pumpT+recT, nstp)]

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-1,# len(times),
        # timprs=times,
        prsity=0.35, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=[perlen/nstp, pumpT/nstp, recT/nstp],
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=500)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=recT/nstp,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()

    return swt 

def cam_with_pumping_and_recovery_homogenous(modelname, Lz, Lx, Ly, nlay, nrow, ncol, qinflow, perlen, dt, nstp, qlay, qrow, qcol, q, pumpT, recT, recQmult, onshore_proportion=0.25):
    
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
        nper=3,
        itmuni=4,
        delr=delr,
        delc=delc,
        laycbd=0,
        top=top,
        botm=botm,
        perlen=[perlen, pumpT, recT], # COULD BE A POTENTIAL PROBLEM
        nstp=nstp,
        steady=[True, False, False]
    )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol*onshore_proportion):ncol] = -1

    strt = np.zeros((nlay, nrow, ncol))

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laytyp[0] = 1
    laywet=np.zeros(50)
    laywet[0] = 1
    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)
    
    oc_spd = {}
    for kper in range(3):
        for kstp in range(nstp):
            oc_spd[(kper, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(
        swt,
        stress_period_data=oc_spd,
        compact=True
    )

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
        for j in range(nrow):
            offshore_boundary_cells.append([0, j, i])
       # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # wel data
    wel_data = {}
    wel_sp1 = []
    wel_sp2 = []
    wel_sp3 = []
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []
    for cell in onshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["WEL"]])
        wel_sp1.append([cell[0], cell[1], cell[2], qinflow[0] / nlay])
        wel_sp2.append([cell[0], cell[1], cell[2], qinflow[1] / nlay])
        wel_sp3.append([cell[0], cell[1], cell[2], qinflow[2] / nlay])
    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])
    wel_data[0] = wel_sp1

    ssm_data[0] = ssm_sp1
    ssm_data[1] = np.append(ssm_sp1, [[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    ssm_data[2] = np.append(ssm_sp1, [[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    chd_data[0] = chd_sp1
    chd_data[1] = chd_sp1
    chd_data[2] = chd_sp1

    wel_sp2 = np.append(wel_sp2, [(qlay, qrow, qcol, -q)], axis=0)
    wel_data[1] = wel_sp2
    wel_sp3 = np.append(wel_sp3, [(qlay, qrow, qcol, q*recQmult)], axis=0)
    wel_data[2] = wel_sp3
         
    wel = flopy.modflow.ModflowWel(
        swt, 
        stress_period_data=wel_data, 
        ipakcb=53
    )  

    chd = flopy.modflow.ModflowChd(
        swt,
        stress_period_data=chd_data,
        ipakcb = ipakcb
    )

    sconc = np.zeros((nlay, nrow, ncol))
    for col in range(int(ncol*onshore_proportion), ncol):
        for row in range(nrow): 
            for lay in range(nlay):
                sconc[lay, row, col] = 35.0

    # times=np.concatenate([t for t in np.linspace(perlen/nstp, perlen, nstp)],
    #                      [t for t in np.linspace(perlen+pumpT/nstp, perlen+pumpT, nstp)],
    #                      [t for t in np.linspace(perlen+pumpT+recT/nstp, perlen+pumpT+recT, nstp)])

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-1,#len(times),
        # timprs=times,
        prsity=0.35, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=[perlen/nstp, pumpT/nstp, recT/nstp],
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=500)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=recT/nstp,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    swt.write_input()

    return swt 

def cam_steady(modelname, Lz, Lx, Ly, nlay, nrow, ncol, perlen, dt, nstp, onshore_proportion, k_type, head=None, qinflow=None):
    
    model_ws = f".\\model_files\\{modelname}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = modelname+"_steady"
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")

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
            itmuni=4,
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=perlen, # COULD BE A POTENTIAL PROBLEM
            nstp=nstp,
            steady=True
        )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol*onshore_proportion):ncol] = -1
    strt = np.zeros((nlay, nrow, ncol))
    if k_type == "heterogenous":
        ibound[:, :, 0] = -1 

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laytyp[0] = 1
    laywet=np.zeros(50)
    laywet[0] = 1

    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet)
    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)

    oc_spd = {} 
    for kstp in range(nstp):
            oc_spd[(0, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
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
        if k_type == "heterogenous":
            ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["BAS6"]])
            chd_sp1.append([cell[0], cell[1], cell[2], head, head])
        else:
            ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["WEL"]])
            wel_sp1.append([cell[0], cell[1], cell[2], qinflow / nlay])

    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])
        
    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1
    wel_data[0] = wel_sp1

    if k_type == "homogenous":
        wel = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb=53)  

    chd = flopy.modflow.ModflowChd(swt, stress_period_data=chd_data, ipakcb = ipakcb)

    sconc = np.zeros((nlay, nrow, ncol))
    for col in range(int(ncol*onshore_proportion), ncol):
        for row in range(nrow): 
            for lay in range(nlay):
                sconc[lay, row, col] = 35.0

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-1,
        prsity=0.35, # can change this
        sconc=sconc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=dt
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=500)

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

def cam_transient(modelname, Lz, Lx, Ly, nlay, nrow, ncol, pumpT, recT, dt, onshore_proportion, k_type, q, qlay, qrow, qcol, recQmult, start_head, start_conc, head=None, qinflow=None):

    model_ws = f".\\model_files\\{modelname}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    modelname = modelname+"_transient"
    swt = flopy.seawat.Seawat(modelname, model_ws=model_ws, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")

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
            nper=2,
            itmuni=4,
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=[pumpT, recT], # COULD BE A POTENTIAL PROBLEM
            nstp=[pumpT/dt, recT/dt],
            steady=[False, False]
        )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, -1] = -1
    ibound[0, :, int(ncol*onshore_proportion):ncol] = -1
    strt = start_head
    if k_type == "heterogenous":
        ibound[:, :, 0] = -1 

    bas = flopy.modflow.ModflowBas(swt, ibound, strt=strt)

    laytyp=np.zeros(50)
    laytyp[0] = 1
    laywet=np.zeros(50)
    laywet[0] = 1

    lpf = flopy.modflow.ModflowLpf(swt, hk=0, vka=0, ipakcb=ipakcb, laytyp=laytyp, laywet=laywet)
    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-5, npcond=0, mxiter=500)

    oc_spd = {} 
    perlens = [pumpT, recT]
    for kper in range(2):
        steps = int(perlens[kper]/dt)
        for kstp in range(steps):
            oc_spd[(kper, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        for j in range(nrow):
            onshore_boundary_cells.append([k, j, 0])
            offshore_boundary_cells.append([k, j, ncol-1])

    for i in range(int(ncol*onshore_proportion), ncol):
        for j in range(nrow):
            offshore_boundary_cells.append([0, j, i])

    # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # wel data
    wel_data = {}
    wel_sp1 = []
    wel_sp2 = []
    # Set up ssm 
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []

    for cell in onshore_boundary_cells:
        if k_type == "heterogenous":
            ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["BAS6"]])
            chd_sp1.append([cell[0], cell[1], cell[2], head, head])
        else:
            ssm_sp1.append([cell[0], cell[1], cell[2], 0.0, itype["WEL"]])
            wel_sp1.append([cell[0], cell[1], cell[2], qinflow[0] / nlay])
            wel_sp2.append([cell[0], cell[1], cell[2], qinflow[1] / nlay])


    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], 0.0, 0.0])
        
    ssm_data[0] = np.append(ssm_sp1,[[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    ssm_data[1] = np.append(ssm_sp1,[[qlay, qrow, qcol, 0.0, itype['WEL']]], axis=0)
    chd_data[0] = chd_sp1
    chd_data[1] = chd_sp1
    wel_sp1.append([qlay, qrow, qcol, -q])
    wel_sp2.append([qlay, qrow, qcol, q*recQmult])
    wel_data[0] = wel_sp1
    wel_data[1] = wel_sp2

    wel = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb=53)  

    chd = flopy.modflow.ModflowChd(swt, stress_period_data=chd_data, ipakcb = ipakcb)

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-1,
        prsity=0.35, # can change this
        sconc=start_conc, # can change this: to each area having different starti,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=dt
    )

    adv = flopy.mt3d.Mt3dAdv(swt, 
        mixelm=0,
        dceps=1.0e-5,
        nplane=1,
        npl=16,
        nph=16,
        npmin=4,
        npmax=32,
        dchmoc=1.0e-3,
        nlsink=1,
        npsink=16,
        percel=0.5)
    # sip = flopy.modflow.ModflowSip(swt)
    # lmt = flopy.modflow.ModflowLmt(swt)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=2, cclose=1e-5)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=500)

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