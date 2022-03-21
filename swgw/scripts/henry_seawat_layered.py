import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
import os

ref_hk = 864.0 # m/day
max_hk = 9000.0 # m/day
modelname = "henry_layered"
Lx = 2.0
Lz = 1.0
nlay = 50
nrow = 1
ncol = 100
delr = Lx / ncol
delc = 1.0
delv = Lz / nlay
henry_top = 1.0
henry_botm = np.linspace(henry_top - delv, 0.0, nlay)
qinflow = 5.702  # m3/day
dmcoef = 0.57024  # m2/day  Could also try 1.62925 as another case of the Henry problem

def build_model():
    ''' Build base model excluding hk '''

    model_ws = os.path.join(os.getcwd(), "henry_layered")
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)
        
    swt = flopy.seawat.Seawat(modelname, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe", model_ws=model_ws)
    print(swt.namefile)

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
        top=henry_top,
        botm=henry_botm,
        perlen=1.5,
        nstp=15,
    )

    # Variables for the BAS package
    ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    ibound[:, :, -1] = -1

    bas = flopy.modflow.ModflowBas(swt, ibound, 0)

    lpf = flopy.modflow.ModflowLpf(swt, ipakcb=ipakcb)

    pcg = flopy.modflow.ModflowPcg(swt, hclose=1.0e-8)

    oc = flopy.modflow.ModflowOc(
        swt,
        stress_period_data={(0, 0): ["save head", "save budget"]},
        compact=True,
    )

    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    wel_data = {}
    ssm_data = {}
    wel_sp1 = []
    ssm_sp1 = []
    for k in range(nlay):
        wel_sp1.append([k, 0, 0, qinflow / nlay])
        ssm_sp1.append([k, 0, 0, 0.0, itype["WEL"]])
        ssm_sp1.append([k, 0, ncol - 1, 35.0, itype["BAS6"]])
    wel_data[0] = wel_sp1
    ssm_data[0] = ssm_sp1
    wel = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb=ipakcb)

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-5,
        prsity=0.35,
        sconc=35.0,
        ifmtcn=0,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=0.001,
    )
    adv = flopy.mt3d.Mt3dAdv(swt, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.0, trpt=1.0, trpv=1.0, dmcoef=dmcoef)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=1, cclose=1e-7)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=1e-3,
    )
    return swt


def change_hk(swt, hk):
    # change hk
    swt.lpf.hk = [ref_hk]*25 +[hk]*12+[ref_hk]*13 
    swt.lpf.vka = [ref_hk]*25 +[hk]*12+[ref_hk]*13

    return swt

def run_model(swt):
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def extract_results(swt, hk):
    ucnobj = bf.UcnFile("henry_layered\\MT3D001.UCN", model=swt)
    times = ucnobj.get_times()
    concentration = ucnobj.get_data(totim=times[-1])
    cbbobj = bf.CellBudgetFile("henry_layered\\henry_layered.cbc")
    times = cbbobj.get_times()
    qx = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
    qy = np.zeros((nlay, nrow, ncol), dtype=float)
    qz = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]
    headobj = bf.HeadFile("henry_layered\\henry_layered.hds")
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    results_ws = os.path.join(os.getcwd(), "henry_layered", "results")
    if not os.path.exists(results_ws): os.makedirs(results_ws)

    np.savetxt(results_ws + f"\\concentration_{hk:08.2f}.csv", np.squeeze(concentration))
    np.savetxt(results_ws + f"\\qx_{hk:08.2f}.csv", np.squeeze(qx))
    np.savetxt(results_ws + f"\\qy_{hk:08.2f}.csv", np.squeeze(qy))
    np.savetxt(results_ws + f"\\qz_{hk:08.2f}.csv", np.squeeze(qz))
    np.savetxt(results_ws + f"\\head_{hk:08.2f}.csv", np.squeeze(head))

def main():
    swt = build_model()

    for hk in np.linspace(ref_hk, max_hk, 2):
        swt = change_hk(swt, hk)
        run_model(swt)
        extract_results(swt, hk)

    ucnobj = bf.UcnFile("henry_layered\\MT3D001.UCN", model=swt)
    times = ucnobj.get_times()
    concentration = ucnobj.get_data(totim=times[-1])

    cbbobj = bf.CellBudgetFile("henry_layered\\henry_layered.cbc")
    times = cbbobj.get_times()
    qx = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
    qy = np.zeros((nlay, nrow, ncol), dtype=float)
    qz = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]

    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(concentration)
    pmv.plot_vector(qx, qy, -qz, color="white", kstep=3, hstep=3)
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Concentrations");
    plt.show()

if __name__ == "__main__":
    main()
