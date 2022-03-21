import os
import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt


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

qinflow = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 5.702)  # m3/day
dmcoef = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 0.57024)  # m2/day  Could also try 1.62925 as another case of the Henry problem
hk = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 864.0)  # m/day
porosity = np.multiply([0.8, 0.9, 1.0, 1.1, 1.2], 0.35)

extent = np.zeros((4, 5))
mixing_zone_size = np.zeros((4, 5))
outflow = np.zeros((4, 5))
inflow= np.zeros((4, 5))

def build_base():
    modelname = "henry"
    swt = flopy.seawat.Seawat(modelname, exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
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

    lpf = flopy.modflow.ModflowLpf(swt, hk=hk[2], vka=hk[2], ipakcb=ipakcb)

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
        wel_sp1.append([k, 0, 0, qinflow[2] / nlay])
        ssm_sp1.append([k, 0, 0, 0.0, itype["WEL"]])
        ssm_sp1.append([k, 0, ncol - 1, 35.0, itype["BAS6"]])
    wel_data[0] = wel_sp1
    ssm_data[0] = ssm_sp1
    wel = flopy.modflow.ModflowWel(swt, stress_period_data=wel_data, ipakcb=ipakcb)

    btn = flopy.mt3d.Mt3dBtn(
        swt,
        nprs=-5,
        prsity=porosity[2],
        sconc=35.0,
        ifmtcn=0,
        chkmas=False,
        nprobs=10,
        nprmas=10,
        dt0=0.001,

    )
    adv = flopy.mt3d.Mt3dAdv(swt, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.0, trpt=1.0, trpv=1.0, dmcoef=dmcoef[2])
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

def change(swt, i, bqinflow=False, bdmcoef=False, bhk=False, bporosity=False):

    
    swt.lpf.hk = hk[i]*bhk + hk[2]*int(not(bhk))
    swt.lpf.vha = hk[i]*bhk + hk[2]*int(not(bhk))
    swt.wel.stress_period_data[0]['flux'] = np.transpose(len(swt.wel.stress_period_data[0])*[(qinflow[i]*bqinflow+qinflow[2]*int(not(bqinflow)))/nlay])
    swt.btn.porosity = porosity[i]*bporosity + porosity[2]*int(not(bporosity))
    flopy.mt3d.Mt3dDsp(swt, al=0.0, trpt=1.0, trpv=1.0, dmcoef=dmcoef[i]*bdmcoef + dmcoef[2]*int(not(bdmcoef)))

    return swt

def run_model(swt):
    swt.write_input()
    success, buff = swt.run_model(silent=True, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def post_proc(swt, par_name, par_value, m, n):
    ucnobj = bf.UcnFile("MT3D001.UCN", model=swt)
    times = ucnobj.get_times()
    concentration = ucnobj.get_data(totim=times[-1])

    cbbobj = bf.CellBudgetFile("henry.cbc")
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
    plt.savefig(f'.\\results\\sensitivity\\concentration_{par_name}_{par_value}.png')

    headobj = bf.HeadFile("henry.hds")
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(head)
    contours = pmv.contour_array(head, colors="white")
    ax.clabel(contours, fmt="%2.2f")
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Heads");
    plt.savefig(f'.\\results\\sensitivity\\head_{par_name}_{par_value}.png')

    # calculate extent
    # find toe --> concentration greater than 0.5 at base
    for i in range(100):
        if concentration[-1][0][i] > 0.5:
            extent[m, n] = 0.02*i + 0.01
            break

    # calculate size of mixing zone
    # --> area between 0.1 and 0.9 salinity
    mix = 0
    for i in range(concentration.shape[0]):
        for j in range(concentration.shape[2]):
            if 3.5 <= concentration[i][0][j] <= 31.5:
                mix += 1
    mixing_zone_size[m, n] = 2*mix/5000

    # calculate outflow
    # --> positive flows on the rhs
    # calculate inflow
    # --> negative flows on the rhs

    end = qx[:, 0, -2]
    for i in range(len(end)):
        if end[i] > 0:
            outflow[m,n] += end[i]
        else:
            inflow[m,n] += end[i]
 
def main():
    # build model
    params = [qinflow, dmcoef, hk, porosity]
    pbool = [False, False, False, False]
    pnames = ['qinflow', 'dmcoef', 'hk', 'porosity']
    swt = build_base()
    # loop through configurations
        # change
        # run 
        # post proc
    for m in range(4):
        pbool[m] = True
        for n in range(5):
            swt = change(swt, n, pbool[0], pbool[1], pbool[2], pbool[3])
            run_model(swt)
            post_proc(swt, pnames[m], params[m][n], m, n)
        pbool[m] = False
    # save csvs
    pass
    np.savetxt(r'.\results\sensitivity\extent.csv', extent, delimiter=',')
    np.savetxt(r'.\results\sensitivity\mixing_zone_size.csv', mixing_zone_size, delimiter=',')
    np.savetxt(r'.\results\sensitivity\outflow.csv', outflow, delimiter=',')
    np.savetxt(r'.\results\sensitivity\inflow.csv', inflow, delimiter=',')
    # plot changes

if __name__ == "__main__":
    main()