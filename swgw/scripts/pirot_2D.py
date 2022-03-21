import numpy as np
import flopy
import numpy as np
import flopy.utils.binaryfile as bf
import matplotlib.pyplot as plt
from scipy.io import loadmat
import os

# parameters

Lx = 800.0
Ly = 10.0
Lz = 25.0
nlay = 50
nrow = 1
ncol = 80
delr = Ly/nrow
delc = Lx/ncol
delv = Lz/nlay
top = 0.0
botm = np.linspace(top - delv, -Lz, nlay)
qinflow = 1e3

def build_model():
    modelname = "pirot2d"
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
        top=top,
        botm=botm,
        perlen=1500000000, # COULD BE A POTENTIAL PROBLEM
        nstp=15
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
        stress_period_data={(0, 14): ["save head", "save budget"]},
        compact=True,
    )

    # find inland boundary cells 
    onshore_boundary_cells = []
    offshore_boundary_cells = []
    for k in range(nlay):
        onshore_boundary_cells.append([k, 0, 0])
        offshore_boundary_cells.append([k, 0, ncol-1])

    for i in range(int(ncol/4), ncol):
        offshore_boundary_cells.append([0, 0, i])

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
        chd_sp1.append([cell[0], cell[1], cell[2], 0.6, 0.6])
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
        dt0=1e-3*1e10
    )

    adv = flopy.mt3d.Mt3dAdv(swt, mixelm=0)
    dsp = flopy.mt3d.Mt3dDsp(swt, al=0.4, trpt=0.1, trpv=0.01, dmcoef=1e-9)
    gcg = flopy.mt3d.Mt3dGcg(swt, iter1=500, mxiter=1, isolve=1, cclose=1e-7)
    ssm = flopy.mt3d.Mt3dSsm(swt, stress_period_data=ssm_data, mxss=10000)

    vdf = flopy.seawat.SeawatVdf(
        swt,
        iwtable=0,
        densemin=0,
        densemax=0,
        denseref=1000.0,
        denseslp=0.7143,
        firstdt=1e-3*1e10,
    )

    fname = r"./MT3D001.UCN"
    if os.path.isfile(fname):
        os.remove(fname)

    return swt 

def load_field():
    fields = loadmat(r'./fields/80x40x50aTest')
    Kh = fields['model_Kh']
    Kv = fields['model_Kv']

    return Kh, Kv

def change_model(swt, Kh, Kv):

    Kh_big = np.zeros((nlay, ncol))
    Kv_big = np.zeros((nlay, ncol))
    # for i in range(4):
    #     for j in range(4):
    #         Kh_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kh, (2, 0, 1))[:,(10*i+j)%40,:]
    #         Kv_big[i*50:i*50+50,j*80:j*80+80] = np.transpose(Kv, (2, 0, 1))[:,(10*i+j)%40,:]

    swt.lpf.hk = np.transpose(Kh, (2, 0, 1))[:,0,:]
    swt.lpf.vka = np.transpose(Kv, (2, 0, 1))[:,0,:]

    # swt.lpf.hk = Kh_big
    # swt.lpf.vka = Kv_big

    return swt


def run_model(swt):
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def load_model():
    swt = flopy.seawat.swt.Seawat.load(r'.\pirot2d.nam',  exe_name=r"C:\Users\ccl124\bin\swt_v4x64.exe")
    return swt

def post_proc(swt, par_name, par_value, m, n):
    ucnobj = bf.UcnFile("MT3D001.UCN", model=swt)
    times = ucnobj.get_times()
    concentration = ucnobj.get_data(totim=times[-1])

    cbbobj = bf.CellBudgetFile("pirot2d.cbc")
    times = cbbobj.get_times()
    qx = cbbobj.get_data(text="flow right face", totim=times[-1])[0]
    #qy = cbbobj.get_data(text="flow front face", totim=times[-1])[0]
    qy = np.zeros((nlay, nrow, ncol), dtype=float)
    qz = cbbobj.get_data(text="flow lower face", totim=times[-1])[0]

    fig = plt.figure(figsize=(12,9))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(concentration)
    pmv.plot_vector(qx, qy, -qz, color="white", kstep=3, hstep=3)
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Concentrations");
    plt.savefig(f'.\\results\\pirot2d\\concentration_{par_name}_{par_value}.png')

    headobj = bf.HeadFile("pirot2d.hds")
    times = headobj.get_times()
    head = headobj.get_data(totim=times[-1])

    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    ax.set_aspect(10)
    pmv = flopy.plot.PlotCrossSection(model=swt, ax=ax, line={"row": 0})
    arr = pmv.plot_array(head)
    contours = pmv.contour_array(head, colors="white")
    ax.clabel(contours, fmt="%2.2f")
    plt.colorbar(arr, shrink=0.5, ax=ax)
    ax.set_title("Simulated Heads");
    plt.savefig(f'.\\results\\pirot2d\\head_{par_name}_{par_value}.png')

def run_pirot():

    swt = build_model()

    # swt = load_model()
    Kh, Kv = load_field()
    swt = change_model(swt, Kh, Kv)

    swt.write_input()

    run_model(swt)

    #swt = load_model()
    post_proc(swt, 'none', 'yet', 1, 1)
    pass

if __name__ == "__main__":
    run_pirot()

