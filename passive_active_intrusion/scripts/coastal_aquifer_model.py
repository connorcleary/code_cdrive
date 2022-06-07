import numpy as np
import flopy
import os
import pickle
import flopy.utils.binaryfile as bf

class ModelParameters:
    """
        Class to store model parameters

        Inputs:
            name: The name of the model/scenario/realization
            Lx: Length of the aquifer model [m]
            Ly: Alongshore dimension - set as 1 for the 2D case [m]
            Lz: Depth of the model [m]
            offshore_proportion: portion of the aquifer which is submarine 
            sea_level: datum for sea level (gives space for water table mound) [m]
            ncol: number of columns
            nrow: number of rows
            nlay: number of layers
            K: Horizontal hydraulic conductivity [m/day]
            anis: ratio between K and vertical hydraulic conductivity
            sy: specific yield 
            ss: specific storage [m^-1]
            n: porosity
            alpha_L: longitudinal dispersivity
            alpha_anisT: ratio between longitunal and transverse dispersivity
            alpha_anisV: ratio between longitudinal and vertical dispersivity
            diff: molecular diffusion coefficient
            perlen: simulation period length [days]
            dt: timestep [days]
            h_b: inland boundary head
            W_net: recharge [m\day]
            rho_f: fresh water density [kg/m^3]
            rho_s: salt water density [kg/m^3]
            exe_path: path to the seawat executable
    """
    
    def __init__(self, name, Lx=200, Ly=1, Lz=5.5, offshore_proportion=0.025, 
                sea_level=5, ncol=400, nrow=1, nlay=110, K=10, anis=1, sy=0.24, 
                ss=1e-5, n=0.3, alpha_L=1, alpha_anisT=0.1, alpha_anisV=0.1, 
                diff=8.64e-5, perlen=1e9, dt=1e6, h_b=0, W_net=0.00285,
                rho_f=1000, rho_s=1025, exe_path=r"C:\Users\ccl124\bin\swt_v4x64.exe"):

        self.name=name
        self.Lx=Lx
        self.Ly=Ly
        self.Lz=Lz
        self.offshore_proportion=offshore_proportion
        self.sea_level=sea_level
        self.ncol=ncol
        self.nrow=nrow
        self.nlay=nlay
        self.K=K
        self.anis=anis 
        self.sy=sy
        self.ss=ss
        self.n=n
        self.alpha_L=alpha_L
        self.alpha_anisT=alpha_anisT
        self.alpha_anisV=alpha_anisV
        self.diff=diff
        self.perlen=perlen
        self.dt=dt
        self.h_b=h_b
        self.W_net=W_net
        self.rho_f=rho_f
        self.rho_s=rho_s
        self.exe_path=exe_path
        self.save_parameters()

    def save_parameters(self):
        """
            Save object
        """
        model_ws = f".\\model_files\\{self.name}"
        if not os.path.exists(model_ws):
            os.makedirs(model_ws)

        f = open(f"{model_ws}\\pars", 'ab')
        pickle.dump(self, f)
        f.close()

    def load_parameters(self, name):

        model_ws = f".\\model_files\\{self.name}"
        f = open(f"{model_ws}\\pars", 'ab')
        self = pickle.load(f)


def build_steady_model(pars):
    '''
        A function to build a coastal aquifer model.

        Input: 
            pars: parameters object

        Outputs:
            swt: groundwater flow and transport model
    '''

    
        
    # create model workspace
    model_ws = f".\\model_files\\{pars.name}"
    if not os.path.exists(model_ws):
        os.makedirs(model_ws)

    # create base seawat model
    swt = flopy.seawat.Seawat(pars.name, model_ws=model_ws, exe_name=pars.exe_path)

    # calculate cell dimension
    delr = pars.Lx/pars.ncol
    delc = pars.Ly/pars.nrow
    delv = pars.Lz/pars.nlay

    # define top and bottoms
    top = pars.Lz
    botm = np.linspace(top-delv, 0, pars.nlay)

    # something I've copied
    ipakcb = 53

    # define discretization package
    dis = flopy.modflow.ModflowDis(
            model=swt,
            nlay=pars.nlay,
            nrow=pars.nrow,
            ncol=pars.ncol,
            nper=1,
            itmuni=4, # four for days
            delr=delr,
            delc=delc,
            laycbd=0,
            top=top,
            botm=botm,
            perlen=pars.perlen,
            nstp=pars.perlen/pars.dt,
            steady=True
        )

    # define cell groups
    inactive_cells = []
    offshore_boundary_cells = []
    onshore_boundary_cells = []
    #  surface_boundary_cells = []

    # add inactive cells
    for i in range(int(pars.ncol*pars.offshore_proportion)):
        for j in range(pars.nrow):
            for k in range(0, int((pars.Lz-pars.sea_level)/delv)-1):
                inactive_cells.append([k, j, i])
        
    # add cells on ends of domain
    for k in range(pars.nlay):
        for j in range(pars.nrow):
            if k >= (pars.Lz-pars.sea_level)/delv: # if the cell is below sea level
                offshore_boundary_cells.append([k, j, 0])

            if k >= (pars.Lz-pars.h_b)/delv:
                onshore_boundary_cells.append([k, j, pars.ncol-1])

    # add the seafloor
    for i in range(int(pars.ncol*pars.offshore_proportion)):
         for j in range(pars.nrow):
            offshore_boundary_cells.append([int((pars.Lz-pars.sea_level)/delv), j, i])

    # # add the recharge surface 
    # for i in range(int(offshore_proportion*ncol), ncol):
    #     for j in range(nrow):
    #         surface_boundary_cells.append([0, j, i])

    # create ibound array
    ibound = np.ones((pars.nlay, pars.nrow, pars.ncol), dtype=np.int32)
    for cell in inactive_cells:
        ibound[cell[0], cell[1], cell[2]] = 0
    for cell in onshore_boundary_cells+offshore_boundary_cells:
        ibound[cell[0], cell[1], cell[2]] = -1

    # define starting heads
    strt = pars.sea_level*np.ones((pars.nlay, pars.nrow, pars.ncol))

    # create basic package
    bas = flopy.modflow.ModflowBas(
            model=swt, 
            ibound=ibound, 
            strt=strt
        )

    # define layer types and wetting, this one I'm not sure about
    laytyp=np.zeros(pars.nlay)
    laytyp[0] = 1
    laywet=np.zeros(pars.nlay)
    laywet[0] = 1

    # create layer property flow package
    lpf = flopy.modflow.ModflowLpf(
            swt, 
            hk=pars.K, 
            vka=pars.anis, 
            ipakcb=ipakcb, 
            laytyp=laytyp, 
            laywet=laywet,
            ss=pars.ss, # not sure about these ones
            sy=pars.sy
        )

    # create solver package
    pcg = flopy.modflow.ModflowPcg(
            swt, 
            hclose=1.0e-5, 
            npcond=1, 
            mxiter=500
        )
    

    oc_spd = {} 
    for kstp in range(int(pars.perlen/pars.dt)):
            oc_spd[(0, kstp)] = ["save head", "save budget"]

    oc = flopy.modflow.ModflowOc(swt, stress_period_data=oc_spd, compact=True)

    # set up constant head stress period data
    chd_data = {}
    chd_sp1 = []
    # Set up source sink data
    itype = flopy.mt3d.Mt3dSsm.itype_dict()
    ssm_data = {}
    ssm_sp1 = []

    # define onshore boundary data
    for cell in onshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], pars.h_b, pars.h_b])

    # define offshore boundary data
    for cell in offshore_boundary_cells:
        ssm_sp1.append([cell[0], cell[1], cell[2], 35.0, itype["BAS6"]])
        chd_sp1.append([cell[0], cell[1], cell[2], pars.sea_level, pars.sea_level])
        
    ssm_data[0] = ssm_sp1
    chd_data[0] = chd_sp1

    # define constant head package
    chd = flopy.modflow.ModflowChd(
            model=swt, 
            stress_period_data=chd_data, 
            ipakcb = ipakcb
        )

    # create recharge package
    rch = flopy.modflow.ModflowRch(
            model=swt,
            rech=pars.W_net
        )

    # set starting concentrations
    sconc = 35*np.ones((pars.nlay, pars.nrow, pars.ncol))

    # define basic transport package
    btn = flopy.mt3d.Mt3dBtn(
            swt,
            nprs=-1,
            prsity=pars.n,
            sconc=sconc,
            chkmas=False,
            nprobs=10,
            nprmas=10,
            dt0=pars.dt
        )

    # define advection package
    adv = flopy.mt3d.Mt3dAdv(
            model=swt, 
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
            percel=0.5
        )

    # define dispersion package
    dsp = flopy.mt3d.Mt3dDsp(
            swt, 
            al=pars.alpha_L, 
            trpt=pars.alpha_anisT, 
            trpv=pars.alpha_anisV, 
            dmcoef=pars.diff,
        )

    # define transport solver
    gcg = flopy.mt3d.Mt3dGcg(
            model=swt, 
            iter1=500, 
            mxiter=1, 
            isolve=2, 
            cclose=1e-5
        )

    # find number of sinks and sources
    mxss = int(np.ceil(2*pars.nlay*pars.nrow + 
                        pars.nrow*pars.ncol*pars.offshore_proportion+1 +
                        pars.nrow*pars.ncol))

    # define source sink mixing package
    ssm = flopy.mt3d.Mt3dSsm(
            model=swt,
            stress_period_data=ssm_data, 
            mxss=mxss
        )

    vdf = flopy.seawat.SeawatVdf(
            swt,
            iwtable=0,
            densemin=0,
            densemax=0,
            denseref=1000.0,
            denseslp=0.7143,
            firstdt=pars.dt,
        )
    # write input
    swt.write_input() 

    return swt


def run_model(swt):
    """
        A function to run the seawat model

        Inputs: 
            swt: model object
        Outputs:
            None
    """
    swt.write_input()
    success, buff = swt.run_model(silent=False, report=True)
    if not success:
        raise Exception("SEAWAT did not terminate normally.")

def extract_results(pars, only_final=True):
    """
        Open model results from binary files

        Inputs:
            swt: model object
            only_final: flag to only save the final timestep
        Outputs:
            head: head matrix [nstp, nlay, nrow, ncol]
            qx: longitudinal flux matrix [nstp, nlay, nrow, ncol]
            qy: transverse flux matrix matrix [nstp, nlay, nrow, ncol]
            qz: vertical flux matrix matrix [nstp, nlay, nrow, ncol]
            concentration: concentration matrix [nstp, nlay, nrow, ncol]
    """
    name = pars.name
    model_ws = f".\\model_files\\{name}"

    if only_final:
        nstp = 1
    else:
        nstp = int(pars.perlen/pars.dt)

    head = np.zeros(pars.nstep, pars.nlay, pars.nrow, pars.ncol)
    qx = np.zeros(pars.nstep, pars.nlay, pars.nrow, pars.ncol)
    qy = np.zeros(pars.nstep, pars.nlay, pars.nrow, pars.ncol)
    qz = np.zeros(pars.nstep, pars.nlay, pars.nrow, pars.ncol)
    concentration = np.zeros(pars.nstep, pars.nlay, pars.nrow, pars.ncol)

    # open binary files
    ucnobj = bf.UcnFile(os.path.join(model_ws, "MT3D001.UCN"))
    cbbobj = bf.CellBudgetFile(os.path.join(model_ws, f'{name}.cbc'))
    headobj = bf.HeadFile(os.path.join(model_ws, f'{name}.hds'))


pars = ModelParameters("test")
swt = build_steady_model(pars)
run_model(swt)


