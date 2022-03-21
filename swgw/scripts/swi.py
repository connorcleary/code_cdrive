import os
import sys
import matplotlib.pyplot as plt
import flopy
import numpy as np

mf6exe = r"C:\Users\ccl124\bin\mf6.exe"
figure_size = (6, 4)

length_units = "m"
time_units = "seconds"

nper = 1  # Number of periods
nstp = 500  # Number of time steps
perlen = 60*60*24*365  # Simulation time length ($s$)
   
xb = 3500 # [L] metres 
z0 = 20 
hydraulic_conductivity = 10 * 1/(60*60*24) # convert m/d to m/s
Wnet = 20 * 1/1000 * 1/(60*60*24*365) # convert mm/yr to m/s

Lx = xb + 1000 
Ly = z0 + 10

nlay = Ly*2  # Number of layers
nrow = 1  # Number of rows
ncol = int(Lx/50)  # Number of columns
delr = 50  # Column width ($m$)
delc = 50 # Row width ($m$)
delv = 0.5  # Layer thickness
top = 10  # Top of the model ($m$)

botm = [top - k * delv for k in range(1, nlay + 1)]

initial_concentration = 35.0  # Initial concentration (unitless)
porosity = 0.3  # porosity (unitless)
diffusion_coefficient = 10e-9  # diffusion coefficient ($m^2/s$)

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-6, 0.97

shore_cells = [[[0, 0, int(ncol/2)]]]
for lay in range(1, int(nlay/2)): 
    shore_cells.append([[lay, 0, col] for col in range(shore_cells[-1][-1][-1]+1, int(ncol/2+(ncol*(lay+1))/nlay))])

def build_model(sim_folder):

    ws = os.path.join(os.getcwd(), 'data')
    if not os.path.exists(ws):
        os.makedirs(ws)

    print("Building model...{}".format(sim_folder))
    name = "flow"
    sim_ws = os.path.join(ws, sim_folder)
    sim = flopy.mf6.MFSimulation(
        sim_name=name, sim_ws=sim_ws, exe_name=mf6exe
    )

    tdis_ds = ((perlen, nstp, 1.0),)
    flopy.mf6.ModflowTdis(
        sim, nper=nper, perioddata=tdis_ds, time_units=time_units
    )
    gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)
    ims = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwf.name),
    )
    sim.register_ims_package(ims, [gwf.name])

    idomain = np.ones((nlay,nrow,ncol), dtype=int)
    for lay, cells in enumerate(shore_cells): 
        idomain[lay][0][cells[-1][-1]+1:ncol] = [0]*(ncol - cells[-1][-1] - 1)
        # get the axis right

    flopy.mf6.ModflowGwfdis(
        gwf,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    ) 

    flopy.mf6.ModflowGwfnpf(
        gwf,
        save_specific_discharge=True,
        icelltype=0,
        k=hydraulic_conductivity,
    )

    # drn_cells_shore = [cell for lay in shore_cells[:int(nlay/4)+1] for cell in lay]
    # drn_spd = [[(cell[0], cell[1], cell[2]), 10-cell[0]*delv, 10000, -0.25] 
    #                 for cell in drn_cells_shore]
    # flopy.mf6.ModflowGwfdrn(
    #     gwf,
    #     pname="DRN-1",
    #     auxiliary=["ddrn"],
    #     auxdepthname="ddrn",
    #     print_input=False,
    #     print_flows=False,
    #     mover=True,
    #     stress_period_data=drn_spd,  # wel_spd established in the MVR setup
    #     boundnames=False,
    #     save_flows=True
    # )

    flopy.mf6.ModflowGwfic(gwf, strt=initial_concentration)
    pd = [(0, 0.7, 0.0, "trans", "concentration")]

    ghb_cells_shore = [cell for lay in shore_cells[int(nlay/4)+1:] for cell in lay]

    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)
    ghbcond = hydraulic_conductivity * delv * delc / (0.5 * delr)
    ghbspd = [[(cell[0], cell[1], cell[2]), 10, ghbcond, 35.0] for cell in ghb_cells_shore]
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        pname="GHB-1",
        auxiliary="CONCENTRATION",
    )

    rch_cells_top = [[0,0,col] for col in range(int(ncol/2))]
    rch_cells_shore = [cell for lay in shore_cells[:int(nlay/4)+1] for cell in lay]
    rch_cells = np.concatenate((rch_cells_top, rch_cells_shore))
    # rch_cells = rch_cells_top
    rch_spd = [[(cell[0], cell[1], cell[2]), Wnet, 0.0] for cell in rch_cells] 
    flopy.mf6.ModflowGwfrch(
        gwf, 
        stress_period_data=rch_spd,
        pname="RCH-1",
        auxiliary="CONCENTRATION",
    )

    head_filerecord = "{}.hds".format(name)
    budget_filerecord = "{}.bud".format(name)
    flopy.mf6.ModflowGwfoc(
        gwf,
        head_filerecord=head_filerecord,
        budget_filerecord=budget_filerecord,
        saverecord=[("HEAD", "ALL"), ("BUDGET", "ALL")],
    )

    gwt = flopy.mf6.ModflowGwt(sim, modelname="trans")
    imsgwt = flopy.mf6.ModflowIms(
        sim,
        print_option="ALL",
        outer_dvclose=hclose,
        outer_maximum=nouter,
        under_relaxation="NONE",
        inner_maximum=ninner,
        inner_dvclose=hclose,
        rcloserecord=rclose,
        linear_acceleration="BICGSTAB",
        scaling_method="NONE",
        reordering_method="NONE",
        relaxation_factor=relax,
        filename="{}.ims".format(gwt.name),
    )
    sim.register_ims_package(imsgwt, [gwt.name])
    flopy.mf6.ModflowGwtdis(
        gwt,
        length_units=length_units,
        nlay=nlay,
        nrow=nrow,
        ncol=ncol,
        delr=delr,
        delc=delc,
        top=top,
        botm=botm,
    )
    flopy.mf6.ModflowGwtmst(gwt, porosity=porosity)
    flopy.mf6.ModflowGwtic(gwt, strt=initial_concentration)
    flopy.mf6.ModflowGwtadv(gwt, scheme="UPSTREAM")
    flopy.mf6.ModflowGwtdsp(gwt, xt3d_off=True, diffc=diffusion_coefficient)
    sourcerecarray = [
        ("GHB-1", "AUX", "CONCENTRATION"),
        ("RCH-1", "AUX", "CONCENTRATION")
    ]
    flopy.mf6.ModflowGwtssm(gwt, sources=sourcerecarray)
    flopy.mf6.ModflowGwtoc(
        gwt,
        budget_filerecord="{}.cbc".format(gwt.name),
        concentration_filerecord="{}.ucn".format(gwt.name),
        concentrationprintrecord=[
            ("COLUMNS", 10, "WIDTH", 15, "DIGITS", 6, "GENERAL")
        ],
        saverecord=[("CONCENTRATION", "ALL")],
        printrecord=[("CONCENTRATION", "LAST"), ("BUDGET", "LAST")],
    )
    flopy.mf6.ModflowGwfgwt(
        sim, exgtype="GWF6-GWT6", exgmnamea=gwf.name, exgmnameb=gwt.name
    )

    fig, ax = plt.subplots(1, 1, figsize=(9, 3), constrained_layout=True)
    # first subplot
    ax.set_title("Cross section")
    modelmap = flopy.plot.PlotCrossSection(
        model=gwf,
        ax=ax,
        line={"row": 0}
    )

    inactive = modelmap.plot_inactive(idomain)
    rch = modelmap.plot_bc(name="RCH-1", color='r')
    ghb = modelmap.plot_bc(name="GHB-1", color='c')
    # drn = modelmap.plot_bc(name="DRN-1", color='y')
    plt.show()

    return sim

def run_model(sim, silent=True):
    success = True
    success, buff = sim.run_simulation(silent=silent)
    if not success:
        print(buff)
    return success

def plot_conc(sim, name):
    sim_name = name 
    gwf = sim.get_model("flow")
    gwt = sim.get_model("trans")

    fig = plt.figure(figsize=figure_size)
    fig.tight_layout()

    # get MODFLOW 6 concentration
    conc = gwt.output.concentration().get_data()

    ax = fig.add_subplot(1, 1, 1, aspect="equal")
    pxs = flopy.plot.PlotCrossSection(model=gwf, ax=ax, line={"row": 0})
    pxs.plot_array(conc, cmap="jet")
    levels = [35 * f for f in [0.01, 0.1, 0.5, 0.9, 0.99]]
    cs = pxs.contour_array(
        conc, levels=levels, colors="w", linewidths=1.0, linestyles="-"
    )
    ax.set_xlabel("x position (m)")
    ax.set_ylabel("z position (m)")
    plt.clabel(cs, fmt="%4.2f", fontsize=5)
    ax.set_aspect(100)
    # save figure

    fpth = os.path.join(
        ".", "figures", "{}-conc{}".format(sim_name, "1")
    )
    fig.savefig(fpth)
    return

def main():
    name = "swi_test"
    sim = build_model(name)
    sim.write_simulation()
    success = run_model(sim, silent=False)
    plot_conc(sim, name)

if __name__ == "__main__":
    main()

