import os
import sys
import matplotlib.pyplot as plt
import flopy
import numpy as np

mf6exe = r"C:\Users\ccl124\bin\mf6.exe"
figure_size = (6, 4)

ws = os.path.join(os.getcwd(), 'data')
if not os.path.exists(ws):
    os.makedirs(ws)

parameters = {
    "ex-gwt-henry-a": {"inflow": 5.7024,},
    "ex-gwt-henry-b": {"inflow": 2.851,},
}



parameter_units = {
    "inflow": "$m^3/d$",
}

length_units = "cm"
time_units = "seconds"

nper = 1  # Number of periods
nstp = 500  # Number of time steps
perlen = 0.5  # Simulation time length ($d$)
nlay = 40  # Number of layers
nrow = 40  # Number of rows
ncol = 80  # Number of columns
system_length = 2.0  # Length of system ($m$)
delr = 0.025  # Column width ($m$)
delc = 0.025 # Row width ($m$)
delv = 0.025  # Layer thickness
top = 0.1  # Top of the model ($m$)
hydraulic_conductivity = 864.0  # Hydraulic conductivity ($m d^{-1}$)
initial_concentration = 35.0  # Initial concentration (unitless)
porosity = 0.35  # porosity (unitless)
diffusion_coefficient = 0.57024  # diffusion coefficient ($m^2/d$)

botm = [top - k * delv for k in range(1, nlay + 1)]

nouter, ninner = 100, 300
hclose, rclose, relax = 1e-10, 1e-6, 0.97

def build_model(sim_folder, inflow):
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
    flopy.mf6.ModflowGwfic(gwf, strt=initial_concentration)
    pd = [(0, 0.7, 0.0, "trans", "concentration")]
    flopy.mf6.ModflowGwfbuy(gwf, packagedata=pd)
    ghbcond = hydraulic_conductivity * delv * delc / (0.5 * delr)
    ghbspd = [[(lay, row, ncol - 1), top, ghbcond, 35.0] for row in range(nrow) for lay in range(nlay)]
    flopy.mf6.ModflowGwfghb(
        gwf,
        stress_period_data=ghbspd,
        pname="GHB-1",
        auxiliary="CONCENTRATION",
    )

    welspd = [[(lay, row, 0), inflow / nlay / nrow, 0.0] for row in range(nrow) for lay in range(nlay)]
    flopy.mf6.ModflowGwfwel(
        gwf,
        stress_period_data=welspd,
        pname="WEL-1",
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
        ("WEL-1", "AUX", "CONCENTRATION"),
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
    return sim

def write_model(sim, silent=True):
    sim.write_simulation(silent=silent)
    return

def run_model(sim, silent=True):
    success = True
    success, buff = sim.run_simulation(silent=silent)
    if not success:
        print(buff)
    return success

def plot_conc(sim, idx):
    sim_name = list(parameters.keys())[idx]
    sim_ws = os.path.join(ws, sim_name)
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
    ax.set_aspect('equal')
    # save figure

    fpth = os.path.join(
        ".", "figures", "{}-conc{}".format(sim_name, "1")
    )
    plt.show()
    fig.savefig(fpth)
    return

def plot_results(sim, idx):

    plot_conc(sim, idx)
    return

def scenario(idx, silent=False):
    key = list(parameters.keys())[idx]
    parameter_dict = parameters[key]
    sim = build_model(key, **parameter_dict)
    write_model(sim, silent=silent)
    success = run_model(sim, silent=silent)
    if success:
        plot_results(sim, idx)

def main():
    scenario(0)

if __name__ == "__main__":
    main()