#!/usr/bin/env python
import flopy

import numpy as np
import matplotlib.pyplot as plt


def get_chd_stress_period():
    for i in range(60):
        yield ((0, i, 0), 3820)
    for i in range(60):
        yield ((0, i, 39), 3824)


ws = './models/2_simple_model'
name = '2_simple_model'

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='modflow-mf6')

tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

dis = flopy.mf6.ModflowGwfdis(gwf,
                              nlay=10,
                              nrow=60,
                              ncol=40,
                              delc=50,
                              delr=50,
                              top=3832,
                              botm=np.linspace(3832, 3600, 11)[1:])

initial_head = np.ones((10, 60, 40)) * 3832
for gp, head in get_chd_stress_period():
    initial_head[gp] = head
ic = flopy.mf6.ModflowGwfic(gwf, strt=initial_head)
ic = flopy.mf6.ModflowGwfic(gwf)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0055)
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              k=4.0,
                              save_specific_discharge=True)
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_chd_stress_period()))
budget_file = name + '.bud'
head_file = name + '.hds'
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord=budget_file,
                            head_filerecord=head_file,
                            saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')])
sim.write_simulation()
sim.run_simulation()

head_arr = gwf.output.head().get_data()
bud = gwf.output.budget()



spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)

pmv = flopy.plot.PlotMapView(gwf)
pmv.plot_array(head_arr)
pmv.plot_grid(colors='white', linewidths=0.3)
pmv.contour_array(head_arr, linewidths=1., c_label=True, cmap='Wistia')
# flopy.plot.styles.graph_legend()
pmv.plot_vector(qx, qy, normalize=True, color="white")
plt.savefig("./images/01_00_plan.png")

plt.show()


def plot_x_section(**kwargs):
    fig, ax = plt.subplots(1, 1, figsize=(9, 3), constrained_layout=True)
    # first subplot
    title_text = "; ".join((f'{k}={v}' for k, v in kwargs.items()))
    ax.set_title(f"X-Section ({title_text})")
    modelmap = flopy.plot.PlotCrossSection(
        model=gwf,
        ax=ax,
        line=kwargs,
    )
    pa = modelmap.plot_array(head_arr, vmin=3600, vmax=3832)
    quadmesh = modelmap.plot_bc("CHD")
    linecollection = modelmap.plot_grid(lw=0.2, color="white")
    minor_contours = modelmap.contour_array(
        head_arr,
        levels=np.arange(3600, 3832, .1),
        linewidths=0.2,
        colors='black'
    )
    contours = modelmap.contour_array(
        head_arr,
        levels=np.arange(3600, 3832, .5),
        linewidths=0.8,
        colors='black'
    )
    ax.clabel(contours, fmt="%2.1f")
    pv = modelmap.plot_vector(qx, qy, qz,
                              headwidth=3, headlength=4, width=2e-3,
                              pivot='mid', minshaft=2, hstep=4, scale=2,
                              color='blue')
    filename = "_".join((f'{k}-{v}' for k, v in kwargs.items()))
    plt.savefig(f"./images/01_01_{filename}.png")
    plt.show()


plot_x_section(row=20)
