import flopy

import numpy as np
import matplotlib.pyplot as plt

import flopy
import matplotlib.pyplot as plt

ws = './models/model-2'
name = 'model-2'

sim = flopy.mf6.MFSimulation(sim_name=name, sim_ws=ws, exe_name='modflow-mf6')

tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

dis = flopy.mf6.ModflowGwfdis(gwf,
                              nrow=60,
                              ncol=40,
                              delc=50,
                              delr=50,
                              top=3832,
                              botm=np.linspace(3832, 3600, 10))

ic = flopy.mf6.ModflowGwfic(gwf)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=0.0055)
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              k=4.0,
                              save_specific_discharge=True)
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=(
        [((0, i, 0), 3820) for i in range(60)] +
        [((0, i, 39), 3824) for i in range(60)]))
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
plt.savefig(f"{ws}/plot.png")

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
    pa = modelmap.plot_array(head_arr, vmin=3832, vmax=3600)
    quadmesh = modelmap.plot_bc("CHD")
    linecollection = modelmap.plot_grid(lw=0.5, color="blue")
    contours = modelmap.contour_array(
        head_arr,
        levels=np.arange(3815, 3830, .2),
    )
    ax.clabel(contours, fmt="%2.1f")
    plt.colorbar(pa, shrink=0.5, ax=ax)
    plt.show()


plot_x_section(row=20)
