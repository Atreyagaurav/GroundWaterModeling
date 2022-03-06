import flopy
import matplotlib.pyplot as plt
import numpy as np

from shapely import geometry
from collections import namedtuple

Rectangle = namedtuple("Rectangle", "x y w h")

# Simulation parameters
X0 = 0
XN = 7800
NC = 78
ΔX = 7800/NC

Y0 = 0
YN = 4000
NR = 40
ΔY = 4000/NR

NLay = 4
Top = 15
Height = 220
Bottom = Top-Height


grid_points = np.mgrid[X0:XN+1:ΔX, Y0:YN+1:ΔY].reshape(2, -1).T
x_grids = np.linspace(X0, XN+1, NC)


def rect_2_poly(rect):
    "Convert custom Rectangle data to shapely Polygon."
    return geometry.Polygon([
        (rect.x, rect.y),
        (rect.x + rect.w, rect.y),
        (rect.x + rect.w, rect.y + rect.h),
        (rect.x, rect.y + rect.h),
    ])


domain = rect_2_poly(Rectangle(X0, Y0, XN, YN))
river = rect_2_poly(Rectangle(300, 0, 850, YN))
stream = rect_2_poly(Rectangle(6150, 0, 270, YN))
well = geometry.Point((5720, 2000))


def get_stress_period():
    "gives the stress_period_data on the grid_points"
    for i, gp in enumerate(grid_points):
        col = i // (NR+1)           # might have to swap these two.
        row = i % (NR+1)
        pt = geometry.Point(gp[0], gp[1])
        if river.contains(pt):
            yield (0, row, col), 10
        elif stream.contains(pt):
            yield (0, row, col), 10.5


sp = list(get_stress_period())

ipoints = np.ones((NLay, NR, NC))
for i, _ in sp:
    ipoints[i] = -1
x = [l[0][2] for l in sp]
y = [l[0][1] for l in sp]


plt.scatter(x, y)
plt.show()

ws = './models/model-2'
name = 'model-2'

sim = flopy.mf6.MFSimulation(sim_name=name,
                             sim_ws=ws,
                             exe_name='modflow-mf6')

tdis = flopy.mf6.ModflowTdis(sim)
ims = flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

bot = Top - np.linspace(Height / NLay, Height, NLay)
dis = flopy.mf6.ModflowGwfdis(gwf,
                              # nogrb=True,
                              # idomain=ipoints,
                              nlay=NLay,
                              nrow=NR,
                              ncol=NC,
                              delc=ΔX,
                              delr=ΔY,
                              top=Top,
                              botm=bot)

ic = flopy.mf6.ModflowGwfic(gwf)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=2)
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              # icelltype=[1, 0, 0, 0],
                              k=[30, 3, 150, 150],
                              k33=[3, .01, 15, 15],
                              save_specific_discharge=True)
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_stress_period())
)
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


def plot_plan(layer=0):
    pmv = flopy.plot.PlotMapView(gwf)
    pmv.plot_array(head_arr[layer])
    pmv.plot_grid(colors='white', linewidths=0.3)
    pmv.contour_array(head_arr[layer],
                      linewidths=1.,
                      c_label=True,
                      cmap='Wistia')
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
    pa = modelmap.plot_array(head_arr, vmin=Bottom, vmax=Top)
    quadmesh = modelmap.plot_bc("CHD")
    linecollection = modelmap.plot_grid(lw=0.5, color="blue")
    contours = modelmap.contour_array(
        head_arr,
        levels=[100*i for i in range(100)]
    )
    ax.clabel(contours, fmt="%2.1f")
    plt.colorbar(pa, shrink=0.5, ax=ax)
    plt.show()


plot_plan(layer=3)
plot_x_section(row=20)
