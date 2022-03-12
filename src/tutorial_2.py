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


Top = 15
Height = 220
Bottom = Top-Height

# geo layers
lyr_thickness = [50, 20, 150]
lyr_subdivisions = [1, 4, 10]
lyr_k_hz = [30.0, 3.0, 150.0]
lyr_k_vt = [3.0, 0.01, 15.0]

# computational layers
NLay = sum(lyr_subdivisions)
lookup_table = np.concatenate(
    list(np.ones(s, dtype=int)*i for i, s in
         enumerate(lyr_subdivisions)))

thickness = np.zeros(NLay)
k_hz = np.ones(NLay)
k_vt = np.ones(NLay)
bot = np.ones(NLay)

for lay in range(NLay):
    geo_lay = lookup_table[lay]
    thickness[lay] = lyr_thickness[geo_lay]/lyr_subdivisions[geo_lay]
    k_hz[lay] = lyr_k_hz[geo_lay]
    k_vt[lay] = lyr_k_vt[geo_lay]
    bot[lay] = Top-sum(thickness)


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
river_top = 10
river_height = 25
river_bottom = river_top - river_height
stream = rect_2_poly(Rectangle(6150, 0, 270, YN))
stream_top = 10.5
stream_height = 3.3
stream_bottom = stream_top - stream_height
well = geometry.Point((5720, 2000))


def get_layers(top=Top, bottom=Bottom):
    if bottom >= bot[0]:
        yield 0
    else:
        yield from (i for i, b in enumerate(bot) if b > bottom)


def get_grid_points(shape, layers=None):
    if not layers:
        layers = [0]
    else:
        layers = list(layers)
    for i, gp in enumerate(grid_points):
        col = i // (NR+1)           # might have to swap these two.
        row = i % (NR+1)
        pt = geometry.Point(gp[0], gp[1])
        if shape.contains(pt):
            # cellid, stage, cond, rbot, aux, boundname
            for j in layers:
                yield (j, row, col)


def get_riv_stress_period():
    "gives the stress_period_data on the grid_points for river grids."

    riv_layers = get_layers(bottom=river_bottom)
    for grid_pt in get_grid_points(river, layers=riv_layers):
        # cellid, stage, cond, rbot, aux, boundname
        yield grid_pt, 10, 50, 10-25

    stream_layers = get_layers(bottom=stream_bottom)
    for grid_pt in get_grid_points(stream, layers=stream_layers):
        yield grid_pt, 10.5, 1, 10.5-3.3


def get_chd_stress_period():
    "gives the stress_period_data on the grid_points for constant head points."
    # river grid points
    riv_layers = get_layers(bottom=river_bottom)
    for grid_pt in get_grid_points(river, layers=riv_layers):
        # cellid, head
        yield grid_pt, 10

    # stream grid points
    stream_layers = get_layers(bottom=stream_bottom)
    for grid_pt in get_grid_points(stream, layers=stream_layers):
        yield grid_pt, 10.5


# sp = list(get_riv_stress_period())

# ipoints = np.ones((NLay, NR, NC))
# for i, _ in sp:
#     ipoints[i] = -1
# x = [l[0][2] for l in sp]
# y = [l[0][1] for l in sp]


# plt.scatter(x, y)
# plt.show()

ws = './models/model-2'
name = 'model-2'

sim = flopy.mf6.MFSimulation(sim_name=name,
                             sim_ws=ws,
                             exe_name='modflow-mf6')

tdis = flopy.mf6.ModflowTdis(sim,
                             time_units='days')
ims = flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

dis = flopy.mf6.ModflowGwfdis(gwf,
                              # nogrb=True,
                              # idomain=ipoints,
                              length_units='FEET',
                              nlay=NLay,
                              nrow=NR,
                              ncol=NC,
                              delc=ΔX,
                              delr=ΔY,
                              top=Top,
                              botm=bot)

ic = flopy.mf6.ModflowGwfic(gwf)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=1/365)
rivers = flopy.mf6.ModflowGwfriv(
    gwf,
    stress_period_data=list(get_riv_stress_period()))
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              # icelltype=[1, 0, 0, 0],
                              icelltype=1,
                              k=k_hz,
                              k33=k_vt,
                              save_specific_discharge=True)
chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_chd_stress_period()))
budget_file = name + '.bud'
head_file = name + '.hds'
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord=budget_file,
                            head_filerecord=head_file,
                            saverecord=[('HEAD', 'ALL'),
                                        ('BUDGET', 'ALL')])
sim.write_simulation()
result, _ = sim.run_simulation()

if not result:
    print("Error in Simulation")
    exit(1)

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
    linecollection = modelmap.plot_grid(lw=0.5, color="white")
    contours = modelmap.contour_array(
        head_arr,
        levels=np.arange(0, 25, .1),
        lw=0.5,
        colors='black'
    )
    ax.clabel(contours, fmt="%2.1f")
    # plt.colorbar(pa, shrink=0.5, ax=ax)
    plt.show()

plot_plan(layer=1)
plot_x_section(row=0)
