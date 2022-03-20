#!/usr/bin/env python
import flopy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp2d

from shapely import geometry


# Utils functions
def csv_2_shape(csvfile, shape=geometry.Polygon):
    df = pd.read_csv(csvfile)
    return shape([
        (row.x, row.y) for i, row in df.iterrows()])


# Simulation parameters
X0 = 0
XN = 3000
NC = 3000//30
ΔX = XN/NC

Y0 = 0
YN = 2230
NR = 2230//30
ΔY = YN/NR


Top = 3                         # m
Height = 100*0.3048
Bottom = Top-Height

recharge = 17/365/12

# geo layers
geolyr_thickness = [Height]
geolyr_subdivisions = [4]


SECOND_WELL_ON = False

xy_grid_points = np.mgrid[X0:XN:ΔX, YN:Y0:-ΔY].reshape(2, -1).T
x_grids = np.linspace(X0, XN+1, NC)


domain = geometry.box(X0, Y0, XN, YN)
lake = csv_2_shape('./data/4_lake.csv', geometry.Polygon)
lake_top = 0
lake_height = Height
lake_bottom = lake_top - lake_height
lake_conductance = 50 * ΔX * ΔY  # leakance to conductance
river_df = pd.read_csv('./data/4_river.csv')
river = geometry.LineString([
        (row.x, row.y) for i, row in river_df.iterrows()])
get_river_head = interp2d(river_df.x,
                           river_df.y,
                           river_df.h,
                           kind='linear')
river_top = 0
river_height = Height
river_bottom = river_top - river_height
river_conductance = 50 * ΔX * ΔY  # leakance to conductance

well1 = geometry.Point((2373.8920225624497, 1438.255033557047))
# 2373.8920225624497,1438.255033557047
# 1871.071716357776,1030.7494407158838

well1_top = 0
well1_bottom = Bottom
well1_rate = -200 * 5.451  # GPM → m³/day


def get_layers(top=Top, bottom=Bottom):
    all_layers = [(i, b) for i, b in enumerate(bot) if b < top]
    b = top
    for i, b in all_layers:
        if b > bottom:
            yield i, top, b
        else:
            break
        top = b
    if b <= bottom:
        yield i, top, bottom


def get_grid_points(shape, /, xy_grid_points, layers=None):
    if not layers:
        layers = [0]
    else:
        layers = list(layers)

    grid_pts = enumerate(map(geometry.Point, xy_grid_points))
    grid_boxes = enumerate(map(lambda x: geometry.box(
            x[0]-ΔX/2, x[1]-ΔY/2, x[0]+ΔX/2, x[1]+ΔY/2),
                            xy_grid_points))

    if isinstance(shape, geometry.Polygon):
        points = filter(lambda gp: shape.contains(gp[1]), grid_pts)
    elif isinstance(shape, geometry.Point):
        nearest = min(grid_pts, key=lambda gp: shape.distance(gp[1]))
        points = [nearest]
    elif isinstance(shape, geometry.LineString):
        points = filter(lambda gp: shape.intersects(gp[1]), grid_boxes)

    for i, _ in points:
        col = i // (NR)
        row = i % (NR)
        for j in layers:
            yield (j, row, col), xy_grid_points[i]


# computational layers
NLay = sum(geolyr_subdivisions)
lookup_table = np.concatenate(
    list(np.ones(s, dtype=int)*i for i, s in
         enumerate(geolyr_subdivisions)))

lyr_k_hz = [12.8]
lyr_k_vt = [12.8]


thickness = np.zeros(NLay)
k_hz = [0 for i in range(NLay)]
k_vt = [0 for i in range(NLay)]
bot = np.ones(NLay)

for lay in range(NLay):
    geo_lay = lookup_table[lay]
    thickness[lay] = geolyr_thickness[geo_lay]/geolyr_subdivisions[geo_lay]
    k_hz[lay] = lyr_k_hz[geo_lay]
    k_vt[lay] = lyr_k_vt[geo_lay]
    bot[lay] = Top-sum(thickness)


def get_riv_stress_period():
    "gives the stress_period_data on the grid_points for river grids."
    for grid_pt, pt in get_grid_points(river, xy_grid_points=xy_grid_points):
        # cellid, stage, cond, rbot, aux, boundname
        stage = get_river_head(pt[0], pt[1])[0]
        rbot = stage-1
        yield ((0, grid_pt[1], grid_pt[2]), stage, 2, rbot)


def get_chd_stress_period():
    "gives the stress_period_data on the grid_points for constant head points."
    layers_tuple = list(get_layers(top=lake_top, bottom=lake_bottom))
    for grid_pt, _ in get_grid_points(lake, xy_grid_points=xy_grid_points):
        for lay, thk, bottom in layers_tuple:
            # cellid, head
            yield ((lay, grid_pt[1], grid_pt[2]), lake_top)

    for grid_pt, pt in get_grid_points(river, xy_grid_points=xy_grid_points):
        # cellid, head
        stage = get_river_head(pt[0], pt[1])[0]
        yield ((0, grid_pt[1], grid_pt[2]), stage)


def get_well_stress_period():
    # temp fix
    well1_layers = [l[0] for l in get_layers(well1_top, well1_bottom)]
    well_pts = get_grid_points(well1, xy_grid_points=xy_grid_points,
                               layers=well1_layers)
    rate = well1_rate/len(well1_layers)
    return {0: [(wpt, rate) for wpt, _ in well_pts]}


# sp = list(get_chd_stress_period())

# x = [l[0][2] for l in sp]+[0]
# y = [l[0][1] for l in sp]+[0]


# plt.scatter(x, y)
# plt.show()


# MODELING STARTS FROM HERE:
ws = './models/4_calibration'
name = '4_calibration'

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

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=recharge)
# rivers = flopy.mf6.ModflowGwfriv(
#     gwf,
#     stress_period_data=list(get_riv_stress_period()))
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              # icelltype=[1, 0, 0, 0],
                              icelltype=1,
                              k=k_hz,
                              k33=k_vt,
                              save_specific_discharge=True)

# EXample to modify the k values after it is defined.
# k_values = npf.k.get_data()
# kv_values = npf.k33.get_data()
# for lay in layers_2nd:
#     k_values[lay] = k_2nd_layer
#     kv_values[lay] = kv_2nd_layer
# npf.k.set_data(k_values)
# npf.k33.set_data(kv_values)

chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_chd_stress_period()))
rivers = flopy.mf6.ModflowGwfriv(
    gwf,
    stress_period_data=list(get_riv_stress_period()))

wells = flopy.mf6.ModflowGwfwel(
    gwf,
    stress_period_data=get_well_stress_period())

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
    fig, ax = plt.subplots(1, 1, figsize=(9, 3), constrained_layout=True)
    ax.set_title(f'Layer-{layer}')
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    pmv.plot_array(head_arr[layer])
    pmv.plot_grid(colors='white', linewidths=0.3)
    contours = pmv.contour_array(head_arr[layer],
                                 levels=np.arange(0, 100, 1),
                                 linewidths=1.,
                                 colors='black')
    ax.clabel(contours, fmt="%.0f")
    pmv.contour_array(head_arr[layer],
                      levels=np.arange(0, 100, .2),
                      linewidths=.4,
                      colors='black')
    # flopy.plot.styles.graph_legend()
    pmv.plot_vector(qx[layer, :, :], qy[layer, :, :],
                    headwidth=3, headlength=4, width=1.4e-3, scale=20,
                    normalize=False, istep=10, jstep=10, color="white")
    shps= pmv.plot_shapefile('./data/4_river',
                             edgecolor='red',
                             linewidth=2)
    # ax.plot([x for x, y in river.coords],
    #         [y for x, y in river.coords],
    #         color='red',
    #         linewidth=2)
    plt.savefig(f"./images/03_00_plan_layer-{layer}.png")
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
    k_values = npf.k.get_data()
    pa = modelmap.plot_array(k_values, alpha=0.6)
    quadmesh = modelmap.plot_bc("CHD")
    linecollection = modelmap.plot_grid(lw=0.2, color="white")
    minor_contours = modelmap.contour_array(
        head_arr,
        levels=np.arange(0, 100, .2),
        linewidths=0.4,
        colors='black'
    )
    contours = modelmap.contour_array(
        head_arr,
        head=head_arr,
        levels=np.arange(0, 100, 1),
        linewidths=0.8,
        colors='black'
    )
    ax.clabel(contours, fmt="%.0f")
    pv = modelmap.plot_vector(qx, qy, qz,
                              headwidth=3, headlength=4, width=1e-3,
                              pivot='mid', minshaft=2, hstep=4,
                              scale=10, linewidths=0.1,
                              color='blue')
    # plt.colorbar(pa, shrink=0.5, ax=ax)
    filename = "_".join((f'{k}-{v}' for k, v in kwargs.items()))
    plt.savefig(f"./images/03_01_{filename}.png")
    plt.show()


plot_plan(layer=NLay//2)

well_gp,_ = next(get_grid_points(well1, xy_grid_points=xy_grid_points))
plot_x_section(row=well_gp[1])
