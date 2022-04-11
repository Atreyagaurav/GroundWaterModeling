#!/usr/bin/env python
import math
import flopy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
NC = 100
ΔX = XN/NC

Y0 = 0
YN = 2230
NR = 75
ΔY = YN/NR


Top = 3                         # m
Height = 100*0.3048
Bottom = Top-Height

# geo layers
geolyr_thickness = [Height]
geolyr_subdivisions = [50]


SECOND_WELL_ON = False

xy_grid_points = np.mgrid[X0:XN:ΔX, YN:Y0:-ΔY].reshape(2, -1).T
x_grids = np.linspace(X0, XN+1, NC)


domain = geometry.box(X0, Y0, XN, YN)
lake = csv_2_shape('./data/4_lake.csv', geometry.Polygon)
lake_top = 0
lake_height = Height
lake_bottom = lake_top - lake_height
lake_bed_thickness = 1                               # m
lake_conductance = 5 #* ΔX * ΔY / lake_bed_thickness  # leakance to conductance how?

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
river_width = 1
riverbed_thickness = 1

well1 = geometry.Point((2373.8920225624497, 1438.255033557047))
well1_top = 0
well1_bottom = Bottom
well1_rate = -200 * 5.451  # GPM → m³/day

well2 = geometry.Point((1871.071716357776, 1030.7494407158838))
well2_top = 0
well2_bottom = Bottom
well2_rate = -250 * 5.451  # GPM → m³/day


# CALIBRATION PARAMETERS
Kh = 4.48
riv_cond = .01
# between 10-18 inch/year
Rch = 18                        # inch/year
rech = Rch * 0.0254 / 365                  # m/day


def get_riv_conductance(intersect_length):
    "Give conductance based on river intersection on grid."
    return (riv_cond *
            intersect_length * river_width / (ΔX*ΔY)  # factor of area covered
            / riverbed_thickness)


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
        points = filter(lambda gp: shape.contains(gp[1]), grid_boxes)
        points = map(lambda gp: (gp[0], shape.intersection(gp[1]).area), points)
    elif isinstance(shape, geometry.Point):
        nearest = min(grid_pts, key=lambda gp: shape.distance(gp[1]))
        points = [(nearest[0], nearest[1].area)]
    elif isinstance(shape, geometry.LineString):
        points = filter(lambda gp: shape.intersects(gp[1]), grid_boxes)
        points = map(lambda gp: (gp[0], shape.intersection(gp[1]).length), points)

    for i, insec in points:
        col = i // (NR)
        row = i % (NR)
        for j in layers:
            yield (j, row, col), xy_grid_points[i], insec



# Calib data
calib_wells = pd.read_csv("./data/4_wells.csv")
calib_wells_grid_pts = list(calib_wells.apply(
                lambda row: next(get_grid_points(
                    geometry.Point(row.x, row.y),
                    xy_grid_points=xy_grid_points))[0], axis=1))

# computational layers
NLay = sum(geolyr_subdivisions)
lookup_table = np.concatenate(
    list(np.ones(s, dtype=int)*i for i, s in
         enumerate(geolyr_subdivisions)))

lyr_k_hz = [Kh]
lyr_k_vt = [Kh]


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
    for grid_pt, pt, length in get_grid_points(river, xy_grid_points=xy_grid_points):
        # cellid, stage, cond, rbot, aux, boundname
        stage = get_river_head(pt[0], pt[1])[0]
        rbot = stage-1
        lyrs = get_layers(stage, rbot)
        for l, t, b in lyrs:
            yield ((l, grid_pt[1], grid_pt[2]), stage,
                   get_riv_conductance(length), b)


def get_chd_stress_period():
    "gives the stress_period_data on the grid_points for constant head points."
    layers_tuple = list(get_layers(top=lake_top, bottom=lake_bottom))
    for grid_pt, _, _ in get_grid_points(lake, xy_grid_points=xy_grid_points):
        for lay, thk, bottom in layers_tuple:
            # cellid, head
            yield ((lay, grid_pt[1], grid_pt[2]), lake_top)

    for grid_pt, pt, _ in get_grid_points(river, xy_grid_points=xy_grid_points):
        # cellid, head
        stage = get_river_head(pt[0], pt[1])[0]
        rbot = stage-1
        lyrs = get_layers(stage, rbot)
        for l, t, b in lyrs:
            yield ((l, grid_pt[1], grid_pt[2]), stage)


def get_well_stress_period():
    # temp fix
    well1_layers = [l[0] for l in get_layers(well1_top, well1_bottom)]
    well1_pts = get_grid_points(well1, xy_grid_points=xy_grid_points,
                               layers=well1_layers)
    rate1 = well1_rate/len(well1_layers)
    spd = [(wpt, rate1) for wpt, _, _ in well1_pts]
    if SECOND_WELL_ON:
        well2_layers = [l[0] for l in get_layers(well2_top, well2_bottom)]
        well2_pts = get_grid_points(well2, xy_grid_points=xy_grid_points,
                                   layers=well2_layers)
        rate2 = well2_rate/len(well2_layers)
        spd += [(wpt, rate2) for wpt, _, _ in well2_pts]
    return {0: spd}



#  TO plot the heads
sp = list(get_chd_stress_period())

x = [l[0][2] for l in sp]+[0]
y = [l[0][1] for l in sp]+[0]
c = [l[1] for l in sp] + [None]

plt.scatter(x, y, c=c)
plt.xlim(left=0, right=NC)
plt.ylim(bottom=NR, top=0)
plt.colorbar()
plt.show()


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
                              length_units='METERS',
                              nlay=NLay,
                              nrow=NR,
                              ncol=NC,
                              delc=ΔX,
                              delr=ΔY,
                              top=Top,
                              botm=bot)

initial_head = np.ones((NLay, NR, NC)) * Top
for gp, head in get_chd_stress_period():
    initial_head[gp] = head

ic = flopy.mf6.ModflowGwfic(gwf, strt=initial_head)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=rech)


k_vt_new = np.ones(shape=(NLay, NR, NC))*Kh

# for gp, _ in get_chd_stress_period():
#     k_vt_new[gp] = lake_conductance
for gp, _, cond, _ in get_riv_stress_period():
    k_vt_new[gp] = cond

npf = flopy.mf6.ModflowGwfnpf(gwf,
                              icelltype=1,
                              k=k_hz,
                              k33=k_vt_new,
                              save_specific_discharge=True)

chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_chd_stress_period()))

rivers = flopy.mf6.ModflowGwfriv(
    gwf,
    stress_period_data=list(get_riv_stress_period()))

# flopy.mf6.ModflowGwfdrn(
#     gwf,
#     stress_period_data=list(get_drn_stress_period()))


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

watertable = flopy.utils.postprocessing.get_water_table(head_arr, -1e30)



spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)


model_heads = map(lambda x: watertable[(x[1], x[2])], calib_wells_grid_pts)


calib_wells.loc[:, 'model_h'] = pd.Series(model_heads)
calib_wells.loc[:, 'err'] = calib_wells.model_h - calib_wells.h
calib_wells.loc[:, 'sq_err'] = calib_wells.err * calib_wells.err
calib_wells.loc[:, 'pt_size'] = calib_wells.err.map(lambda x: abs(x))
calib_wells.loc[:, 'pt_color'] = calib_wells.err.map(lambda x: 'red' if x>0 else 'blue')

rmse = math.sqrt(calib_wells.sq_err.sum())
nse = 1 - calib_wells.sq_err.sum()/(calib_wells.h - calib_wells.h.mean()).map(lambda x: x**2).sum()

print(f'Rch={Rch} inch/year; K={Kh} m/day; RK={riv_cond} ; RMSE={rmse}; NSE={nse}')


# plots start here
gs = gs = gridspec.GridSpec(1, 5)
fig = plt.figure(constrained_layout=True)
ax1 = fig.add_subplot(gs[0, :4])


ax1.set_title('Water Table')
pmv = flopy.plot.PlotMapView(gwf, ax=ax1)
pmv.plot_array(watertable)
pmv.plot_grid(colors='white', linewidths=0.3)
# pmv.plot_bc('CHD')
contours = pmv.contour_array(watertable,
                             levels=np.arange(0, 100, 1),
                             linewidths=1.,
                             colors='black')
ax1.clabel(contours, fmt="%.0f")
pmv.contour_array(watertable,
                  levels=np.arange(0, 100, .2),
                  linewidths=.4,
                  colors='black')
shps = pmv.plot_shapefile('./data/4_river',
                         edgecolor='red',
                         linewidth=2)

ax1.scatter(calib_wells.x, calib_wells.y,
           s=calib_wells.pt_size,
           c=calib_wells.pt_color)


ax2 = fig.add_subplot(gs[0, 4])
ax2.scatter(calib_wells.h, calib_wells.model_h, c=calib_wells.pt_color)
max_h = max(calib_wells.h.max(), calib_wells.model_h.max())
plt.plot([0, max_h], [0, max_h])
fig.tight_layout()
plt.show()
