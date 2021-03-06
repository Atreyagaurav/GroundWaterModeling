import flopy
import matplotlib.pyplot as plt
import numpy as np

from shapely import geometry

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

geolyr_thickness = [50, 20, 150]
geolyr_subdivisions = [10, 4, 30]

WELL_ON = True

xy_grid_points = np.mgrid[X0:XN:ΔX, Y0:YN:ΔY].reshape(2, -1).T
x_grids = np.linspace(X0, XN, NC)

domain = geometry.box(X0, Y0, XN, YN)
river = geometry.box(300, 0, 300+850, YN)
river_top = 10
river_height = 25
river_bottom = river_top - river_height
river_conductance = 50 * ΔX * ΔY  # leakance to conductance
stream = geometry.box(6150, 0, 6150+270, YN)
stream_top = 10.5
stream_height = 3.3
stream_conductance = 1 * ΔX * ΔY  # leakance to conductance
stream_bottom = stream_top - stream_height
well = geometry.Point((5720, 2000))
well_top = -150
well_bottom = -160
well_rate = -200 * 0.1336801*60*24  # GPM → ft³/day

def get_layers(top=Top, bottom=Bottom):
    all_layers = [(i, b) for i, b in enumerate(bot) if b < top]
    for i, b in all_layers:
        if b > bottom:
            yield i, top, b
        else:
            break
        top = b
    if b <= bottom:
        yield i, top, bottom

def get_grid_points(shape, layers=None):
    if not layers:
        layers = [0]
    else:
        layers = list(layers)
    for i, gp in enumerate(xy_grid_points):
        col = i // (NR)           # might have to swap these two.
        row = i % (NR)
        pt = geometry.Point(gp[0], gp[1])
        if shape.contains(pt):
            # layer, row, col
            for j in layers:
                yield (j, row, col)

# computational layers
NLay = sum(geolyr_subdivisions)
lookup_table = np.concatenate(
    list(np.ones(s, dtype=int)*i for i, s in
         enumerate(geolyr_subdivisions)))

# hetereogeiniety in 2nd geolayer
k_2nd_layer = np.ones(shape=(NR, NC))*3.0
kv_2nd_layer = np.ones(shape=(NR, NC))*.01
for cell in get_grid_points(river):
    k_2nd_layer[cell[1], cell[2]] = 30.0
    kv_2nd_layer[cell[1], cell[2]] = 3.0

lyr_k_hz = [30.0,
            k_2nd_layer,
            150.0]
lyr_k_vt = [3.0,
            kv_2nd_layer,
            15.0]

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

    layers_tuple = list(get_layers(top=river_top, bottom=river_bottom))
    for grid_pt in get_grid_points(river):
        for lay, thk, bottom in layers_tuple:
            # cellid, stage, cond, rbot, aux, boundname
            yield ((lay, grid_pt[1], grid_pt[2]),
                   thk, river_conductance, bottom)
    layers_tuple = list(get_layers(top=stream_top, bottom=stream_bottom))
    for grid_pt in get_grid_points(stream):
        for lay, thk, bottom in layers_tuple:
            yield ((lay, grid_pt[1], grid_pt[2]),
                   thk, stream_conductance, bottom)

def get_chd_stress_period():
    "gives the stress_period_data on the grid_points for constant head points."
    # river grid points
    layers_tuple = list(get_layers(top=river_top, bottom=river_bottom))
    for grid_pt in get_grid_points(river):
        for lay, thk, bottom in layers_tuple:
            # cellid, head
            yield ((lay, grid_pt[1], grid_pt[2]), 10)

    # stream grid points
    layers_tuple = list(get_layers(top=stream_top, bottom=stream_bottom))
    for grid_pt in get_grid_points(stream):
        for lay, thk, bottom in layers_tuple:
            yield ((lay, grid_pt[1], grid_pt[2]), 10.5)

_gps = map(geometry.Point, xy_grid_points)
_well_gp = min(enumerate(_gps), key=lambda x: well.distance(x[1]))
well_row = _well_gp[0] % (NR)
well_col = _well_gp[0] // (NR)

def get_well_stress_period():
    well_layers = list(get_layers(well_top, well_bottom))
    return {0: [((i, well_row, well_col),
                 well_rate/len(well_layers)) for i, _, _ in
                well_layers]}

sp = list(get_chd_stress_period())

ipoints = np.ones((NLay, NR, NC))
for i, _ in sp:
    ipoints[i] = -1
x = [l[0][2] for l in sp]
y = [l[0][1] for l in sp]
c = [l[1] for l in sp]

plt.scatter(x, y, c=c)
plt.colorbar()
plt.show()

ws = './models/3_water_withdrawal_controversy'
name = '3_water_wd'

sim = flopy.mf6.MFSimulation(sim_name=name,
                             sim_ws=ws,
                             exe_name='modflow-mf6')

tdis = flopy.mf6.ModflowTdis(sim,
                             time_units='days')
ims = flopy.mf6.ModflowIms(sim)
gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

dis = flopy.mf6.ModflowGwfdis(gwf,
                              length_units='FEET',
                              nlay=NLay,
                              nrow=NR,
                              ncol=NC,
                              delc=ΔX,
                              delr=ΔY,
                              top=Top,
                              botm=bot)

initial_head = np.ones((NLay, NR, NC)) * Top
ic = flopy.mf6.ModflowGwfic(gwf, strt=initial_head)

recharge = flopy.mf6.ModflowGwfrcha(gwf, recharge=1/365)
rivers = flopy.mf6.ModflowGwfriv(
    gwf,
    stress_period_data=list(get_riv_stress_period()))
npf = flopy.mf6.ModflowGwfnpf(gwf,
                              icelltype=1,
                              k=k_hz,
                              k33=k_vt,
                              save_specific_discharge=True)

# EXample to modify the k values after it is defined.
# k_values = npf.k.get_data()
# kv_values = npf.k33.get_data()
# layers_2nd = [i for i, v in enumerate(lookup_table) if v == 1]
# for lay in layers_2nd:
#     k_values[lay] = k_2nd_layer
#     kv_values[lay] = kv_2nd_layer
# npf.k.set_data(k_values)
# npf.k33.set_data(kv_values)

chd = flopy.mf6.ModflowGwfchd(
    gwf,
    stress_period_data=list(get_chd_stress_period()))

if WELL_ON:
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

chd_bud = bud.get_data(text='CHD')


spdis = bud.get_data(text='DATA-SPDIS')[0]
qx, qy, qz = flopy.utils.postprocessing.get_specific_discharge(spdis, gwf)
watertable = flopy.utils.postprocessing.get_water_table(head_arr, -1e30)

plt.imshow(watertable)
plt.colorbar()
plt.show()

def plot_plan(ext='pdf', layer=0):
    fig, ax = plt.subplots(1, 1, figsize=(9, 3), constrained_layout=True)
    ax.set_title(f'Layer-{layer}')
    pmv = flopy.plot.PlotMapView(gwf, ax=ax)
    pmv.plot_array(head_arr[layer])
    pmv.plot_grid(colors='white', linewidths=0.3)
    pmv.contour_array(head_arr[layer],
                      linewidths=1.,
                      cmap='Wistia')
    # flopy.plot.styles.graph_legend()
    pmv.plot_vector(qx[layer, :, :], qy[layer, :, :],
                    normalize=False, istep=2, jstep=2, color="white")
    filename = f"./images/3_plan_layer-{layer}.{ext}"
    plt.savefig(filename)
    plt.show()
    return filename

def plot_x_section(ext='pdf', **kwargs):
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
        levels=np.arange(0, 25, .2),
        linewidths=0.4,
        colors='black'
    )
    contours = modelmap.contour_array(
        head_arr,
        head=head_arr,
        levels=np.arange(0, 25, 1),
        linewidths=0.8,
        colors='black'
    )
    ax.clabel(contours, fmt="%.0f")
    pv = modelmap.plot_vector(qx, qy, qz,
                              headwidth=3, headlength=4, width=2e-3,
                              pivot='mid', minshaft=2, hstep=4,
                              scale=3,
                              color='blue')
    # plt.colorbar(pa, shrink=0.5, ax=ax)
    filename = "_".join((f'{k}-{v}' for k, v in kwargs.items()))
    saveas = f"./images/3_section_{filename}.{ext}"
    plt.savefig(saveas)
    plt.show()
    return saveas

plot_plan(ext='png', layer=geolyr_subdivisions[0]-1)

plot_plan(ext='png', layer=sum(geolyr_subdivisions[:2])-1)

lyr_index = sum(map(lambda b: b > well_bottom, bot))
plot_plan(ext='png', layer=lyr_index)

plot_x_section(ext='png', row=20)

plot_x_section(ext='png', column=60)

zones = np.ones((NLay, NR, NC), dtype=int)
for p in get_grid_points(stream, layers=[0, 1]):
    zones[p] = 2

bm = gwf.output.zonebudget(zones)

bm.change_model_name(name)
bm.change_model_ws(ws)

bm.write_input()
bm.run_model(exe_name='modflow-zbud6')

bm.get_budget()
