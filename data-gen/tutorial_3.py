from shapely import geometry
import geopandas as gpd
import pandas as pd
from collections import namedtuple

Rectangle = namedtuple("Rectangle", "x y w h")


def rect_2_poly(rect):
    return geometry.Polygon([
        (rect.x, rect.y),
        (rect.x + rect.w, rect.y),
        (rect.x + rect.w, rect.y + rect.h),
        (rect.x, rect.y + rect.h),
    ])


shapes = dict(
    domain=rect_2_poly(Rectangle(0, 0, 7800, 4000)),
    river=rect_2_poly(Rectangle(300, 0, 850, 4000)),
    stream=rect_2_poly(Rectangle(6150, 0, 270, 4000)),
    well=geometry.Point((5720, 2000))
)


for name, geom in shapes.items():
    df = pd.DataFrame(dict(id=[1]))
    gdf = gpd.GeoDataFrame(data=df, geometry=[geom])
    gdf.to_file(f"./data/tutorial-3-{name}.shp")
