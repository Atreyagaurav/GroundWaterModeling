from shapely import geometry
from collections import namedtuple

import pandas as pd

Rectangle = namedtuple("Rectangle", "x y w h")


def rect_2_poly(rect):
    "Convert custom Rectangle data to shapely Polygon."
    return geometry.Polygon([
        (rect.x, rect.y),
        (rect.x + rect.w, rect.y),
        (rect.x + rect.w, rect.y + rect.h),
        (rect.x, rect.y + rect.h),
    ])


def csv_2_poly(csvfile):
    df = pd.read_csv(csvfile)
    return geometry.Polygon([
        (row.x, row.y) for i, row in df.iterrows()])
