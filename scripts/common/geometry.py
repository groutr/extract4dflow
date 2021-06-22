from shapely.geometry import MultiPoint, Polygon
from scipy.spatial import cKDTree


def clip_point_to_roi(vertices, points):
    """Clip points to a polygon region of interest

    Arguments:
      vertices: Ordered array of vertices that describe polygon
      points: an Nx2 array of points

    Returns:
      (list): Boolean mask of size N. True if point lies inside polygon
        False otherwise.
    """
    rv = []
    polygon = Polygon(vertices)
    for pt in MultiPoint(points):
        rv.append(polygon.contains(pt))
    return rv

def kd_nearest_neighbor(vertices, points):
    tree = cKDTree(vertices)
    return tree.query(points)