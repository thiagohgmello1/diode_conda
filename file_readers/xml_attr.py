import skgeom


def create_points(attributes: list):
    polygons = list()
    for polygon in attributes:
        moveto_points = moveto(polygon.get('M'))


def moveto(points_str: str, *args):
    points_str = points_str.split(' ')[1:-1]
    return [points.split(',') for points in points_str]


def create_points2(points_list):
    points = list()
    for point in points_list:
        point = [float(p) for p in point]
        points.append(skgeom.Point2(point[0], point[1]))
    return points


func_dict = {'M': moveto}
