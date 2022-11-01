import skgeom


def moveto(points: str):
    points = points.split(' ')[1:-1]
    return [points.split(',') for points in points]


def lineto(next_point: str, *args):
    next_point = next_point.split(' ')[1:-1]
    return [points.split(',') for points in next_point]


def horizontal(new_x_coord, previous_point):
    previous_point[0] = new_x_coord
    return [previous_point]


def vertical(new_y_coord, previous_point):
    previous_point[1] = new_y_coord
    return [previous_point]


def create_points(points_list):
    points = list()
    for point in points_list:
        point = [float(p) for p in point]
        points.append(skgeom.Point2(point[0], point[1]))
    return points


func_dict = {'M': moveto, 'L': lineto, 'H': horizontal, 'V': vertical}
