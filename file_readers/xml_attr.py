def moveto(points: str):
    """
    Create a list of string from moveto SVG pattern element

    :param points: moveto points
    :return: list of points as strings
    """
    points = points.split(' ')[1:-1]
    return [points.split(',') for points in points]


def lineto(next_point: str, *args):
    """
    Create a list of string from lineto SVG pattern element

    :param next_point: polygon vertex
    :param args: None
    :return: list of points as strings
    """
    next_point = next_point.split(' ')[1:-1]
    return [points.split(',') for points in next_point]


def horizontal(new_x_coord, previous_point):
    """
    Create a list of string from horizontal SVG pattern element

    :param new_x_coord: new x coordinate
    :param previous_point: previous polygon vertex
    :return: list of points as strings
    """
    previous_point[0] = new_x_coord
    return [previous_point]


def vertical(new_y_coord, previous_point):
    """
    Create a list of string from vertical SVG pattern element

    :param new_y_coord: new y coordinate
    :param previous_point: previous polygon vertex
    :return: list of points as strings
    """
    previous_point[1] = new_y_coord
    return [previous_point]


func_dict = {'M': moveto, 'L': lineto, 'H': horizontal, 'V': vertical}
