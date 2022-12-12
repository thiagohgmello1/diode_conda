import numpy as np
from skgeom import Vector2, Point2, Segment2


def dot_prod(vec_1: Vector2, vec_2: Vector2) -> float:
    """
    Calculate dot product between two vectors

    :param vec_1: vector 1
    :param vec_2: vector 2
    :return: dot product
    """
    x = vec_1.x() * vec_2.x()
    y = vec_1.y() * vec_2.y()
    return float(x + y)


def norm(vec: Vector2) -> float:
    """
    Calculate vector norm

    :param vec: vector
    :return: vector norm
    """
    return np.sqrt(float(vec.squared_length()))


def angle_between(vec_1: Vector2, vec_2: Vector2) -> float:
    """
    Calculate the angle between two vectors

    :param vec_1: vector 1
    :param vec_2: vector 2
    :return: radian angle between vectors 1 and 2
    """
    norm_1 = norm(vec_1)
    norm_2 = norm(vec_2)
    return np.arccos(dot_prod(vec_1, vec_2) / (norm_1 * norm_2))


def mirror(vec: Vector2, normal_vec: Vector2) -> Vector2:
    """
    Mirror vector according specific boundary

    :param vec: vector to be mirrored
    :param normal_vec: boundary normal vector (if exists)
    :return: mirrored vector, if applicable
    """
    if normal_vec:
        return vec - 2 * (dot_prod(vec, normal_vec)) / (norm(normal_vec) ** 2) * normal_vec
    else:
        return vec


def calc_normal(segment: Segment2, point: Point2):
    """
    Calculate segment normal vector

    :param segment: segment to define normal vector
    :param point: point to direct normal vector
    :return: normal vector (or None, if there is no segment)
    """
    if not segment:
        return None
    supporting_line = segment.supporting_line()
    perpendicular_line = supporting_line.perpendicular(point)
    return perpendicular_line.to_vector()


def vec_to_point(vector: Vector2) -> Point2:
    """
    Convert vector to point

    :param vector: vector to be converted
    :return: vector as a point
    """
    return Point2(vector.x(), vector.y())


def equal(point_1: Point2, point_2: Point2):
    """
    Check if two points are equal

    :param point_1: point 1
    :param point_2: point 2
    :return: boolean indicating if points are equal
    """
    x_coord = np.isclose(point_1.x().__float__(), point_2.x().__float__(), atol=0)
    y_coord = np.isclose(point_1.y().__float__(), point_2.y().__float__(), atol=0)
    return x_coord and y_coord


def create_segments(points_list: list[Point2]) -> list:
    """
    Create sgments from points

    :param points_list: list of points to be converted in segments
    :return: list of segments
    """
    segments_list = list()
    for pos in range(len(points_list) - 1):
        segments_list.append(Segment2(points_list[pos], points_list[pos + 1]))
    return segments_list


def calc_distance_between(segment: Segment2, point: Point2):
    line_p_1 = segment[0]
    line_p_2 = segment[1]
    numerator = np.linalg.norm(
        float(
            (line_p_2.x() - line_p_1.x()) * (line_p_1.y() - point.y())
            - (line_p_1.x() - point.x()) * (line_p_2.y() - line_p_1.y())
        )
    )
    denominator = np.sqrt(float((line_p_2.x() - line_p_1.x())) ** 2 + float((line_p_2.y() - line_p_1.y())) ** 2)
    return numerator / denominator
