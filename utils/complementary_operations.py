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


def point_to_vec(point: Point2) -> Vector2:
    """
    Convert point to vector

    :param point: point to be converted
    :return: point as a vector
    """
    return Vector2(float(point.x()), float(point.y()))


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


def round_vec_coord(vec_to_round: Vector2, number_digits: int):
    x_rounded = round(float(vec_to_round.x()), number_digits)
    y_rounded = round(float(vec_to_round.y()), number_digits)
    return Vector2(x_rounded, y_rounded)


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
    point_1 = segment[0]
    point_2 = segment[1]
    z_numerator = (point - point_1)
    z_numerator = np.sqrt(float(z_numerator.x()) ** 2 + float(z_numerator.y()) ** 2)
    z_denominator = (point_2 - point_1)
    z_denominator = np.sqrt(float(z_denominator.x()) ** 2 + float(z_denominator.y()) ** 2)
    z = z_numerator / z_denominator
    if z <= 1:
        numerator = abs(
            float(
                (point_2.x() - point_1.x()) * (point_1.y() - point.y())
                - (point_1.x() - point.x()) * (point_2.y() - point_1.y())
            )
        )
        denominator = np.sqrt(float((point_2.x() - point_1.x())) ** 2 + float((point_2.y() - point_1.y())) ** 2)
        return numerator / denominator
    else:
        dist_1 = point - point_1
        dist_1 = np.sqrt(float(dist_1.x()) ** 2 + float(dist_1.y()) ** 2)
        dist_2 = point - point_2
        dist_2 = np.sqrt(float(dist_2.x()) ** 2 + float(dist_2.y()) ** 2)
        return min(dist_1, dist_2)
