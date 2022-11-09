import numpy as np
from skgeom import Vector2, Point2, Segment2


def dot_prod(vec_1: Vector2, vec_2: Vector2):
    x = vec_1.x() * vec_2.x()
    y = vec_1.y() * vec_2.y()
    return float(x + y)


def norm(vec: Vector2):
    return np.sqrt(float(vec.squared_length()))


def angle_between(vec_1: Vector2, vec_2: Vector2):
    norm_1 = norm(vec_1)
    norm_2 = norm(vec_2)
    return np.arccos(dot_prod(vec_1, vec_2) / (norm_1 * norm_2))


def mirror(vec: Vector2, normal_vec: Vector2):
    if normal_vec:
        return vec - 2 * (dot_prod(vec, normal_vec)) / (norm(normal_vec) ** 2) * normal_vec
    else:
        return vec


def calc_normal(segment: Segment2, point: Point2):
    if not segment:
        return None
    supporting_line = segment.supporting_line()
    perpendicular_line = supporting_line.perpendicular(point)
    return perpendicular_line.to_vector()


def vec_to_point(vector: Vector2):
    return Point2(vector.x(), vector.y())


def equal(point_1: Point2, point_2: Point2):
    x_coord = np.isclose(point_1.x().__float__(), point_2.x().__float__())
    y_coord = np.isclose(point_1.y().__float__(), point_2.y().__float__())
    return x_coord and y_coord


def create_segments(points_list: list[Point2]):
    segments_list = list()
    for pos in range(len(points_list) - 1):
        segments_list.append(Segment2(points_list[pos], points_list[pos + 1]))
    return segments_list


if __name__ == '__main__':
    vec1 = Vector2(-1, -1)
    vec2 = Vector2(1, 0)
    seg = Segment2(Point2(0, 0), Point2(2, 0))
    p = Point2(1, 1)
    normal = calc_normal(seg, p)
    dot_prod(vec1, vec2)
