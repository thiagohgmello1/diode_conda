import numpy as np
import skgeom as sg


def dot_prod(vec_1: sg.Vector2, vec_2: sg.Vector2):
    x = vec_1.x() * vec_2.x()
    y = vec_1.y() * vec_2.y()
    return float(x + y)


def norm(vec: sg.Vector2):
    return np.sqrt(float(vec.squared_length()))


def angle_between(vec_1: sg.Vector2, vec_2: sg.Vector2):
    norm_1 = norm(vec_1)
    norm_2 = norm(vec_2)
    return np.arccos(dot_prod(vec_1, vec_2) / (norm_1 * norm_2))


def mirror(vec: sg.Vector2, normal_vec: sg.Vector2):
    return vec - 2 * (dot_prod(vec, normal_vec)) / (norm(normal_vec) ** 2) * normal_vec


def calc_normal(segment: sg.Segment2, point: sg.Point2):
    supporting_line = segment.supporting_line()
    perpendicular_line = supporting_line.perpendicular(point)
    return perpendicular_line.direction()


if __name__ == '__main__':
    vec1 = sg.Vector2(-1, -1)
    vec2 = sg.Vector2(1, 0)
    seg = sg.Segment2(sg.Point2(0, 0), sg.Point2(2, 0))
    p = sg.Point2(1, 1)
    normal = calc_normal(seg, p)
    dot_prod(vec1, vec2)
