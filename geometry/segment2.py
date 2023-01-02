import numpy as np
from geometry.point2 import Point2
from geometry.vector2 import Vector2


class Segment2:
    def __repr__(self):
        return f'Segment2({self.p0}, {self.p1})'


    def __str__(self):
        return f'Segment2({self.p0}, {self.p1})'


    def __ne__(self, other):
        return not np.array_equal(self.array, other.array)


    def __eq__(self, other):
        return np.array_equal(self.array, other.array)


    def __iter__(self):
        return self.points


    def __truediv__(self, other):
        return self.array / other.array


    def __copy__(self):
        return Segment2(self.p0.__copy__(), self.p1.__copy__())


    def __init__(self, p_0: Point2, p_1: Point2):
        self.p0 = p_0
        self.p1 = p_1
        self.points = [self.p0, self.p1]
        self.array = np.array([p_0.array, p_1.array])


    def to_vector(self):
        return Vector2(self.p1 - self.p0)


    def length(self):
        return (self.p1 - self.p0).length()


    def perpendicular_vec(self, guide_vec: Vector2):
        p0 = Vector2(self.p0)
        seg_vec = self.to_vector().normalize()
        proj_vec = p0 + seg_vec * ((guide_vec - p0).dot_prod(seg_vec))
        return (guide_vec - proj_vec).normalize()


    def is_collinear(self, point: Point2):
        ab = point - self.p0
        bc = self.p1 - point
        cross_prod = float(np.cross(ab.array, bc.array))
        if np.isclose(cross_prod, 0):
            return True
        else:
            return False


    def angle_between(self, point: Point2):
        a = self.p0 - point
        b = self.p1 - point
        dot_prod = np.dot(a.array, b.array)
        return np.arccos(dot_prod / (a.length() * b.length()))


    def angle_between2(self, point: Point2):
        a = self.p0 - point
        b = self.p1 - point
        ang = np.arctan2(a.array[1], a.array[0]) - np.arctan2(b.array[1], b.array[0])
        while ang > np.pi:
            ang -= 2 * np.pi
        while ang < -1 * np.pi:
            ang += 2 * np.pi
        return ang


    def intersection(self, other):
        x0 = self.p0.x
        y0 = self.p0.y
        x1 = self.p1.x
        y1 = self.p1.y

        x2 = other.p0.x
        y2 = other.p0.y
        x3 = other.p1.x
        y3 = other.p1.y

        num_x = (x0 * y1 - y0 * x1) * (x2 - x3) - (x0 - x1) * (x2 * y3 - y2 * x3)
        den_x = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3)
        num_y = (x0 * y1 - y0 * x1) * (y2 - y3) - (y0 - y1) * (x2 * y3 - y2 * x3)
        den_y = (x0 - x1) * (y2 - y3) - (y0 - y1) * (x2 - x3)

        if den_x == 0 or den_y == 0:
            return Point2(np.inf, np.inf)

        point_x = num_x / den_x
        point_y = num_y / den_y

        return Point2(point_x, point_y)


    def flip(self):
        return Segment2(self.p1, self.p0)


    def contains(self, point: Point2):
        if np.inf in point.array:
            return False
        collinear = self.is_collinear(point)
        dist1 = (point - self.p0).length()
        dist2 = (point - self.p1).length()
        if collinear and all([dist1, dist2] < self.length()):
            return True
        else:
            return False


    def distance(self, point: Point2):
        x0 = point.x
        y0 = point.y
        x1 = self.p0.x
        y1 = self.p0.y
        x2 = self.p1.x
        y2 = self.p1.y
        num = (x2 - x1) * (y1 - y0) - (x1 - x0) * (y2 - y1)
        den = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        return abs(num) / den


if __name__ == '__main__':
    point_0 = Point2(0, 0)
    point_1 = Point2(4, 4)
    point_2 = Point2(1, 4)
    vec = Vector2(4, 2)
    seg1 = Segment2(point_0, point_1)
    print(seg1)
    print(seg1.to_vector())
    print(seg1.length())
    print(seg1.perpendicular_vec(vec))
    print(seg1.is_collinear(point_2))
    print(seg1.angle_between(point_2) * 180 / np.pi)
    print(seg1.angle_between2(point_2) * 180 / np.pi)
