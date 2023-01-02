import numpy as np
from geometry.point2 import Point2
from geometry.segment2 import Segment2


class Polygon2:
    def __repr__(self):
        return f'Polygon2({[vertex for vertex in self.vertices]})'


    def __len__(self):
        return len(self.vertices)


    def __add__(self, other):
        polygon_points = list()
        polygon_points.extend(self.intersection_points(other))
        for point in other.vertices:
            if not self.is_inside(point) and point not in polygon_points:
                polygon_points.append(point)
        for point in self.vertices:
            if not other.is_inside(point) and point not in polygon_points:
                polygon_points.append(point)
        print('ei')


    def __init__(self, points: list[Point2]):
        self.vertices = np.array(points)
        self.edges = self._edges()
        self.array = np.array([point.array for point in points])
        self._set_orientation()


    def _edges(self):
        edges = list()
        for pos in range(len(self.vertices) - 1):
            edge = Segment2(self.vertices[pos], self.vertices[pos + 1])
            edges.append(edge)
        edge = Segment2(self.vertices[-1], self.vertices[0])
        edges.append(edge)
        return np.array(edges)


    def get_point(self, pos):
        return self.vertices[pos]


    def area(self):
        cross_sum = 0
        for pos in range(self.array.shape[0] - 1):
            cross_sum += np.cross(self.array[pos], self.array[pos + 1])
        cross_sum += np.cross(self.array[-1], self.array[0])
        return abs(cross_sum) / 2


    def orientation(self):
        ori_sum = 0
        for pos in range(self.array.shape[0] - 1):
            ori_sum += (self.array[pos + 1][0] - self.array[pos][0]) * (self.array[pos + 1][1] + self.array[pos][1])
        ori_sum += (self.array[0][0] - self.array[-1][0]) * (self.array[0][1] + self.array[-1][1])
        if ori_sum > 0:
            return -1
        else:
            return 1


    def reverse_orientation(self):
        self.vertices = np.ascontiguousarray(self.vertices[::-1])
        for pos in range(len(self.edges)):
            self.edges[pos] = self.edges[pos].flip()
        self.edges = np.ascontiguousarray(self.edges[::-1])
        self.array = np.ascontiguousarray(self.array[::-1])


    def _set_orientation(self):
        if self.orientation() == -1:
            self.reverse_orientation()


    def bbox(self):
        x_min = min(self.array[:, 0])
        y_min = min(self.array[:, 1])
        x_max = max(self.array[:, 0])
        y_max = max(self.array[:, 1])
        p0 = Point2(x_min, y_min)
        p1 = Point2(x_max, y_min)
        p2 = Point2(x_max, y_max)
        p3 = Point2(x_min, y_max)
        return Polygon2([p0, p1, p2, p3])


    def is_inside(self, point: Point2):
        total_angle = 0
        for edge in self.edges:
            total_angle += edge.angle_between2(point)
        if np.isclose(abs(total_angle), 2 * np.pi):
            return True
        elif np.isclose(abs(total_angle), 0):
            for edge in self.edges:
                if edge.contains(point):
                    return True
            return False
        else:
            return None


    def _intersection_points(self, segment: Segment2):
        points = list()
        for edge in self.edges:
            intersect_point = segment.intersection(edge)
            if edge.contains(intersect_point) and segment.contains(intersect_point):
                points.append([edge, intersect_point])
        return points


    def intersection_points(self, other):
        points = list()
        for edge in self.edges:
            for other_edge in other.edges:
                intersect_point = edge.intersection(other_edge)
                if edge.contains(intersect_point) and other_edge.contains(intersect_point) or intersect_point in other_edge:
                    points.append(intersect_point)
        return points


if __name__ == '__main__':
    A = Point2(1, 3.5)
    B = Point2(3.5, 2.5)
    C = Point2(6.5, 4)
    D = Point2(2.5, 1)
    E = Point2(6, 1)
    F = Point2(1, 0)
    G = Point2(3.5, 2.5)
    H = Point2(5, 0)
    I = Point2(5, 4)
    J = Point2(5, 6)
    K = Point2(3, 6)
    p = Point2(5, 2)
    # pol1 = Polygon2([A, F, E, C, B, D])  # CCW (1)
    pol1 = Polygon2([A, D, B, C, E, F])  # CW (-1)
    pol2 = Polygon2([G, H, I, J, K])  # CCW (1)
    print(pol1.intersection_points(pol2))
    # print(pol1 + pol2)
    print(pol1.area())
    print(pol1.orientation())
    pol1.reverse_orientation()
    print(pol1.orientation())
    print(pol1.bbox())
    print(len(pol1))
    print(pol1.is_inside(p))
