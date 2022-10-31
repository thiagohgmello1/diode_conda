import skgeom
import numpy as np

from skgeom.draw import draw
from skgeom import boolean_set
from file_readers.xml_reader import XMLReader


class Topology:
    def __init__(self, topologies):
        self.topologies: list[skgeom.Polygon] = topologies


    @classmethod
    def from_file(cls, file_name: str):
        topologies = XMLReader(file_name)
        return cls(topologies.get_geometries())


    @classmethod
    def from_points(cls, points: list[list]):
        topologies = list()
        for geometry in points:
            topologies.append(skgeom.Polygon(cls.create_points(geometry)))
        return cls(topologies)


    def get_segments(self, pos=None):
        segments = list()
        if pos:
            segments.append(self.topologies[pos].edges)
        else:
            for polygon in self.topologies:
                segments.append(polygon.edges)
        return [edge for segment in segments for edge in segment]


    def union(self):
        boolean_set.join(*self.topologies)


    @staticmethod
    def cross_point(line_0: np.array, line_1: np.array):
        den = (line_0[0, 0] - line_0[1, 0]) * (line_1[0, 1] - line_1[1, 1]) - \
              (line_0[0, 1] - line_0[1, 1]) * (line_1[0, 0] - line_1[1, 0])
        if den == 0:
            return np.array([]), False
        x_cross_num = (line_0[0, 0] * line_0[1, 1] - line_0[1, 0] * line_0[0, 1]) * (line_1[0, 0] - line_1[1, 0]) - \
                      (line_1[0, 0] * line_1[1, 1] - line_1[1, 0] * line_1[0, 1]) * (line_0[0, 0] - line_0[1, 0])
        y_cross_num = (line_0[0, 0] * line_0[1, 1] - line_0[1, 0] * line_0[0, 1]) * (line_1[0, 1] - line_1[1, 1]) - \
                      (line_1[0, 0] * line_1[1, 1] - line_1[1, 0] * line_1[0, 1]) * (line_0[0, 1] - line_0[1, 1])
        x_cross = x_cross_num / den
        y_cross = y_cross_num / den

        return np.array([x_cross, y_cross]), True


    @staticmethod
    def create_points(points_list):
        points = list()
        for point in points_list:
            point = [float(p) for p in point]
            points.append(skgeom.Point2(point[0], point[1]))
        return points


if __name__ == "__main__":
    A = [0, 0]
    B = [0, 1]
    C = [1, 0]
    D = [1, 1]
    # pol1 = Topology.from_points([[A, B, D, C]])
    # segments1 = pol1.get_segments()
    pol2 = Topology.from_file('../tests/test.svg')
    segments2 = pol2.get_segments()
    print('ei')
