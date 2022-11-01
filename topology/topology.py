import skgeom as sg
import numpy as np

from skgeom.draw import draw
from file_readers.xml_reader import XMLReader


class Topology:
    def __init__(self, topologies: list[sg.Polygon], time_scale: float, space_scale: float):
        topologies = self.set_orientation(topologies)
        self.topologies: sg.PolygonSet = sg.PolygonSet(topologies)
        self.time_scale = time_scale
        self.space_scale = space_scale


    @classmethod
    def from_file(cls, file_name: str, time_scale: float, space_scale: float):
        topologies = XMLReader(file_name)
        return cls(topologies.get_geometries(), time_scale, space_scale)


    @classmethod
    def from_points(cls, points: list[list], time_scale: float, space_scale: float):
        topologies = list()
        for geometry in points:
            topologies.append(sg.Polygon(cls.create_points(geometry)))
        return cls(topologies, time_scale, space_scale)


    def is_inside_polygon_set(self, point: sg.Point2):
        is_inside = self.topologies.locate(point)
        if is_inside:
            return True
        else:
            return False


    def diff_polygon(self, polygon: sg.Polygon):
        self.topologies = self.topologies.difference(polygon)


    def get_segments(self, pos=None):
        segments = list()
        if pos:
            segments.append(self.topologies[pos].edges)
        else:
            for polygon in self.topologies.polygons:
                segments.append(polygon.edges)
        return [edge for segment in segments for edge in segment]


    @staticmethod
    def set_orientation(polygons):
        for polygon in polygons:
            if polygon.orientation() == -1:
                polygon.reverse_orientation()
        return polygons


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
            points.append(sg.Point2(point[0], point[1]))
        return points


if __name__ == "__main__":
    A = [100, 80]
    B = [100, 100]
    C = [120, 100]
    D = [120, 80]
    # pol1 = Topology.from_points([[A, B, D, C]], 1e-9, 1e-9)
    # segments1 = pol1.get_segments()
    test_point = sg.Point2(0, 0)
    pol2 = Topology.from_file('../tests/test2.svg', 1e-9, 1e-9)
    pol = sg.Polygon([sg.Point2(*D), sg.Point2(*C), sg.Point2(*B), sg.Point2(*A)])
    print(pol2.is_inside_polygon_set(test_point))
    # segments2 = pol2.get_segments()
    print('ei')
