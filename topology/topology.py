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
        self.boundaries = dict()
        self.get_boundaries_points()
        self.segments = dict()
        self.get_segments()


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


    def contains(self, point: sg.Point2):
        is_inside = self.topologies.locate(point)
        if is_inside:
            return True
        else:
            return False


    def get_boundaries_points(self):
        external_boundaries = list()
        internal_boundaries = list()
        for polygon in self.topologies.polygons:
            external_boundaries.extend([polygon.outer_boundary()])
            internal_boundaries.extend(list(polygon.holes))
        self.boundaries['internal'] = internal_boundaries
        self.boundaries['external'] = external_boundaries


    def diff_polygon(self, polygon: sg.Polygon):
        self.topologies = self.topologies.difference(polygon)
        self.get_boundaries_points()


    def union_polygon(self, polygon: sg.Polygon):
        self.topologies = self.topologies.union(polygon)
        self.get_boundaries_points()


    def get_segments(self):
        def extract_segments(pols):
            segments = list()
            for pol in pols:
                segments.extend(list(pol.edges))
            return segments

        internal_segments = list()
        external_segments = list()
        for region, polygons in self.boundaries.items():
            if region == 'internal':
                internal_segments.extend(extract_segments(polygons))
            else:
                external_segments.extend(extract_segments(polygons))
        self.segments['internal'] = internal_segments
        self.segments['external'] = external_segments


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
    print(pol2.contains(test_point))
    pol2.diff_polygon(pol)
    pol2.get_boundaries_points()
    # segments2 = pol2.get_segments()
    print('ei')
    seg = sg.Segment2(sg.Point2(*D), sg.Point2(*C))
    seg.supporting_line()
