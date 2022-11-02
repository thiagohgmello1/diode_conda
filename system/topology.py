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
        self.get_boundaries_polygons()
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


    def get_boundaries_polygons(self):
        external_boundaries = list()
        internal_boundaries = list()
        for polygon in self.topologies.polygons:
            external_boundaries.extend([polygon.outer_boundary()])
            internal_boundaries.extend(list(polygon.holes))
        self.boundaries['internal'] = internal_boundaries
        self.boundaries['external'] = external_boundaries


    def diff_polygon(self, polygon: sg.Polygon):
        polygon = self.set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.difference(pol)
        self.get_boundaries_polygons()
        self.get_segments()


    def union_polygon(self, polygon: sg.Polygon):
        polygon = self.set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.union(pol)
        self.get_boundaries_polygons()
        self.get_segments()


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


    def intersection_points(self, segment_to_compare: sg.Segment2):
        intersection_points = list()
        for segments in self.segments.values():
            for segment in segments:
                intersection_point = sg.intersection(segment, segment_to_compare)
                if intersection_point:
                    intersection_points.append([intersection_point, segment])
        return intersection_points


    @staticmethod
    def set_orientation(polygons: list[sg.Polygon]):
        for polygon in polygons:
            if polygon.orientation() == -1:
                polygon.reverse_orientation()
        return polygons


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
    test_point = sg.Point2(0, 0)
    test_segment = sg.Segment2(sg.Point2(110, 80), sg.Point2(80, 70))
    pol2 = Topology.from_file('../tests/test2.svg', 1e-9, 1e-9)
    pol = sg.Polygon([sg.Point2(*D), sg.Point2(*C), sg.Point2(*B), sg.Point2(*A)])
    print(pol2.contains(test_point))
