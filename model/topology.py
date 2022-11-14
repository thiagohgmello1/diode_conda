import skgeom as sg

from utils.complementary_operations import equal, vec_to_point
from skgeom.draw import draw
from file_readers.xml_reader import XMLReader


class Topology:
    def __init__(self, topologies: list[sg.Polygon], scale: float):
        topologies = self.set_orientation(topologies)
        self.bbox = None
        self.topologies: sg.PolygonSet = sg.PolygonSet(topologies)
        self.boundaries = dict()
        self.get_boundaries_polygons()
        self.segments = dict()
        self.get_segments()
        self.scale = scale


    @classmethod
    def from_file(cls, file_name: str, scale: float):
        """
        Create topology from file

        :param file_name: file name
        :param scale: scale dimension
        :return: class instantiation
        """
        topologies = XMLReader(file_name)
        return cls(topologies.get_geometries(), scale)


    @classmethod
    def from_points(cls, points: list[list], scale: float):
        """
        Create topology from specified points

        :param points: geometry vertices
        :param scale: scale dimension
        :return: class instantiation
        """
        topologies = list()
        for geometry in points:
            topologies.append(sg.Polygon(cls.create_points(geometry)))
        return cls(topologies, scale)


    def contains(self, point: sg.Point2):
        """
        Check if point is inside geometry

        :param point: point to be checked
        :return: boolean indicating if point is or is not inside geometry
        """
        is_inside = self.topologies.locate(point)
        if is_inside:
            return True
        else:
            return False


    def get_boundaries_polygons(self):
        """
        Get all polygons boundaries and divide them into internal and external boundaries

        :return: None
        """
        external_boundaries = list()
        internal_boundaries = list()
        for polygon in self.topologies.polygons:
            external_boundaries.extend([polygon.outer_boundary()])
            internal_boundaries.extend(list(polygon.holes))
        self.boundaries['internal'] = internal_boundaries
        self.boundaries['external'] = external_boundaries
        self.bbox = self.boundaries['external'][0].bbox()


    def diff_polygon(self, polygon: sg.Polygon):
        """
        Subtract polygon from topology

        :param polygon: polygon to be subtracted
        :return: None
        """
        polygon = self.set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.difference(pol)
        self.get_boundaries_polygons()
        self.get_segments()


    def union_polygon(self, polygon: sg.Polygon):
        """
        Add polygon to topology

        :param polygon: polygon to be added
        :return: None
        """
        polygon = self.set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.union(pol)
        self.get_boundaries_polygons()
        self.get_segments()


    def get_segments(self):
        """
        Get all geometry segments

        :return: None
        """
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


    def intersection_points(self, segment_to_compare: sg.Segment2, actual_particle_vec: sg.Vector2) -> list:
        """
        Define all possible intersection points

        :param segment_to_compare: line segment eventually travelled by particle
        :param actual_particle_vec: actual particle position vector
        :return: all possible intersection points
        """
        actual_pos = vec_to_point(actual_particle_vec)
        intersection_points = list()
        for segments in self.segments.values():
            for segment in segments:
                intersection_point = sg.intersection(segment, segment_to_compare)
                if intersection_point and not equal(intersection_point, actual_pos):
                    intersection_points.append([intersection_point, segment])
        return intersection_points


    @staticmethod
    def set_orientation(polygons: list[sg.Polygon]) -> list:
        """
        Set polygons orientation
        :param polygons: polygons to be standardized
        :return: list of standardized polygons
        """
        for polygon in polygons:
            if polygon.orientation() == -1:
                polygon.reverse_orientation()
        return polygons


    @staticmethod
    def create_points(points_list) -> list:
        """
        Create points from list of floats
        :param points_list: list of float points
        :return: list of Point2
        """
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
    polygon_1 = Topology.from_file('../tests/test2.svg', 1e-9)
    polygon_2 = sg.Polygon([sg.Point2(*D), sg.Point2(*C), sg.Point2(*B), sg.Point2(*A)])
    polygon_2.bbox()
    print(polygon_1.contains(test_point))
