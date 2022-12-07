import numpy as np
import skgeom as sg
import matplotlib.pyplot as plt

from skgeom.draw import draw
from matplotlib.widgets import Button
from file_readers.xml_reader import XMLReader
from matplotlib.backend_bases import MouseButton
from utils.complementary_operations import equal, vec_to_point, calc_distance_between


class Topology:
    def __init__(self, topologies: list[sg.Polygon], scale: float):
        """
        Create simulated topology

        :param topologies: list of polygon that form desired topology to be simulated
        :param scale: length scale in meters (ex.: 10e-9 for nanometer)
        """
        topologies = self._set_orientation(topologies)
        self.bbox = None
        self.topologies: sg.PolygonSet = sg.PolygonSet(topologies)
        self.boundaries = dict()
        self._get_boundaries_polygons()
        self.area = self.calc_topology_area()
        self.segments = dict()
        self._get_segments()
        self.scale = scale
        self.current_computing_elements = {'direct': list(), 'reverse': list()}
        self._get_current_computing_elements()
        print(self.current_computing_elements)

    @classmethod
    def from_file(cls, file_name: str, scale: float):
        """
        Create topology from file

        :param file_name: file name
        :param scale: scale dimension (ex.: 1e-6; 1e-9)
        :return: class instantiation
        """
        topologies = XMLReader(file_name, scale)
        return cls(topologies.get_geometries(), scale)

    @classmethod
    def from_points(cls, points: list[list], scale: float):
        """
        Create topology from specified points

        :param points: geometry vertices
        :param scale: scale dimension (ex.: 1e-6; 1e-9)
        :return: class instantiation
        """
        topologies = list()
        for geometry in points:
            topologies.append(sg.Polygon(cls._create_points(geometry, scale)))
        return cls(topologies, scale)

    def calc_topology_area(self):
        area = {'external': 0, 'internal': 0}
        for boundary_type, polygons in self.boundaries.items():
            for polygon in polygons:
                area[boundary_type] += polygon.area()
        return float(area['external'] - area['internal'])

    def contains(self, point: sg.Point2) -> bool:
        """
        Check if point is inside geometry

        :param point: point to be checked
        :return: boolean indicating if point is or is not inside geometry
        """
        return self.topologies.locate(point)

    def _get_boundaries_polygons(self):
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
        polygon = self._set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.difference(pol)
        self._get_boundaries_polygons()
        self._get_segments()

    def union_polygon(self, polygon: sg.Polygon):
        """
        Add polygon to topology

        :param polygon: polygon to be added
        :return: None
        """
        polygon = self._set_orientation([polygon])
        for pol in polygon:
            self.topologies = self.topologies.union(pol)
        self._get_boundaries_polygons()
        self._get_segments()

    def _get_segments(self):
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

    def intersection_points(self, traveled_path: sg.Segment2) -> list:
        """
        Define all possible intersection points

        :param traveled_path: line segment eventually travelled by particle
        :return: all possible intersection points
        """
        actual_pos = vec_to_point(traveled_path[0])
        intersection_points = list()
        for segments in self.segments.values():
            for segment in segments:
                intersection_point = sg.intersection(segment, traveled_path)
                if intersection_point and not equal(intersection_point, actual_pos, self.scale):
                    intersection_points.append([intersection_point, segment])
        return intersection_points

    def get_closer_segment(self, selected_point: sg.Point2):
        distance = np.inf
        selected_segment = None
        for segments in self.segments.values():
            for segment in segments:
                dist = calc_distance_between(segment, selected_point)
                if dist < distance:
                    distance = dist
                    selected_segment = segment
        return selected_segment

    def _get_current_computing_elements(self):
        fig = plt.figure(num='Current elements choice')
        draw(self.topologies)
        ax_select_segments = fig.add_axes(rect=[0.77, 0.9, 0.1, 0.05])
        select_segments = Button(ax_select_segments, 'Current')
        select_segments.on_clicked(self._on_click_segments)
        plt.show()

    def _on_click_segments(self, event):
        if event.button is MouseButton.LEFT:
            self.binding_id = plt.connect('button_press_event', self._on_click_event)

    def _on_click_event(self, event):
        if event.button is MouseButton.LEFT and event.xdata and event.ydata:
            point = sg.Point2(event.xdata, event.ydata)
            segment = self.get_closer_segment(point)
            if segment not in self.current_computing_elements['direct']:
                self.current_computing_elements['direct'].append(segment)
        elif event.button is MouseButton.RIGHT and event.xdata and event.ydata:
            point = sg.Point2(event.xdata, event.ydata)
            segment = self.get_closer_segment(point)
            if segment not in self.current_computing_elements['reverse']:
                self.current_computing_elements['reverse'].append(segment)

    @staticmethod
    def _set_orientation(polygons: list[sg.Polygon]) -> list:
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
    def _create_points(points_list, scale) -> list:
        """
        Create points from list of floats

        :param points_list: list of float points
        :return: list of Point2
        """
        points = list()
        for point in points_list:
            point = [float(p) * scale for p in point]
            points.append(sg.Point2(point[0], point[1]))
        return points
