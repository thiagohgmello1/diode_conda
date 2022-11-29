import re

from xml.dom import minidom
from skgeom.draw import draw
from skgeom import Polygon, Point2
from file_readers.xml_attr import func_dict

SVG_D_ATTR = r'(M|L|H|V|C|S|Q|T|A|Z)'


class XMLReader:
    def __init__(self, file_name: str, scale: float):
        self.file_name = file_name
        self.min_pos = list()
        self.scale = scale


    def get_general_polygons(self) -> list[Polygon]:
        """
        Create general polygons

        :return: created polygons
        """
        polygons = list()
        polygons_attributes = self.get_polygons_attributes()
        self.get_geometry_limits(polygons_attributes)
        for polygon in polygons_attributes:
            polygons.append(Polygon(self.create_points(polygon)))
        return polygons


    def get_polygons_attributes(self) -> list:
        """
        Get parameters from SVG file

        :return: polygons attributes as a list of strings. Return contains parameters as SVG PATH pattern
        """
        document = minidom.parse(self.file_name)
        paths = [
            re.split(SVG_D_ATTR, path.getAttribute('d')) for path in document.getElementsByTagName('path')
        ]
        rectangles = self.get_rectangles()
        if len(rectangles):
            paths.append(rectangles)
        polygons_attributes = self.create_attributes(paths)
        document.unlink()
        return polygons_attributes


    def get_rectangles(self) -> list:
        """
        Get rectangles from SVG file (if exist)

        :return: list of rectangles parameters
        """
        rectangles = list()
        document = minidom.parse(self.file_name)
        for rectangle in document.getElementsByTagName('rect'):
            x_0 = float(rectangle.getAttribute('x'))
            y_0 = float(rectangle.getAttribute('y'))
            x_1 = x_0 + float(rectangle.getAttribute('width'))
            y_1 = y_0 + float(rectangle.getAttribute('height'))
            value = ['', 'M', f' {x_0},{y_0} ', 'V', f'{y_1}', 'H', f'{x_1}', 'V', f'{y_0}', 'Z', '']
            rectangles.extend(value)

        document.unlink()
        return rectangles


    def get_geometries(self) -> list[Polygon]:
        """
        Get all geometry polygons

        :return: list of topology polygons
        """
        polygons = self.get_general_polygons()
        return polygons


    def get_geometry_limits(self, list_of_polygons):
        """
        Get geometry limits to set all elements in (0,0) origin

        :param list_of_polygons: list of polygons to extract vertex positions
        :return: None
        """
        def get_min(elements, pos):
            return min(
                [min(polygon, key=lambda x: float(x[pos]))[pos] for polygon in elements], key=lambda x: float(x)
            )
        self.min_pos.append(float(get_min(list_of_polygons, 0)))
        self.min_pos.append(float(get_min(list_of_polygons, 1)))


    def create_points(self, points_list) -> list[Point2]:
        """
        Create scaled points to polygons
        :param points_list: list of extracted points
        :return: list of created points
        """
        points = list()
        for point in points_list:
            point = [(float(point[pos]) - self.min_pos[pos]) * self.scale for pos in range(2)]
            points.append(Point2(point[0], point[1]))
        return points


    @staticmethod
    def create_attributes(paths) -> list:
        """
        Create attributes from list of SVG patterns

        :param paths: SVG patterns
        :return: list of extracted attributes
        """
        polygons = list()
        attributes = [list(filter(None, path))[:-1] for path in paths]
        attributes = [tuple(zip(attribute[::2], attribute[1::2])) for attribute in attributes]
        for polygon in attributes:
            params = list()
            for attribute, value in polygon:
                if attribute == 'M':
                    params.extend(func_dict[attribute](value))
                else:
                    params.extend(func_dict[attribute](value, params[-1].copy()))
            polygons.append(params)
        return polygons
