from xml.dom import minidom
import skgeom
import re
from file_readers.xml_attr import create_points
from skgeom.draw import draw


class XMLReader:
    def __init__(self, file_name: str):
        self.file_name = file_name


    def get_general_polygons(self):
        polygons = list()
        document = minidom.parse(self.file_name)
        paths_str = [path.getAttribute('d').split(' ') for path in document.getElementsByTagName('path')]
        polygons_attributes = self.get_polygons_attributes()
        create_points(polygons_attributes)
        for path_str in paths_str:
            path = [path.split(',') for path in path_str[1:-1]]
            polygons.append(skgeom.Polygon(self.create_points2(path)))
        document.unlink()
        return polygons


    def get_polygons_attributes(self):
        document = minidom.parse(self.file_name)
        paths = [
            re.split(r'(M|L|H|V|C|S|Q|T|A|Z)', path.getAttribute('d')) for path in document.getElementsByTagName('path')
        ]

        polygons_attributes = self.create_attributes(paths)
        return polygons_attributes


    @staticmethod
    def create_attributes(paths):
        attributes = [list(filter(None, path))[:-1] for path in paths]
        attributes = [dict(zip(attribute[::2], attribute[1::2])) for attribute in attributes]
        return attributes


    def get_rectangles(self) -> list:
        rectangles = list()
        document = minidom.parse(self.file_name)
        for rectangle in document.getElementsByTagName('rect'):
            x_0 = float(rectangle.getAttribute('x'))
            y_0 = float(rectangle.getAttribute('y'))
            p_0 = skgeom.Point2(x_0, y_0)
            x_1 = x_0 + float(rectangle.getAttribute('width'))
            y_1 = y_0 + float(rectangle.getAttribute('height'))
            p_1 = skgeom.Point2(x_1, y_1)
            rectangles.append(skgeom.Polygon(p_0, p_1))
        document.unlink()
        return rectangles


    def get_geometries(self):
        polygons = self.get_general_polygons()
        rectangles = self.get_rectangles()
        polygons.extend(rectangles)
        return polygons


    @staticmethod
    def create_points2(points_list):
        points = list()
        for point in points_list:
            point = [float(p) for p in point]
            points.append(skgeom.Point2(point[0], point[1]))
        return points
