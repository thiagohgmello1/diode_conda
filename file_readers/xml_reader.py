from xml.dom import minidom
import skgeom as sg
import re
from file_readers.xml_attr import func_dict
from skgeom.draw import draw

SVG_D_ATTR = r'(M|L|H|V|C|S|Q|T|A|Z)'


class XMLReader:
    def __init__(self, file_name: str):
        self.file_name = file_name


    def get_general_polygons(self):
        polygons = list()
        polygons_attributes = self.get_polygons_attributes()
        for polygon in polygons_attributes:
            polygons.append(sg.Polygon(self.create_points(polygon)))
        return polygons


    def get_polygons_attributes(self):
        document = minidom.parse(self.file_name)
        paths = [
            re.split(SVG_D_ATTR, path.getAttribute('d')) for path in document.getElementsByTagName('path')
        ]

        polygons_attributes = self.create_attributes(paths)
        document.unlink()
        return polygons_attributes


    def get_rectangles(self) -> list:
        rectangles = list()
        document = minidom.parse(self.file_name)
        for rectangle in document.getElementsByTagName('rect'):
            x_0 = float(rectangle.getAttribute('x'))
            y_0 = float(rectangle.getAttribute('y'))
            x_1 = x_0 + float(rectangle.getAttribute('width'))
            y_1 = y_0 + float(rectangle.getAttribute('height'))
            p_00 = sg.Point2(x_0, y_0)
            p_01 = sg.Point2(x_0, y_1)
            p_10 = sg.Point2(x_1, y_0)
            p_11 = sg.Point2(x_1, y_1)
            rectangles.append(sg.Polygon([p_00, p_01, p_11, p_10]))
        document.unlink()
        return rectangles


    def get_geometries(self):
        polygons = self.get_general_polygons()
        rectangles = self.get_rectangles()
        polygons.extend(rectangles)
        return polygons


    @staticmethod
    def create_attributes(paths):
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


    @staticmethod
    def create_points(points_list):
        points = list()
        for point in points_list:
            point = [float(p) for p in point]
            points.append(sg.Point2(point[0], point[1]))
        return points
