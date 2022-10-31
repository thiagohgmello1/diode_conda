from file_readers.xml_reader import XMLReader
from skgeom.draw import draw

svg_file = XMLReader('tests/test.svg')

rectangles = svg_file.get_rectangles()
polygons = svg_file.get_general_polygons()
print('ei')
