from xml.dom import minidom
import numpy as np
import skgeom


class XMLReader:
    def __init__(self, file_name: str):
        self.file_name = file_name


    def get_paths(self):
        paths = dict()
        document = minidom.parse(self.file_name)
        paths_str = [path.getAttribute('d').split(' ') for path in document.getElementsByTagName('path')]
        paths_id = [structure.getAttribute('id') for structure in document.getElementsByTagName('path')]
        for path_id, path_str in zip(paths_id, paths_str):
            paths[path_id] = np.array([path.split(',') for path in path_str[1:-1]]).astype(float)

        document.unlink()
        return paths


    def get_rectangles(self):
        paths = dict()
        document = minidom.parse(self.file_name)
        rectangles = document.getElementsByTagName('rect')
        rectangles_id = [structure.getAttribute('id') for structure in document.getElementsByTagName('path')]

