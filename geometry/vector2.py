import numpy as np
from geometry.point2 import Point2


class Vector2:
    def __repr__(self):
        return f'Vector2({self.x}, {self.y})'


    def __add__(self, other):
        return Vector2(self.array + other.array)


    def __sub__(self, other):
        return Vector2(self.array - other.array)


    def __mul__(self, other):
        return Vector2(self.array * other)


    def __eq__(self, other):
        return all(self.array == other.array)


    def __truediv__(self, other):
        return Vector2(self.array / other)


    def __init__(self, *args):
        if len(args) == 2 and all(isinstance(arg, float) or isinstance(arg, int) for arg in args):
            self.x = args[0]
            self.y = args[1]
        elif len(args) == 1 and isinstance(args[0], Point2):
            self.x = args[0].x
            self.y = args[0].y
        elif isinstance(args[0], np.ndarray):
            self.x = args[0][0]
            self.y = args[0][1]
        else:
            raise AttributeError('Bad attribute. It must be numerical or Point2 or np.ndarray')
        self.array = np.array([self.x, self.y])


    def dot_prod(self, other):
        return np.dot(self.array, other.array)


    def cross_prod(self, other):
        return float(np.cross(self.array, other.array))


    def length(self):
        return np.linalg.norm(self.array)


    def rotate(self, ang: float):
        cos_ang = np.cos(ang)
        sin_ang = np.sin(ang)
        rotate_matrix = np.array([[cos_ang, -1 * sin_ang], [sin_ang, cos_ang]])
        return Vector2(np.dot(rotate_matrix, self.array))


    def normalize(self):
        norm = self.length()
        return Vector2(self.array / norm)


    def to_point(self):
        return Point2(self.array[0], self.array[1])


    def perpendicular(self):
        return self.rotate(np.pi / 2)


if __name__ == '__main__':
    p0 = Point2(1, 1)
    p1 = Point2(1, 2)
    vec_point1 = Vector2(p0)
    vec_point2 = Vector2(p1)
    vec_ar1 = Vector2(np.array([4, 5]))
    vec_ar2 = Vector2(np.array([4, 5]))
    vec_float = Vector2(0, 1)
    print(vec_point1 + vec_point2)
    print(vec_point1 - vec_point2)
    print(vec_point1 * 2)
    print(vec_point1 == vec_ar1)
    print(vec_ar1 == vec_ar2)
    print(vec_ar1 / 2)
    print(vec_float.dot_prod(vec_ar1))
    print(vec_float.cross_prod(vec_ar1))
    print(vec_ar1.cross_prod(vec_float))
    print(vec_float.length())
    print(vec_float.rotate(np.pi / 2))
    print(vec_ar1.normalize())
    print(vec_ar1.to_point())
