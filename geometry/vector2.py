import numpy as np
from geometry.point2 import Point2


class Vector2:
    def __repr__(self):
        return f'{self.array}'


    def __add__(self, other):
        x = self.x + other.x
        y = self.y + other.y
        return Vector2(x, y)


    def __sub__(self, other):
        x = self.x - other.x
        y = self.y - other.y
        return Vector2(x, y)


    def __mul__(self, other):
        x = self.x * other
        y = self.y * other
        return Vector2(x, y)


    def __eq__(self, other):
        return self.x == other.x and self.y == other.y


    def __truediv__(self, other: float):
        return Vector2(self.x / other, self.y / other)


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
        return np.cross(self.array, other.array)


    def norm2(self):
        return np.linalg.norm(self.array)


if __name__ == '__main__':
    p0 = Point2(1, 1)
    p1 = Point2(1, 2)
    vec_point1 = Vector2(p0)
    vec_point2 = Vector2(p1)
    vec_ar1 = Vector2(np.array([4, 5]))
    vec_ar2 = Vector2(np.array([4, 5]))
    vec_float = Vector2(1, 1)
    print(vec_point1 + vec_point2)
    print(vec_point1 - vec_point2)
    print(vec_point1 * 2)
    print(vec_point1 == vec_ar1)
    print(vec_ar1 == vec_ar2)
    print(vec_ar1 / 2)
    print(Vector2.dot_prod(vec_float, vec_ar1))
    print(Vector2.cross_prod(vec_float, vec_ar1))
    print(Vector2.norm2(vec_float))
