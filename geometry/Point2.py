import numpy as np


class Point2:
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y
        self.array = np.array([self.x, self.y])


    def __sub__(self, other):
        return self.array - other.array


    def __add__(self, other):
        return self.array + other.array


    def __repr__(self):
        return f'{self.array}'


if __name__ == '__main__':
    a = Point2(1, 2)
    b = Point2(1, 3)
    print(a - b)
    print(a + b)
