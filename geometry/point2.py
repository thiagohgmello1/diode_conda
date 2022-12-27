import numpy as np


class Point2:
    def __repr__(self):
        return f'{self.array}'


    def __sub__(self, other):
        """
        Subtract two 2D points
        :param other: point
        :return: subtracted point
        """
        sub = self.array - other.array
        return Point2(sub[0], sub[1])


    def __add__(self, other):
        """
        Add two 2D points
        :param other: point
        :return: added point
        """
        add = self.array + other.array
        return Point2(add[0], add[1])


    def __eq__(self, other):
        """
        Check if two points have the same coordinates
        :param other:
        :return: boolean
        """
        return self.x == other.x and self.y == other.y


    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y
        self.array = np.array([self.x, self.y])


if __name__ == '__main__':
    a = Point2(1, 2)
    b = Point2(1, 3)
    c = Point2(1, 2)
    print(a)
    print(a == b)
    print(a == c)
    a -= b
    print(a)
    print(a + b)
