import numpy as np
import point2


class Segment2:
    def __init__(self, p_0: point2, p_1: point2):
        self.init_point = p_0
        self.end_point = p_1
