from scipy.integrate import quad, dblquad
from scipy.constants import c, elementary_charge
from numpy import cos, sin, pi, tan, arctan, deg2rad, sqrt, sign, linspace


class Analytical:
    def __init__(self, a, h, alpha, v0, voltage, mobility, m_max=1000):
        self.a = a
        self.h = h
        self.alpha = alpha
        self.d = self.a + self.h / tan(self.alpha)
        self.v0 = v0
        self.voltage = voltage
        self.mobility = mobility
        self.m_max = m_max
        self.x_am = list()
        self.y_am = list()
        self.calc_m_max()


    def calc_m_max(self):
        m = 2
        x_list = list()
        y_list = list()
        x_list.append(None)
        y_list.append(None)
        x_am1 = -self.h
        y_am1 = self.d
        x_list.append(x_am1)
        y_list.append(y_am1)

        x_am = x_am1 - 2 * self.d * sin((m - 1) * (2 * self.alpha - pi))
        y_am = y_am1 + 2 * self.d * cos((m - 1) * (2 * self.alpha - pi))
        x_list.append(x_am)
        y_list.append(y_am)

        while (x_am1 < 0 and x_am < 0) and m < self.m_max:
            x_am1 = x_am
            y_am1 = y_am
            m += 1
            x_am = x_am1 - 2 * self.d * sin((m - 1) * (2 * self.alpha - pi))
            y_am = y_am1 + 2 * self.d * cos((m - 1) * (2 * self.alpha - pi))
            x_list.append(x_am)
            y_list.append(y_am)

        self.m_max = m
        self.x_am = x_list
        self.y_am = y_list


    def y_min(self, m, x_am):
        if m < self.m_max or (m == self.m_max and self.m_max % 2 == 0):
            return -self.d
        else:
            # return self.d + x_am / abs(sin(2 * m * self.alpha))
            return -self.d


    def y_max(self, m, x_am):
        if m < self.m_max or (m == self.m_max and self.m_max % 2 == 1):
            return self.d
        else:
            # return -self.d - x_am / abs(sin(2 * m * self.alpha))
            return self.d


    def theta_min(self, m, x_am, y_am):
        if m == 0:
            return lambda y: arctan((-y - self.a) / self.h)
        else:
            return lambda y: (-1) ** m * (m * (pi - 2 * self.alpha) +
                                          arctan((-abs(sin(2 * m * self.alpha)) / tan(2 * m * self.alpha) *
                                                  (self.d + (-1) ** m * y) + y_am + (-1) ** m * self.a) /
                                                 (x_am + abs(sin(2 * m * self.alpha)) * (self.d + (-1) ** m * y))
                                                 )
                                          )


    def theta_max(self, m, x_am, y_am):
        if m == 0:
            return lambda y: arctan((-y + self.a) / self.h)
        else:
            return lambda y: (-1) ** m * (m * (pi - 2 * self.alpha) +
                                          arctan((-abs(sin(2 * m * self.alpha)) / tan(2 * m * self.alpha) *
                                                  (self.d + (-1) ** m * y) + y_am - (-1) ** m * self.a) /
                                                 (x_am + abs(sin(2 * m * self.alpha)) * (self.d + (-1) ** m * y))
                                                 )
                                          )


    def calc_f1(self):
        f1_list = list()
        for m in range(self.m_max + 1):
            x_am = self.x_am[m]
            y_am = self.y_am[m]
            y_min = self.y_min(m, x_am)
            y_max = self.y_max(m, x_am)
            theta_min = self.theta_min(m, x_am, y_am)
            theta_max = self.theta_max(m, x_am, y_am)
            f1 = dblquad(lambda theta, y: 1, y_min, y_max, theta_min, theta_max)[0]
            f1_list.append(f1)
        f1 = sum(f1_list) / (pi * self.d)
        print(f1)


    def calc_f1_test(self, num_elem=11):
        f1_list = list()
        for m in range(self.m_max + 1):
            x_am = self.x_am[m]
            y_am = self.y_am[m]
            y_min = self.y_min(m, x_am)
            y_max = self.y_max(m, x_am)
            y_values = linspace(y_min, y_max, num_elem)
            delta_y = (y_max - y_min) / num_elem
            f1_inter = list()
            theta_min = self.theta_min(m, x_am, y_am)
            theta_max = self.theta_max(m, x_am, y_am)
            for y in y_values:
                f1_inter.append((theta_max(y) - theta_min(y)) * delta_y)
            # f1_list_test = quad(lambda y: theta_max(y) - theta_min(y), y_min, y_max)[0]
            f1_list.append(sum(f1_inter))
        f1 = sum(f1_list) / (2 * pi * self.d)
        print(f1)


if __name__ == '__main__':
    analytical = Analytical(a=1, h=1, alpha=deg2rad(40), v0=1, voltage=0, mobility=1, m_max=1000)
    analytical.calc_f1()
    # analytical.calc_f1_test()
