from scipy.integrate import quad, dblquad
from scipy.constants import c, elementary_charge
from numpy import cos, sin, pi, tan, arctan, deg2rad, sqrt


class Analytical:
    def __init__(
            self,
            fermi_vel: float,
            neck_height: float,
            diode_length: float,
            angle: float,
            voltage: float,
            carrier_concentration: float,
            mobility: float,
            max_reflections: int = 1000
    ):
        self.angle = angle
        self.fermi_vel = fermi_vel
        self.reflection_counter = 0
        self.half_neck = neck_height / 2
        self.length = diode_length
        self.e_field = voltage / self.length
        self.max_reflections = max_reflections
        self.drift_vel = mobility * self.e_field
        self.vel_ratio = self.drift_vel / self.fermi_vel
        self.carrier_concentration = carrier_concentration
        self.half_shoulder = self.half_neck + self.length / tan(self.angle)


    def vt_theta_y(self):
        return lambda theta, y: self.fermi_vel * sqrt(1 + self.vel_ratio ** 2 + 2 * self.vel_ratio * cos(theta))


    def sin_theta(self):
        return lambda theta: sin(theta) / sqrt(1 + self.vel_ratio ** 2 + 2 * self.vel_ratio * cos(theta))


    def mean_vel(self):
        f = lambda theta: self.fermi_vel * sqrt(1 + self.vel_ratio ** 2 + 2 * self.vel_ratio * cos(theta)) * cos(theta)
        return 1 / (2 * pi) * quad(f, -pi / 2, pi / 2)[0]


    def calc_next_coordinates(self, x_previous, y_previous):
        if self.reflection_counter <= 1:
            return -self.length, self.half_shoulder
        else:
            x_next = x_previous - 2 * self.half_shoulder * sin((self.reflection_counter - 1) * (2 * self.angle - pi))
            y_next = y_previous + 2 * self.half_shoulder * cos((self.reflection_counter - 1) * (2 * self.angle - pi))
            return x_next, y_next


    def y_min(self, x_m, max_reflection):
        if max_reflection and (self.reflection_counter % 2 == 1):
            return self.half_shoulder + x_m / abs(sin(2 * self.reflection_counter * self.angle))
        else:
            return -1 * self.half_shoulder


    def y_max(self, x_m, max_reflection):
        if max_reflection and (self.reflection_counter % 2 == 0):
            return -1 * self.half_shoulder - x_m / abs(sin(2 * self.reflection_counter * self.angle))
        else:
            return self.half_shoulder


    @staticmethod
    def func_1(y, reflection_counter, angle, half_neck, half_shoulder, x_m, y_m):
        atan_num = (-1 * abs(sin(2 * reflection_counter * angle)) * (half_shoulder + ((-1) ** reflection_counter) * y) /
                    tan(2 * reflection_counter * angle) + y_m + ((-1) ** reflection_counter) * half_neck)
        atan_den = (x_m + abs(sin(2 * reflection_counter * angle)) * (half_shoulder + ((-1) ** reflection_counter) * y))
        if atan_den == 0:
            return ((-1) ** reflection_counter) * (reflection_counter * (pi - 2 * angle) + pi / 2)
        else:
            return ((-1) ** reflection_counter) * (reflection_counter * (pi - 2 * angle) + arctan(atan_num / atan_den))

    def theta_min(self, x_m, y_m):
        if self.reflection_counter == 0:
            return lambda y: arctan((-y - self.half_neck) / self.length)
        else:
            return lambda y: ((-1) ** self.reflection_counter) * (
                    self.reflection_counter * (pi - 2 * self.angle)
                    + arctan((-1 * abs(sin(2 * self.reflection_counter * self.angle)) *
                              (self.half_shoulder + ((-1) ** self.reflection_counter) * y) /
                              tan(2 * self.reflection_counter * self.angle) + y_m + ((-1) ** self.reflection_counter) *
                              self.half_neck) / (x_m + abs(sin(2 * self.reflection_counter * self.angle)) *
                                                 (self.half_shoulder + ((-1) ** self.reflection_counter) * y))
                             )
            )


    def theta_max(self, x_m, y_m):
        if self.reflection_counter == 0:
            return lambda y: arctan((self.half_neck - y) / self.length)
        else:
            return lambda y: ((-1) ** self.reflection_counter) * (
                    self.reflection_counter * (pi - 2 * self.angle)
                    + arctan((-1 * abs(sin(2 * self.reflection_counter * self.angle)) *
                              (self.half_shoulder + ((-1) ** self.reflection_counter) * y) /
                              tan(2 * self.reflection_counter * self.angle) + y_m - ((-1) ** self.reflection_counter) *
                              self.half_neck) / (x_m + abs(sin(2 * self.reflection_counter * self.angle)) *
                                                 (self.half_shoulder + ((-1) ** self.reflection_counter) * y))
                             )
            )


    def m_probability(self, theta_min, theta_max, y_min, y_max):
        num = 2 * dblquad(self.vt_theta_y(), y_min, y_max, theta_min, theta_max)[0]
        den = dblquad(self.vt_theta_y(), -self.half_shoulder, self.half_shoulder, -pi / 2, pi / 2)[0]
        if num / den < 0:
            print('Menor que zero')
        return num / den


    def m_probability_zb(self, theta_min, theta_max, y_min, y_max):
        num = 1 / (pi * self.half_shoulder) * dblquad(lambda theta, y: 1, y_min, y_max, theta_min, theta_max)[0]
        if num:
            print('Menor que zero')
        return num


    def calc_direct_current(self, normalize=False):
        total_direct_current = 2 * elementary_charge * self.carrier_concentration * self.fermi_vel * self.half_shoulder / pi
        x_m1, y_m1 = self.calc_next_coordinates(None, None)
        max_reflection = False
        stop_condition = False
        transport_probability = 0

        while not stop_condition:
            stop_condition = max_reflection
            x_m = x_m1
            y_m = y_m1
            y_min = self.y_min(x_m, max_reflection)
            y_max = self.y_max(x_m, max_reflection)
            t_min = self.theta_min(x_m, y_m)
            t_max = self.theta_max(x_m, y_m)
            # transport_probability += self.m_probability(t_min, t_max, y_min, y_max)
            transport_probability += self.m_probability_zb(t_min, t_max, y_min, y_max)
            self.reflection_counter += 1
            x_m1, y_m1 = self.calc_next_coordinates(x_m, y_m)
            if self.reflection_counter > self.max_reflections:
                break
            if x_m < 0 and x_m1 > 0:
                max_reflection = True

        self.reflection_counter = 0
        print(f'Transport probability: {transport_probability}')
        direct_current = total_direct_current * transport_probability
        normalize_factor = (2 * elementary_charge * self.carrier_concentration * self.fermi_vel * self.half_neck / pi)
        if normalize:
            return direct_current / normalize_factor
        else:
            return direct_current


    def calc_reverse_current(self, normalize=False):
        reverse_current = elementary_charge * self.carrier_concentration * self.half_neck * 2 * self.mean_vel()
        normalize_factor = (2 * elementary_charge * self.carrier_concentration * self.fermi_vel * self.half_neck / pi)
        if normalize:
            return reverse_current / normalize_factor
        else:
            return reverse_current


    def calc_total_current(self, normalize=False):
        total_current = self.calc_direct_current() - self.calc_reverse_current()
        normalize_factor = (2 * elementary_charge * self.carrier_concentration * self.fermi_vel * self.half_neck / pi)
        if normalize:
            return total_current / normalize_factor
        else:
            return total_current


if __name__ == '__main__':
    mob = 4
    volt = 0
    v0 = 1
    carrier_c = 1
    neck_h = 10
    length = 5
    alpha = deg2rad(45)
    shoulder = neck_h + length / tan(alpha)
    e_field = volt / length
    drift_vel = mob * e_field
    g = drift_vel / v0
    analitic = Analytical(v0, neck_h, length, alpha, volt, carrier_c, mob, 1000)
    direct_curr = analitic.calc_direct_current(True)
    reverse_curr = analitic.calc_reverse_current(True)
    total_curr = analitic.calc_total_current(True)
    print(f'Direct current: {direct_curr}')
    print(f'Reverse current: {reverse_curr}')
    print(f'Total current: {total_curr}')
