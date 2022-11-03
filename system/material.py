

class Material:
    def __init__(self, relax_time, fermi_velocity, mass_scale=1, permittivity=1, permeability=1):
        self.relax_time = relax_time
        self.fermi_velocity = fermi_velocity
        self.mass_scale = mass_scale
        self.permittivity = permittivity
        self.permeability = permeability
