

class Material:
    def __init__(self, relax_time, scalar_fermi_velocity, mass_scale=1, permittivity=1, permeability=1):
        self.relax_time = relax_time
        self.scalar_fermi_velocity = scalar_fermi_velocity
        self.mass_scale = mass_scale
        self.permittivity = permittivity
        self.permeability = permeability
