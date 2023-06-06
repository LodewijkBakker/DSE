class Material:
    def __init__(self, params=None):
        if params is not None:
            self.thermal_conductivity = params['thermal_conductivity']
            self.thermal_capacitance = params['thermal_capacitance']
            self.absorptivity = params['absorptivity']
            self.emissivity = params['emissivity']

    def aluminium(self) -> 'Material':
        return Material({'thermal_conductivity': 239, 'thermal_capacitance': 897, 'absorptivity': 0.25, 'emissivity': 0.15})

    def OSR(self) -> 'Material':
        return Material({'thermal_conductivity': 4000, 'thermal_capacitance': 897, 'absorptivity': 0.06, 'emissivity': 0.88})

    def battery(self) -> 'Material':
        return Material({'thermal_conductivity': 20, 'thermal_capacitance': 1000, 'absorptivity': 0.14, 'emissivity': 0.035})

    def solar_panel(self) -> 'Material':
        return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900, 'absorptivity': 0.42, 'emissivity': 0.52})

    def coated_solar_panel(self) -> 'Material':
        return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900, 'absorptivity': 0.42, 'emissivity': 0.83})

    def MLI(self) -> 'Material':
        return Material({'thermal_conductivity': 5, 'thermal_capacitance': 1009, 'absorptivity': 0.14, 'emissivity': 0.035})

    def white_coating(self) -> 'Material':
        return Material({'thermal_conductivity': 210, 'thermal_capacitance': 897, 'absorptivity': 0.27, 'emissivity': 0.83})

    def PCB(self) -> 'Material':
        return Material({'thermal_conductivity': 18, 'thermal_capacitance': 1150, 'absorptivity': 0.14, 'emissivity': 0.035})
