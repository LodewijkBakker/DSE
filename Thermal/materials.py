class Material:
    def __init__(self, params=None):
        if params is not None:
            self.thermal_conductivity = params['thermal_conductivity']
            self.thermal_capacitance = params['thermal_capacitance']
            self.absorptivity = params['absorptivity']
            self.emissivity = params['emissivity']

    def aluminium(self) -> 'Material':
        return Material({'thermal_conductivity': 210, 'thermal_capacitance': 897, 'absorptivity': 0.1, 'emissivity': 0.07})

    def anodized_aluminium(self) -> 'Material':
        return Material({'thermal_conductivity': 210, 'thermal_capacitance': 897, 'absorptivity': 0.86, 'emissivity': 0.86})

    def carbon_fibre(self) -> 'Material':
        return Material({'thermal_conductivity': 200, 'thermal_capacitance': 800, 'absorptivity': 0.11, 'emissivity': 0.5})

    def battery(self) -> 'Material':
        return Material({'thermal_conductivity': 20, 'thermal_capacitance': 1000, 'absorptivity': 0, 'emissivity': 0})

    def solar_panel(self) -> 'Material':
        return Material({'thermal_conductivity': 150, 'thermal_capacitance': 712, 'absorptivity': 0.8, 'emissivity': 0.8})
