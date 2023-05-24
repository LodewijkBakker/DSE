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
        return Material({'thermal_conductivity': 20, 'thermal_capacitance': 1000, 'absorptivity': 0.1, 'emissivity': 0.1})

    def solar_panel(self) -> 'Material':    # https://apps.dtic.mil/sti/pdfs/AD1170386.pdf
        return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900, 'absorptivity': 0.42, 'emissivity': 0.52})

    def kapton2mm(self) -> 'Material':
        return Material({'thermal_conductivity': 210, 'thermal_capacitance': 1170, 'absorptivity': 0.23, 'emissivity': 0.86})

    def white_coating(self) -> 'Material':
        return Material({'thermal_conductivity': 210, 'thermal_capacitance': 897, 'absorptivity': 0.06, 'emissivity': 0.88})

    def PCB(self) -> 'Material':
        return Material({'thermal_conductivity': 18, 'thermal_capacitance': 1150, 'absorptivity': 0.1, 'emissivity': 0.1})
