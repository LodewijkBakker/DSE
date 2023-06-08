class Material:
    def __init__(self, params=None):
        if params is not None:
            self.thermal_conductivity = params['thermal_conductivity']
            self.thermal_capacitance = params['thermal_capacitance']

    def aluminium(self) -> 'Material':
        return Material({'thermal_conductivity': 239, 'thermal_capacitance': 897})

    def OSR(self) -> 'Material':
        return Material({'thermal_conductivity': 4000, 'thermal_capacitance': 897})

    def battery(self) -> 'Material':
        return Material({'thermal_conductivity': 20, 'thermal_capacitance': 1000})

    def solar_panel(self) -> 'Material':
        return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900})

    def coated_solar_panel(self) -> 'Material':
        return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900})

    def MLI(self) -> 'Material':
        return Material({'thermal_conductivity': 5, 'thermal_capacitance': 1009})

    def PCB(self) -> 'Material':
        return Material({'thermal_conductivity': 18, 'thermal_capacitance': 1150})


class Coating:
    def __init__(self, params=None):
        if params is not None:
            self.absorptivity = params['absorptivity']
            self.emissivity = params['emissivity']

    def white_paint(self) -> 'Coating':
        return Coating({'absorptivity': 0.27, 'emissivity': 0.83})

    def black_paint(self) -> 'Coating':
        return Coating({'absorptivity': 0.94, 'emissivity': 0.94})

    def bare_Al(self) -> 'Coating':
        return Coating({'absorptivity': 0.25, 'emissivity': 0.15})

    def MLI(self) -> 'Coating':
        return Coating({'absorptivity': 0.14, 'emissivity': 0.035})

    def bare_solar_panel(self) -> 'Coating':
        return Coating({'absorptivity': 0.42, 'emissivity': 0.52})

    def coated_solar_panel(self) -> 'Coating':
        return Coating({'absorptivity': 0.42, 'emissivity': 0.83})

    def OSR(self) -> 'Coating':
        return Coating({'absorptivity': 0.06, 'emissivity': 0.80})

