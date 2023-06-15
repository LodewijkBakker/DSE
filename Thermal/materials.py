class Material:
    def __init__(self, params=None):
        if params is not None:
            self.thermal_conductivity = params['thermal_conductivity']
            self.thermal_capacitance = params['thermal_capacitance']


class Coating:
    def __init__(self, params=None):
        if params is not None:
            self.absorptivity = params['absorptivity']
            self.emissivity = params['emissivity']

# MATERIALS
def aluminium() -> 'Material':
    return Material({'thermal_conductivity': 239, 'thermal_capacitance': 897})

def OSR_mat() -> 'Material':
    return Material({'thermal_conductivity': 1000, 'thermal_capacitance': 897})

def battery() -> 'Material':
    return Material({'thermal_conductivity': 20, 'thermal_capacitance': 1000})

def solar_panel() -> 'Material':
    return Material({'thermal_conductivity': 168, 'thermal_capacitance': 900})

def MLI_mat() -> 'Material':
    return Material({'thermal_conductivity': 5, 'thermal_capacitance': 1009})

def PCB() -> 'Material':
    return Material({'thermal_conductivity': 18, 'thermal_capacitance': 1150})

def hastealloy():
    return Material({'thermal_conductivity': 10.5, 'thermal_capacitance': 427})

def low_conductance_tape():
    return Material({'thermal_conductivity': 3, 'thermal_capacitance': 1000})



# COATINGS
def white_paint() -> 'Coating':
    return Coating({'absorptivity': 0.17, 'emissivity': 0.92})

def black_paint() -> 'Coating':
    return Coating({'absorptivity': 0.94, 'emissivity': 0.94})

def bare_Al() -> 'Coating':
    return Coating({'absorptivity': 0.25, 'emissivity': 0.15})

def MLI_coat() -> 'Coating':
    return Coating({'absorptivity': 0.14, 'emissivity': 0.035})

def bare_solar_panel() -> 'Coating':
    return Coating({'absorptivity': 0.42, 'emissivity': 0.52})

def coated_solar_panel() -> 'Coating':
    return Coating({'absorptivity': (0.06 + 0.42)/2, 'emissivity': (0.52 + 0.88)/2})

def OSR_coat() -> 'Coating':
    return Coating({'absorptivity': 0.06, 'emissivity': 0.9})

def thermal_tape_1() -> 'Coating':
    return Coating({'absorptivity': 0.2, 'emissivity': 0.3})

