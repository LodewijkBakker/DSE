from urllib.request import urlopen
from skyfield.api import load

def Aqua_TLE():
    text = [i.decode('UTF-8') for i in urlopen("https://celestrak.org/NORAD/elements/active.txt").readlines()]
    for i in range(0, len(text), 3):
        if (text[i].split()[0] == "AQUA"):
            l1 = text[i + 1]
            l2 = text[i + 2]
            break
    return l1, l2
l1, l2 = Aqua_TLE()
sat = Satrec.twoline2rv(l1, l2)
params = {"inc": float(l2[8:16].strip()),  # Inclination
          # Right Ascension of Ascending Node
          "RAAN": float(l2[17:25].strip()),
          "ecc": float("0."+l2[26:33].strip()),  # Eccentricity
          "AoP": float(l2[34:42].strip()),  # Argument of perigee
          "M": float(l2[43:51].strip()),  # Mean anomally
          "n": float(l2[52:64].strip()),  # Mean motion
          "Rn": int(l2[64:69].strip())}  # Revolutions at epoch
sgp4(l1,)
ts = load.timescale()     # Timescale object
t = ts.now()   # Time object
RA_hours = t.gmst - west_longitude_degrees/15.