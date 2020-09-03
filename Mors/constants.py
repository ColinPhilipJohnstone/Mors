
"""Module for holding some basic constants."""

# Imports for Mors modules
import Mors.stellarevo as SE

# Time units
sec = 1.0                   # seconds
minute = 60.0*sec             # 
hour = 60.0*minute        # 
day = 24.0*hour           # 
year = 365.24*day
Myr = 1.0e6*year
Gyr = 1.0e3*Myr

# Distance units in cm
A = 1.0e-8                # Angstrom
nm = 1.0e-7               # nanometer
micron = 1.0e-4           # micrometer
mm = 0.1                  # milimeter
cm = 1.0                  # centimeter
km = 1.0e5                # kilometer
m = 1.0e3                 # meter
Rearth = 6.371e8          # Radius of Earth
Rjup = 6.9911e9           # Radius of Jupiter
Rsun = 6.957e10           # cm

# Area units
cm2 = 1.0                 # cm^2

# Mass units in g
g = 1.0                   # gram
kg = 1.0e3                # kilogram
Mearth = 5.972e27         # Mass of Earth
Mjup = 1.89813e30         # Mass of Juputer
Msun = 1.99e33            # Mass of Sun
Mproton = 1.6726219e-24   # Mass of proton
Melec = 9.10938356e-28    # Mass of electron

# Thermal units
K = 1.0                   # Kelvin
erg = 1.0                 # erg
ergs_cm_3 = 1.0           # erg s^-1 cm^-3

# Densities
cm_3 = 1.0                # cm^-3

# Speeds
kms_ = 1.0e-5             # km s^-1 in cm s^-1

# Radiation units
ergs_cm_2 = 1.0           # erg s^-1 cm^-2

# Some other things
Msunyr_ = Msun / year      # g s^-1

# Other solar quantities
LbolSun = 3.828e33        # erg s^-1
AgeSun = 4567.0           # Myr
OmegaSun = 2.67e-6        # rad s^-1 - our assumed solar rotation
ProtSun = 2.0*3.142 / ( OmegaSun ) / 86400.0    # days - corresponding solar rotation period
#tauConvSun = SE.tauConv(1.0,AgeSun)             # days - the Sun's convective turnover time using stellar evo models
tauConvSun = 28.436   # days - the Sun's convective turnover time using stellar evo models of Spada et al. (2013)
RoSun = ProtSun / tauConvSun                    # Sun's current Rossby number using stellar evo models

# Standard constants
GravConstant = 6.674e-8   # cm^3 g^-1 s^-2
kB = 1.38064852e-16       # erg K^-1
Pi = 3.14159265359        # you know what this is
