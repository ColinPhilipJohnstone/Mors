
# The basic star class that people should be using
from .star import Star

# Parameters
from .parameters import paramsDefault , PrintParams , NewParams

# Useful function for stellar evolution stuff
from .stellarevo import StellarEvoModel , Value , LoadTrack

# Stellar evolution basic quantities
from .stellarevo import ( Rstar , Lbol , Teff , Itotal , Icore , Ienv , Mcore , Menv , Rcore , tauConv ,
                         dItotaldt , dIcoredt , dIenvdt , dIenvdt , dMcoredt , dRcoredt )

# Functions from physical model
#from .physicalmodel import Bdip , Mdot , OmegaBreak

# Physical quantities and default parameters from the physical model
#from .physicalmodel import OmegaSun 

# Rotational evolution stuff
from .rotevo import EvolveRotation , EvolveRotationStep

