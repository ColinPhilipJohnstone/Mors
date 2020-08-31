
"""Module for holding default parameters and related functions."""

import copy

#==================================================================================================================

# The following dictionary hold default values for all parameters in the code that the user might want to 
# change in order to determine the behavior of the code.

# Make empty dictionary
paramsDefault = {}

# Time integration parameters
paramsDefault['TimeIntegrationMethod'] = 'RungeKutta4'    # options are 'ForwardEuler', 'RungeKutta4',...
paramsDefault['nStepMax'] = 10**6                         # maximum number of timesteps to allow before error
paramsDefault['DeltaDesired'] = 1.0e-5                    # desired Delta value to use in numerical solvers such as Runge-Kutta-Fehlberg and Rosenbrock

# Physical processes to include
paramsDefault['CoreEnvelopeDecoupling'] = True      # should core envelope decoupling be included
paramsDefault['MomentInertiaChangeTorque'] = True   # should the code include effects of moment of inertia changes on rotation (do not turn this off)
paramsDefault['CoreGrowthTorque'] = True            # should the code include effects of core-growth on core and envelope rotation (do not turn this off)
paramsDefault['WindTorque'] = True                  # should the code include wind spin-down torque
paramsDefault['CoreEnvelopeTorque'] = True          # should the code include core-envelope coupling torque
paramsDefault['DiskLocking'] = True                 # should disk-locking be used

# Core envelope decoupling parameters
paramsDefault['IcoreThresholdCE'] = 0.01        # value of Icore/Itotal below which core-envelope decoupling is not used
paramsDefault['MstarThresholdCE'] = 0.35        # Msun - value of Mstar below which core-envelope decoupling is not used

# Disk locking parameters
paramsDefault['aDiskLock'] = 13.5   # a in tDiskLock = a * Omega^b
paramsDefault['bDiskLock'] = -0.5   # b in tDiskLock = a * Omega^b
paramsDefault['ageDLmax'] = 15.0    # Myr - maximum disk locking age

# Wind parameters
paramsDefault['Kwind'] = 11.0                 # additional constant in torque formula
paramsDefault['BdipSun'] = 1.35               # G - estimate for rotation models from Johnstone et al. (2015)
paramsDefault['aBdip'] = -1.32                # from Vidotto et al. (2014)
paramsDefault['RoSatBdip'] = 0.0605           # saturation Rossby number for dipole field strength from Johnstone et al. (2020) using Spada et al. (2013) models
paramsDefault['MdotSun'] = 1.4e-14            # Msun yr^-1 - solar mass loss rate
paramsDefault['aMdot'] = -1.7591              # power law index on Ro in mass loss rate equation from Johnstone et al. (2020)
paramsDefault['bMdot'] = 0.6494               # power law index on Mstar in mass loss rate equation from Johnstone et al. (2020)
paramsDefault['RoSatMdot'] = 0.0605           # saturation Rossby number for mass loss rate from Johnstone et al. (2020) using Spada et al. (2013) models
paramsDefault['BreakupMdotIncrease'] = True   # should we increase Mdot as approaching breakup
paramsDefault['fracBreakThreshold'] = 0.1     # fraction of breakup to start increasing Mdot

# Core-envelope coupling stuff
paramsDefault['aCoreEnvelope'] = 25.6015        # a in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
paramsDefault['bCoreEnvelope'] = -3.4817E-002   # b in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
paramsDefault['cCoreEnvelope'] = -0.4476        # c in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
paramsDefault['timeCEmin'] = 0.0                # minimum value of core-envelope coupling timescale

#paramsDefault[''] =   # 
#paramsDefault[''] =   # 
#paramsDefault[''] =   # 
#paramsDefault[''] =   # 
#paramsDefault[''] =   # 


#==================================================================================================================

def PrintParams(params=paramsDefault):
  
  """Takes parameter dictionary, prints values to screen."""
  
  # Loop over parameters and print each
  for param in params:
    print(param , " = " , params[param] )
  
  return

#==================================================================================================================

def NewParams(**kwargs):
  
  """Takes list of parameters, returns new parameter dictionary using defaults for non-specified parameters."""
  
  # Make new parameter dictionary based on paramsDefault
  params = copy.deepcopy(paramsDefault)
  
  # Change parameters based on user specifications
  for param in kwargs:
    params[param] = kwargs[param]
  
  return params

#==================================================================================================================

def SetParams(**kwargs):
  
  """
  Takes set of keyword parameters, returns parameter dictionary.
  """
  
  return params

#==================================================================================================================
