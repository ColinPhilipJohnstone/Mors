
"""Module for holding default parameters and related functions."""

# Imports for standard stuff needed here
import copy
import numpy as np

# Make empty dictionary to hold default values for all parameters in the code. Some of these
# are parameters that the user might want to change in order to determine the behavior of the code.
paramsDefault = {}

def SetDefaultParameters(paramsDefault):
    
    """Sets up the default parameters."""
    
    # Time integration solver parameters
    paramsDefault['TimeIntegrationMethod'] = 'Rosenbrock'     # options are 'ForwardEuler', 'RungeKutta4', 'RungeKuttaFehlberg', 'Rosenbrock'
    paramsDefault['nStepMax'] = 10**6                         # maximum number of timesteps to allow before error
    paramsDefault['DeltaDesired'] = 1.0e-5                    # desired Delta value to use in numerical solvers such as Runge-Kutta-Fehlberg and Rosenbrock
    paramsDefault['AgeMinDefault'] = 1.0                      # Myr - default age to start evolutionary tracks
    paramsDefault['AgeMaxDefault'] = 5000.0                   # Myr - default age to end evolutionary tracks
    paramsDefault['deltaJac'] = 1.0e-5                        # pertubation to use when calculating Jacobian for Rosenbrock solver
    paramsDefault['CoefficientsRB'] = _CoefficientsRB()       # all coefficents needed for the Rosenbrock solver
    paramsDefault['dAgeMin'] = 1.0e-5                         # Myr - minimum timestep to allow
    paramsDefault['dAgeMax'] = 5000.0                         # Myr - maximum timestep to allow
    
    # Physical processes to include
    paramsDefault['CoreEnvelopeDecoupling'] = True            # should core envelope decoupling be included
    paramsDefault['MomentInertiaChangeTorque'] = True         # should the code include effects of moment of inertia changes on rotation (do not turn this off)
    paramsDefault['CoreGrowthTorque'] = True                  # should the code include effects of core-growth on core and envelope rotation (do not turn this off)
    paramsDefault['WindTorque'] = True                        # should the code include wind spin-down torque
    paramsDefault['CoreEnvelopeTorque'] = True                # should the code include core-envelope coupling torque
    paramsDefault['DiskLocking'] = True                       # should disk-locking be used
    
    # Core envelope decoupling parameters
    paramsDefault['IcoreThresholdCE'] = 0.01                  # value of Icore/Itotal below which core-envelope decoupling is not used
    paramsDefault['MstarThresholdCE'] = 0.35                  # Msun - value of Mstar below which core-envelope decoupling is not used
    
    # Disk locking parameters
    paramsDefault['aDiskLock'] = 13.5                         # a in tDiskLock = a * Omega^b
    paramsDefault['bDiskLock'] = -0.5                         # b in tDiskLock = a * Omega^b
    paramsDefault['ageDLmax'] = 15.0                          # Myr - maximum disk locking age
    
    # Wind parameters
    paramsDefault['Kwind'] = 11.0                             # additional constant in torque formula
    paramsDefault['BdipSun'] = 1.35                           # G - estimate for rotation models from Johnstone et al. (2015)
    paramsDefault['aBdip'] = -1.32                            # from Vidotto et al. (2014)
    paramsDefault['RoSatBdip'] = 0.0605                       # saturation Rossby number for dipole field strength from Johnstone et al. (2020) using Spada et al. (2013) models
    paramsDefault['MdotSun'] = 1.4e-14                        # Msun yr^-1 - solar mass loss rate
    paramsDefault['aMdot'] = -1.7591                          # power law index on Ro in mass loss rate equation from Johnstone et al. (2020)
    paramsDefault['bMdot'] = 0.6494                           # power law index on Mstar in mass loss rate equation from Johnstone et al. (2020)
    paramsDefault['RoSatMdot'] = 0.0605                       # saturation Rossby number for mass loss rate from Johnstone et al. (2020) using Spada et al. (2013) models
    paramsDefault['BreakupMdotIncrease'] = True               # should we increase Mdot as approaching breakup
    paramsDefault['fracBreakThreshold'] = 0.1                 # fraction of breakup to start increasing Mdot
    
    # Core-envelope coupling stuff
    paramsDefault['aCoreEnvelope'] = 25.6015                  # a in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
    paramsDefault['bCoreEnvelope'] = -3.4817E-002             # b in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
    paramsDefault['cCoreEnvelope'] = -0.4476                  # c in tCoreEnvelope = a * (|OmegaCore-OmegaEnv|)^b * Mstar^d
    paramsDefault['timeCEmin'] = 0.0                          # minimum value of core-envelope coupling timescale
    
    # Output parameters
    paramsDefault['ExtendedTracks'] = False                   # should the code return evolutionary tracks for all quantities 
    
    # XUV emission parameters
    paramsDefault['RoSatXray'] = 0.0605                       # saturation Rossby number for X-ray emission from Johnstone et al. (2020) using models of Spada et al. (2013)
    paramsDefault['RxSatXray'] = 5.135e-4                     # saturation Rx
    paramsDefault['beta1Xray'] = -0.135                       # beta from Rx = C * Ro^beta for the saturated regime
    paramsDefault['beta2Xray'] = -1.889                       # beta from Rx = C * Ro^beta for the unsaturated regime
    paramsDefault['C1Xray'] = paramsDefault['RxSatXray'] / paramsDefault['RoSatXray']**paramsDefault['beta1Xray'] # C from Rx = C * Ro^beta for the saturated regime
    paramsDefault['C2Xray'] = paramsDefault['RxSatXray'] / paramsDefault['RoSatXray']**paramsDefault['beta2Xray'] # C from Rx = C * Ro^beta for the unsaturated regime
    paramsDefault['sigmaXray'] = 0.359                        # dex - standard deviation of normal PDF for X-ray scatter
    
    # Rotation fitting parameters
    paramsDefault['AgeMinFit'] = 1.0                          # Myr - starting age for fitting rotation track
    paramsDefault['Omega0FitMin'] = 0.1                       # OmegaSun - minimum starting rotation rate to consider when fitting initial rotation
    paramsDefault['Omega0FitMax'] = 50.0                      # OmegaSun - maximum starting rotation rate to consider when fitting initial rotation
    paramsDefault['nStepMaxFit'] = 1000                       # maximum number of steps to take when fitting rotation rate
    paramsDefault['toleranceFit'] = 1.0e-5                    # tolerance for fitting initial rotation rate
    
    # Other
    paramsDefault['dMstarPer'] = 0.1                          # Msun - half width of bins for working out percentiles of distribution
    paramsDefault['AgeHZ'] = 5000.0                           # Myr - age to use for calculating habitable zone boundaries
    
    return paramsDefault

def _CoefficientsRB():
    
    """Returns dictionary with basic coefficients for Rosenbrock solver."""

    # Using here values from RODAS3 of Sandu et al. (1997)
    
    # Set s value (i.e. number of steps)
    s = 4
    
    # Set alpha_ij
    alpha = np.zeros((s,s))
    alpha[2,0] = 1.0
    alpha[3,0] = 3.0/4.0
    alpha[3,1] = -1.0/4.0
    alpha[3,2] = 1.0/2.0
    
    # Set gamma_ij
    gamma = np.zeros((s,s))
    gamma[0,0] = 1.0/2.0
    gamma[1,0] = 1.0
    gamma[2,0] = -1.0/4.0
    gamma[3,0] = 1.0/12.0
    gamma[1,1] = 1.0/2.0
    gamma[2,1] = -1.0/4.0
    gamma[3,1] = 1.0/12.0
    gamma[2,2] = 1.0/2.0
    gamma[3,2] = -2.0/3.0
    gamma[3,3] = 1.0/2.0

    # Set b_i
    b = np.zeros(s)
    b[0] = 5.0/6.0
    b[1] = -1.0/6.0
    b[2] = -1.0/6.0
    b[3] = 1.0/2.0
    
    # Set b2_i
    b2 = np.zeros(s)
    b2[0] = 3.0/4.0
    b2[1] = -1.0/4.0
    b2[2] = 1.0/2.0
    b2[3] = 0.0
    
    # Set the order of the scheme
    order = 3
    
    # Put these into a dictionary for returning
    CoefficientsRB = {}
    CoefficientsRB['s'] = s
    CoefficientsRB['alpha'] = alpha
    CoefficientsRB['gamma'] = gamma
    CoefficientsRB['b'] = b
    CoefficientsRB['b2'] = b2
    CoefficientsRB['order'] = order
    
    return CoefficientsRB

def PrintParams(params=paramsDefault):
    
    """Takes parameter dictionary, prints values to screen."""
    
    # Loop over parameters and print each
    for param in params:
        print(param , " = " , params[param] )
    
    return

def NewParams(**kwargs):
    
    """Takes list of parameters, returns new parameter dictionary using defaults for non-specified parameters."""
    
    # Make new parameter dictionary based on paramsDefault
    params = copy.deepcopy(paramsDefault)
    
    # Change parameters based on user specifications
    for param in kwargs:
        params[param] = kwargs[param]
    
    return params

# Setup the default parameter dictionary
paramsDefault = SetDefaultParameters(paramsDefault)
