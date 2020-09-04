
"""Module for holding functions for calculating physical quantities for the physical rotation and activity model."""

# Imports for standard stuff needed here
import numpy as np
import sys as sys
import copy as copy

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.stellarevo as SE
import Mors.constants as const
import Mors.parameters as params

#====================================================================================================================

def dOmegadt(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns rates of change of rotation.
  
  This is the main function to be used by the solvers to get the rates of change of the core and envelope 
  rotation rates given basic stellar parameters. The actual calculations of these are done in RotationQuantities() 
  which itself also returns a very large number of other things used during the calculation of the rates 
  of change.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  OmegaEnv : float
      Envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaCore : float
      Core rotation rate in OmegaSun (=2.67e-6 rad/s).
  params : dict , optional
      Dictionary holding model parameters. 
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  dOmegaEnvdt : float
      Rate of change of envelope rotation rate in OmegaSun Myr^-1.
  dOmegaCoredt : float
      Rate of change of core rotation rate in OmegaSun Myr^-1.
  
  """
  
  # Make sure parameters are set
  if Mstar is None: 
    misc._PrintErrorKill("argument Mstar must be set in call to function")
  if Age is None: 
    misc._PrintErrorKill("argument Age must be set in call to function")
  if OmegaEnv is None: 
    misc._PrintErrorKill("argument OmegaEnv must be set in call to function")
  if OmegaCore is None: 
    misc._PrintErrorKill("argument OmegaCore must be set in call to function")
  
  # If an instance of the StarEvo class is not input, load it with defaults
  if StarEvo is None:
    StarEvo = SE.StarEvo()
  
  # Get state of system 
  StarState = RotationQuantities(Mstar=Mstar,Age=Age,OmegaEnv=OmegaEnv,OmegaCore=OmegaCore,params=params,StarEvo=StarEvo)
  
  return ( StarState['dOmegaEnvdt'] , StarState['dOmegaCoredt'] )

#====================================================================================================================

def RotationQuantities(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns rotation quantities including rates of change of rotation.
  
  This is the main function to be used to get the rates of change of the core and envelope rotation 
  rates given basic stellar parameters. The actual calculations of these are done in GetState() which
  itself also returns a very large number of other things used during the calculation of the rates 
  of change.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  OmegaEnv : float
      Envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaCore : float
      Core rotation rate in OmegaSun (=2.67e-6 rad/s).
  params : dict , optional
      Dictionary holding model parameters. 
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  StarState : dict
      Set of rotation quantities.
  
  """
  
  # Make sure parameters are set
  if Mstar is None: 
    misc._PrintErrorKill("argument Mstar must be set in call to function")
  if Age is None: 
    misc._PrintErrorKill("argument Age must be set in call to function")
  if OmegaEnv is None: 
    misc._PrintErrorKill("argument OmegaEnv must be set in call to function")
  if OmegaCore is None: 
    misc._PrintErrorKill("argument OmegaCore must be set in call to function")
  
  #--------------------------------------------------------
  
  # Setup stellar parameter dictionary
  # This will be used to hold all of the 
  
  # If an instance of the StarEvo class is not input, load it with defaults
  if StarEvo is None:
    StarEvo = SE.StarEvo()
  
  # Start dictionary
  StarState = {}
  StarStateUnits = {}
  
  # Add basic input parameters 
  StarState['Mstar'] = Mstar
  StarState['Age'] = Age
  StarState['OmegaEnv'] = OmegaEnv
  StarState['OmegaCore'] = OmegaCore
  
  # Get rotation period from surface rotation in days
  StarState['Prot'] = _Prot(OmegaEnv)
  
  # Radius in cm
  StarState['Rstar'] = StarEvo.Rstar(Mstar,Age)
  
  # Convective turnover time in days
  StarState['tauConv'] = StarEvo.tauConv(Mstar,Age)
  
  # Rossby number
  StarState['Ro'] = _Ro(StarState)
  
  # Get moments of interia in g cm^2
  StarState['Itotal'] = StarEvo.Itotal(Mstar,Age)
  StarState['Ienv'] = StarEvo.Ienv(Mstar,Age)
  StarState['Icore'] = StarEvo.Icore(Mstar,Age)
  
  # Rates of change of moment of inertia in g cm^2 s^-1
  StarState['dItotaldt'] = StarEvo.dItotaldt(Mstar,Age) / const.Myr
  StarState['dIenvdt'] = StarEvo.dIenvdt(Mstar,Age) / const.Myr
  StarState['dIcoredt'] = StarEvo.dIcoredt(Mstar,Age) / const.Myr
  
  # Rcore in cm
  StarState['Rcore'] = StarEvo.Rcore(Mstar,Age) * const.Rsun
  
  # dMcoredt in g s^-1
  StarState['dMcoredt'] = StarEvo.dMcoredt(Mstar,Age) * const.Msun / const.Myr
  
  # Mass loss rate in g s^-1
  StarState['Mdot'] = _Mdot(StarState,params=params) * const.Msunyr_
  
  # Get dipole field strength in G
  StarState['Bdip'] = _Bdip(StarState,params=params)
  
  # Get surface escape velocity in cm s^-1
  StarState['vEsc'] = _vEsc(Mstar,StarState['Rstar'])
  
  #--------------------------------------------------------
  
  # Determine if doing core-envelope decoupling
  CoreEnvelopeDecoupling = _shouldCoreEnvelopeDecoupling( StarState , params=params )
  
  #--------------------------------------------------------
  
  # Add moment of interia change torque
  if params['MomentInertiaChangeTorque']:
    StarState['torqueEnvMom'] , StarState['torqueCoreMom'] = _torqueMoment( StarState , CoreEnvelopeDecoupling )
  else:
    StarState['torqueEnvMom'] = StarState['torqueCoreMom'] = 0.0
  
  # Add core-growth torque
  if params['CoreGrowthTorque'] and CoreEnvelopeDecoupling:
    StarState['torqueEnvCG'] , StarState['torqueCoreCG'] = _torqueCoreGrowth( StarState )
  else:
    StarState['torqueEnvCG'] = StarState['torqueCoreCG'] = 0.0
  
  # Add wind spin-down torque
  if params['WindTorque']:
    StarState['torqueEnvWind'] = _torqueWind( StarState , params=params )
  else:
    StarState['torqueEnvWind'] = 0.0
  
  # Add wind spin-down torque
  if params['CoreEnvelopeTorque'] and CoreEnvelopeDecoupling:
    StarState['torqueEnvCE'] , StarState['torqueCoreCE'] = _torqueCoreEnvelope( StarState , params=params )
  else:
    StarState['torqueEnvCE'] = StarState['torqueCoreCE'] = 0.0
  
  # Add disk locking torque
  if params['DiskLocking']:
    StarState['torqueEnvDL'] = _torqueDiskLocking( StarState , CoreEnvelopeDecoupling , params=params )
  else:
    StarState['torqueEnvDL']  = 0.0
  
  #--------------------------------------------------------
  # Envelope rotation rate of change
  
  # Get total torque on envelope
  StarState['torqueEnv'] = StarState['torqueEnvMom'] + StarState['torqueEnvCG'] + StarState['torqueEnvWind'] + StarState['torqueEnvCE'] + StarState['torqueEnvDL']
  
  # Get rates of change of rotation in rad s^-2
  StarState['dOmegaEnvdt'] = StarState['torqueEnv'] / StarState['Ienv']
  
  # Get envelope rates of change in OmegaSun Myr^-1
  StarState['dOmegaEnvdt'] *= const.Myr / const.OmegaSun
  
  #--------------------------------------------------------
  # Core rotation rate of change
  
  # Have core spin-up with envelope if not developed yet
  if CoreEnvelopeDecoupling:
    
    # Get total torque on core
    StarState['torqueCore'] = StarState['torqueCoreMom'] + StarState['torqueCoreCG'] + StarState['torqueCoreCE']
    
    # Get rates of change of rotation in rad s^-2
    StarState['dOmegaCoredt'] = StarState['torqueCore'] / StarState['Icore']
    
    # Get envelope rates of change in OmegaSun Myr^-1
    StarState['dOmegaCoredt'] *= const.Myr / const.OmegaSun
    
  else:
    
    # There should be no torque on the core and it should spin up or down with the envelope
    StarState['torqueCore'] = 0.0
    StarState['dOmegaCoredt'] = StarState['dOmegaEnvdt']
  
  #--------------------------------------------------------
  
  return StarState

#====================================================================================================================

def ExtendedQuantities(StarState=None,Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns extended set of physical quantities not calculated by RotationQuantities.
  
  This is the main function to be used to get things like stellar XUV emission. The user can either input as arguments the 
  dictionary StarState, which is returned by RotationQuantities, or the basic stellar parameters Mstar, Age, OmegaEnv, and
  OmegaCore (final one not necessary). If StarState is not input, this function will first call RotationQuantities to get it.
  Either way, the return will be the same dictionary with additional parameters added to it.
  
  Parameters
  ----------
  StarState : dict , optional
      Set of quantities already calculated by RotationQuantities().
  Mstar : float , optional
      Mass of star in Msun.
  Age : float , optional
      Age of star.
  OmegaEnv : float 
      Envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaCore : float 
      Core rotation rate in OmegaSun (=2.67e-6 rad/s).
  params : dict , optional
      Dictionary holding model parameters. 
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  StarState : dict
      Set of extended quantities.
  
  """
  
  # If an instance of the StarEvo class is not input, load it with defaults
  if StarEvo is None:
    StarEvo = SE.StarEvo()
    
  # Get StarState if it was not already set
  if StarState is None: 
    
    # Make sure required parameters were set
    if Mstar is None: 
      misc._PrintErrorKill("argument Mstar must be set in call to function")
    if Age is None: 
      misc._PrintErrorKill("argument Age must be set in call to function")
    if OmegaEnv is None: 
      misc._PrintErrorKill("argument OmegaEnv must be set in call to function")
    if OmegaCore is None: 
      misc._PrintErrorKill("argument OmegaCore must be set in call to function")
    
    # Get rotation quantities
    StarState = RotationQuantities(Mstar=Mstar,Age=Age,OmegaEnv=OmegaEnv,OmegaCore=OmegaCore,params=params,StarEvo=StarEvo)
  
  # Get Lbol in erg s^-1
  StarState['Lbol'] = StarEvo.Lbol(Mstar,Age) * const.LbolSun
  
  # Get effective temperature in K
  StarState['Teff'] = StarEvo.Teff(Mstar,Age)
  
  # Get Lx in erg s^-1, Fx in erg s^-1 cm^-2, and Rx
  StarState['Lx'] , StarState['Fx'] , StarState['Rx'] = _Xray(StarState,params=params)
  
  # Get Tcor in MK
  StarState['Tcor'] = _Tcor(StarState,params=params)
  
  # Get Leuv1 in erg s^-1, Feuv1 in erg s^-1 cm^-2, and Reuv1 (10-36 nm)
  StarState['Leuv1'] , StarState['Feuv1'] , StarState['Reuv1'] = _EUV1(StarState,params=params)
  
  # Get Leuv2 in erg s^-1, Feuv2 in erg s^-1 cm^-2, and Reuv2 (36-92 nm)
  StarState['Leuv2'] , StarState['Feuv2'] , StarState['Reuv2'] = _EUV2(StarState,params=params)
  
  # Get Leuv in erg s^-1, Feuv in erg s^-1 cm^-2, and Reuv
  StarState['Leuv'] , StarState['Feuv'] , StarState['Reuv'] = _EUV(StarState,params=params)
  
  # Get Lly in erg s^-1, Fly in erg s^-1 cm^-2, and Rly
  StarState['Lly'] , StarState['Fly'] , StarState['Rly'] = _Lymanalpha(StarState,params=params)
  
  # Get habitable zone orbital distance
  
  # Get FxHZ, FeuvHZ, and Fla
  
  # Get rate of flares above given energy
  
  
  
  return StarState 

#====================================================================================================================

def QuantitiesUnits(StarState=None):
  
  """
  Returns dictionary of units for each physical quantity in calculated by RotationQuantities() and ExtendedQuantities().
  
  For output purposes, this is essential since it gives the user a way to know the units of the quantities output. It is 
  recommended to input also 'StarState' which will cause the function to also check that all quantities in StarState do
  in fact have units and will cause the function to only return units for quantities that are included, so for instance if
  StarState was only calculated by RotationQuantities() and not by ExtendedQuantities(), then only units for the rotation
  quantities will be included.
  
  Parameters
  ----------
  StarState : dict , optional
      Set of quantities already calculated by RotationQuantities() and possibly ExtendedQuantities().
  
  Returns
  ----------
  StarStateUnits : dict
      Dictionary of strings holding units for quantities.
  
  """
  
  # Create empty dictionary
  StarStateUnits = {}
  
  # Get initial dictionary of units
  StarStateUnits['Mstar'] = 'Msun'
  StarStateUnits['dAge'] = 'Myr'
  StarStateUnits['Age'] = 'Myr'
  StarStateUnits['OmegaEnv'] = 'OmegaSun =2.67e-6 rad s^-1'
  StarStateUnits['OmegaCore'] = 'OmegaSun =2.67e-6 rad s^-1'  
  StarStateUnits['Prot'] = 'days'
  StarStateUnits['Rstar'] = 'Rsun'
  StarStateUnits['tauConv'] = 'days'
  StarStateUnits['Ro'] = ''
  StarStateUnits['Itotal'] = 'g cm^2'
  StarStateUnits['Ienv'] = 'g cm^2'
  StarStateUnits['Icore'] = 'g cm^2'  
  StarStateUnits['dItotaldt'] = 'g cm^2 s^-1'
  StarStateUnits['dIenvdt'] = 'g cm^2 s^-1'
  StarStateUnits['dIcoredt'] = 'g cm^2 s^-1'
  StarStateUnits['Rcore'] = 'Rsun'
  StarStateUnits['dMcoredt'] = 'g s^-1'
  StarStateUnits['Mdot'] = 'g s^-1'
  StarStateUnits['Bdip'] = 'G'
  StarStateUnits['vEsc'] = 'cm s^-1'
  StarStateUnits['torqueEnvMom'] = 'erg'
  StarStateUnits['torqueCoreMom'] = 'erg'
  StarStateUnits['torqueEnvCG'] = 'erg'
  StarStateUnits['torqueCoreCG'] = 'erg'
  StarStateUnits['torqueEnvWind'] = 'erg'
  StarStateUnits['torqueEnvCE'] = 'erg'
  StarStateUnits['torqueCoreCE'] = 'erg'
  StarStateUnits['torqueEnvDL'] = 'erg'
  StarStateUnits['torqueEnv'] = 'erg'
  StarStateUnits['dOmegaEnvdt'] = 'OmegaSun Myr^-1'
  StarStateUnits['torqueCore'] = 'erg'
  StarStateUnits['dOmegaCoredt'] = 'OmegaSun Myr^-1'
  StarStateUnits['Lbol'] = 'erg s^-1'
  StarStateUnits['Teff'] = 'K'
  StarStateUnits['Lx'] = 'erg s^-1'
  StarStateUnits['Fx'] = 'erg s^-1 cm^-2'
  StarStateUnits['Rx'] = ''
  StarStateUnits['Tcor'] = 'MK'
  StarStateUnits['Leuv1'] = 'erg s^-1'
  StarStateUnits['Feuv1'] = 'erg s^-1 cm^-2'
  StarStateUnits['Reuv1'] = ''
  StarStateUnits['Leuv2'] = 'erg s^-1'
  StarStateUnits['Feuv2'] = 'erg s^-1 cm^-2'
  StarStateUnits['Reuv2'] = ''
  StarStateUnits['Leuv'] = 'erg s^-1'
  StarStateUnits['Feuv'] = 'erg s^-1 cm^-2'
  StarStateUnits['Reuv'] = ''  
  StarStateUnits['Lly'] = 'erg s^-1'
  StarStateUnits['Fly'] = 'erg s^-1 cm^-2'
  StarStateUnits['Rly'] = ''
  
  # If StarState was set in call to function then only include units for quantities in StarState and
  # make sure all quantities in StarState have units given (even if the quantity is dimensionless)
  if not StarState is None:
    
    # Make deep copy of StarStateUnits
    StarStateUnitsTemp = copy.deepcopy(StarStateUnits)
    
    # Start new StarStateUnits
    StarStateUnits = {}
    
    # Loop over quantities in StarState and add each one to StarStateUnits
    for quantity in StarState:
      
      # Make sure it is included in list of units above
      if not ( quantity in StarStateUnitsTemp ):
        misc._PrintErrorKill("units for "+quantity+" not included in units list")
      
      # Add to dictionary
      StarStateUnits[quantity] = StarStateUnitsTemp[quantity]
  
  return StarStateUnits

#====================================================================================================================

def Lxuv(Mstar=None,Age=None,Omega=None,OmegaEnv=None,Prot=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns X-ray, EUV, and Ly-alpha luminosities in erg s^-1.
  
  The user can supply a mass, age, and surface rotation rate for a star and it returns the XUV luminosities as 
  a dictionary. The rotation rate can be either the rotation velocity as a multiple of the solar rotation 
  rate (2.67e-6 rad s^-1) using either Omega or OmegaEnv keyword arguments, or as a rotation period in days using 
  the Prot keyword argument, The user can optionally also give a parameter dictionary and an instance of the 
  StarEvo class, though this is not necessary and if these are not specified then the defaults will be used.
  
  Parameters
  ----------
  Band : str
      Wavelength band to calculate.
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  Omega : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaEnv : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  Prot : float , optional
      Surface rotation period in days.
  params : dict , optional
      Dictionary holding model parameters.
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  LxuvDict : dict
      Dictionary of luminositie in erg s^-1.
  
  """
  
  # If an instance of the StarEvo class is not input, load it with defaults
  if StarEvo is None:
    StarEvo = SE.StarEvo()
  
  # Make sure mass and age are specified
  if Mstar is None:
    misc._PrintErrorKill("Mstar must be set in call to function")
  if Age is None:
    misc._PrintErrorKill("Age must be set in call to function")
  
  # Make sure rotation is set
  if ( Omega is None ) and ( OmegaEnv is None ) and ( Prot is None ):
    misc._PrintErrorKill("Omega, OmegaEnv, or Prot must be set in call to function")
  
  # Make sure only one rotation rate is set (add up number of set parameters)
  nSet = 0
  if not Omega is None:
    nSet += 1
  if not OmegaEnv is None:
    nSet += 1
  if not Prot is None:
    nSet += 1
  if not ( nSet == 1 ):
    misc._PrintErrorKill("can only set one of Omega, OmegaEnv, and Prot in call to function")
  
  # Get Omega not set
  if not OmegaEnv is None:
    Omega = OmegaEnv
  if not Prot is None:
    Omega = _Omega(Prot)
  
  # The function _Xray needs Ro, Rstar, and Lbol so make a StarState dictionary with these
  StarState = {}
  StarState['OmegaEnv'] = Omega
  StarState['Rstar'] = StarEvo.Rstar(Mstar,Age)
  StarState['Prot'] = _Prot(Omega)
  StarState['tauConv'] = StarEvo.tauConv(Mstar,Age)
  StarState['Ro'] = _Ro(StarState)
  StarState['Lbol'] = StarEvo.Lbol(Mstar,Age) * const.LbolSun
  
  # Get X-ray
  StarState['Lx'] , StarState['Fx'] , StarState['Rx'] = _Xray(StarState,params=params)
  
  # Get EUV
  StarState['Leuv1'] , StarState['Feuv1'] , StarState['Reuv1'] = _EUV1(StarState,params=params)
  StarState['Leuv2'] , StarState['Feuv2'] , StarState['Reuv2'] = _EUV2(StarState,params=params)
  StarState['Leuv'] , StarState['Feuv'] , StarState['Reuv'] = _EUV(StarState,params=params)
  
  # Get Ly-alpha 
  StarState['Lly'] , StarState['Fly'] , StarState['Rly'] = _Lymanalpha(StarState,params=params)
  
  # Make dictionary with luminosities
  LxuvDict = {}
  LxuvDict['Lxuv'] = StarState['Lx'] + StarState['Leuv']
  LxuvDict['Lx'] = StarState['Lx']
  LxuvDict['Leuv1'] = StarState['Leuv1']
  LxuvDict['Leuv2'] = StarState['Leuv2']
  LxuvDict['Leuv'] = StarState['Leuv']
  LxuvDict['Lly'] = StarState['Lly']
  
  return LxuvDict

#====================================================================================================================

def Lx(Mstar=None,Age=None,Omega=None,OmegaEnv=None,Prot=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns X-ray luminosity in erg s^-1.
  
  The user can supply a mass, age, and surface rotation rate for a star and it returns the X-ray luminosity. The
  rotation rate can be either the rotation velocity as a multiple of the solar rotation rate (2.67e-6 rad s^-1) 
  using either Omega or OmegaEnv keyword arguments, or as a rotation period in days using the Prot keyword 
  argument, The user can optionally also give a parameter dictionary and an instance of the StarEvo class, though 
  this is not necessary and if these are not specified then the defaults will be used.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  Omega : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaEnv : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  Prot : float , optional
      Surface rotation period in days.
  params : dict , optional
      Dictionary holding model parameters.
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  Lx : float
      X-ray luminosity in erg s^-1.
  
  """
  
  # Get dictionary with all luminosities
  LxuvDict = Lxuv(Mstar=Mstar,Age=Age,Omega=Omega,OmegaEnv=OmegaEnv,Prot=Prot,params=params,StarEvo=StarEvo)
  
  return LxuvDict['Lx']

#====================================================================================================================

def _Xray(StarState,params=params.paramsDefault):
  
  """Takes star state, returns Lx in erg s^-1, Fx in erg s^-1 cm^-1, and Rx."""
  
  # Get Rx depending on Ro and regime
  if ( StarState['Ro'] < params['RoSatXray'] ):
    Rx = params['C1Xray'] * StarState['Ro']**params['beta1Xray']
  else:
    Rx = params['C2Xray'] * StarState['Ro']**params['beta2Xray']
  
  # Get Lx from that
  Lx = StarState['Lbol'] * Rx
  
  # Get Fx
  Fx = Lx / ( 4.0 * const.Pi * (StarState['Rstar']*const.Rsun)**2.0 )
  
  return Lx , Fx , Rx

#====================================================================================================================

def _Tcor(StarState,params=params.paramsDefault):
  
  """Takes star state, returns average coronal temperature in MK."""
  
  # Get Tcor depending on Fx using equations from Johnstone & Guedel (2015)
  Tcor = 0.11 * StarState['Fx']**0.26
  
  return Tcor

#====================================================================================================================

def Leuv(Mstar=None,Age=None,Omega=None,OmegaEnv=None,Prot=None,band=0,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns EUV luminosity in erg s^-1.
  
  The user can supply a mass, age, and surface rotation rate for a star and it returns the EUV luminosity. 
  By default, this will be the luminosity in the 10-92 nm band, but using the band keyword argument, the 
  user can set specify either the 10-36 nm range with 'band=1' or the 36-92 nm range for 'band=2'. The default
  10-92 nm range is set with 'band=0'. The rotation rate can be either the rotation velocity as a multiple 
  of the solar rotation rate (2.67e-6 rad s^-1) using either Omega or OmegaEnv keyword arguments, or as a 
  rotation period in days using the Prot keyword argument, The user can optionally also give a parameter 
  dictionary and an instance of the StarEvo class, though this is not necessary and if these are not specified 
  then the defaults will be used.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  Omega : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaEnv : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  Prot : float , optional
      Surface rotation period in days.
  band : float , optional
      Which wavelength band to return (=0 for 10-92 nm; =1 for 10-32 nm; =2 for 32-92 nm).
  params : dict , optional
      Dictionary holding model parameters.
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  Leuv : float
      EUV luminosity in erg s^-1.
  
  """
  
  # Get dictionary with all luminosities
  LxuvDict = Lxuv(Mstar=Mstar,Age=Age,Omega=Omega,OmegaEnv=OmegaEnv,Prot=Prot,params=params,StarEvo=StarEvo)
  
  # Return depending on the desired band
  if ( band == 0 ):
    return LxuvDict['Leuv']
  elif ( band == 1 ):
    return LxuvDict['Leuv1']
  elif ( band == 2 ):
    return LxuvDict['Leuv2']
  else:
    misc._PrintErrorKill("invalid band set in call to function")
    
  return -1

#====================================================================================================================

def _EUV1(StarState,params=params.paramsDefault):
  
  """Takes star state, returns Leuv1 in erg s^-1, Feuv in erg s^-1 cm^-1, and Reuv (10-36 nm)."""
  
  # Get Feuv from Fx using equation from Johnstone et al. (2020)
  Feuv1 = 10.0**( 2.04 + 0.681 * np.log10(StarState['Fx']) )
  
  # Get Leuv1 from this
  Leuv1 = Feuv1 * ( 4.0 * const.Pi * (StarState['Rstar']*const.Rsun)**2.0 )
  
  # Get Reuv1
  Reuv1 = Leuv1 / StarState['Lbol'] 
  
  return Leuv1 , Feuv1 , Reuv1

#====================================================================================================================

def _EUV2(StarState,params=params.paramsDefault):
  
  """Takes star state, returns Leuv2 in erg s^-1, Feuv2 in erg s^-1 cm^-1, and Reuv2 (36-92 nm)."""
  
  # Get Feuv2 from Feuv1 using equation from Johnstone et al. (2020)
  Feuv2 = 10.0**( -0.341 + 0.92 * np.log10(StarState['Feuv1']) )
  
  # Get Leuv1 from this
  Leuv2 = Feuv2 * ( 4.0 * const.Pi * (StarState['Rstar']*const.Rsun)**2.0 )
  
  # Get Reuv2
  Reuv2 = Leuv2 / StarState['Lbol']
  
  return Leuv2 , Feuv2 , Reuv2

#====================================================================================================================

def _EUV(StarState,params=params.paramsDefault):
  
  """Takes star state, returns Leuv in erg s^-1, Feuv in erg s^-1 cm^-1, and Reuv (10-92 nm)."""
  
  # Just add quantities 
  Leuv = StarState['Leuv1'] + StarState['Leuv2']
  Feuv = StarState['Feuv1'] + StarState['Feuv2']
  Reuv = StarState['Reuv1'] + StarState['Reuv2']
  
  return Leuv , Feuv , Reuv

#====================================================================================================================

def Lly(Mstar=None,Age=None,Omega=None,OmegaEnv=None,Prot=None,params=params.paramsDefault,StarEvo=None):
  
  """
  Takes basic stellar parameters, returns Ly-alpha luminosity in erg s^-1.
  
  The user can supply a mass, age, and surface rotation rate for a star and it returns the Ly-alpha luminosity. The
  rotation rate can be either the rotation velocity as a multiple of the solar rotation rate (2.67e-6 rad s^-1) 
  using either Omega or OmegaEnv keyword arguments, or as a rotation period in days using the Prot keyword 
  argument, The user can optionally also give a parameter dictionary and an instance of the StarEvo class, though 
  this is not necessary and if these are not specified then the defaults will be used.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float
      Age of star.
  Omega : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaEnv : float , optional
      Surface rotation rate in OmegaSun (=2.67e-6 rad/s).
  Prot : float , optional
      Surface rotation period in days.
  params : dict , optional
      Dictionary holding model parameters.
  StarEvo : Mors.stellarevo.StarEvo , optional
      Instance of StarEvo class holding stellar evolution model data.
  
  Returns
  ----------
  Lly : float
      Lyman-alpha luminosity in erg s^-1.
  
  """
  
  # Get dictionary with all luminosities
  LxuvDict = Lxuv(Mstar=Mstar,Age=Age,Omega=Omega,OmegaEnv=OmegaEnv,Prot=Prot,params=params,StarEvo=StarEvo)
  
  return LxuvDict['Lly']

#====================================================================================================================

def _Lymanalpha(StarState,params=params.paramsDefault):
  
  """Takes star state, returns Lly in erg s^-1, Fly in erg s^-1 cm^-1, and Rly (10-92 nm)."""
  
  # Get Fly using from Fx using equation from Johnstone et al. (2020)
  Fly = 10.0**( 3.97 + 0.375 * np.log10(StarState['Fx']) )
  
  # Get Lly
  Lly = Fly * ( 4.0 * const.Pi * (StarState['Rstar']*const.Rsun)**2.0 )
  
  # Get Rly
  Rly = Lly / StarState['Lbol']
  
  return Lly , Fly , Rly

#====================================================================================================================

def _torqueDiskLocking(StarState,CoreEnvelopeDecoupling,params=params.paramsDefault):
  
  """Takes stellar parameters, returns disk-locking torque."""
  
  # As it is currently setup, this will only work if the other torques are calculated first
  
  # Get disk locking age
  ageDL = params['aDiskLock'] * StarState['OmegaEnv']**params['bDiskLock']
  ageDL = min( ageDL , params['ageDLmax'] )
  
  # Return nothing if after disk locking age
  if ( StarState['Age'] > ageDL ):
    return 0.0
  
  # Get torque on envelope
  if ( CoreEnvelopeDecoupling ):
    torqueEnvDL = - ( StarState['torqueEnvWind'] + StarState['torqueEnvCE'] + StarState['torqueEnvCG'] + StarState['torqueEnvMom'] )
  else:
    torqueEnvDL = - ( StarState['torqueEnvWind'] + StarState['torqueEnvMom'] )
  
  return torqueEnvDL

#====================================================================================================================

def _torqueCoreEnvelope(StarState,params=params.paramsDefault):
  
  """Takes stellar parameters, returns core-envelope torque."""
    
  # Get deltaJ
  deltaJ = StarState['Icore'] * StarState['Ienv'] / ( StarState['Icore'] + StarState['Ienv'] ) * ( StarState['OmegaCore'] - StarState['OmegaEnv'] ) * const.OmegaSun
  
  # Get core-envelope coupling timescale in Myr
  timeCE = params['aCoreEnvelope'] * StarState['OmegaEnv']**params['bCoreEnvelope'] * StarState['Mstar']**params['cCoreEnvelope']
  timeCE = max( timeCE , params['timeCEmin'] )
  
  # Convert timescale to seconds
  timeCE *= const.Myr
  
  # Get torques
  torqueEnvMoment = deltaJ / timeCE
  torqueCoreMoment = - torqueEnvMoment
  
  return ( torqueEnvMoment , torqueCoreMoment )

#====================================================================================================================

def _torqueMoment(StarState,CoreEnvelopeDecoupling):
  
  """Takes stellar parameters, returns moment of inertia change torque in erg."""
  
  # Note that the torque calculated in this function is not actually a real physical torque
  # but it is convenient to write it in the equations and in the code as if it was. 
  
  # Get torque in erg depending on if core-envelope decoupling is being done
  if ( CoreEnvelopeDecoupling ):
    torqueEnvMoment = - StarState['dIenvdt'] * StarState['OmegaEnv']*const.OmegaSun
    torqueCoreMoment = - StarState['dIcoredt'] * StarState['OmegaCore']*const.OmegaSun
  else:
    torqueEnvMoment = - StarState['dItotaldt'] * StarState['OmegaEnv']*const.OmegaSun
    torqueCoreMoment = 0.0
    
  return ( torqueEnvMoment , torqueCoreMoment )

#====================================================================================================================

def _torqueCoreGrowth(StarState):
  
  """Takes stellar parameters, returns core-growth torque in erg."""
  
  # Get torque in erg
  torqueEnvCG = - 2.0/3.0 * StarState['Rcore']**2.0 * StarState['OmegaEnv']*const.OmegaSun * StarState['dMcoredt']
  torqueCoreCG = - torqueEnvCG
  
  return ( torqueEnvCG , torqueCoreCG )

#====================================================================================================================

def _torqueWind(StarState,params=params.paramsDefault):
  
  # Some constants from Matt et al. (2012)
  K1Matt = 1.3
  K2Matt = 0.0506
  mMatt = 0.2177
  
  # Calculate torque
  torqueWind = K1Matt**2.0 * StarState['Bdip']**(4.0*mMatt) * StarState['Mdot']**(1.0-2.0*mMatt) * (StarState['Rstar']*const.Rsun)**(4.0*mMatt+2.0) * StarState['OmegaEnv'] * const.OmegaSun / ( K2Matt**2.0 * StarState['vEsc']**2.0 + (StarState['OmegaEnv']*const.OmegaSun)**2.0 * (StarState['Rstar']*const.Rsun)**2.0 )**mMatt
  
  # Multiply by extra factor and make negitive so it is a spin-down torque
  torqueWind = - params['Kwind'] * torqueWind
  
  return torqueWind

#====================================================================================================================

def _Omega(Prot):
  """Takes rotation period in days, returns angular velocity in OmegaSun."""
  return 2.0*const.Pi / ( Prot*const.day ) / const.OmegaSun

#====================================================================================================================

def _Prot(Omega):
  """Takes angular velocity in OmegaSun, returns rotation period in days."""
  return 2.0*const.Pi / ( Omega*const.OmegaSun ) / const.day

#====================================================================================================================

def _Ro(StarState):
  
  """Takes state of star, returns Rossby number."""
  
  # Get Rossby 
  Ro = StarState['Prot'] / StarState['tauConv']
  
  return Ro

#====================================================================================================================

def _Bdip(StarState,params=params.paramsDefault):
  
  """Takes state of star, returns dipole field strength."""
  
  # Use Rossby number to get dipole field, depending if saturated or not
  if ( StarState['Ro'] >= params['RoSatBdip'] ):
    BdipStar = params['BdipSun'] * ( StarState['Ro'] / const.RoSun )**params['aBdip']
  else:
    BdipStar = params['BdipSun'] * ( params['RoSatBdip'] / const.RoSun )**params['aBdip']
  
  return BdipStar

#====================================================================================================================

def _Mdot(StarState,params=params.paramsDefault):
  
  # Use Rossby number to get mass loss rate
  if ( StarState['Ro'] >= params['RoSatMdot'] ):
    MdotWind = params['MdotSun'] * StarState['Rstar']**2.0 * (StarState['Ro']/const.RoSun)**params['aMdot'] * StarState['Mstar']**params['bMdot']
  else:
    MdotWind = params['MdotSun'] * StarState['Rstar']**2.0 * (params['RoSatMdot']/const.RoSun)**params['aMdot'] * StarState['Mstar']**params['bMdot']
  
  # Multiply Mdot by additional mutliplicative factor for rapid rotators if desired
  if params['BreakupMdotIncrease']:
    MdotWind *= _MdotFactor(StarState['Mstar'],StarState['Rstar'],StarState['OmegaEnv'])
  
  return MdotWind

#====================================================================================================================

def _MdotFactor(Mstar,Rstar,OmegaSurface,params=params.paramsDefault):
  
  """Takes stellar mass, radius, and rotation rate, returns factor to increase wind mass loss rate."""
  
  # Initially assume factor is unity
  factor = 1.0
  
  # Rotation fraction of breakup
  fracBreak = OmegaSurface / OmegaBreak(Mstar,Rstar)
  
  # Work out if close to break up 
  if ( fracBreak > params['fracBreakThreshold'] ):
    
    # The equation to use should depend on the fraction of break-up that the star is at
    # basically I (C.P.Johnstone) did two fits and the better fit is the first one here
    # but this is unstable when fracBreak>1.009 so at that point, the code switches to 
    # the other fit. This should not actually happen since stars should not get that 
    # rapidly rotating. 
    if ( fracBreak < 1.009 ):
    
      # For the relation below, do not allow frac larger than 1.009 since that
      # gives us NaNs 
      fracBreak = min( fracBreak , 1.009 )
      
      # Get the factor using a better fit to the factor-frac relation based on Johnstone (2017)
      factor = 0.93 * ( 1.01 - fracBreak )**(-0.43) * np.exp( 0.31 * fracBreak**7.5 )
      
    else:
      
      #Get the multiplicative factor (ORIGINAL NOT SO GOOD EQUATION)
      factor = (5.94401988e+03) * fracBreak**7.0 + (-2.12411847e+04) * fracBreak**6.0 + (3.07933216e+04) * fracBreak**5.0 + (-2.32693067e+04) * fracBreak**4.0 + (9.80417648e+03) * fracBreak**3.0 + (-2.27923380e+03) * fracBreak**2.0 + (2.68629816e+02) * fracBreak**1.0 + (-1.13054867e+01) * fracBreak**0.0
      
  return factor

#====================================================================================================================

def OmegaBreak(Mstar,Rstar):
  
  """
  Takes stellar mass and radius in solar units, returns breakup rotation rate in OmegaSun.
    
  Parameters
  ----------
  Mstar : float 
      Mass of star in Msun.
  Rstar : float 
      Radius of star in Rsun.
  
  Returns
  ----------
  OmegaBreak : float
      Breakup angular velocity in OmegaSun.
  
  """
  
  # Assuming here that Rpole = Rstar  
  return (2.0/3.0)**(3.0/2.0) * ( const.GravConstant * Mstar*const.Msun )**0.5 / ( Rstar*const.Rsun )**(3.0/2.0) / const.OmegaSun

#====================================================================================================================

def _vEsc(Mstar,Rstar):
  
  """
  Takes stellar mass and radius in solar units, returns surface escape velocity in cm s^-1.
    
  Parameters
  ----------
  Mstar : float 
      Mass of star in Msun.
  Rstar : float 
      Radius of star in Rsun.
  
  Returns
  ----------
  vEsc : float
      Surface escape velocity in cm s^-1.
  
  """
  
  return ( 2.0 * const.GravConstant * Mstar*const.Msun / (Rstar*const.Rsun) )**0.5

#====================================================================================================================

def _shouldCoreEnvelopeDecoupling(StarState,params=params.paramsDefault):
  
  """Takes state of star, determines if core-envelope decoupling should be used."""
  
  # For this, there are several things that need to be fulfilled:-
  #     1. The parameter CoreEnvelopeDecoupling must be set to True
  #     2. The ratio of the core to total stellar moment of inertia must be above the parameter IcoreThresholdCE
  #     3. The stellar mass must be above the parameter MstarThresholdCE
  
  # Test condition 1
  if not ( params['CoreEnvelopeDecoupling'] ):
    return False
  
  # Test condition 2
  if ( StarState['Icore']/StarState['Itotal'] < params['IcoreThresholdCE'] ):
    return False
  
  # Test condition 3
  if ( StarState['Mstar'] < params['MstarThresholdCE'] ):
    return False
  
  # If got to here, then core-envelope decoupling should happen
  return True

#====================================================================================================================



























#====================================================================================================================

#def Mdot(RoStar=None,Mstar=None,Rstar=None,Omega=None,Prot=None,Age=AgeDefault,params=params.paramsDefault):
  
  #"""
  #Takes stellar rotation related parameters, returns wind mass loss rate.
  
  #The user can use this function to get the mass loss rate given either the star's Rossby number or its
  #mass, rotation rate, and age. If RoStar is given, the other stellar parameters will not be considered. If RoStar is
  #not given, the mass and either Omega or Prot must be set. The age is set to 5000 Myr by default but can also
  #be specified. For example, to get the Sun's Mdot from our model, the user can use
  
  #>>> mors.Mdot(Mstar=1.0,Omega=1.0,age=4567.0)
  
  #In addition, all of the parameters of the model, such as the dipole field strength of the Sun, the slope of 
  #the Bdip--Ro relation in the unsaturated regime, and the saturation Ro value can be input, though generally
  #this is not recommended. If the radius is not set, it will be calculated based on the mass and age.
  
  #Parameters
  #----------
  #RoStar : float , optional
      #Rossby number.
  #Mstar : float , optional
      #Mass of star in Msun.
  #Rstar : float , optional
      #Radius of star in Rsun.
  #Omega : float , optional
      #Rotation rate in OmegaSun (=2.67e-6 rad/s).
  #Prot : float , optional
      #Rotation period in days (cannot input both Prot and Omega).
  #Age : float , optional
      #Age of star in Myr (default is 5000 Myr).
  
  #Returns
  #----------
  #MdotWind : float
      #Stellar wind mass loss rate in Msun yr^-1.
  
  #"""
  
  #print("Doing nothing yet.")
  
  #return





#====================================================================================================================

#def Bdip(RoStar=None,Mstar=None,Omega=None,Prot=None,Age=AgeDefault,BdipSun=BdipSun,RoSun=RoSun,aBdip=aBdip,RoSatBdip=RoSatBdip):
  
  #"""
  #Takes stellar rotation related parameters, returns empirical dipole field strength.
  
  #The user can use this function to get the dipole field strength given either its Rossby number or given its
  #mass, rotation rate, and age. If RoStar is given, the other stellar parameters will not be considered. If RoStar is
  #not given, the mass and either Omega or Prot must be set. The age is set to 5000 Myr by default but can also
  #be specified. For example, to get the Sun's Bdip from our model, the user can use
  
  #>>> mors.Bdip(Mstar=1.0,Omega=1.0,age=4567.0)
  
  #In addition, all of the parameters of the model, such as the dipole field strength of the Sun, the slope of 
  #the Bdip--Ro relation in the unsaturated regime, and the saturation Ro value can be input, though generally
  #this is not recommended.
  
  #Parameters
  #----------
  #RoStar : float , optional
      #Rossby number.
  #Mstar : float , optional
      #Mass of star in Msun.
  #Omega : float , optional
      #Rotation rate in OmegaSun (=2.67e-6 rad/s).
  #Prot : float , optional
      #Rotation period in days (cannot input both Prot and Omega).
  #Age : float , optional
      #Age of star in Myr (default is 5000 Myr).
  #BdipSun : float , optional
      #Solar dipole field strength in G. 
  #RoSun : float , optional
      #Solar Rossby number.
  #aBdip : float , optional
      #Power law index a in Bdip ~ Ro^a for the unsaturated regime.
  #RoSatBdip : float , optional
      #Saturation Rossby number for dipole field strength.
  
  #Returns
  #----------
  #BdipStar : float
      #Dipole field strength in G.
  
  #"""
  
  ## Get Rossby number if it was not set as argument
  #if RoStar is None:
    #RoStar = Ro(Mstar=Mstar,Omega=Omega,Prot=Prot,Age=Age)
    
  ## Use Rossby number to get dipole field, depending if saturated or not
  #if ( RoStar >= RoSatBdip ):
    #BdipStar = BdipSun * ( RoStar / RoSun )**aBdip
  #else:
    #BdipStar = BdipSat
  
  #return BdipStar





#====================================================================================================================

#def Ro(Mstar=None,Omega=None,Prot=None,Age=AgeDefault):
  
  #"""
  #Takes stellar properties, returns Rossby number.
  
  #The user can use this to get the Rossby number for a star with a given mass, rotation rate, and 
  #age. The mass and either Omega or Prot must be set. The age is set to 5000 Myr by default. For example
  #the user can use
  
  #>>> mors.Ro(Mstar=1.0,Omega=1.0,Age=4567.0)
  
  #to get the Rossby number of the Sun as calculated assuming our assumption for the Sun's rotation rate
  #and the stellar evolution models we use.
  
  #Parameters
  #----------
  #Mstar : float
      #Mass of star in Msun.
  #Omega : float , optional
      #Rotation rate in OmegaSun (=2.67e-6 rad/s).
  #Prot : float , optional
      #Rotation period in days (cannot input both Prot and Omega).
  #Age : float , optional
      #Age of star in Myr (default is 5000 Myr).
  
  #Returns
  #----------
  #Ro : float
      #Rossby number.
  
  #"""
  
  ## Make sure Mstar is set
  #if Mstar is None:
    #misc._PrintErrorKill("if Ro not given as argument, must input stellar mass")
  
  ## Make sure either Omega or Prot (and not both) are set
  #if ( ( Omega is None ) and ( Prot is None ) ):
    #misc._PrintErrorKill("if Ro not given as argument, must input either Omega or Prot")
  #if ( ( not Omega is None ) and ( not Prot is None ) ):
    #misc._PrintErrorKill("cannot give both Omega and Prot as arguments")
  
  ## Get Prot if not set, get it from Omega
  #if Prot is None: 
    #Prot = 2.0 * const.Pi / ( Omega*const.OmegaSun ) / const.day
  
  ## Get convective turnover time
  #tauConv = SE.tauConv(Mstar,Age)
  
  ## Get Rossby number
  #Ro = Prot / tauConv
  
  #return Ro


