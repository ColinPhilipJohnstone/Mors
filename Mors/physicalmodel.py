
"""Module for holding functions for calculating physical quantities for the physical rotation and activity model."""

# Imports for standard stuff needed here
import numpy as np
import sys as sys

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.stellarevo as SE
import Mors.constants as const
import Mors.parameters as params

#====================================================================================================================

def dOmegadt(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,params=params.paramsDefault):
  
  """
  Takes basic stellar parameters, returns rates of change of rotation.
  
  This is the main function to be used to get the rates of change of the core and envelope rotation 
  rates given basic stellar parameters. The actual calculations of these are done in GetState() which
  itself also returns a very large number of other things used during the calculation of the rates 
  of change.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Age : float , optional
      Age of star.
  OmegaEnv : float 
      Envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
  OmegaCore : float 
      Core rotation rate in OmegaSun (=2.67e-6 rad/s).
  
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
  
  # Get state of system 
  StarState = GetState(Mstar=Mstar,Age=Age,OmegaEnv=OmegaEnv,OmegaCore=OmegaCore,params=params)
  
  return ( StarState['dOmegaEnvdt'] , StarState['dOmegaCoredt'] )

#====================================================================================================================

def GetState(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,params=params.paramsDefault):
  
  #--------------------------------------------------------
  
  # Setup stellar parameter dictionary
  # This will be used to hold all of the 
  
  # Start dictionary
  StarState = {}
  
  # Add basic input parameters 
  StarState['Mstar'] = Mstar
  StarState['Age'] = Age
  StarState['OmegaEnv'] = OmegaEnv
  StarState['OmegaCore'] = OmegaCore
  
  # Radius in cm
  StarState['Rstar'] = SE.Rstar(Mstar,Age)
  
  # Convective turnover time in days
  StarState['tauConv'] = SE.tauConv(Mstar,Age)
  
  # Rossby number
  StarState['Ro'] = _Ro(StarState)
  
  # Get moments of interia in g cm^2
  StarState['Itotal'] = SE.Itotal(Mstar,Age)
  StarState['Ienv'] = SE.Ienv(Mstar,Age)
  StarState['Icore'] = SE.Icore(Mstar,Age)
  
  # Rates of change of moment of inertia in g cm^2 s^-1
  StarState['dItotaldt'] = SE.dItotaldt(Mstar,Age) / const.Myr
  StarState['dIenvdt'] = SE.dIenvdt(Mstar,Age) / const.Myr
  StarState['dIcoredt'] = SE.dIcoredt(Mstar,Age) / const.Myr
  
  # Rcore in cm
  StarState['Rcore'] = SE.Rcore(Mstar,Age) * const.Rsun
  
  # dMcoredt in g s^-1 Msun Myr^-1
  StarState['dMcoredt'] = SE.dMcoredt(Mstar,Age) * const.Msun / const.Myr
  
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

def _Ro(StarState):
  
  """Takes state of star, returns Rossby number."""
  
  # Get rotation period from surface rotation in days
  Prot = 2.0 * const.Pi / ( StarState['OmegaEnv']*const.OmegaSun ) / const.day
  
  # Get Rossby 
  Ro = Prot / StarState['tauConv']
  
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


