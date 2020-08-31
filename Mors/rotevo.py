
"""Module for holding functions used to evovle the rotation of stars."""

# Imports for standard stuff needed here
import numpy as np

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.physicalmodel as phys
import Mors.parameters as params

#====================================================================================================================

AgeMinDefault = 1.0       # when to start evolution (Myr)
AgeMaxDefault = 5000.0    # when to end evolution (Myr)

#====================================================================================================================

def EvolveRotation(Mstar=None,Omega0=None,OmegaEnv0=None,OmegaCore0=None,AgeMin=AgeMinDefault,AgeMax=AgeMaxDefault,params=params.paramsDefault):
  
  """
  Takes stellar mass and rotation rate, evolves rotation. 
  
  This is the main function for rotational evolution calculations. It will take some basic stellar parameters and
  evolve rotation between two ages.
  
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  Omega0 : float 
      Starting rotation rate in OmegaSun (=2.67e-6 rad/s).
  AgeMin : float , optional
      Starting age (default = 1 Myr).
  AgeMax : float , optional
      Final age (default = 5 Gyr).
  params : dict , optional
      Parameters to determine behavior of code (default given in parameters.py).
  
  Returns
  ----------
  Tracks : dict
      Dictionary with evolutionary tracks and related information. 
  
  """
  
  # Make sure stellar mass is set
  if Mstar is None: 
    misc._PrintErrorKill("argument Mstar must be set in call to function")
  
  # Make sure initial rotation is set correctly
  if Omega0 is None: 
    if ( ( OmegaEnv0 is None ) or ( OmegaCore0 is None ) ):
      misc._PrintErrorKill("if Omega0 is not set then both OmegaEnv0 and OmegaCore0 must be set in call to function")
  else:
    if not ( ( OmegaEnv0 is None ) and ( OmegaCore0 is None ) ):
      misc._PrintErrorKill("if Omega0 is set then neither OmegaEnv0 nor OmegaCore0 can be set in call to function")
  
  # Starting values for time
  iAge = 0
  Age = AgeMin
  
  # Starting rotation rates
  if Omega0 is None:
    OmegaEnv = OmegaEnv0
    OmegaCore = OmegaCore0
  else:
    OmegaEnv = Omega0
    OmegaCore = Omega0
  
  # Make dictionary for holding tracks 
  Tracks = { 
    'Age' : np.array([]) ,
    'dAge' : np.array([]) ,
    'OmegaEnv' : np.array([]) ,
    'OmegaCore' : np.array([]) ,
    }
  
  # Add initial values to tracks
  Tracks['Age'] = np.append( Tracks['Age'] , Age )
  Tracks['dAge'] = np.append( Tracks['dAge'] , 0.0 )
  Tracks['OmegaEnv'] = np.append( Tracks['OmegaEnv'] , OmegaEnv )
  Tracks['OmegaCore'] = np.append( Tracks['OmegaCore'] , OmegaCore )
  
  # Get starting timestep in Myr
  dAge = 0.5
  
  # Start looping, end only when end time is reached (additional 0.999 factor to stop code getting stuck)
  while ( Age < 0.999*AgeMax ):
    
    print(Age,dAge,OmegaEnv,OmegaCore)
    
    # Get maximum timestep in Myr
    dAgeMax = AgeMax - Age
    
    # Make sure timestep is not too long
    dAge = min( dAge , dAgeMax )
    
    # Do timestep
    dAge , dAgeNew , OmegaEnv , OmegaCore = EvolveRotationStep( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , dAgeMax=dAgeMax , dAge=dAge , params=params )
    
    # Update time
    iAge += 1
    Age += dAge
    
    # Get new age step
    dAge = dAgeNew
    
    # Add time to tracks
    Tracks['Age'] = np.append( Tracks['Age'] , Age )
    Tracks['OmegaEnv'] = np.append( Tracks['OmegaEnv'] , OmegaEnv )
    Tracks['OmegaCore'] = np.append( Tracks['OmegaCore'] , OmegaCore )
    
    # Determine if too many steps taken
    if ( iAge >= params['nStepMax'] ):
      misc._PrintErrorKill("too many timesteps taken")
      
  print(Age,dAge,OmegaEnv,OmegaCore)
  
  return Tracks

#====================================================================================================================

def EvolveRotationStep(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAgeMax=None,dAge=None,params=params.paramsDefault):
  
  """
  Evolves rotation by a single step.
  
  This is the main function to be used to do timestep updates of the core and envelope rotation rates
  for a star given its basic parameters. The five arguments Mstar, Age, OmegaEnv, OmegaCore, 
  and either dAgeMax or dAge must be set in the call to this function. For certain solvers, the timestep 
  is determined automatically and this value is returned along with the updated envelope and
  core rotation rates
  
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
  dAgeMax : float , optional
      Maximum timestep to take in Myr.
  dAgeMax : float , optional
      Timestep to take in Myr.
  
  Returns
  ----------
  dAge : numpy.ndarray
      Age step in Myr.
  OmegaEnv : numpy.ndarray
      Envelope rotation rates in OmegaSun (=2.67e-6 rad/s).
  OmegaCore : numpy.ndarray
      Core rotation rates in OmegaSun (=2.67e-6 rad/s).
  
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
  if ( ( dAgeMax is None ) and ( dAge is None ) ): 
    misc._PrintErrorKill("either dAge or dAgeMax argument must be set in call to function")
  
  # Do step based on method used
  if ( params['TimeIntegrationMethod'] == 'ForwardEuler' ):
    dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepFE( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , dAge=dAge , params=params )
  elif ( params['TimeIntegrationMethod'] == 'RungeKutta4' ):
    dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepRK4( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , dAge=dAge , params=params )
  elif ( params['TimeIntegrationMethod'] == 'RungeKuttaFehlberg' ):
    dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepRKF( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , dAge=dAge , dAgeMax=dAgeMax , params=params )
  else:
    misc._PrintErrorKill("invalid value of TimeIntegrationMethod in parameters")

  
  return ( dAge , dAgeNew , OmegaEnv , OmegaCore )

#====================================================================================================================

def _EvolveRotationStepFE(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAge=None,params=params.paramsDefault):
  
  """Takes basic stellar parameters, evolves by timestep using forward Euler method."""
  
  # Get rates of change
  dOmegaEnvdt , dOmegaCoredt = phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params )
  
  # Do update
  OmegaEnv += dAge * dOmegaEnvdt
  OmegaCore += dAge * dOmegaCoredt
  
  return ( dAge , dAge , OmegaEnv , OmegaCore )

#====================================================================================================================

def _EvolveRotationStepRK4(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAge=None,params=params.paramsDefault):
  
  """Takes basic stellar parameters, evolves by timestep using classical Runge-Kutta method."""
  
  # Get rates of change
  k1Env , k1Core = phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params )
  k2Env , k2Core = phys.dOmegadt( Mstar=Mstar , Age=Age+0.5*dAge , OmegaEnv=OmegaEnv+0.5*dAge*k1Env , OmegaCore=OmegaCore+0.5*dAge*k1Core , params=params )
  k3Env , k3Core = phys.dOmegadt( Mstar=Mstar , Age=Age+0.5*dAge , OmegaEnv=OmegaEnv+0.5*dAge*k2Env , OmegaCore=OmegaCore+0.5*dAge*k2Core , params=params )
  k4Env , k4Core = phys.dOmegadt( Mstar=Mstar , Age=Age+dAge , OmegaEnv=OmegaEnv+dAge*k3Env , OmegaCore=OmegaCore+dAge*k3Core , params=params )
    
  # Do update
  OmegaEnv += dAge * ( k1Env + 2.0*k2Env + 2.0*k3Env + k4Env ) / 6.0
  OmegaCore += dAge * ( k1Core + 2.0*k2Core + 2.0*k3Core + k4Core ) / 6.0
  
  return ( dAge , dAge , OmegaEnv , OmegaCore )

#====================================================================================================================

def _EvolveRotationStepRKF(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAgeMax=None,dAge=None,params=params.paramsDefault):

  """Takes basic stellar parameters, evolves by timestep using Runge-Kutta-Fehlberg method."""
  
  # Start with dAge as timestep
  dAgeNew = dAge
  
  # Setup array to hold integration values
  X = np.array([OmegaEnv,OmegaCore])
  
  # Start iterating until got accuracy desired
  while True:
    
    # Set dAge
    dAge = dAgeNew
    
    # k1
    k1 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=X[0] , OmegaCore=X[1] , params=params ) )
    k1Env , k1Core = phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params )
    
    # k2
    Age2 = Age + (1.0/5.0)*dAge
    X2 = X + (1.0/5.0)*dAge*k1
    k2 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age2 , OmegaEnv=X2[0] , OmegaCore=X2[1] , params=params ) )
    #k2Env , k2Core = phys.dOmegadt( Mstar=Mstar , Age=Age+(1.0/5.0)*dAge , OmegaEnv=OmegaEnv+(1.0/5.0)*dAge*k1Env , OmegaCore=OmegaCore+(1.0/5.0)*dAge*k2Core , params=params )
    #dXdt2(:) = get_dXdt( t+(1.0/5.0)*dt , X(:)+(1.0/5.0)*dt*dXdt1(:) )
    
    # k3
    Age3 = Age + (3.0/10.0)*dAge
    X3 = X + (3.0/40.0)*dAge*k1 + (9.0/40.0)*dAge*k2
    k3 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age3, OmegaEnv=X3[0] , OmegaCore=X3[1] , params=params ) )
    #dXdt3(:) = get_dXdt( t+(3.0/10.0)*dt , X(:)+(3.0/40.0)*dt*dXdt1(:)+(9.0/40.0)*dt*dXdt2(:) )
    
    # k4
    Age4 = Age + (3.0/5.0)*dAge
    X4 = X + (3.0/10.0)*dAge*k1 - (9.0/10.0)*dAge*k2 + (6.0/5.0)*dAge*k3
    k4 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age4, OmegaEnv=X4[0] , OmegaCore=X4[1] , params=params ) )
    #dXdt4(:) = get_dXdt( t+(3.0/5.0)*dt , X(:)+(3.0/10.0)*dt*dXdt1(:)-(9.0/10.0)*dt*dXdt2(:)+(6.0/5.0)*dt*dXdt3(:) )
    
    # k5
    Age5 = Age + dAge
    X5 = X - (11.0/54.0)*dAge*k1 + (5.0/2.0)*dAge*k2 - (70.0/27.0)*dAge*k3 + (35.0/27.0)*dAge*k4
    k5 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age5, OmegaEnv=X5[0] , OmegaCore=X5[1] , params=params ) )
    #dXdt5(:) = get_dXdt( t+dt , X(:)-(11.0/54.0)*dt*dXdt1(:)+(5.0/2.0)*dt*dXdt2(:)-(70.0/27.0)*dt*dXdt3(:)+(35.0/27.0)*dt*dXdt4(:) )
    
    # k6
    Age6 = Age + (7.0/8.0)*dAge
    X6 = X + (1631.0/55296.0)*dAge*k1 + (175.0/512.0)*dAge*k2 + (575.0/13824.0)*dAge*k3 + (44275.0/110592.0)*dAge*k4 + (253.0/4096.0)*dAge*k5
    k6 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age6, OmegaEnv=X6[0] , OmegaCore=X6[1] , params=params ) )
    #dXdt6(:) = get_dXdt( t+(7.0/8.0)*dt , X(:)+(1631.0/55296.0)*dt*dXdt1(:)+(175.0/512.0)*dt*dXdt2(:)+(575.0/13824.0)*dt*dXdt3(:)+(44275.0/110592.0)*dt*dXdt4(:)+(253.0/4096.0)*dt*dXdt5(:) )
    
    # Calculate candidate update
    Xnew = X + (37.0/378.0)*dAge*k1 + (250.0/621.0)*dAge*k3 + (125.0/594.0)*dAge*k4 + (512.0/1771.0)*dAge*k6 
    
    # Calculate forth order update
    Xforth = X + (2825.0/27648.0)*dAge*k1 + (18575.0/48384.0)*dAge*k3 + (13525.0/55296.0)*dAge*k4 + (277.0/14336.0)*dAge*k5 + (1.0/4.0)*dAge*k6
    
    # Get error in this timestep for each species and grid point
    Delta =  Xnew - Xforth
        
    # Get the factor by which to change dt
    # Need to take care here if Delta=0 since this can happen if dX/dt=0 over the entire timestep
    if ( np.min(np.abs(Delta)) == 0.0 ):
      dAgeFactor = 1.5
    else:
      dAgeFactor = np.min( np.abs( params['DeltaDesired']*X / Delta ) )
    
    # Get new dt
    dAgeNew = 0.9 * dAge * dAgeFactor**0.2
    
    # Make sure dt is not too big
    dAgeNew = min( dAgeNew , dAgeMax )
    
    # Check for stopping of iteration
    if ( dAgeFactor > 1.0 ):
      break
    
  # Save new estimate
  OmegaEnv = Xnew[0]
  OmegaCore = Xnew[1]
  
  return dAge , dAgeNew , OmegaEnv , OmegaCore
    
#====================================================================================================================

