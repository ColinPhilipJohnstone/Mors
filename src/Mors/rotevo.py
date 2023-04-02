
"""Module for holding functions used to evovle the rotation of stars."""

# Imports for standard stuff needed here
import numpy as np
import math
import copy

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.physicalmodel as phys
import Mors.parameters as params
import Mors.stellarevo as SE
import sys

AgeMinDefault = 1.0       # when to start evolution (Myr)
AgeMaxDefault = 5000.0    # when to end evolution (Myr)

def FitRotation(Mstar=None,Age=None,Omega=None,AgeMin=None,params=params.paramsDefault,StarEvo=None):
    """Takes stellar mass, age, and surface rotation rate, returns initial rotation rate of star. 
    
    This function can be used if the age and surface rotation rate of the star is known and the user wants to get
    the corresponding initial (usually 1 Myr) rotation rate for the star. The star's mass, age, and rotation rates 
    must be input. The starting age to get the initial rotation rate for can be set, and by default this is set to
    AgeMinFit in the parameters and it is not recommended to change this. In addition, model parameters and an instance
    of the StarEvo class can be input. The return value of this function should be the initial rotation rate but if 
    the input Omega is below the minimum that can be fit, it will be -1, if it is above the maximum that can be 
    fit, it will be -2, and if it is in the range but no solution could be found it will be -3.
    
    Parameters
    ----------
    Mstar : float
        Mass of star in Msun.
    Age : float
        Age of the star in Myr. 
    Omega : float 
        Rotation rate in OmegaSun (=2.67e-6 rad/s).
    AgeMin : float , optional
        Starting age (default = 1 Myr).
    params : dict , optional
        Parameters to determine behavior of code (default given in parameters.py).
    StarEvo : Mors.stellarevo.StarEvo
        Instance of StarEvo class holding stellar evolution models to use.
    
    Returns
    ----------
    Omega0 : dict
        Starting rotation rate in OmegaSun (=2.67e-6 rad/s) or negative if cannot get initial rotation.
    
    """
    
    #---------------------------------------------------------------------
        
    # Make sure stellar mass is set
    if Mstar is None: 
        misc._PrintErrorKill("argument Mstar must be set in call to function")
    
    # Make sure initial rotation is set correctly
    if Omega is None: 
        misc._PrintErrorKill("argument Omega must be set in call to function")
    
    # Make sure AgeMin is set
    if AgeMin is None:
        AgeMin = params['AgeMinFit']
    
    # If an instance of the StarEvo class is not input, load it with defaults
    if StarEvo is None:
        StarEvo = SE.StarEvo()
    
    # Make sure this exact stellar mass is loaded (will do nothing if it already is)
    StarEvo.LoadTrack(Mstar)
        
    #---------------------------------------------------------------------
    
    # First thing to be checked is if the model can actually fit this rotation rate at this age since
    # we have defined minimum and maximum initial rotation rates in the parameters, given by
    
    # Set starting Omega0Min and Omega0Max for fitting
    Omega0Min = params['Omega0FitMin']
    Omega0Max = params['Omega0FitMax']
    
    # Get minimum rotation rate at Age corresponding to Omega0Min
    Tracks = EvolveRotation(Mstar=Mstar,Omega0=Omega0Min,AgeMin=params['AgeMinFit'],AgeMax=Age,params=params,StarEvo=StarEvo)
    OmegaMin = Tracks['OmegaEnv'][-1]
    
    # Get maximum rotation rate at Age corresponding to Omega0Max
    Tracks = EvolveRotation(Mstar=Mstar,Omega0=Omega0Max,AgeMin=params['AgeMinFit'],AgeMax=Age,params=params,StarEvo=StarEvo)
    OmegaMax = Tracks['OmegaEnv'][-1]
    
    # Check that Omega is between OmegaMin and OmegaMax
    if ( Omega < OmegaMin ):
        return -1
    if ( Omega > OmegaMax ):
        return -2
    
    #---------------------------------------------------------------------
    
    # Initially assume not found solution
    found = False
        
    # Start looping until solution has been found
    for iStep in range(0,params['nStepMaxFit']):
            
        # Get middle Omega0
        Omega0Mid = 0.5 * ( Omega0Min + Omega0Max )
        
        # Get rotation at mid spot
        Tracks = EvolveRotation(Mstar=Mstar,Omega0=Omega0Mid,AgeMin=params['AgeMinFit'],AgeMax=Age,params=params,StarEvo=StarEvo)
        OmegaMid = Tracks['OmegaEnv'][-1]
        
        # Check if OmegaMid is close enough to the desired answer
        if ( abs(OmegaMid/Omega-1.0) < params['toleranceFit'] ):
            found = True
            Omega0 = Omega0Mid
            break
        
        # Check if Omega0Min and Omega0Max are is close enough together (closest they will get)
        if ( abs(Omega0Min/Omega0Max-1.0) < params['toleranceFit'] ):
            found = True
            Omega0 = Omega0Mid
            break
    
        # Change either min or max Omega0
        # Note that we already know that the solution is between the min and max
        if ( Omega < OmegaMid ):
            Omega0Max = Omega0Mid
            OmegaMax = OmegaMid
        else:
            Omega0Min = Omega0Mid
            OmegaMin = OmegaMid
    
    # Return -3 if not found result
    if not found:
        return -3
  
    #---------------------------------------------------------------------
    
    return Omega0

def EvolveRotation(Mstar=None,Omega0=None,OmegaEnv0=None,OmegaCore0=None,AgeMin=None,AgeMax=None,AgesOut=None,params=params.paramsDefault,StarEvo=None):
    """Takes stellar mass and rotation rate, evolves rotation. 
    
    This is the main function for rotational evolution calculations. It will take some basic stellar parameters and
    evolve rotation between two ages.
    
    Parameters
    ----------
    Mstar : float
        Mass of star in Msun.
    Omega0 : float , optional
        Starting rotation rate in OmegaSun (=2.67e-6 rad/s).
    OmegaEnv0 : float , optional
        Starting envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
    OmegaCore0 : float , optional
        Starting envelope rotation rate in OmegaSun (=2.67e-6 rad/s).
    AgeMin : float , optional
        Starting age (default = 1 Myr).
    AgeMax : float , optional
        Final age (default = 5 Gyr).
    AgesOut : float or numpy.ndarray , optional
        Ages to output data for in Myr.
    params : dict , optional
        Parameters to determine behavior of code (default given in parameters.py).
    StarEvo : Mors.stellarevo.StarEvo
        Instance of StarEvo class holding stellar evolution models to use.
    
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
        
    # Make sure AgeMin is set
    if AgeMin is None:
        AgeMin = params['AgeMinDefault']
        
    # Make sure AgeMax is set
    if AgeMax is None:
        AgeMax = params['AgeMaxDefault']
        
    # Make editable and np array form of AgesOut and use this from now on
    if not AgesOut is None:
        if isinstance(AgesOut,(float,int)):
            AgesOut2 = np.array([AgesOut])
        else:
            AgesOut2 = copy.deepcopy(AgesOut)
    else:
        AgesOut2 = None
    
    # If AgesOut has been set, use the final age from that as ending age
    if not AgesOut2 is None:
        AgeMax = AgesOut2[-1]
    
    # If an instance of the StarEvo class is not input, load it with defaults
    if StarEvo is None:
        StarEvo = SE.StarEvo()
    
    # Make sure this exact stellar mass is loaded (will do nothing if it already is)
    StarEvo.LoadTrack(Mstar)
    
    # Starting values for time
    iAge = 0
    Age = AgeMin
    
    # Starting timestep in Myr
    dAge = 0.5
    
    # If AgesOut was set, make sure no ages are below AgeMin
    if not AgesOut2 is None:
        AgesOut2 = np.delete( AgesOut2 , np.where(AgesOut2<AgeMin) )
        
    # Starting rotation rates
    if Omega0 is None:
        OmegaEnv = OmegaEnv0
        OmegaCore = OmegaCore0
    else:
        OmegaEnv = Omega0
        OmegaCore = Omega0
        
    # Make dictionary for holding tracks and add initial values
    Tracks = _CreateTracks(Mstar,Age,dAge,OmegaEnv,OmegaCore,params,StarEvo)
        
    # Start looping, end only when end time is reached (additional 0.999 factor to stop code getting stuck)
    while ( Age < 0.999*AgeMax ):
        
        # Make sure timestep is not too long
        dAge , dAgeMax = _dAgeCalc(dAge,Age,AgeMax,AgesOut2,params)
        
        # Do timestep
        dAge , dAgeNew , OmegaEnv , OmegaCore = EvolveRotationStep( Mstar=Mstar , Age=Age , 
                                                                OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , 
                                                                dAgeMax=dAgeMax , dAge=dAge , 
                                                                params=params , StarEvo=StarEvo )
        
        # Update time
        iAge += 1
        Age += dAge
        
        # Work out if should add this age to output
        AgesOut2 , shouldAppend = _shouldAppend(Age,AgesOut2)
            
        # Add current state to tracks if should output
        if shouldAppend:
            Tracks = _AppendTracks(Tracks,Mstar,Age,dAge,OmegaEnv,OmegaCore,params,StarEvo)
        
        # Check for bad data
        _CheckBadData(Tracks)
        
        # Get new age step
        dAge = dAgeNew
        
        # Determine if too many steps taken
        if ( iAge >= params['nStepMax'] ):
            misc._PrintErrorKill("too many timesteps taken")
        
    # In Tracks, replace Mstar array with simple Mstar float
    Tracks['Mstar'] = Mstar
    Tracks['nAge'] = len(Tracks['Age'])
    
    return Tracks

def _dAgeCalc(dAge,Age,AgeMax,AgesOut,params):
    """Takes information about age in evolutionary calculation, returns age step to use."""
    
    # Get maximum timestep in Myr
    dAgeMax = AgeMax - Age
    
    # If AgesOut has been set, work out maximum from that
    if not AgesOut is None:
            
        # If AgesOut is just a number (float or int) then it is easy
        if isinstance(AgesOut,(float,int)):
            
            # Get new maximum age
            dAgeMax = AgesOut - Age
            
        else:
            
            # Get index of next age in AgesOut after current Age
            index = misc._getIndexGT(AgesOut,Age)
            
            # Get new maximum age
            dAgeMax = AgesOut[index] - Age
                
    # Make sure timestep is not too long
    dAge = min( dAge , dAgeMax )
    
    # Make sure dAge is within range in parameter file
    dAge = max( dAge , params['dAgeMin'] )
    dAge = min( dAge , params['dAgeMax'] )
    
    return dAge , dAgeMax

def _shouldAppend(Age,AgesOut):
    """Takes information about age in evolutionary calculation, returns if age should be output."""
    
    # If AgesOut is not set, then should always output
    if AgesOut is None:
        return AgesOut , True
    
    # If AgesOut is set, work out if this is one of the output ages
    if not AgesOut is None:
            
        # If AgesOut is just a number (float or int) then it is easy
        if isinstance(AgesOut,(float,int)):
            
            # Check if close enough to AgesOut to output
            if ( abs(Age/float(AgesOut))-1.0 < 1.0e-6 ):
                return AgesOut , True
            else:
                return AgesOut , False
            
        else:
            
            # Get smallest distance between age in AgesOut and current Age
            indexMin = np.argmin( np.abs(AgesOut-Age) )
            deltaAgeMin = np.min( np.abs(AgesOut-Age) )
            
            # Since we are removing this output age each time, indexMin should always be zero
            if not indexMin == 0:
                misc._PrintErrorKill("indexMin is not zero, probably because the input AgesOut is not in ascending order")
            
            # If age is close enough to element in AgesOut then remove this element and return True
            if deltaAgeMin < 1.0e-3:
                AgesOut = np.delete(AgesOut,indexMin)
                return AgesOut , True
            
            else:
                return AgesOut , False
  
    return AgesOut , False

def EvolveRotationStep(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAgeMax=None,dAge=None,params=params.paramsDefault,StarEvo=None):
    """Evolves rotation by a single step.
    
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
        
    # If an instance of the StarEvo class is not input, load it with defaults
    if StarEvo is None:
        StarEvo = SE.StarEvo()
        
    # Do step based on method used
    if ( params['TimeIntegrationMethod'] == 'ForwardEuler' ):
        
        dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepFE( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , 
                                                                    OmegaCore=OmegaCore , dAge=dAge , params=params , StarEvo=StarEvo )
    
    elif ( params['TimeIntegrationMethod'] == 'RungeKutta4' ):
        
        dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepRK4( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , 
                                                                        OmegaCore=OmegaCore , dAge=dAge , params=params , StarEvo=StarEvo )
    
    elif ( params['TimeIntegrationMethod'] == 'RungeKuttaFehlberg' ):
        
        dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepRKF( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , 
                                                                        OmegaCore=OmegaCore , dAge=dAge , dAgeMax=dAgeMax ,
                                                                        params=params , StarEvo=StarEvo )
        
    elif ( params['TimeIntegrationMethod'] == 'Rosenbrock' ):
        
        dAge , dAgeNew , OmegaEnv , OmegaCore =  _EvolveRotationStepRB( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , 
                                                                        OmegaCore=OmegaCore , dAge=dAge , dAgeMax=dAgeMax ,
                                                                        params=params , StarEvo=StarEvo )
    
    else:
        misc._PrintErrorKill("invalid value of TimeIntegrationMethod in parameters")
    
    return ( dAge , dAgeNew , OmegaEnv , OmegaCore )

def _EvolveRotationStepFE(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAge=None,params=params.paramsDefault,StarEvo=None):
    """Takes basic stellar parameters, evolves by timestep using forward Euler method."""
    
    # Get rates of change
    dOmegaEnvdt , dOmegaCoredt = phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params , StarEvo=StarEvo )
    
    # Do update
    OmegaEnv += dAge * dOmegaEnvdt
    OmegaCore += dAge * dOmegaCoredt
    
    return ( dAge , dAge , OmegaEnv , OmegaCore )

def _EvolveRotationStepRK4(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAge=None,params=params.paramsDefault,StarEvo=None):
    """Takes basic stellar parameters, evolves by timestep using classical Runge-Kutta method."""
    
    # Get rates of change
    k1Env , k1Core = phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params , StarEvo=StarEvo )
    k2Env , k2Core = phys.dOmegadt( Mstar=Mstar , Age=Age+0.5*dAge , OmegaEnv=OmegaEnv+0.5*dAge*k1Env , OmegaCore=OmegaCore+0.5*dAge*k1Core , params=params , StarEvo=StarEvo )
    k3Env , k3Core = phys.dOmegadt( Mstar=Mstar , Age=Age+0.5*dAge , OmegaEnv=OmegaEnv+0.5*dAge*k2Env , OmegaCore=OmegaCore+0.5*dAge*k2Core , params=params , StarEvo=StarEvo )
    k4Env , k4Core = phys.dOmegadt( Mstar=Mstar , Age=Age+dAge , OmegaEnv=OmegaEnv+dAge*k3Env , OmegaCore=OmegaCore+dAge*k3Core , params=params , StarEvo=StarEvo )
        
    # Do update
    OmegaEnv += dAge * ( k1Env + 2.0*k2Env + 2.0*k3Env + k4Env ) / 6.0
    OmegaCore += dAge * ( k1Core + 2.0*k2Core + 2.0*k3Core + k4Core ) / 6.0
    
    return ( dAge , dAge , OmegaEnv , OmegaCore )

def _EvolveRotationStepRKF(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAgeMax=None,dAge=None,params=params.paramsDefault,StarEvo=None):
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
        k1 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age , OmegaEnv=X[0] , OmegaCore=X[1] , params=params , StarEvo=StarEvo ) )
        
        # k2
        Age2 = Age + (1.0/5.0)*dAge
        X2 = X + (1.0/5.0)*dAge*k1
        k2 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age2 , OmegaEnv=X2[0] , OmegaCore=X2[1] , params=params , StarEvo=StarEvo ) )
        
        # k3
        Age3 = Age + (3.0/10.0)*dAge
        X3 = X + (3.0/40.0)*dAge*k1 + (9.0/40.0)*dAge*k2
        k3 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age3, OmegaEnv=X3[0] , OmegaCore=X3[1] , params=params , StarEvo=StarEvo ) )
        
        # k4
        Age4 = Age + (3.0/5.0)*dAge
        X4 = X + (3.0/10.0)*dAge*k1 - (9.0/10.0)*dAge*k2 + (6.0/5.0)*dAge*k3
        k4 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age4, OmegaEnv=X4[0] , OmegaCore=X4[1] , params=params , StarEvo=StarEvo ) )
        
        # k5
        Age5 = Age + dAge
        X5 = X - (11.0/54.0)*dAge*k1 + (5.0/2.0)*dAge*k2 - (70.0/27.0)*dAge*k3 + (35.0/27.0)*dAge*k4
        k5 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age5, OmegaEnv=X5[0] , OmegaCore=X5[1] , params=params , StarEvo=StarEvo ) )
        
        # k6
        Age6 = Age + (7.0/8.0)*dAge
        X6 = X + (1631.0/55296.0)*dAge*k1 + (175.0/512.0)*dAge*k2 + (575.0/13824.0)*dAge*k3 + (44275.0/110592.0)*dAge*k4 + (253.0/4096.0)*dAge*k5
        k6 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age6, OmegaEnv=X6[0] , OmegaCore=X6[1] , params=params , StarEvo=StarEvo ) )
        
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

def _EvolveRotationStepRB(Mstar=None,Age=None,OmegaEnv=None,OmegaCore=None,dAgeMax=None,dAge=None,params=params.paramsDefault,StarEvo=None):
    """Takes basic stellar parameters, evolves by timestep using Runge-Kutta-Fehlberg method."""
    
    # Get coefficients for solver
    CoefficientsRB = params['CoefficientsRB']
    
    # Start with dAge as timestep
    dAgeNew = dAge
    
    # Setup array to hold integration values
    X = np.array([OmegaEnv,OmegaCore])
    nVar = len(X)
    
    # Get Jacobian and assume it is a constant over timestep
    Jac = _JacobianRB( Mstar , Age+dAge , X , nVar , params , StarEvo )
  
    # Start iterating until got accuracy desired
    while True:
        
        # Set dAge
        dAge = dAgeNew
        
        # Calculate the k coefficients
        kCoeff = _kCoeffRB( Mstar , Age , dAge , X , Jac , nVar , CoefficientsRB , params , StarEvo )
        
        # Get updated values
        Xnew = copy.deepcopy(X)
        for i in range(0,CoefficientsRB['s']):
            Xnew += CoefficientsRB['b'][i] * kCoeff[:,i]
    
        # Get less accurate updated values
        Xnew2 = copy.deepcopy(X)
        for i in range(0,CoefficientsRB['s']):
            Xnew2 += CoefficientsRB['b2'][i] * kCoeff[:,i]
        
        # Get error in this timestep for each species and grid point
        Error = np.sqrt( np.sum( ( ( Xnew - Xnew2 ) / params['DeltaDesired'] )**2.0 ) / float(nVar) ) 
        
        # Get the factor by which to change dt (if Error is zero then just increase dAge by 50%)
        if not ( Error == 0.0 ):
            dAgeFactor = min( 10.0 , max( 0.1 , 0.9/Error**(1.0/CoefficientsRB['order']) ) )
        else:
            dAgeFactor = 1.5
        
        # get new dt
        dAgeNew = 0.99 * dAge * dAgeFactor
        
        # Make sure dt is not too big
        dAgeNew = min( dAgeNew , dAgeMax )
        
        # Check for stopping of iteration
        if ( dAgeFactor > 1.0 ):
            break
            
    # Save new estimate
    OmegaEnv = Xnew[0]
    OmegaCore = Xnew[1]
    
    return dAge , dAgeNew , OmegaEnv , OmegaCore

def _JacobianRB(Mstar,Age,X,nVar,params,StarEvo):
    """Calculates Jacobian d(dOmega/dt)/dOmega terms for Rosenbrock solver."""
    
    # Make array
    Jac = np.zeros((nVar,nVar))
    
    # Loop over elements and fill in Jacobian
    for iVar in range(0,nVar):
        
        # Set perturbed X values to X
        X1 = copy.deepcopy(X)
        X2 = copy.deepcopy(X)
        
        # Perturb this variable
        X1[iVar] = X1[iVar] - params['deltaJac'] * X1[iVar]
        X2[iVar] = X2[iVar] + params['deltaJac'] * X2[iVar]
        
        # Get rates based on pertubed quantities
        dXdt1 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age, OmegaEnv=X1[0] , OmegaCore=X1[1] , params=params , StarEvo=StarEvo ) )
        dXdt2 = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age, OmegaEnv=X2[0] , OmegaCore=X2[1] , params=params , StarEvo=StarEvo ) )
        
        # Fill in this column of Jacobian
        Jac[:,iVar] = ( dXdt2[:] - dXdt1[:] ) / ( X2[iVar] - X1[iVar] )
        
    return Jac

def _kCoeffRB(Mstar,Age,dAge,X,Jac,nVar,CoefficientsRB,params,StarEvo):
    """Calculates kCoeff for Rosenbrock solver."""
    
    # Make array for holding result
    kCoeff = np.zeros((nVar,CoefficientsRB['s']))
    
    # Get the function involving the Jacobian, given by (I - dt * gamma_ii * J)
    # note: since all gamma_ii coefficients are equal, this function only needs to
    #       be calculated once and can then be used for all values of i
    # however, it must be recalculated if the timestep is changed
    
    # Set the function to -dt*J for all elements
    JacFunc = - dAge * CoefficientsRB['gamma'][0,0] * Jac
    
    # Loop over diagonal elements and add 1 to each
    for iVar in range(0,nVar):
        JacFunc[iVar,iVar] += 1.0
        
    # Get k1
    dXdt = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age+dAge, OmegaEnv=X[0] , OmegaCore=X[1] , params=params , StarEvo=StarEvo ) )
    RateFunction = dAge * dXdt
    kCoeff[:,0] = _GaussianElimination( JacFunc , RateFunction )
        
    # Get k2,...ks
    for i in range(1,CoefficientsRB['s']):
        
        # Get intermediate values
        Xmid = copy.deepcopy(X)
        for i2 in range(0,i):
            Xmid += CoefficientsRB['alpha'][i,i2] * kCoeff[:,i2]
            
        # Get intermediate rates
        dXdt = np.array( phys.dOmegadt( Mstar=Mstar , Age=Age+dAge, OmegaEnv=Xmid[0] , OmegaCore=Xmid[1] , params=params , StarEvo=StarEvo ) )
            
        # Get the rate function
        RateFunction = dAge * dXdt
        for i2 in range(0,i):
            for iVar in range(0,nVar):
                RateFunction[iVar] += dAge * CoefficientsRB['gamma'][i,i2] * np.sum( Jac[iVar,:] * kCoeff[:,i2] )
            
        # Get k value for this step
        kCoeff[:,i] = _GaussianElimination( JacFunc , RateFunction )
  
    return kCoeff

def _GaussianElimination(Ain,B):
    """Solves the equation A*X = B using Gaussian elimination."""
    
    # This is meant to be used when getting the k coefficients in the Rosenbrock
    # method. This subroutine is able to be relatively simple because the NxN 
    # matrix A is given by I-dt*gamma*J, where I is the identity matrix. It is 
    # highly unlikely that any of the diagonal terms in A are zero, since that
    # would only really be the case if a diagonal term in J is exactly 1.0/(dt*gamma),
    # which would be quite a coincidence. If any diagonal terms in A are zero,
    # this subroutine will give a nonsense result. I have added an if statement
    # below to catch these errors, so the user will at least know what is going on.
    
    #---------------------------------------
    # GET N AND ALLOCATE X AND A
    
    N = len(B)
    X = np.zeros(N)
    A = np.zeros((N,N))
    
    # ---------------------------------------
    # COPY A and B
    
    A[:,:] = Ain[:,:]
    X[:] = B[:]
    
    #---------------------------------------
    # GET IN ROW ECHELON FORMAT
        
    # Loop over rows
    for i in range(0,N):
    
        # Make sure diagonal term not zero
        if ( A[i,i] == 0.0 ):
            misc._PrintErrorKill("diagonal term in zero")
        
        # Divide row by value in diagonal element
        X[i] = X[i] / A[i,i]
        A[i,i:N] = A[i,i:N] / A[i,i]
        
        # Get rid of zeros in columns below this diagonal
        # unless doing final row
        if ( i < N-1 ):
            for k in range(i+1,N):
                if ( A[k,i] != 0.0 ):
                    X[k] += - X[i] * A[k,i]
                    A[k,:] += - A[i,:] * A[k,i]
        
    #---------------------------------------
    # GET IN REDUCED ROW ECHELON FORMAT
    
    for i in range(N-2,-1,-1):
        for j in range(N-1,i,-1):
            if ( A[i,j] != 0.0 ):
                X[i] += - A[i,j] * X[j]
                A[i,:] += - A[i,j] * A[j,:]
    
    #---------------------------------------
    
    return X

def _CreateTracks(Mstar,Age,dAge,OmegaEnv,OmegaCore,params,StarEvo):
    """Takes simulation parameters, returns dictionary with empty arrays for evolutionary tracks."""
    
    # Create empty dictionary
    Tracks = {}
    
    # Add basic tracks that are always included
    Tracks['Age'] = np.array([])
    Tracks['dAge'] = np.array([])
    Tracks['OmegaEnv'] = np.array([])
    Tracks['OmegaCore'] = np.array([])
    
    # If doing extended tracks, add these
    if params['ExtendedTracks']:
  
        # Add extended activity quantities to dictionary
        # note: do not need to call phys.RotationQuantities since without StarState being passed as an argument, it will be called 
        #       by ExtendedQuantities() anyway.
        StarState = phys.ExtendedQuantities( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params , StarEvo=StarEvo )
        
        # Add empty arrays for each of these values 
        for quantity in StarState:
            Tracks[quantity] = np.array([])
  
    else:
        
        # Just set StarState to None since it is not needed
        StarState = None
    
    # Add initial values to tracks
    Tracks = _AppendTracks(Tracks,Mstar,Age,dAge,OmegaEnv,OmegaCore,params,StarEvo,StarState=StarState)
    
    return Tracks
  
def _AppendTracks(Tracks,Mstar,Age,dAge,OmegaEnv,OmegaCore,params,StarEvo,StarState=None):
    """Takes state of evolutionary simulation, appends state onto arrays in Tracks dictionary."""
    
    # Basic quantities that are always included
    Tracks['Age'] = np.append( Tracks['Age'] , Age )
    Tracks['dAge'] = np.append( Tracks['dAge'] , 0.0 )
    Tracks['OmegaEnv'] = np.append( Tracks['OmegaEnv'] , OmegaEnv )
    Tracks['OmegaCore'] = np.append( Tracks['OmegaCore'] , OmegaCore )
    
    # If doing extended tracks, calculate and add them
    if params['ExtendedTracks']:
        
        # Calculate extended quantities if not input as keyword argument
        if StarState is None:
            
            # note: do not need to call phys.RotationQuantities since without StarState being passed as an argument, it will be called 
            #       by ExtendedQuantities() anyway.
            StarState = phys.ExtendedQuantities( Mstar=Mstar , Age=Age , OmegaEnv=OmegaEnv , OmegaCore=OmegaCore , params=params , StarEvo=StarEvo )
            
        # Append these quantities, but only if they are not already present
        for quantity in StarState:
            
            # Do not do quantity if it was already added above
            if not ( quantity in ['Age','dAge','OmegaEnv','OmegaCore'] ):
                Tracks[quantity] = np.append( Tracks[quantity] , StarState[quantity] )
        
    return Tracks

def _CheckBadData(Tracks):
    """Takes current set of tracks, checks last entries for bad data."""
    
    # Originally assume no error and sould continue
    error = False
    
    # Make sure rotation rates are positive
    if ( Tracks['OmegaEnv'][-1] <= 0.0 ):
        misc._PrintError("negative value of OmegaEnv in rotational evolution calculation")
        error = True
    if ( Tracks['OmegaCore'][-1] <= 0.0 ):
        misc._PrintError("negative value of OmegaCore in rotational evolution calculation")
        error = True
        
    # Loop over elements in Tracks and check each one for NaN or Inf
    for track in Tracks:
        if math.isnan(Tracks[track][-1]):
            misc._PrintError("NaN found for "+track+" in rotational evolution calculation")
            error = True
        if math.isinf(Tracks[track][-1]):
            misc._PrintError("inf found for "+track+" in rotational evolution calculation")
            error = True
    
    # If there was an error, kill the code
    if error:
        misc._PrintErrorKill("ending simulation due to bad data found")
    
    return
