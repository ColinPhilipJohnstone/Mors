
"""Module holding the Cluster class and related functions."""

# Imports for standard stuff needed here
import sys
import inspect
import numpy as np
import os
import pickle

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.stellarevo as SE
import Mors.parameters as params
import Mors.rotevo as RE
import Mors.physicalmodel as phys
import Mors.star as star

class Cluster:
  """
  A class for star objects that hold all information about a star. 
  
  Attributes
  ------------
  
  Methods
  ------------
  
  """
    
  def __init__(self,Mstar=None,Age=None,Omega=None,OmegaEnv=None,OmegaCore=None,AgesOut=None,starEvoDir=None,evoModels=None,params=params.paramsDefault,verbose=False):
    """
    Initialises instance of Cluster class.
    
    This is the main function that is run when creating an instance of the Cluster class and it sets up all
    the things needed including calculating evolutionary tracks for each of the stars star. The function requires 
    that an array of stellar masses (in Msun) and initial rotation rates (in OmegaSun=2.67e-6 rad s^-1) are input.
    using the Mstar and Omega keyword arguments. Alternatively, OmegaEnv and OmegaCore can be set, in which case
    Omega does not need to be specified. If the argument Age (in Myr) is also specified, then the code will find 
    the evolutionary tracks for each star that passes through this rotation rate at this age, otherwise if Age is 
    not set then it will calculate evolutionary tracks assuming this Omega as the initial (1 Myr) rotation rate. 
    The user should not specify OmegaCore and Age simultaneously, and if Age is set then either Omega or OmegaEnv 
    can be used to specify the surface rotation rate. Age can also be set as an array giving different ages for 
    each star.
    
    """
    
    # Make sure Mstar was set
    if Mstar is None:
      misc._PrintErrorKill("Mstar keyword argument not set")
    
    # Make sure Omega, OmegaEnv, and OmegaCore are well set
    # (if Omega is set, OmegaEnv and OmegaCore will both be set to that value)
    Omega , OmegaEnv , OmegaCore = _CheckInputRotation(Age,Omega,OmegaEnv,OmegaCore)
    
    # Make sure arrays all same length
    if not ( len(Mstar) == len(Omega) ):
      misc._PrintErrorKill("Mstar and Omega have different lengths")
    if not ( len(Mstar) == len(OmegaEnv) ):
      misc._PrintErrorKill("Mstar and OmegaEnv have different lengths")
    if not ( len(Mstar) == len(OmegaCore) ):
      misc._PrintErrorKill("Mstar and OmegaCore have different lengths")
    
    # Set number of stars
    self.nStars = len(Mstar)
    
    # Set parameters
    self.params = params
    
    # Set the ExtendedTracks parameter to True so we get all parameters
    self.params['ExtendedTracks'] = True
    
    # Set Mstar arrays
    self.Mstar = Mstar
    
    # Set AgesOut 
    self.AgesOut = AgesOut
    
    # Set stellar evo model
    self.starEvoDir = starEvoDir
    self.evoModels = evoModels
    
    # Get evolutionary tracks 
    self._LoadEvoTracks(Age,OmegaEnv,OmegaCore,verbose)
    
    # Get HZ boundaries
    self.aOrbHZ = phys.aOrbHZ(Mstar=self.Mstar,params=self.params)
    
    return
    
  def _LoadEvoTracks(self,Age,OmegaEnv0,OmegaCore0,verbose):
    """Loads rotation and activity tracks for each star."""
    
    # Start dictionary holding stars
    self.stars = []
    
    # Loop over each star and load each one
    for iStar in range(0,len(self.Mstar)):
      
      # Print to screen if verbose was set by user
      if verbose:
        #print(f"\rLoading star "+str(iStar)+" out of "+str(self.nStars),end="")
        print("Loading star "+str(iStar)+" out of "+str(self.nStars)+" in cluster",end="\r")
      
      # Get parameters for this star in right type
      MstarStar = float(self.Mstar[iStar])
      OmegaEnvStar = float(OmegaEnv0[iStar])
      OmegaCoreStar = float(OmegaCore0[iStar])
      
      # Get Age depending on if it is an array or a float/int
      if Age is None:
        AgeIn = None
      elif isinstance(Age,(float,int)):
        AgeIn = Age
      else:
        AgeIn = Age[iStar]
      
      # Create instance of star class for this star
      starTemp = star.Star(Mstar=MstarStar,Age=AgeIn,OmegaEnv=OmegaEnvStar,OmegaCore=OmegaCoreStar,
                           AgesOut=self.AgesOut,starEvoDir=self.starEvoDir,evoModels=self.evoModels,params=self.params)
      
      # Append dictionary
      self.stars.append(starTemp)
      
      # Also make this star an attribute of the class
      setattr( self , "star"+str(iStar) , starTemp )
    
    # Make new line if printing status to screen
    if verbose:
      print("")
    
    # Make functions for each quantity that return this quantity at a given age as attributes of class
    self._setupQuantityFunctions()
    
    return
    
  def _setupQuantityFunctions(self):
    """Makes functions for each quantity that return this quantity at a given age as attributes of class."""
    
    # A description of how the code below works is given in the star class function with the same name
    # as this function. In this case, the code does it based on the quantities held in the first star. 
    for track in self.stars[0].Tracks:
      
      exec( "def "+track+"(self,Age):\n  return self.Values(Age=Age,Quantity='"+track+"')" )
      exec( "setattr( self.__class__ , '"+track+"' , "+track+" )" )
    
    return
  
  def Values(self,Age=None,Quantity=None):
    """
    Takes age in Myr and a string with name of quantity to output, returns value of that quantity at specified age for all stars.
    
    
    Parameters
    ----------
    Age : float
        Age of cluster.
    Quantity : str
        String holding name of parameter to get.
    
    Returns
    ----------
    values : numpy.ndarray
        Set of extended quantities.
    
    """
    
    # Make sure input parameters are set
    if Age is None:
      misc._PrintErrorKill("keyword parameter Age not set in call to function")
    if Quantity is None:
      misc._PrintErrorKill("keyword parameter Quantity not set in call to function")
    
    # Make array to hold values
    values = np.zeros(self.nStars)
    
    # Loop over stars and get values for each
    for iStar in range(0,self.nStars):
      values[iStar] = self.stars[iStar].Value(Age=Age,Quantity=Quantity)
    
    return values
    
  def PrintStars(self):
    """Prints list of stars in cluster to screen."""
    
    # Header
    print("The following is a list of masses for stars in this cluster.")
    
    # Loop over each star and print basic parameters
    for iStar in range(0,self.nStars):
      print("   "+str(iStar)+". "+str(self.Mstar[iStar])+" Msun")
    
    # Total number of stars
    print("Number of stars in cluster = "+str(self.nStars))
    
    return
  
  def Save(self,filename='cluster.pickle'):
    """Takes filename (default is 'cluster.pickle'), saves cluster to this file using pickle."""
    
    with open(filename,'wb') as f:
      pickle.dump(self,f)
    
    return
  
  def save(self,filename='cluster.pickle'):
    """Same as Save()."""
    self.Save(filename=filename)
    return    
    
  def Percentile(self,Mstar=None,Age=None,Omega=None,Prot=None,percentile=None):
    """
    Gets rotation rate of percentile or percentile of rotation rate in the rotation distribution at given age.
    
    This function can be used for two purposes
      1. to determine the percentile in a rotation distribution of a star given its mass, rotation rate, and age
      2. to determine the rotation rate of a star in a rotation distribution given its mass, percentile, and age
    In the first case, the user should specify the rotation rate using either the Omega or Prot keyword arguments (given 
    in OmegaSun and days respectively) and the mass. In the second case, the user should specify the percentile using the
    percentile keyword argument (in the range 0 to 100) and the mass. The mass should be specified in Msun using the Mstar
    keyword argument and the age in Myr using the Age keyword argument. These should only be given as floats. 
    
    Parameters
    ----------
    Mstar : float
        Stellar mass in Msun.
    Age : float
        Age of star in Myr.
    Omega : float , optional
        Rotation rate of star in OmegaSun.
    Prot : float , optional
        Rotation period of star in days.
    percentile : float , optional
        percentile in distribution (between 0 and 100).
    
    Returns
    ----------
    result : float
        Either rotation rate for percentile or percentile for rotation rate.
    
    """
    
    # Make sure the age is specified
    if Age is None:
      misc._PrintErrorKill("keyword parameter Age not set in call to function")
    
    # If percentile was set to a string, get float version
    if isinstance(percentile,str):
      if ( percentile == 'slow' ):
        percentile = 5.0
      elif ( percentile == 'medium' ):
        percentile = 50.0
      elif ( percentile == 'fast' ):
        percentile = 95.0
      else:
        misc._PrintErrorKill( "invalid percentile string (options are 'slow', 'medium', or 'fast')" )
    
    # Get the rotation distribution at this age
    OmegaDist = self.Values( Age=Age , Quantity='OmegaEnv' )
    
    # Get the result
    result = star.Percentile( Mstar=Mstar , Omega=Omega , Prot=Prot , percentile=percentile , 
                            MstarDist=self.Mstar , OmegaDist=OmegaDist , params=self.params )
    
    return result
  
  def ActivityLifetime(self,Quantity=None,Threshold=None,AgeMax=None):
    """
    Takes threshold value, returns ages at which each star last drops below this threshold.
    
    This function can be used to determine when each star's emission crosses a given threshold value for a few
    activity quantities. These are Lx, Fx, Rx, and FxHZ for X-rays, and similar values for EUV1, EUV2, EUV, 
    XUV, and Ly-alpha. If the star crosses the threshold (from above it to below it) multiple times, this 
    will find the final time it will cross the threshold. If the user wants to set a maximum age so that the 
    code only looks for crossings of the threshold below this age then this can be done using the AgeMax
    keyword argument.
        
    Parameters
    ----------
    Quantity : str
        Gives which quantity to consider (e.g. 'Lx').
    Threshold : float or str
        Value for threshold in units of quantity or string of 'sat'.
    AgeMax : float , optional
        End age of track to consider in Myr.
    
    Returns
    ----------
    AgeActive : float
        Activity lifetime in Myr.
    
    """
    
    # Make sure Quantity is set
    if Quantity is None:
      misc._PrintErrorKill("Quantity not set in call to function")
      
    # Make sure Quantity is string
    if not isinstance(Quantity,str):
      misc._PrintErrorKill("Quantity must be string")
    
    # Make array
    AgeActive = np.zeros(self.nStars)
    
    # Loop over stars and get each
    for iStar in range(self.nStars):
      AgeActive[iStar] = self.stars[iStar].ActivityLifetime(Quantity=Quantity,Threshold=Threshold,AgeMax=AgeMax)
    
    return AgeActive
    
  def IntegrateEmission(self,AgeMin=None,AgeMax=None,Band=None,aOrb=None):
    """
    Takes age range, returns integrated emission in band within that range for each star.
    
    This code can be used to the luminosities of the stars between two ages. This can be applied to any wavelength band
    and the result is a total energy emitted in this time in erg. If the user also specifies an orbital distance using 
    the aOrb keyword argument, the code integrates the flux at this obital distance returns a fluence in erg cm^-2. The
    user can specify aOrb as a string to get the fluences at various habitable zone boundaries (using the HZ calculated 
    at the age defined in params when creating this cluster). Options are 'RecentVenus', 'RunawayGreenhouse', 'MoistGreenhouse', 
    'MaximumGreenhouse', 'EarlyMars', and 'HZ'.
        
    Parameters
    ----------
    AgeMin : float
        Start of time period to integrate in Myr.
    AgeMax : float
        End of time period to integrate in Myr.
    Band : str
        Gives which wavelength band to consider (options are 'XUV', 'Xray', 'EUV1', 'EUV2', 'EUV', 'Lyman', 'bol').
    aOrb : float or str, optional
        Orbital distance to get fluence at in AU or string identifying HZ boundary.
    
    Returns
    ----------
    Energy : numpy.ndarray
        Integrated luminosity or flux in erg or erg cm^-2.
    
    """
    
    # Make array
    Energy = np.zeros(self.nStars)
    
    # Loop over stars and get each
    for iStar in range(self.nStars):
      Energy[iStar] = self.stars[iStar].IntegrateEmission(AgeMin=AgeMin,AgeMax=AgeMax,Band=Band,aOrb=aOrb)
    
    return Energy

def _CheckInputRotation(Age,Omega,OmegaEnv,OmegaCore):
  """Takes input rotation, checks if values are setup correctly."""
  
  # Make sure if Age is set that OmegaCore is not set
  if ( not Age is None ) and ( not OmegaCore is None ):
    misc._PrintErrorKill( "cannot set both Age and OmegaCore as arguments of Cluster" )
  
  # Make sure at least one rotation rate is set
  if ( Omega is None ):
    
    # Since Omega is not set, make sure both OmegaEnv and OmegaCore are set
    if ( OmegaEnv is None ) or ( OmegaCore is None ):
      misc._PrintErrorKill( "must set either Omega or both OmegaEnv and OmegaCore as arugment of Star" )
  
  else:
    
    # Since Omega is set, make sure both OmegaEnv and OmegaCore are not set
    if not ( ( OmegaEnv is None ) and ( OmegaCore is None ) ):
      misc._PrintErrorKill( "cannot set OmegaEnv and OmegaCore when Omega is set as arugment of Star" )
    
    # Set both OmegaEnv and OmegaCore to Omega
    OmegaEnv = Omega
    OmegaCore = Omega
  
  
  return Omega , OmegaEnv , OmegaCore
