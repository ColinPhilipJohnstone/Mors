
"""Module holding the Cluster class and related functions."""

#====================================================================================================================

# Imports for standard stuff needed here
import sys
import inspect
import numpy as np

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.stellarevo as SE
import Mors.parameters as params
import Mors.rotevo as RE
import Mors.physicalmodel as phys
import Mors.star as star

#====================================================================================================================

class Cluster:
  
  """
  A class for star objects that hold all information about a star. 
  
  Attributes
  ------------
  
  Methods
  ------------
  
  """
  
  #---------------------------------------------------------------------------------------
  
  def __init__(self,Mstar=None,Age=None,Omega=None,OmegaEnv=None,OmegaCore=None,AgesOut=None,starEvoDir=None,evoModels=None,params=params.paramsDefault):
    
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
    
    # Set parameters
    self.params = params
    
    # Set the ExtendedTracks parameter to True so we get all parameters
    self.params['ExtendedTracks'] = True
    
    # Set Mstar and rotation arrays
    self.Mstar = Mstar
    self.Omega = Omega
    self.OmegaEnv = OmegaEnv
    self.OmegaCore = OmegaCore
    
    # Set AgesOut 
    self.AgesOut = AgesOut
    
    # Set stellar evo model
    self.starEvoDir = starEvoDir
    self.evoModels = evoModels
    
    # Get evolutionary tracks 
    self._LoadEvoTracks(Age)
    
    return
  
  
  #---------------------------------------------------------------------------------------
  
  def _LoadEvoTracks(self,Age):
    
    """Loads rotation and activity tracks for each star."""
    
    # Start dictionary holding stars
    self.stars = []
    
    # Loop over each star and load each one
    for iStar in range(0,len(self.Mstar)):
      
      # Get parameters for this star in right type
      Mstar = float(self.Mstar[iStar])
      Omega = float(self.Omega[iStar])
      OmegaEnv = float(self.OmegaEnv[iStar])
      OmegaCore = float(self.OmegaCore[iStar])
      
      # Get Age depending on if it is an array or a float/int
      if Age is None:
        AgeIn = None
      elif ( type(Age) == float ) or ( type(Age) == int ):
        AgeIn = Age
      else:
        AgeIn = Age[iStar]
      
      # Create instance of star class for this star
      starTemp = star.Star(Mstar=Mstar,Age=AgeIn,OmegaEnv=OmegaEnv,OmegaCore=OmegaCore,AgesOut=self.AgesOut,starEvoDir=self.starEvoDir,evoModels=self.evoModels,params=self.params)
      
      # Append dictionary
      self.stars.append(starTemp)
      
      # Also make this an attribute of the class
      setattr( self , "star"+str(iStar) , starTemp )
    
    return
  
  #---------------------------------------------------------------------------------------
  
  
  #---------------------------------------------------------------------------------------
  
#====================================================================================================================

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

#====================================================================================================================
