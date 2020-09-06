
"""Module holding the Cluster class and related functions."""

#====================================================================================================================

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
    
    return
  
  #---------------------------------------------------------------------------------------
  
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
      elif ( type(Age) == float ) or ( type(Age) == int ):
        AgeIn = Age
      else:
        AgeIn = Age[iStar]
      
      # Create instance of star class for this star
      starTemp = star.Star(Mstar=MstarStar,Age=AgeIn,OmegaEnv=OmegaEnvStar,OmegaCore=OmegaCoreStar,AgesOut=self.AgesOut,starEvoDir=self.starEvoDir,evoModels=self.evoModels,params=self.params)
      
      # Append dictionary
      self.stars.append(starTemp)
      
      # Also make this star an attribute of the class
      setattr( self , "star"+str(iStar) , starTemp )
    
    # Make new line if printing status to screen
    if verbose:
      print("")
    
    # Make functions for each quantity that return this quantity at a given age as attributes of class
    # A description of how the code below works is given in the star class function with the same name
    # as this function. In this case, the code does it based on the quantities held in the first star. 
    for track in self.stars[0].Tracks:
      
      exec( "def "+track+"(self,Age):\n  return self.Values(Age=Age,Quantity='"+track+"')" )
      exec( "setattr( self.__class__ , '"+track+"' , "+track+" )" )
    
    return
  
  #---------------------------------------------------------------------------------------
  
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
  
  #---------------------------------------------------------------------------------------
  
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
  
  #---------------------------------------------------------------------------------------
  
  def Save(self,filename='cluster.pickle'):
    
    """Takes filename (default is 'cluster.pickle'), saves cluster to this file using pickle."""
    
    with open(filename,'wb') as f:
      pickle.dump(self,f)
    
    return
  
  def save(self,filename='cluster.pickle'):
    """Same as Save()."""
    self.Save(filename=filename)
    return
  
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

def LoadModelCluster():
  
  """Reads the model cluster used in Johnstone et al. (2020)."""
  
  # Get directory where package is installed
  packageDir = misc._GetPackageDirectory()
  
  # Set filename including path
  filename = packageDir + "ModelDistribution.dat"
  
  # Read lines of file
  with open(filename) as f:
      content = f.readlines()
  
  # Set number of header lines
  nHeader = 1
  
  # Get number of steps in file
  nStars = len(content) - nHeader
  
  # Make arrays
  Mstar = np.zeros(nStars)
  Omega = np.zeros(nStars)
  
  # Loop over lines and fill arrays
  iStar = 0
  for line in content[nHeader:len(content)]:
    data = line.split()
    Mstar[iStar] = data[0]
    Omega[iStar] = data[1]
    iStar += 1
  
  return Mstar , Omega

#====================================================================================================================
