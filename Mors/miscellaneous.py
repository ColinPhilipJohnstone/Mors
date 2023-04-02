
"""Module for holding miscellaneous functions that can be used by many other modules."""

# Imports for standard stuff needed here
import inspect
import sys
import numpy as np
import pickle
import copy

# Imports for Mors modules
import Mors.constants as const

def _PrintErrorKill(errorString):
  """Takes string with error message, prints message to screen and end code execution."""
  
  # Get info about the calling function
  frame = inspect.stack()[1]
  filename = frame.filename
  function = frame.function
  
  # From filename, just get relative to kompot/ or kvass/ directory
  temp = filename.split('/')
  filename = temp[-1]
  
  # Print error to screen
  print( "Error in "+function+"() in file "+filename+": "+errorString )
  
  # Stop the code
  sys.tracebacklimit = 1
  raise ValueError()
  
  return

def _PrintError(errorString):
  """Takes string with error message, prints message to screen."""

  # Get info about the calling function
  frame = inspect.stack()[1]
  filename = frame.filename
  function = frame.function
  
  # From filename, just get relative to kompot/ directory
  temp = filename.split('/')
  filename = temp[-1]
  
  # Print error to screen
  print( "Error in "+function+"() in file "+filename+": "+errorString )
  
  return

def _GetPackageDirectory():
  """Gets path to main directory where the code is installed."""
  
  # Get information about function
  frame = inspect.stack()[0]
  
  # Get name of this file
  filename = frame.filename
  
  # Split by / and then reconstruct string removing anything after last /
  temp = filename.split("/")
  path = ''
  for item in temp[0:-1]:
    path += item + '/'
  
  return path

def _convertFloatArray(Xin):
  """Takes a value and returns either as float or numpy.ndarray of floats."""
  
  # If input is float or int, return float version
  if isinstance(Xin,(float,int)):
    return float(Xin)
  
  # If input is a list, numpy array version (will only work if elements can be converted into floats)
  if isinstance(Xin,list):
    X = np.zeros(len(Xin))
    for i in range(0,len(Xin)):
      X[i] = float(Xin[i])
    return X
  
  # If input is a numpy object then do stuff
  if ( type(Xin).__module__ == 'numpy' ):
    
    # First see if it can be converted into a float (i.e. a 1-element array)
    isNumber = True
    try:
      X = float(Xin)
    except:
      isNumber = False
    
    # If it is just a single value, return float version
    if ( isNumber ):
      return float(X)
    
    # Otherwise, turn it into an array
    X = np.zeros(len(Xin))
    for i in range(0,len(Xin)):
      X[i] = float(Xin[i])
    return X
  
  return

def Load(filename):
  """Takes filename of saved star or cluster and loads object."""
  
  # Load the object
  with open(filename,'rb') as f:
    obj = pickle.load(f)
  
  # Setup the quantity functions (which for some reason do not work on objects loaded from pickle)
  obj._setupQuantityFunctions()
  
  return obj

def ModelCluster():
  """Reads the model cluster used in Johnstone et al. (2020)."""
  
  # Get directory where package is installed
  packageDir = _GetPackageDirectory()
  
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

def ActivityLifetime(Age=None,Track=None,Threshold=None,AgeMax=None):
  """
  Takes evolutionary track for parameter, calculates when value drops below threshold.

  This function can be used to determine when a star's emission crosses a given threshold value for a few
  activity quantities. These are Lx, Fx, Rx, and FxHZ for X-rays, and similar values for EUV1, EUV2, EUV, 
  XUV, and Ly-alpha. If the star crosses the threshold (from above it to below it) multiple times, this 
  will find the final time it will cross the threshold. If the user wants to set a maximum age so that the 
  code only looks for crossings of the threshold below this age then this can be done using the AgeMax
  keyword argument. If the track is always below this threshold value then the function returns 0.0.
      
  Parameters
  ----------
  Age : numpy.ndarray
      Age array for evolutionary track.
  Track : numpy.ndarray
      Value array for evolutionary track.
  Threshold : float
      Value for threshold in units of quantity.
  AgeMax : float , optional
      End age of track to consider in Myr.
  
  Returns
  ----------
  AgeActive : float
      Activity lifetime in Myr.
  
  """
  
  # Make sure parameters set
  if Age is None:
    misc._PrintErrorKill("required argument Age not set")
  if Track is None:
    misc._PrintErrorKill("required argument Track not set")
  if Threshold is None:
    misc._PrintErrorKill("required argument Threshold not set")
  
  # Make sure Age and Track same length
  if not ( len(Age) == len(Track) ):
    misc._PrintErrorKill("Age and Track are different lengths")
  
  # Change Age and Track to only be ages below AgeMax, if AgeMax is set
  if not AgeMax is None:
    includeAges = np.where( Age <= AgeMax )
    Age = Age[includeAges]
    Track = Track[includeAges]
  
  # See if final value below threshold, otherwise return final age
  if ( Track[-1] > Threshold ):
    return Age[-1]
  
  # Start at end of array and loop backwards
  for iAge in range(len(Age)-1,0,-1):
    
    # Make sure next one back is not equal to threshold
    if ( Track[iAge-1] == Threshold ):
      return Age[iAge-1]
    
    # Check if this age bin is when it drops past the threshold 
    if ( ( Track[iAge] < Threshold ) and ( Track[iAge-1] > Threshold ) ):
      
      # Do interpoaltion to get when it crosses 
      # Assume here Age = m*Track + c
      mInterp = ( Age[iAge] - Age[iAge-1] ) / ( Track[iAge] - Track[iAge-1] )
      cInterp = Age[iAge] - mInterp * Track[iAge]
      return mInterp*Threshold + cInterp
      
  return 0.0

def IntegrateEmission(AgeMin=None,AgeMax=None,Age=None,Luminosity=None,aOrb=None):
  """
  Takes evolutionary track for parameter, calculates when value drops below threshold.
  
  This code can be used to integrate a star's luminosity between two ages. This can be applied to Lx, Leuv, or any
  other luminosity and the result is a total energy emitted in this time in erg. If the user also specifies an orbital 
  distance using the aOrb keyword argument, the code integrates the flux at this obital distance returns a fluence 
  in erg cm^-2.
      
  Parameters
  ----------
  AgeMin : float
      Start of time period to integrate in Myr.
  AgeMax : float
      End of time period to integrate in Myr.
  Age : numpy.ndarray
      Age array for evolutionary track.
  Luminosity : numpy.ndarray
      Value array for evolutionary track.
  aOrb : float , optional
      Orbital distance to get fluence at in AU.
  
  Returns
  ----------
  Energy : float
      Integrated luminosity of flux in erg or erg cm^-2.
  
  """
  
  # Make sure parameters set
  if AgeMin is None:
    misc._PrintErrorKill("required argument AgeMin not set")
  if AgeMax is None:
    misc._PrintErrorKill("required argument AgeMax not set")
  if Age is None:
    misc._PrintErrorKill("required argument Age not set")
  if Luminosity is None:
    misc._PrintErrorKill("required argument Luminosity not set")
  
  # Make sure Age and Luminosity same length
  if not ( len(Age) == len(Luminosity) ):
    misc._PrintErrorKill("Age and Luminosity are different lengths")
  
  # Make sure AgeMin and AgeMax is in track
  if not ( ( AgeMin >= Age[0] ) and ( AgeMin <= Age[-1] ) ):
    misc._PrintErrorKill("AgeMin not in range of evolutionary track")
  if not ( ( AgeMax >= Age[0] ) and ( AgeMax <= Age[-1] ) ):
    misc._PrintErrorKill("AgeMax not in range of evolutionary track")
  
  # Get track to actually integrate
  Track = copy.deepcopy(Luminosity)
  if not aOrb is None:
    Track *= 1.0 / ( 4.0 * const.Pi * (aOrb*const.AU)**2.0 )
  
  # Get index of first age bin with age above AgeMin
  indexMin = SE._getIndexGT(Age,AgeMin)
  
  # Get index of final age bin with age below AgeMax
  indexMax = SE._getIndexLT(Age,AgeMax)
  
  # Initially take no energy then add up energy from bins
  Energy = 0.0
  
  # Add energy from first bin if needed (simple assumption for luminosity constant ove age bin)
  if ( AgeMin < Age[indexMin] ):
    Energy += ( Age[indexMin] - AgeMin ) * Track[indexMin]
  
  # Intergate over all full bins (only if there is something to integrate)
  if not ( indexMin == indexMax ):
    for iAge in range(indexMin,indexMax):
      
      # Do integration using trapezoidal rule (not the chained version)
      Energy += 0.5 * ( Age[iAge+1] - Age[iAge] ) * ( Track[iAge] + Track[iAge+1] )
  
  # Add energy from final bin if needed (simple assumption for luminosity constant ove age bin)
  if ( AgeMax > Age[indexMax] ):
    Energy += ( AgeMax - Age[indexMax] ) * Track[indexMax]
  
  # So far, units of Age were in Myr when integrating, so include also Myr to s conversion
  Energy *= const.Myr
  
  return Energy

def _getIndexLTordered(Xarray,X):
  """Takes min to max ordered array of values and a value, returns index of closest element in array smaller than value."""
  
  # It is assumed that X is between the min and max of Xarray and it is assumed
  # that Xarray is in ascending order
    
  # Initial guesses for i
  i1 = 0
  i2 = len(Xarray)-1
  
  # Check if right value
  if ( Xarray[i1] == X ):
    return i1
  if ( Xarray[i2] == X ):
    return i2
  
  # Start iterating
  for iIter in range(0,len(Xarray)):
    
    # Get iMid
    iMid = int(0.5*(i1+i2))
    
    # Check if iMid is right value
    if ( ( Xarray[iMid] <= X ) and ( Xarray[iMid+1] > X ) ):
      return iMid
    
    # Work out if answer is between i1 and iMid or iMid and i2
    if ( Xarray[iMid] > X ):
      i2 = iMid
    else:
      i1 = iMid
  
  # It should not get here
  misc._PrintErrorKill("did not find index")
  
  return iMid

def _getIndexLT(Xarray,X):
  """Takes array of values and a value, returns index of closest element in array smaller than value."""
  
  # It is assumed that X is greater than the min of Xarray, but check this
  if ( X < np.min(Xarray) ):
    misc._PrintErrorKill("X is less than minimum of Xarray")
  
  # Get smaller value by taking difference of X and all Xarray elements
  # then removing all values smaller than 0.0 (i.e. those larger than X)
  deltaX = X - Xarray
  deltaX[np.where(deltaX<0.0)] = np.max(deltaX)*1.1+0.1
  
  # Get index
  index = np.argmin(deltaX)
  
  return index

def _getIndexGT(Xarray,X):
  """Takes array of values and a value, returns index of closest element in array larger than value."""
  
  # It is assumed that X is less than the max of Xarray, but check this
  if ( X > np.max(Xarray) ):
    misc._PrintErrorKill("X is more than maximum of Xarray")
  
  # Get smaller value by taking difference of X and all Xarray elements
  # then removing all values smaller than 0.0 (i.e. those larger than X)
  deltaX = Xarray - X
  deltaX[np.where(deltaX<0.0)] = np.max(deltaX)*1.1+0.1
  
  # Get index
  index = np.argmin(deltaX)
  
  return index
