
"""Module for holding miscellaneous functions that can be used by many other modules."""

import inspect
import sys
import numpy as np
import pickle

#==================================================================================================================

def _PrintErrorKill(errorString):
  
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
  sys.exit()
  
  return

#==================================================================================================================

def _PrintError(errorString):
  
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

#==================================================================================================================

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

#==================================================================================================================

def _convertFloatArray(Xin):
  
  """Takes a value and returns either as float or numpy.ndarray of floats."""
  
  # If input is float or int, return float version
  if ( ( type(Xin) == float ) or ( type(Xin) == int ) ):
    return float(Xin)
  
  # If input is a list, numpy array version (will only work if elements can be converted into floats)
  if ( type(Xin) == list ):
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

#==================================================================================================================

def Load(filename):
  
  """Takes filename of saved star or cluster and loads object."""
  
  with open(filename,'rb') as f:
    obj = pickle.load(f)
  
  return obj

#==================================================================================================================

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

#====================================================================================================================

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

#====================================================================================================================


