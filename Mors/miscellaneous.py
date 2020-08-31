
"""Module for holding miscellaneous functions that can be used by many other modules."""

import inspect
import sys
import numpy as np

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
