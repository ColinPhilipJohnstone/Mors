
"""Module for loading the stellar evolution tracks and retrieving basic stellar properties."""

# Imports for standard stuff needed here
import numpy as np
import pickle
import os
import math
import sys
import copy
import inspect

# Imports for Mors modules
import Mors.miscellaneous as misc

#====================================================================================================================

# Parameters for stellar evolution models

# Directory for stellar evolution models
starEvoDirDefault = os.getenv('STELLARMODELS')

# Set which set of models to use, i.e. which X, Z, and A values (set to X0p70952_Z0p01631_A1p875 for closest to solar)
evoModelsDefault = "X0p70952_Z0p01631_A1p875"

# Initially set default model data to None so we know it has not been set
ModelDataDefault = None

# Number of decimal places to represent masses and ages when doing interpolations
nDecValue = 7

#====================================================================================================================

class StarEvo:
  
  """
  Class to hold complete information about stellar evolution models.
  
  An instance of this class will hold evolutionary tracks for all masses and all quantities of interest. For
  now, only models from Spada et al. (2013) can be loaded. The main set of information held in this class is
  in ModelData, which is a dictionary of dictionaries. Each dictionary in this dictionary corresponds to a 
  single stellar mass and holds evolutionary tracks for a bunch of important quantities.

  Attributes
  ------------
  ModelData : dict
      This is a dictionary of dictionaries, holding evolution tracks for each stellar mass.
  
  Methods
  ------------
  
  """
  
  #---------------------------------------------------------------------------------------
  
  def __init__(self,starEvoDir=starEvoDirDefault,evoModels=evoModelsDefault):
    
    """Initialises instance of StarEvo class."""
    
    # If starEvoDir and evoModels are None, use defaults
    if starEvoDir is None:
      starEvoDir = starEvoDirDefault
    if evoModels is None:
      evoModels = evoModelsDefault
    
    # Load the models
    self.ModelData = _LoadModels( starEvoDir=starEvoDir , evoModels=evoModels )
    
    return
  
  #---------------------------------------------------------------------------------------
  
  def LoadTrack(self,Mstar,ClearData=False):
    
    """
    Takes stellar mass, loads evolutionary track for a specific mass into the model data.
    
    This can be useful if the user wants to load a specific evolutionary track into a model data for a specific mass
    since this will make getting values for this mass faster since no interpolation between mass bins will be needed.
      
    Parameters
    ----------
    Mstar : float
        Mass of star in Msun.
    ClearData : bool , optional
        If set to True, all other stellar masses will be removed from the ModelData dictionary.
    
    Returns
    ----------
    None
        None
    
    """
    
    self.ModelData = LoadTrack(Mstar,ModelData=self.ModelData,ClearData=ClearData)
    
    return
   
  #---------------------------------------------------------------------------------------
  
  def Value(Mstar,Age,ParamString):
    
    """
    Takes stellar mass, age, and a parameter string, returns values corresponding to named parameter.
    
    The set of models should have already been loaded. With this function, the user can ask for a value
    of one of the parameters for a specific stellar mass and age. All three of these can be input as multiple
    values and an array of values will be returned if this is the case. For example, if MstarIn is input as 
    a 1D array, the function will return a 1D array giving the value for each of these masses. ParamString
    can be input as a list of strings and values for each parameter in that list will be returned in a 1D 
    array. If two are given as arrays or lists, then a 2D array will be returned. If all three then a 3D 
    array with dimensions len(MstarIn)xlen(AgeIn)xlen(ParamString) will be returned.
    
    Parameters
    ----------
    MstarIn : float or int or numpy.ndarray
        Mass of star in Msun.
    AgeIn : float or int or numpy.ndarray
        Age in Myr.
    ParamString : str
        String holding name of parameter to get value for.
    
    Returns
    ----------
    value : float or numpy.ndarray
        Value of parameter at this mass and age. 
    
    """
    
    # Call Value() outside this class to get the value with this model
    value = Value(Mstar,Age,ParamString,ModelData=self.ModelData)
    
    return
  
  #---------------------------------------------------------------------------------------
  # The following functions are for individual parameters that can be called
  
  def Rstar(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns stellar radius."""
    return Value( Mstar , Age , 'Rstar' , ModelData=self.ModelData )
  
  def Lbol(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns bolometric luminosity."""
    return Value( Mstar , Age , 'Lbol' , ModelData=self.ModelData )
  
  def Teff(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns effective temperature."""
    return Value( Mstar , Age , 'Teff' , ModelData=self.ModelData )
  
  def Itotal(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns total moment of interia."""
    return Value( Mstar , Age , 'Itotal' , ModelData=self.ModelData )
  
  def Icore(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns core moment of interia."""
    return Value( Mstar , Age , 'Icore' , ModelData=self.ModelData )
  
  def Ienv(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns envelope moment of interia."""
    return Value( Mstar , Age , 'Ienv' , ModelData=self.ModelData )
  
  def Mcore(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns core mass."""
    return Value( Mstar , Age , 'Mcore' , ModelData=self.ModelData )
  
  def Menv(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns envelope mass."""
    return Value( Mstar , Age , 'Menv' , ModelData=self.ModelData )
  
  def Rcore(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns core radius."""
    return Value( Mstar , Age , 'Rcore' , ModelData=self.ModelData )
  
  def tauConv(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns convective turnover time."""
    return Value( Mstar , Age , 'tauConv' , ModelData=self.ModelData )
  
  def dItotaldt(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns rate of change of total moment of inertia."""
    return Value( Mstar , Age , 'dItotaldt' , ModelData=self.ModelData )
  
  def dIcoredt(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns rate of change of core moment of inertia."""
    return Value( Mstar , Age , 'dIcoredt' , ModelData=self.ModelData )
  
  def dIenvdt(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns rate of change of envelope moment of inertia."""
    return Value( Mstar , Age , 'dIenvdt' , ModelData=self.ModelData )
  
  def dMcoredt(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns rate of change of core mass."""
    return Value( Mstar , Age , 'dMcoredt' , ModelData=self.ModelData )
  
  def dRcoredt(self,Mstar,Age,ModelData=ModelDataDefault):
    """Takes mass and age, returns rate of change of core radius."""
    return Value( Mstar , Age , 'dRcoredt' , ModelData=self.ModelData )
  
  #---------------------------------------------------------------------------------------

#====================================================================================================================

def _LoadModels(starEvoDir=starEvoDirDefault,evoModels=evoModelsDefault):
  
  """Loads evolutionary tracks as a grid of parameters at each mass and age."""
  
  # Check if should compile new grid of evolutionary models or load previous grid
  if _shouldCompileNew(starEvoDir,evoModels):
    ModelData = _CompileNewGrid(starEvoDir,evoModels)
  else:
    ModelData = _LoadSavedGrid(starEvoDir,evoModels)
  
  return ModelData

#====================================================================================================================

def _shouldCompileNew(starEvoDir,evoModels):
  
  """Takes directory for stellar evo models, returns if new grid needs to be compiled."""
  
  # Initially assume need to compile new
  compileNew = True
  
  # Check if previously compiled models already exist
  if ( os.path.isfile(starEvoDir+evoModels+".pickle") ):
    compileNew = False
  
  return compileNew

#====================================================================================================================

def _CompileNewGrid(starEvoDir,evoModels):
  
  """Compiles new dictionary of stellar evo models."""
  
  # List all masses to load
  MstarAll = np.array([ 0.1 , 0.15 , 0.2 , 0.25 , 0.3 , 0.35 , 0.4 , 0.45 , 0.5 , 0.55 , 0.6 , 0.65 , 0.7 , 0.75 , 0.8 , 0.85 , 0.9 , 0.95 , 1.0 , 1.05 , 1.1 , 1.15 , 1.2 , 1.25 ])
  
  # Mass components of the filenames for these masses
  MstarFilenameMiddle = [ '0p10' , '0p15' , '0p20' , '0p25' , '0p30' , '0p35' , '0p40' , '0p45' , '0p50' , '0p55' , '0p60' , '0p65' , '0p70' , '0p75' , '0p80' , '0p85' , '0p90' , '0p95' , '1p00' , '1p05' , '1p10' , '1p15' , '1p20' , '1p25' ]
  
  # Start empty dictionary
  ModelData = {}
  
  # Add stellar mass array 
  ModelData['MstarAll'] = MstarAll
  
  # Add a list of strings holding each of the parameters
  ModelData['ParamsAll'] = [ 'Mstar' , 'Age' , 'Rstar' , 'Lbol' , 'Teff' , 'Itotal' , 'Icore' , 'Ienv' , 'Mcore' , 'Menv' , 'Rcore' , 
                            'tauConv' , 'dItotaldt' , 'dIcoredt' , 'dIenvdt' , 'dMcoredt' , 'dRcoredt' ]
  
  # Loop over masses and add each one to dictionary
  for iMstar in range(0,len(MstarAll)):
    ModelData[MstarAll[iMstar]] = _ReadEvolutionTrack( starEvoDir , evoModels , MstarAll[iMstar] , MstarFilenameMiddle[iMstar] )
  
  # Save compiled models
  with open(starEvoDir+evoModels+".pickle",'wb') as f:
    pickle.dump(ModelData,f)
    
  
  return ModelData

#====================================================================================================================

def _ReadEvolutionTrack(starEvoDir,evoModels,Mstar,MstarFilenameMiddle):
  
  
  # This function loads the original stellar evolution models from Spada et al. (2013) for an individual 
  # mass bin and puts them into a dictionary.
    
  # Set strings for starting and ending of filenames
  filename_prefix = starEvoDir + evoModels + "/M"
  filename_postfix1 = "_" + evoModels + ".track1"
  filename_postfix2 = "_" + evoModels + ".track2"
  
  # Get names of data files holding evo models
  filename1 = filename_prefix + MstarFilenameMiddle + filename_postfix1
  filename2 = filename_prefix + MstarFilenameMiddle + filename_postfix2
  
  # Read contents of files
  with open(filename1,'r') as f:
    content1 = f.readlines()
  with open(filename2,'r') as f:
    content2 = f.readlines()
  
  # Number of header lines
  nHeader = 1
  
  # Remove any duplicate ages
  content1 = _RemoveDuplicateAges(nHeader,content1)
  content2 = _RemoveDuplicateAges(nHeader,content2)
  
  # Make sure they are the same lengths
  if not ( len(content1) == len(content2) ):
    misc._PrintErrorKill("two evo files not same length")
  
  # Get number of age bins
  nAge = len(content1) - nHeader
  
  # Make arrays to hold desired quantities
  Age = np.zeros(nAge)
  Rstar = np.zeros(nAge)
  Lbol = np.zeros(nAge)
  Teff = np.zeros(nAge)
  Itotal = np.zeros(nAge)
  Icore = np.zeros(nAge)
  Ienv = np.zeros(nAge)
  Mcore = np.zeros(nAge)
  Menv = np.zeros(nAge)
  Rcore = np.zeros(nAge)
  tauConv = np.zeros(nAge)
  dItotaldt = np.zeros(nAge)
  dIcoredt = np.zeros(nAge)
  dIenvdt = np.zeros(nAge)
  dMcoredt = np.zeros(nAge)
  dRcoredt = np.zeros(nAge)
  
  # Read quantities from first file into the arrays (has age, Lbol, and Rstar)
  iAge = 0
  for line in content1[nHeader:len(content1)]:
    
    # Split data into list
    data = line.split()
    
    # Read data 
    Age[iAge] = data[0]
    Lbol[iAge] = data[4]
    Teff[iAge] = data[7]
    Rstar[iAge] = data[5]
    
    # Update index
    iAge += 1
    
  # Read quantities from second file into the arrays (has Itotal, Ienv, Menv, Rcore, and tauConv)
  iAge = 0
  for line in content2[nHeader:len(content2)]:
    
    # Split data into list
    data = line.split()
    
    # Read data 
    Itotal[iAge] = data[7]
    Ienv[iAge] = data[6]
    Menv[iAge] = data[2]
    Rcore[iAge] = data[3]
    tauConv[iAge] = data[4]
    
    # Update i
    iAge += 1
    
  # Change some units and convert from log to lin
  Age[:] = Age[:] * 1000.0 # convert from Gyr into Myr
  Rstar[:] = 10.0**Rstar[:] # convert into linear since Spada has in log
  Lbol[:] = 10.0**Lbol[:] # also this one
  Teff[:] = 10.0**Teff[:] # also this one
    
  # Get some core quantities
  Icore[:] = Itotal[:] - Ienv[:]
  Mcore[:] = ( 1.0 - Menv[:] ) * Mstar
  Rcore[:] = Rcore[:] * Rstar[:]
    
  # Get gradients
  dItotaldt[:] = _CalculateGradient( Age , Itotal )
  dIcoredt[:] = _CalculateGradient( Age , Icore )
  dIenvdt[:] = _CalculateGradient( Age , Ienv )
  dMcoredt[:] = _CalculateGradient( Age , Mcore )
  dRcoredt[:] = _CalculateGradient( Age , Rcore )
  
  # Round mass and ages to nDecValue decimal places as specified at top of file
  Mstar = round( Mstar , nDecValue )
  Age = np.round( Age , decimals=nDecValue )
  
  # Add all quantities to the dictionary
  Data = {}
  Data['Mstar'] = Mstar
  Data['Age'] = Age
  Data['Rstar'] = Rstar
  Data['Lbol'] = Lbol
  Data['Teff'] = Teff
  Data['Itotal'] = Itotal
  Data['Icore'] = Icore
  Data['Ienv'] = Ienv
  Data['Mcore'] = Mcore
  Data['Menv'] = Menv
  Data['Rcore'] = Rcore
  Data['tauConv'] = tauConv
  Data['dItotaldt'] = dItotaldt
  Data['dIcoredt'] = dIcoredt
  Data['dIenvdt'] = dIenvdt
  Data['dMcoredt'] = dMcoredt
  Data['dRcoredt'] = dRcoredt
  #Data[''] = 
  
  return Data

#====================================================================================================================

def _RemoveDuplicateAges(nHeader,content):
  
  """Takes lines from evo file, removes duplicate ages."""
  
  # This is necessary because one of the Spada files has an age repeated a few times.
  
  # Start list
  contentNew = []
  
  # Add header lines
  for line in content[0:nHeader]:
    contentNew.append(line)
  
  # Get age of first line
  data = content[nHeader].split()
  AgeLast = float(data[0])
  
  # Loop over remaning lines and decide which to add (all but duplicate ages)
  for line in content[nHeader+1:len(content)]:
    
    # Get age
    data = line.split()
    Age = float(data[0])
    
    # Add line if not same age (within some tolerance)
    if not ( abs(Age/AgeLast-1.0) < 1.0e-5 ):
      contentNew.append(line)
    
    # Update AgeLast
    AgeLast = Age
  
  return contentNew

#====================================================================================================================

def _CalculateGradient(Age,X):
  
  """Takes an evolutionary track for quantity, returns evolution of rate of change of quantity."""
  
  # Number of age bins
  nAge = len(Age)
  
  # Make array
  dXdt = np.zeros(nAge)
  
  # Set first and last elements
  dXdt[0] = ( X[1] - X[0] ) / ( Age[1] - Age[0] )
  dXdt[nAge-1] = ( X[nAge-1] - X[nAge-2] ) / ( Age[nAge-1] - Age[nAge-2] )
  
  # Do rest of ages
  dXdt[1:nAge-1] = ( X[2:nAge] - X[0:nAge-2] ) / ( Age[2:nAge] - Age[0:nAge-2] )
  
  dAge = np.zeros(nAge)
  dAge[0] = Age[1] - Age[0]
  dAge[1:nAge-1] = Age[2:nAge] - Age[0:nAge-2]
  dAge[nAge-1] = Age[nAge-1] - Age[nAge-2]
  
  return dXdt

#====================================================================================================================

def _LoadSavedGrid(starEvoDir,evoModels):
  
  """Takes filename for stellar evo model, returns grid of models."""
  
  # Simply load data
  with open(starEvoDir+evoModels+".pickle",'rb') as f:
    ModelData = pickle.load(f)
    
  return ModelData

#====================================================================================================================

def LoadTrack(Mstar,ModelData=None,ClearData=False):
  global ModelDataDefault
  
  """
  Takes stellar mass, loads evolutionary track for a specific mass into the model data.
  
  This can be useful if the user wants to load a specific evolutionary track into a model data for a specific mass
  since this will make getting values for this mass faster since no interpolation between mass bins will be needed.
  The user might also want to set ClearData to True if doing the calculation only for a single stellar mass since 
  this will cause the function to remove data for all other masses, which saves memory, though usually this will 
  only be necessary if the user is going to have data for a very large number of stars loaded simultaneously. Note
  that this does not change ModelDataDefault.
    
  Parameters
  ----------
  Mstar : float
      Mass of star in Msun.
  ModelData : dict , optional
      Dictionary of dictionaries holding set of stellar evolution models.
  ClearData : bool , optional
      If set to True, all other stellar masses will be removed from the ModelData dictionary.
  
  Returns
  ----------
  None
      None
  
  """
  
  # Use default model data if nothing was specified
  if ModelData is None:
    ModelData = ModelDataDefault

  # Check if ModelData is not None, otherwise need to load defaults
  if ModelData is None:
    ModelData = _LoadDefaultModelData()
  
  # Only do something if Mstar is not already in ModelData
  if not ( Mstar in ModelData ):
    
    # Get evo tracks for this mass and add it to the dictionary
    Data = _LoadTrack(Mstar,ModelData)
    
    # Add extra term to ModelData
    ModelData[Mstar] = Data
    
    # Add this mass to MstarAll in ModelData
    ModelData['MstarAll'] = np.sort( np.append( ModelData['MstarAll'] , Mstar ) ) 
    
    # If ModelDataDefault is not None, add this also to that
    if not ModelDataDefault is None:
      ModelDataDefault[Mstar] = Data
  
  # If ClearData is set, remove all other tracks from return dictionary, but not ModelDataDefault
  if ClearData:
    
    # Do this by making a new dictionary and copying stuff in
    
    # Make deepcopy of dictionary
    ModelDataOld = copy.deepcopy(ModelData)
    
    # Make new dictionary
    ModelData = {}
    
    # Make new MstarAll array
    ModelData['MstarAll'] = np.array([Mstar])
    
    # Add properties
    ModelData['ParamsAll'] = ModelDataOld['ParamsAll']
    
    # Add tracks for this mass bin
    ModelData[Mstar] = ModelDataOld[Mstar]
  
  return ModelData

#====================================================================================================================

def _LoadTrack(Mstar,ModelData):
  
  """Takes stellar mass and model data dictionary, returns dictionary with track for this mass."""
  
  # Round mass to nDecValue decimal places specified at top of this file
  Mstar = round( Mstar , nDecValue )
  
  # Make sure within mass limit
  _CheckMassLimit( ModelData['MstarAll'] , Mstar )
  
  # Get nearest mass bin below
  iMin = _getIndexLT( ModelData['MstarAll'] , Mstar )
  MstarMin = ModelData['MstarAll'][iMin]
  
  # Get nearest mass bin above
  iMax = _getIndexGT( ModelData['MstarAll'] , Mstar )
  MstarMax = ModelData['MstarAll'][iMax]
  
  # Get dictionaries for masses below and above Mstar
  DataMin = ModelData[MstarMin]
  DataMax = ModelData[MstarMax]
  
  # Get min and max ages to include
  #   for min, this is the maximum starting age of the two tracks
  #   for max, this is the minimum ending age of the two tracks
  AgeMin = np.max( [ DataMin['Age'][0] , DataMax['Age'][0] ] )
  AgeMax = np.min( [ DataMin['Age'][-1] , DataMax['Age'][-1] ] )
  
  # Number of age bins to use, just make this minimum of the number in the two tracks
  nAge = np.min( [ len(DataMin['Age']) , len(DataMax['Age']) ] )
  
  # Make age array, log spacing
  Age = np.logspace( np.log10(AgeMin) , np.log10(AgeMax) , nAge )
  
  # Round ages to nDecValue decimal placed specified at top of this file
  Age = np.round( Age , decimals=nDecValue )
  
  # Start the output dictionary
  Data = {}
  
  # Add mass and age array to dictionary
  Data['Mstar'] = Mstar
  Data['Age'] = Age
  
  # Loop over parameters and get tracks for each (do parameters in DataMin)
  for param in DataMin:
    
    # Skip parameter if already added (so don't do Mstar and Age again)
    if ( param in Data ):
      continue
    
    # Get track
    track = Value( Mstar , Age , param , ModelData=ModelData )
    
    # Add to dictionary
    Data[param] = track
  
  return Data

#====================================================================================================================

def _LoadDefaultModelData():
  global ModelDataDefault
  
  """Loads model data for the default model to be used if no StarEvo class object is used."""
  
  # Loading a default model is generally not necessary, but it would be useful if the user doesn't want to mess
  # around with creating a StarEvo class object and instead just wants to quickly get a few values.
  
  # Load the data 
  ModelDataDefault = _LoadModels()
  
  return ModelDataDefault

#====================================================================================================================

def Value(MstarIn,AgeIn,ParamString,ModelData=ModelDataDefault):
  
  """
  Takes stellar mass, age, and a parameter string, returns values corresponding to named parameter.
  
  The set of models should have already been loaded. With this function, the user can ask for a value
  of one of the parameters for a specific stellar mass and age. All three of these can be input as multiple
  values and an array of values will be returned if this is the case. For example, if MstarIn is input as 
  a 1D array, the function will return a 1D array giving the value for each of these masses. ParamString
  can be input as a list of strings and values for each parameter in that list will be returned in a 1D 
  array. If two are given as arrays or lists, then a 2D array will be returned. If all three then a 3D 
  array with dimensions len(MstarIn)xlen(AgeIn)xlen(ParamString) will be returned.
  
  Parameters
  ----------
  MstarIn : float or int or numpy.ndarray
      Mass of star in Msun.
  AgeIn : float or int or numpy.ndarray
      Age in Myr.
  ParamString : str
      String holding name of parameter to get value for.
  ModelData : dict , optional
      Dictionary of dictionaries holding set of stellar evolution models.
  
  Returns
  ----------
  value : float or numpy.ndarray
      Value of parameter at this mass and age. 
  
  """
  
  # Check if ModelData is not None, otherwise need to load defaults
  if ModelData is None:
    ModelData = _LoadDefaultModelData()
  
  # There are now eight scenarios here
  #   1. Mstar and Age are both floats, so return float.
  #   2. Mstar is float and Age is numpy.ndarray, so return numpy.ndarray of length len(Age).
  #   3. Mstar is numpy.ndarray and Age is float, so return numpy.ndarray of length len(Age)
  #   4. Mstar and Age are both numpy.ndarray, so return numpy.ndarray of shape ( len(Mstar) , len(Age) ).
  # Note: could use type(X).__module__ == 'numpy' to test for numpy but it doesn't work if 
  
  # Get Mstar and Age in correct types
  Mstar = misc._convertFloatArray(MstarIn)
  Age = misc._convertFloatArray(AgeIn)
  
  # Make sure Mstar and Age are correct types
  if Mstar is None:
    misc._PrintErrorKill("argument Mstar has invalid type")
  if Age is None:
    misc._PrintErrorKill("argument Age has invalid type")
  
  # Find if scenario 1 (most likely)
  if ( ( type(Mstar) == float ) and ( type(Age) == float ) and ( type(ParamString) == str ) ):
        
    # Scenario 1 so just get value
    value = _ValueSingle( Mstar , Age , ParamString , ModelData=ModelData )
    
  else:
    
    # In this case, the output will be an array with up to three dimensions, so first make 3D array
    # and then remove dimensions with only one element
    
    # Get number of elements in all three directions
    try:
      nMstar = len(Mstar)
    except:
      nMstar = 1
    
    try:
      nAge = len(Age)
    except:
      nAge = 1
    
    if ( type(ParamString) == list ):
      nParam = len(ParamString)
    else:
      nParam = 1
    
    # Make array
    value = np.zeros((nMstar,nAge,nParam))
    
    # Loop over values 
    for iMstar in range(0,nMstar):
      for iAge in range(0,nAge):
        for iParam in range(0,nParam):
          
          # Get values of Mstar, Age, and ParamString
          if ( type(Mstar) == float ):
            MstarValue = Mstar
          else:
            MstarValue = Mstar[iMstar]
            
          if ( type(Age) == float ):
            AgeValue = Age
          else:
            AgeValue = Age[iAge]
            
          if ( type(ParamString) == str ):
            ParamStringValue = ParamString
          else:
            ParamStringValue = ParamString[iParam]
          
          # Now get the value
          value[iMstar,iAge,iParam] = _ValueSingle( MstarValue , AgeValue , ParamStringValue , ModelData=ModelData )
    
    # Now get rid of unwanted dimensions
    value = np.squeeze(value)
    
  return value

#====================================================================================================================

def _ValueSingle(Mstar,Age,ParamString,ModelData=ModelDataDefault):
  
  """Takes stellar mass and age and parameter string, returns value of that parameter at that mass and age."""
  
  ## Make sure Mstar and age are floats
  #if not ( type(Mstar) == float ):
    #misc._PrintErrorKill("argument Mstar must be float in call to Value")
  #if not ( type(Age) == float ):
    #misc._PrintErrorKill("argument Age must be float in call to Value")
  
  # Make sure ParamString is indeed a string
  if not ( type(ParamString) == str ):
    misc._PrintErrorKill("argument ParamString must be string in call to Value")
  
  # Make sure ParamString corresponds to a valid parameter
  if not ( ParamString in ModelData['ParamsAll'] ):
    misc._PrintErrorKill("parameter '"+ParamString+"' is not valid")
  
  # Check if an evolutionary track for this exact stellar mass is present
  if ( Mstar in ModelData ):
    
    # In this case, an evo track for this exact stellar mass is already loaded and all is needed is an interpolation
    # but first check the age is within the correct age limit
    _CheckAgeLimit( ModelData[Mstar]['Age'] , Age )
    value = _Interpolate1D( ModelData[Mstar]['Age'] , ModelData[Mstar][ParamString] , Age )
    
  else:
    
    # In this case, no evo track for this specific mass is loaded and so a 2D interpolation between values from two
    # separate mass bins is needed
    # Note, it is not assumed that the mass bins are loaded into the dictionary in ascending order
    
    # Make sure within mass limit
    _CheckMassLimit( ModelData['MstarAll'] , Mstar )
    
    # Get nearest mass bin below
    iMin = _getIndexLT( ModelData['MstarAll'] , Mstar )
    MstarMin = ModelData['MstarAll'][iMin]
    
    # Get nearest mass bin above
    iMax = _getIndexGT( ModelData['MstarAll'] , Mstar )
    MstarMax = ModelData['MstarAll'][iMax]
    
    # Check age ranges for both (max first since the upper limit is likely lower)
    _CheckAgeLimit( ModelData[MstarMax]['Age'] , Age )
    _CheckAgeLimit( ModelData[MstarMin]['Age'] , Age )
    
    # Do interpolation
    value = _Interpolate2D( MstarMin , MstarMax , Mstar ,
                            ModelData[MstarMin]['Age'] , ModelData[MstarMax]['Age'] , Age ,
                            ModelData[MstarMin][ParamString] , ModelData[MstarMax][ParamString] )
    
  return value

#====================================================================================================================

def _CheckAgeLimit(AgeArray,Age):
  
  """Takes age track and an age, outputs error and stops code if age is not within limits."""
  
  # Round the age to decimal places determined by nDecValue at top of file
  Age = round( Age , nDecValue )
  
  # Do check
  if not ( AgeArray[0] <= Age <= AgeArray[-1] ):
    misc._PrintErrorKill("input age "+str(Age)+" is not within limits of "+str(AgeArray[0])+" to "+str(AgeArray[-1]))
  
  return


#====================================================================================================================

def _CheckMassLimit(MstarArray,Mstar):
  
  """Takes array of masses and an age, outputs error and stops code if mass is not within limits."""
  
  if not ( np.min(MstarArray) <= Mstar <= np.max(MstarArray) ):
    misc._PrintErrorKill("input stellar mass "+str(Mstar)+" is not within limits of "+str(np.min(MstarArray))+" to "+str(np.max(MstarArray)))
  
  return

#====================================================================================================================

def _Interpolate2D(Z1,Z2,Z,Xarray1,Xarray2,X,Yarray1,Yarray2):
  
  """Takes two sets of 1D arrays for corresponding X and Y values, returns interpolated Y value corresponding to input X."""
  
  # Round the input values to decimal places determined by nDecValue at top of file
  Z = round( Z , nDecValue )
  
  # Do interpolations to get Y at X for both tracks
  Y1 = _Interpolate1D( Xarray1 , Yarray1 , X )
  Y2 = _Interpolate1D( Xarray2 , Yarray2 , X )
  
  # Do linear interpolation between Y1 and Y2 in Z direction
  mInterp = ( Y2 - Y1 ) / ( Z2 - Z1 )
  cInterp = Y1 - mInterp * Z1
  Y = mInterp * Z + cInterp
  
  return Y

#====================================================================================================================

def _Interpolate1D(Xarray,Yarray,X):
  
  """Takes 1D arrays for corresponding X and Y values, returns interpolated Y value corresponding to input X."""
  
  # Note that it is assumed here that Xarray is in ascending order and this won't work if it is not
  
  # Round the input values to decimal places determined by nDecValue at top of file
  X = round( X , nDecValue )  
  
  # Make sure X is in limits
  if not ( Xarray[0] <= X <= Xarray[-1] ):
    misc._PrintErrorKill("input value "+str(X)+" is not within limits of "+str(Xarray[0])+" to "+str(Xarray[-1]))
  
  # Get index of X closest to but smaller than value
  iMin = _getIndexLTordered(Xarray,X)
  iMax = iMin + 1
  
  # Check if X is an element in Xarray
  if ( Xarray[iMin] == X ):
    return Yarray[iMin]
  elif ( Xarray[iMax] == X ):
    return Yarray[iMax]
  else:
    
    # Do linear interpolation
    mInterp = ( Yarray[iMax] - Yarray[iMin] ) / ( Xarray[iMax] - Xarray[iMin] )
    cInterp = Yarray[iMin] - mInterp * Xarray[iMin]
    Y = mInterp * X + cInterp
  
  return Y

#====================================================================================================================

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

#====================================================================================================================

def _getIndexLT(Xarray,X):
  
  """Takes array of values and a value, returns index of closest element in array smaller than value."""
  
  # It is assumed that X is greater than the min of Xarray, but check this
  if ( X < np.min(Xarray) ):
    misc._PrintErrorKill("X is less than minimum of Xarray")
  
  # Get smaller value by taking difference of X and all Xarray elements
  # then removing all values smaller than 0.0 (i.e. those larger than X)
  deltaX = X - Xarray
  deltaX[np.where(deltaX<0.0)] = np.max(deltaX)
  
  # Get index
  index = np.argmin(deltaX)
  
  return index

#====================================================================================================================

def _getIndexGT(Xarray,X):
  
  """Takes array of values and a value, returns index of closest element in array larger than value."""
  
  # It is assumed that X is less than the max of Xarray, but check this
  if ( X > np.max(Xarray) ):
    misc._PrintErrorKill("X is more than maximum of Xarray")
  
  # Get smaller value by taking difference of X and all Xarray elements
  # then removing all values smaller than 0.0 (i.e. those larger than X)
  deltaX = Xarray - X
  deltaX[np.where(deltaX<0.0)] = np.max(deltaX)
  
  # Get index
  index = np.argmin(deltaX)
  
  return index

#====================================================================================================================

# The following functions are for individual parameters that can be called

def Rstar(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns stellar radius."""
  return Value( Mstar , Age , 'Rstar' , ModelData=ModelData )

def Lbol(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns bolometric luminosity."""
  return Value( Mstar , Age , 'Lbol' , ModelData=ModelData )

def Teff(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns effective temperature."""
  return Value( Mstar , Age , 'Teff' , ModelData=ModelData )

def Itotal(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns total moment of interia."""
  return Value( Mstar , Age , 'Itotal' , ModelData=ModelData )

def Icore(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns core moment of interia."""
  return Value( Mstar , Age , 'Icore' , ModelData=ModelData )

def Ienv(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns envelope moment of interia."""
  return Value( Mstar , Age , 'Ienv' , ModelData=ModelData )

def Mcore(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns core mass."""
  return Value( Mstar , Age , 'Mcore' , ModelData=ModelData )

def Menv(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns envelope mass."""
  return Value( Mstar , Age , 'Menv' , ModelData=ModelData )

def Rcore(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns core radius."""
  return Value( Mstar , Age , 'Rcore' , ModelData=ModelData )

def tauConv(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns convective turnover time."""
  return Value( Mstar , Age , 'tauConv' , ModelData=ModelData )

def dItotaldt(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns rate of change of total moment of inertia."""
  return Value( Mstar , Age , 'dItotaldt' , ModelData=ModelData )

def dIcoredt(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns rate of change of core moment of inertia."""
  return Value( Mstar , Age , 'dIcoredt' , ModelData=ModelData )

def dIenvdt(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns rate of change of envelope moment of inertia."""
  return Value( Mstar , Age , 'dIenvdt' , ModelData=ModelData )

def dMcoredt(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns rate of change of core mass."""
  return Value( Mstar , Age , 'dMcoredt' , ModelData=ModelData )

def dRcoredt(Mstar,Age,ModelData=ModelDataDefault):
  """Takes mass and age, returns rate of change of core radius."""
  return Value( Mstar , Age , 'dRcoredt' , ModelData=ModelData )

#====================================================================================================================
