
import StellarEvolution.miscellaneous as misc
#from stellarevo import LoadTracks

# Limits on various input values
MstarMin = 0.05
MstarMax = 1.5


#====================================================================================================================

class Star:
  
  """
  A class for star objects that hold all information about a star. 
  
  Attributes
  ------------
  
  Methods
  ------------
  
  """
  
  #---------------------------------------------------------------------------------------

  def __init__(self,Mstar=None,Omega0=None,percentile=None):
    
    """
    Initialises instance of star class
    """
    
    # Load evo tracks if Mstar is set
    if not Mstar is None:
      self.LoadEvoTracks(Mstar)
    
    
    return
  
  #---------------------------------------------------------------------------------------
  
  def LoadEvoTracks(self,Mstar=None):
    
    """Takes stellar mass, loads stellar evolution tracks."""
    
    # Make sure Mstar has good values
    _CheckMstar(Mstar)
    
    # Set Mstar
    self.Mstar = Mstar
    
    # Load evolutionary tracks for all quantities of interest in form of dictionary
    #self.EvoTracks = 
    
    
    return
  
  #---------------------------------------------------------------------------------------
  
#====================================================================================================================

def _CheckMstar(Mstar):
  
  # Make sure Mstar was set
  if Mstar is None:
    misc._PrintErrorKill("stellar mass not given")
  
  # Make sure it is a float
  if not ( type(Mstar) == float ):
    misc._PrintErrorKill("stellar mass must be given as float")
  
  # Make sure it above minimum
  if ( Mstar < MstarMin ):
    misc._PrintErrorKill( "stellar mass cannot be less than lower limit of "+str(MstarMin)+" Msun" )
  
  # Make sure it above minimum
  if ( Mstar > MstarMax ):
    misc._PrintErrorKill( "stellar mass cannot be greater than upper limit of "+str(MstarMax)+" Msun" )
    
  return

#====================================================================================================================

