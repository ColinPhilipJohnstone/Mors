
"""Module holding the Star class and related functions."""

# Imports for standard stuff needed here
import sys
import inspect
import numpy as np
import pickle

# Imports for Mors modules
import Mors.miscellaneous as misc
import Mors.stellarevo as SE
import Mors.parameters as params
import Mors.rotevo as RE
import Mors.physicalmodel as phys

# Limits on various input values
MstarMin = 0.1
MstarMax = 1.25

class Star:
    """A class for star objects that hold all information about a star. 
    
    Attributes
    ------------
    
    Methods
    ------------
    
    """
    
    #---------------------------------------------------------------------------------------

    def __init__(self,Mstar=None,Age=None,percentile=None,Omega=None,Prot=None,OmegaEnv=None,OmegaCore=None,AgesOut=None,starEvoDir=None,evoModels=None,params=params.paramsDefault):
        """Initialises instance of Star class.
        
        This is the main function that is run when creating an instance of the Star class and it sets up all
        the things needed including calculating evolutionary tracks for the star. The function requires that
        the arguments Mstar (the star's mass in Msun) and Omega (in OmegaSun=2.67e-6 rad s^-1) are specified
        in the call. Alternatively, OmegaEnv and OmegaCore can be set, in which case Omega does not need to 
        be specified. If the argument Age (in Myr) is also specified, then the code will find the evolutionary
        track that passes through this rotation rate at this age, otherwise if Age is not set then it will 
        calculate evolutionary tracks assuming this Omega as the initial (1 Myr) rotation rate. The user should
        not specify OmegaCore and Age simultaneously, and if Age is set then either Omega or OmegaEnv can be
        used to specify the surface rotation rate. The user can also specify the initial rotation rate using 
        
        Parameters
        ----------
        Mstar : float
            Stellar mass in Msun.
        Age : float , optional
            Age in Myr of the input rotation rate (default = 1 Myr).
        Omega : float , optional
            Rotation rate of star at Age in OmegaSun.
        Prot : float , optional
            Rotation period of star at Age in days.
        OmegaEnv : float , optional
            Envelope rotation rate of star at Age in OmegaSun.
        OmegaCore : float , optional
            Core rotation rate of star at Age in OmegaSun.
        percentile : float , optional
            Percentile in 1 Myr distribution (between 0 and 100).
            
        params : dict , optional
            Parameter dictionary used for getting width of bin to use for distribution.
        
        Returns
        ----------
        None
            None
        
        """
        
        # Make sure Mstar has good values
        _CheckInputMstar(Mstar)
        
        # Make sure Omega, OmegaEnv, and OmegaCore are well set
        # (if Omega is set, OmegaEnv and OmegaCore will both be set to that value)
        Omega , OmegaEnv , OmegaCore = _InputRotation(Mstar,Age,Omega,OmegaEnv,OmegaCore,Prot,percentile,params)
        
        # Set parameters
        self.params = params
        
        # Set the ExtendedTracks parameter to True so we get all parameters
        self.params['ExtendedTracks'] = True
        
        # Set Mstar
        self.Mstar = Mstar
        
        # Set stellar evo model
        self.starEvoDir = starEvoDir
        self.evoModels = evoModels
        
        # Load evo tracks if Mstar is set
        self._LoadStarEvo()
        
        # Get evolutionary tracks 
        self._LoadEvoTracks(Age,OmegaEnv,OmegaCore,AgesOut=AgesOut)
        
        # Get HZ boundaries
        self.aOrbHZ = phys.aOrbHZ(Mstar=self.Mstar,params=self.params)
        
        return
    
    def _LoadStarEvo(self,starEvoDir=None,evoModels=None):
        """Loads stellar evolution tracks."""
        
        # Get an instance of the StarEvo class holding evolutionary tracks for all masses
        self.StarEvo = SE.StarEvo(starEvoDir=self.starEvoDir,evoModels=self.evoModels)
        
        # Load evolutionary tracks for this mass only (so delete tracks for other masses after getting this one)
        self.StarEvo.LoadTrack(self.Mstar,ClearData=True)
        
        return
    
    def _LoadEvoTracks(self,Age,OmegaEnv0,OmegaCore0,AgesOut=None):
        """Loads rotation and activity tracks."""
        
        
        # Set starting age to 1 Myr
        self.AgeMin = 1.0
        
        # Get ending age as end of stellar evolution track
        self.AgeMax = self.StarEvo.ModelData[self.Mstar]['Age'][-1]
        
        # If needing to get initial rotation rate then get it
        if not Age is None:
            
            # Get initial rotation
            Omega0 = RE.FitRotation( Mstar=self.Mstar , Age=Age , Omega=OmegaEnv0 , 
                                        AgeMin=self.AgeMin , params=self.params , StarEvo=self.StarEvo )
            
            # Make sure a correct initial rotation rate was found
            if ( Omega0 == -1 ):
                misc._PrintErrorKill("input rotation rate too low for given mass and age")
            if ( Omega0 == -2 ):
                misc._PrintErrorKill("input rotation rate too high for given mass and age")
            if ( Omega0 == -3 ):
                misc._PrintErrorKill("input rotation rate in valid range for mass and age but solver was unable to fit track")
            
            # Set initial rotation rates for input into RE.EvolveRotation()
            OmegaEnv0 = Omega0
            OmegaCore0 = Omega0
            
        # Run evolution model using OmegaEnv0 and OmegaCore0 as initial rotation rate
        self.Tracks = RE.EvolveRotation( Mstar=self.Mstar , OmegaEnv0=OmegaEnv0 , OmegaCore0=OmegaCore0 , 
                                        AgeMin=self.AgeMin , AgeMax=self.AgeMax , AgesOut=AgesOut ,
                                        params=self.params , StarEvo=self.StarEvo )
        
        # Set the percentile of this track in the starting distribution
        self.percentile = Percentile( Mstar=self.Mstar , Omega=OmegaEnv0 )
        
        # Make evolutionary tracks attributes of this class 
        # (this is so the user can get a track using e.g. star.Age instead of star.Tracks['Age'])
        for track in self.Tracks:
            
            # Note: the following two options are basically equivalent
            #self.__dict__[track+'Track'] = self.Tracks[track]
            setattr( self , track+'Track' , self.Tracks[track] )
        
        # Load dictionary holding units for each quantity
        self.Units = phys.QuantitiesUnits()
        
        # Make functions for each quantity that return this quantity at a given age as attributes of class
        self._setupQuantityFunctions()
        
        return
            
    def _setupQuantityFunctions(self):
        """Makes functions for each quantity that return this quantity at a given age as attributes of class."""
        
        # Loop over quantities held in self.Tracks
        for track in self.Tracks:
            
            # The method for this is very strange and I don't like it. Better would likely just to be to 
            # add individual functions for each quantity since I was unable to find a solution that did 
            # not involve using exec() here. This should not provide a security risk since exec() is not
            # being run on any user input data. The following two lines first create a two-line function
            # from a string that calls self.Value() with a specific value for Quantity, so the function
            # Lx() that is created calls Value with Quantity='Lx'. This function is then added to the Star 
            # class using setattr on self.__class__. Technically this means that if the user makes multiple
            # instances of Star, the next two lines will not be necessary.
            exec( "def "+track+"(self,Age):\n  return self.Value(Age=Age,Quantity='"+track+"')" )
            exec( "setattr( self.__class__ , '"+track+"' , "+track+" )" )
            
        return
        
    def Value(self,Age=None,Quantity=None):
        """Takes age in Myr and a string with name of quantity to output, returns value of that quantity at specified age."""
        
        # Do 1D linear interpolation to get value
        value = SE._Interpolate1D( self.Tracks['Age'] , self.Tracks[Quantity] , Age )
        
        return value
    
    def Track(self,Quantity=None):
        """Takes a quantity string, return evolutionary track for that quantity."""
        
        # Make sure Quantity is set
        if Quantity is None:
            misc._PrintErrorKill("Quantity not set in call to function")
        
        # Make sure Quantity is string
        if not isinstance(Quantity,str):
            misc._PrintErrorKill("Quantity must be string")
        
        # Make sure Quantity is a valid quantity with a track
        if not Quantity in self.Tracks:
            misc._PrintErrorKill("no available track for "+Quantity)
        
        return self.Tracks[Quantity]
  
    def PrintAvailableTracks(self):
        """Prints all quantities that have evolutionary tracks available."""
        
        # Print header
        print("")
        print("EVOLUTIONARY TRACKS AVAILABLE FOR:-")
        
        # Loop over each element in Tracks and print them but only if they are numpy arrays
        for track in self.Tracks:
            if ( type(self.Tracks[track]).__module__ == 'numpy' ):
                
                # Start string with small space and name of quantity
                outString = "   "+track
                
                # If the units are given, add that in parentheses
                if not ( self.Units[track] == '' ):
                    outString += " ("+self.Units[track]+")"
                
                # Print the string 
                print(outString)
            
        # A free line at the end to look better
        print("")
        
        return
   
    def Save(self,filename='star.pickle'):
        """Takes filename (default is 'star.pickle'), saves cluster to this file using pickle."""
        with open(filename,'wb') as f:
            pickle.dump(self,f)
        return
  
    def save(self,filename='star.pickle'):
        """Same as Save()."""
        self.Save(filename=filename)
        return
  
    def ActivityLifetime(self,Quantity=None,Threshold=None,AgeMax=None):
        """Takes threshold value, returns age at which a star last drops below this threshold.
        
        This function can be used to determine when a star's emission crosses a given threshold value for a few
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
        
        # Allowed quantities
        QuantitiesAllowed = [ 'Lx' , 'Fx' , 'Rx' , 'FxHZ' , 'Leuv1' , 'Feuv1' , 'Reuv1' , 'Feuv1HZ' , 'Leuv2' , 'Feuv2' , 'Reuv2' , 'Feuv2HZ' , 
                            'Leuv' , 'Feuv' , 'Reuv' , 'FeuvHZ' , 'Lly' , 'Fly' , 'Rly' , 'FlyHZ' ]
        
        # Make sure Quantity is set
        if Quantity is None:
            misc._PrintErrorKill("Quantity not set in call to function")
        
        # Make sure Quantity is string
        if not isinstance(Quantity,str):
            misc._PrintErrorKill("Quantity must be string")
            
        # Check input Quantity is valid
        if not Quantity in QuantitiesAllowed:
            QuantitiesString = ''
            for quantity in QuantitiesAllowed:
                QuantitiesString += quantity + " , "
            misc._PrintErrorKill("invalid quantity "+Quantity+"\n valid options are "+QuantitiesString[0:-3])
        
        # Set tracks
        AgeTrack = self.AgeTrack
        Track = self.Tracks[Quantity]
            
        # If Threshold is 'sat' then normalise track to saturation value and set threshold to unity
        if ( Threshold == 'sat' ):
            
            # Get saturation rotation rate as function of age
            OmegaSatTrack = phys.OmegaSat( Mstar=self.Mstar , Age=self.AgeTrack , param='XUV' , params=self.params , StarEvo=self.StarEvo )
            
            # Make track just OmegaEnv/OmegaSat
            Track = self.OmegaEnvTrack / OmegaSatTrack
            
            # Set threshold to unity
            Threshold = 1.0
        
        # Get result
        AgeActive = misc.ActivityLifetime( Age=AgeTrack , Track=Track , Threshold=Threshold , AgeMax=AgeMax )
        
        return AgeActive
 
    def IntegrateEmission(self,AgeMin=None,AgeMax=None,Band=None,aOrb=None):
        """Takes age range, calculates integrated emission in band within that range.
        
        This code can be used to integrate a star's luminosity between two ages. This can be applied to any wavelength band
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
        aOrb : float or str , optional
            Orbital distance to get fluence at in AU or string identifying HZ boundary.
        
        Returns
        ----------
        Energy : float
            Integrated luminosity or flux in erg or erg cm^-2.
        
        """
        
        # Allowed quantities
        BandsAllowed = [ 'XUV', 'Xray', 'EUV1', 'EUV2', 'EUV', 'Lyman', 'bol' ]
        
        # Make sure Band is set
        if Band is None:
            misc._PrintErrorKill("Band not set in call to function")
        
        # Make sure Band is string
        if not isinstance(Band,str):
            misc._PrintErrorKill("Band must be string")
            
        # Check input Band is valid
        if not Band in BandsAllowed:
            BandsString = ''
            for band in BandsAllowed:
                BandsString += band + " , "
            misc._PrintErrorKill("invalid Band "+Band+"\n valid options are "+BandsString[0:-3])
            
        # Get luminosity track to use
        if Band == 'XUV':
            Luminosity = self.LxTrack + self.LeuvTrack
        elif Band == 'Xray':
            Luminosity = self.LxTrack
        elif Band == 'EUV1':
            Luminosity = self.Leuv1Track
        elif Band == 'EUV2':
            Luminosity = self.Leuv2Track
        elif Band == 'EUV':
            Luminosity = self.LeuvTrack
        elif Band == 'Lyman':
            Luminosity = self.LlyTrack
        elif Band == 'bol':
            Luminosity = self.LbolTrack
        else:
            misc._PrintErrorKill("did not find right luminosity track")
            
        # If aOrb was set to a string, do necessary work to get useable aOrb (given by aOrbUse)
        if isinstance(aOrb,str):
            
            # Valid options
            aOrbAllowed = [ 'RecentVenus' , 'RunawayGreenhouse' , 'MoistGreenhouse' , 'MaximumGreenhouse' , 'EarlyMars' , 'HZ' ]
            
            # Check input is valid
            if not aOrb in aOrbAllowed:
                aOrbString = ''
                for a in aOrbAllowed:
                    aOrbString += a + " , "
                misc._PrintErrorKill("invalid aOrb "+aOrb+"\n valid options are a value in AU or "+aOrbString[0:-3])
            
            # Get orbital distance
            aOrbUse = self.aOrbHZ[aOrb]
            
        else:
            
            # Just set to input value
            aOrbUse = aOrb
            
        # Do calculation
        Energy = misc.IntegrateEmission( AgeMin=AgeMin , AgeMax=AgeMax , Age=self.AgeTrack , Luminosity=Luminosity , aOrb=aOrbUse )
        
        return Energy

def _CheckInputMstar(Mstar):
    """Takes stellar mass, checks if value is a valid input mass."""
    
    # Make sure Mstar was set
    if Mstar is None:
        misc._PrintErrorKill("stellar mass not given")
    
    # Make sure it is a float
    #if not isinstance(Mstar,float):
        #misc._PrintErrorKill("stellar mass must be given as float")
    
    # Make sure it above minimum
    if ( Mstar < MstarMin ):
        misc._PrintErrorKill( "stellar mass cannot be less than lower limit of "+str(MstarMin)+" Msun" )
    
    # Make sure it above minimum
    if ( Mstar > MstarMax ):
        misc._PrintErrorKill( "stellar mass cannot be greater than upper limit of "+str(MstarMax)+" Msun" )
        
    return

def _InputRotation(Mstar,Age,Omega,OmegaEnv,OmegaCore,Prot,percentile,params):
    """Takes input rotation values, checks sets up values correctly."""
    
    # Make sure percentile and rotation are not both set
    if not percentile is None: 
        if not Omega is None :
            misc._PrintErrorKill( "cannot set both percentile and Omega as arguments of Star" )
        if not OmegaEnv is None :
            misc._PrintErrorKill( "cannot set both percentile and OmegaEnv as arguments of Star" )
        if not OmegaCore is None :
            misc._PrintErrorKill( "cannot set both percentile and OmegaCore as arguments of Star" )
        if not Prot is None :
            misc._PrintErrorKill( "cannot set both percentile and Prot as arguments of Star" )
        
    # Make sure percentile and OmegaEnv are not both set
    if ( not percentile is None ) and ( not OmegaEnv is None ):
        misc._PrintErrorKill( "cannot set both percentile and OmegaEnv as arguments of Star" )
    
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
        
    # If percentile was set, get Omega
    if not percentile is None:
        Omega = Percentile(Mstar=Mstar,percentile=percentile,params=params)
    
    # Make sure if Age is set that OmegaCore is not set
    if ( not Age is None ) and ( not OmegaCore is None ):
        misc._PrintErrorKill( "cannot set both Age and OmegaCore as arguments of Star" )
    
    # Make sure not both Omega and Prot were set
    if ( not Omega is None ) and ( not Prot is None ):
        misc._PrintErrorKill( "cannot set both Omega and Prot as arguments of Star" )
    
    # If Prot is set, get Omega
    if not Prot is None:
        Omega = phys._Omega(Prot)
    
    # Make sure at least one rotation rate is set
    if ( Omega is None ):
        
        # If Age is set then make sure OmegaEnv is set, otherwise make sure both OmegaEnv and OmegaCore are set
        if ( not Age is None ): 
            if ( OmegaEnv is None ):
                misc._PrintErrorKill( "if Age is set then must set either Omega or OmegaEnv as argument of Star" )
        elif ( OmegaEnv is None ) or ( OmegaCore is None ):
            misc._PrintErrorKill( "must set either Omega or both OmegaEnv and OmegaCore as argument of Star" )
    
    else:
        
        # Since Omega is set, make sure both OmegaEnv and OmegaCore are not set
        if not ( ( OmegaEnv is None ) and ( OmegaCore is None ) ):
            misc._PrintErrorKill( "cannot set OmegaEnv and OmegaCore when Omega is set as arugment of Star" )
        
        # Set both OmegaEnv and OmegaCore to Omega
        OmegaEnv = Omega
        OmegaCore = Omega
    
    return Omega , OmegaEnv , OmegaCore

def Percentile(Mstar=None,Omega=None,Prot=None,percentile=None,MstarDist=None,OmegaDist=None,ProtDist=None,params=params.paramsDefault):
    """Gets rotation rate of percentile or percentile of rotation rate in the model rotation distribution.
    
    This function can be used for two purposes
        1. to determine the percentile in a rotation distribution of a star given its mass and rotation rate
        2. to determine the rotation rate of a star in a rotation distribution given its mass and percentile
    In the first case, the user should specify the rotation rate using either the Omega or Prot keyword arguments (given 
    in OmegaSun and days respectively) and the mass. In the second case, the user should specify the percentile using the
    percentile keyword argument (in the range 0 to 100) and the mass. The mass should be specified in Msun using the Mstar
    keyword argument. These should only be given as floats. The user can also specify the rotation distribution using the 
    MstarDist and the OmegaDist or ProtDist keyword arguments. If this is not done, then the code will assume the 1 Myr 
    rotation distribution from the empirical model distribution used in Johnstone et al. (2020).
    
    Parameters
    ----------
    Mstar : float
        Stellar mass in Msun.
    Omega : float , optional
        Rotation rate of star in OmegaSun.
    Prot : float , optional
        Rotation period of star in days.
    percentile : float , optional
        percentile in distribution (between 0 and 100).
    MstarDist : numpy.ndarray , optional
        Array of masses of stars in distribution.
    OmegaDist : numpy.ndarray , optional
        Array of rotation rates of stars in distribution in OmegaSun.
    ProtDist : numpy.ndarray , optional
        Array of rotation periods of stars in distribution in days.
    params : dict , optional
        Parameter dictionary used for getting width of bin to use for distribution.
    
    Returns
    ----------
    result : float
        Either rotation rate for percentile or percentile for rotation rate.
    
    """
    
    # Make sure Mstar was set
    if Mstar is None:
        misc._PrintErrorKill( "keyword argument Mstar must be set" )
    
    # If percentile is not set, make sure rotation rate is set
    if percentile is None:
        # Not both
        if not Omega is None:
            if not Prot is None:
                misc._PrintErrorKill( "keyword arguments Omega and Prot cannot both be set" )
        # At least one
        if ( Omega is None ) and ( Prot is None ):
            misc._PrintErrorKill( "must set either Omega, Prot, or percentile keyword arguments" )
        # Set Omega is not set then set it
        if Omega is None:
            Omega = phys._Omega(Prot)
            
    # Setup the distribution
    if MstarDist is None:
        MstarDist , OmegaDist = misc.ModelCluster()
    else:
        
        # Check that one of OmegaDist and ProtDist are set
        if ( OmegaDist is None ) and ( ProtDist is None ):
            misc._PrintErrorKill( "one of keyword arguments OmegaDist and ProtDist must be set" )
        
        # Check that OmegaDist and ProtDist are not both set
        if ( not OmegaDist is None ) and ( not ProtDist is None ):
            misc._PrintErrorKill( "keyword arguments OmegaDist and ProtDist cannot both be set" )
        
        # Setup OmegaDist if ProtDist was set
        if OmegaDist is None:
            OmegaDist = misc._Omega(ProtDist)
        
        # Make sure arrays are same length
        if not ( len(MstarDist) == len(OmegaDist) ):
            misc._PrintErrorKill( "MstarDist and OmegaDist (or ProtDist) must be same length" )
            
    # Work out which one to do
    if not percentile is None:
        
        # percentile was set, so get corresponding rotation rates
        result = _OmegaPercentile(Mstar,percentile,MstarDist,OmegaDist,params)
    
    else:
        
        # Rotation rate was set, so get corresponding percentile
        result = _PerPercentile(Mstar,Omega,MstarDist,OmegaDist,params)
        
    
    return result 

def _OmegaPercentile(Mstar,percentile,MstarDist,OmegaDist,params):
    """Gets rotation rate that corresponds to a given percentile of the rotation distribution."""
    
    ## If percentile was a string, get actual value
    
    # Get part of distribution to include
    includeStars = np.where( np.abs(MstarDist-Mstar) <= params['dMstarPer'] )
    
    # Need at least two stars
    if ( len(includeStars[0]) < 2 ):
        misc._PrintErrorKill( "need at least two stars in mass bin to get percentile" )
    
    # Get percentile Omega
    Omega = np.percentile( OmegaDist[includeStars] , percentile )
    
    return Omega

def _PerPercentile(Mstar,Omega,MstarDist,OmegaDist,params):
    """Gets percentile of the rotation distribution for a given mass and rotation rate."""
    
    # Get part of distribution to include
    includeStars = np.where( np.abs(MstarDist-Mstar) <= params['dMstarPer'] )
    
    # Need at least two stars
    if ( len(includeStars[0]) < 2 ):
        misc._PrintErrorKill( "need at least two stars in mass bin to get percentile" )
    
    # If Omega is less than minimum then make this as the 0th percentile
    if ( Omega <= np.min(OmegaDist[includeStars]) ):
        return 0.0
    
    # If Omega is above maximum then make this as the 100th percentile
    if ( Omega >= np.max(OmegaDist[includeStars]) ):
        return 100.0
    
    # Use bisector to invert numpy.percentile function
    
    # Initial limits
    perMin = 0.0
    perMax = 100.0
    
    # Omegas for limits
    OmegaMin = np.percentile( OmegaDist[includeStars] , perMin )
    OmegaMax = np.percentile( OmegaDist[includeStars] , perMax )
    
    # Tolerance and maximum number of iterations
    tol = 1.0e-5
    nIterMax = 10000
        
    # Start iterating
    found = False
    for iIter in range(nIterMax):
        
        # Get mid value and percentle
        perMid = 0.5 * ( perMin + perMax )
        OmegaMid = np.percentile( OmegaDist[includeStars] , perMid )
        
        # See if OmegaMid is equal to or close enough to answer
        if ( abs(OmegaMid/Omega-1.0) < tol ):
            return perMid
        
        # See which limit to change
        if ( OmegaMid > Omega ):
            perMax = perMid
            OmegaMax = OmegaMid
        else:
            perMin = perMid
            OmegaMin = OmegaMid
        
    # Make sure it was found
    if not found:
        misc._PrintErrorKill( "did not find percentile "+str(per)+" for Omega of "+str(Omega) )
    
    return per
