
TBD - better descriptions of various parameters

-------------------------------------------------------------------------------------------------------------

MODEL FOR ROTATION OF STARS (MORS)
Author: Colin P. Johnstone 

-------------------------------------------------------------------------------------------------------------

This code solves the stellar rotation and XUV evolution model presented in Johnstone et al. (2020). The package can be used to calculate evolutionary tracks for stellar rotaton and X-ray, EUV, and Ly-alpha emission for stars with masses between 0.1 and 1.25 Msun and has additional functionality such as allowing the user to get basic stellar parameters such as stellar radius and luminosity as functions of mass and age using the stellar evolution models of Spada et al. (2013). When publishing results that were calculated using this code, both the Johnstone et al. (2020) paper and Spada et al. (2013) should be cited.


-------------------------------------------------------------------------------------------------------------

CONTENTS
1. INSTALLATION
2. EVOLUTIONARY CALCULATIONS 
3. PERCENTILES AND MODEL DISTRIBUTION
4. SETTING SIMULATION PARAMETERS
5. ROTATION AND ACTIVITY QUANTITES
6. STELLAR EVOLUTION QUANTITIES
7. CLUSTER EVOLUTION CALCULATIONS


-------------------------------------------------------------------------------------------------------------

1. INSTALLATION

1. Prerequisites: The code requires only that an up-to-date version of Python is installed with numpy and matplotlib included. 

2. Obtaining and installing the code: Make a clone of the repository on GitHub somewhere on your computer using

$ git clone https://github.com/ColinPhilipJohnstone/Mors

 The code can be installed from the command line inside the main directory, where the setup.py file is kept, using 

$ git install . 

If this does not work, it should be possible to just include the main directory and all its contents in the working directory of your code if you wish to do it that way.

3. Stellar evolution tracks: The code requires also the set of stellar evolution models from this paper

https://ui.adsabs.harvard.edu/abs/2013ApJ...776...87S/abstract

These models can be downloaded from here

http://www.astro.yale.edu/demarque/yyiso.html

Download and extract their package and put it somewhere safe. After this, it is necessary to tell the code where it can find the stellar evolution tracks. This can be done by setting the environmental variable 'STELLARMODELS', which on Ubuntu can be set using

$ export STELLAREVODIR=...

Where ... should be replaced with the path to the main directory holding the stellar models (i.e. the directory holding directories with names such as "X0p70952_Z0p01631_A1p875") and should end with a /. To make this permanent on Ubuntu, use

$ gedit ~/.profile

and add the export command to the bottom of the file. You will probably have to logout and login again for this to work. Alternatively, when creating a star object in your Python script, you can specify the path to this directory using the starEvoDir keyword

>>> import Mors as mors 
>>> myStar = mors.StarEvo(starEvoDir=...)

where ... can be given as the path relative to the current directory. When this is done, no environmental variable needs to be set. 

-------------------------------------------------------------------------------------------------------------

2. EVOLUTIONARY CALCULATIONS

The main way that the user is meant to interact with the code is through the Star class. The user can create an instance of the star class, specifying only the mass and star's initial (1 Myr) rotation rate using the Mstar and Omega0 keyword arguments. This can be done for a star with a mass of 1 Msun and an initial rotation of 10x the modern Sun using

>>> star = mors.Star(Mstar=1.0,Omega=10.0)

The user can instead set the initial rotation rate using the Prot keyword argument giving the surface rotation period in days.

>>> star = mors.Star(Mstar=1.0,Prot=2.7)

Alternatively, the user can specify starting values for both the core and envelope rotation rates using the OmegaEnv and OmegaCore arguments, though this is usually not recommended.

If the user instead specifies an age the Age and Omega keyword argument, the code will look for the track that passes through the specified surface rotation rate at the specified age. The surface rotation rate can be specified either using Omega or OmegaEnv.

>>> star = mors.Star(Mstar=1.0,Age=100.0,Omega=50.0)

In this case, the code will search for a rotation track that has a surface rotation rate that is 50x the current Sun's at an age of 100 Myr. It is not recommended to do this when the age of the star is so large that the rotation rates of stars with different initial rotation rates have converged (i.e. Age should only be specified if it is low enough that there is still a large spread in rotation rates for that mass at that age). In not all cases will it actually be able to find a physically realistic rotation rate, for example if the rotation rate specified is unreaslistically large.

The code will calculate evolutionary tracks for all rotation and activity quantities available. These include surface rotation rate, OmegaEnv, and X-ray luminosity, Lx. To see a list of quantities with calculated evolutionary tracks, use the PrintAvailableTracks() attribute of the Star class.

>>> star.PrintAvailableTracks()

This gives units for all quantities. This can also be used to get units for each of the available quantities. Tracks for each of the quantities are held in numpy arrays and can be accessed using the Tracks dictionary, which is an attribute of the Star class. For example, the user can plot evolutionary tracks for Lx using

>>> plt.plot( star.Tracks['Age'] , star.Tracks['Lx'] )

Alternatively, numpy arrays for each quantity are attributes of the Star class with the name of the quantity followed by 'Track', so the above can be replaced with

>>> plt.plot( star.AgeTrack , star.LxTrack )

Units are held in the dictionary Units, which is also an attribute of the Star class. For example, to see the units of Lx, the user can use

>>> print( star.Units['Lx'] )

The value of one of these quantities at a given age can be output using the Value() attribute of the Star class giving age the Myr and a string with the name of the desired quantity. These can be given either as positional or as keyword arguments. For example, to print Lx at 150 Myr to the screen you can use

>>> print( star.Value(150.0,'Lx') )

or

>>> print( star.Value(Age=150.0,Quantity='Lx') )

More simply, the user can use the function with the name of the desired quantity, so the two above lines could be replaced with

>>> print( star.Lx(150.0) )

Instances of the Star class can be saved using the Save (or save) functions,

>>> star.Save()

By default, the data will be saved into a file called 'star.pickle' but using the user can specify the name o the file using the filename argument in the call to Save(). Saved stars can be loaded using the Load() function.

>>> star = mors.Load("star.pickle")

-------------------------------------------------------------------------------------------------------------

3. PERCENTILES AND MODEL DISTRIBUTION

It is possible when creating an instance of the Star class to specify the initial rotation as a percentile of the model distribution from Johnstone et al. (2020). This model distribution is composed of measured rotation distributions from several clusters with ages of ~150 Myr evolved back to 1 Myr. To get the masses and initial rotation rates of the model cluster, use the function ModelCluster().

>>> Mstar , Omega = ModelCluster()

If this model cluster is used in any research, please cite the papers listed in Table 1 of Johnstone et al. (2020) for the 150 Myr age bin as the original sources for the rotation measurements (note that the function ModelCluster does not return these measurements, but 1 Myr rotation rates for these stars derived from their measurments using the rotational evolution model of Johnstone et al. 2020). To evolve this cluster, use the Cluster class as described below.

When creating an instance of the Star class, use the keyword argument Percentile to set the initial rotation rate to a percentile of this distribution. This can either be given as a float or int between 0 and 100 or as a string, with the three options being 'slow', 'medium', and 'fast', corresponding to the 5th, 50th, and 95th percentiles of the rotation distribution. For example

>>> star = mors.Star(Mstar=1.0,Percentile=5.0)

This creates a solar mass star with an initial rotation rate equal to the 5th percentile of the rotation distribution at 1 Myr, as determined from our model cluster. This is equivalent to

>>> star = mors.Star(Mstar=1.0,Percentile='slow')

To calculate this percentile, the parameter dMstarPer is used and can be set if the user sets up their own parameters as discussed above. This is set by default to 0.1, meaning that all stars within 0.1 Msun of the specified Mstar will be considered.





-------------------------------------------------------------------------------------------------------------

4. SETTING SIMULATION PARAMETERS

In addition to just the masses and initial rotation rates, the basic behavior of the code depends on a large number of parameters all of which have default values and do not need to be changed by the user. These default values cause the code to run the evolutionary model for rotation and XUV emission described in Johnstone et al. (2020) and generally do not need to be changed by the user. However, if the user wishes, these parameters can be changed.

In general, simulation parameters are held in a dictionary called 'params' within the code, and the user can specify this parameter dictionary when creating a instance of the Star class. When this is not done, the code will use the paramsDefault dictionary created in parameters.py in the source code. To use user specified set of parameters, the user must first generate a parameter dictionary using the NewParams() function.

>>> myParams = mors.NewParams()

This will create a dictionary that is identical to paramsDefault that the user can edit as they want. For example, to change the parameter param1 to 1.5 and param2 to 2.5, use

>>> myParams = mors.NewParams()
>>> myParams['param1'] = 1.5
>>> myParams['param2'] = 2.5

Alternatively, and probably preferrably, this can also be done in the initial call to NewParams using keyword arguments.

>>> myParams = mors.NewParams(param1=1.5,param2=2.5)

This user should then input this as a keyword argument when creating an instance of the Star class.

>>> star = mors.Star(Mstar=1.0,Omega=10.0,params=myParams)

To see a complete list of all parameters that should be set, the user can use the PrintParams() function.

>>> mors.PrintParams()

Or can simply look into the parameters.py file in the main directory where Mors is installed.

The user can set the keyword parameter AgesOut when creating a new Star object to a number or a numpy array and this will cause the code to only include these ages and the starting age in the evolutionary tracks. The ages should be given in Myr. For example

>>> star = mors.Star(Mstar=1.0,Omega=10.0,AgesOut=100.0)

will cause the evolutionary tracks to only contain the starting age of 1 Myr and the specified age of 100 Myr. Similarly

>>> star = mors.Star(Mstar=1.0,Omega=10.0,AgesOut=np.array([100.0,200.0,300.0,400.0]))

will cause the evolutionary tracks to contain 1 Myr, and the specified 100, 200, 300, and 400 Myr ages. The simulations will also end at the oldest year in the array. This array should be in ascending order. This is not recommended if the user wants to extract quantities at arbitrary ages along evolutionary tracks, for example using star.Lx(Age), since these functions interpolate between the values and if there are not enough age bins in the evolutionary tracks then these results can be inaccurate. Note that having two many age bins in the AgesOut array, e.g. with AgesOut=np.linspace(1.0,5000.0,10000), will cause the calculations of the evolutionary tracks to be very slow.

-------------------------------------------------------------------------------------------------------------

5. ROTATION AND ACTIVITY QUANTITES

The code gives the user the ability to use the basic functions for calculating many activity related quantities as functions os stellar properties. For example, the user can calculate XUV luminosities from stellar mass, age, and rotation rates using the Lxuv function. For example, 

>>> Lxuv = mors.Lxuv(Mstar=1.0,Age=5000.0)

This gives a dictionary holding X-ray, EUV, and Lyman-alpha luminosities for a 5000 Myr old solar mass star rotating 10 times faster than the Sun. The surface rotation rate can be specified as a rotation velocity using the Omega or OmegaEnv arguments or as a rotation period in days using the Prot argument, e.g.

>>> Lxuv = mors.Lxuv(Mstar=1.0,Age=5000.0,Prot=1.0)

This is similar to above, but for a star with a rotation period of 1 day.

The dictionary returned contains the following luminosities, all in erg s^-1

Lxuv - 0.517 to 92 nm
Lx - 0.517 to 12.4 nm (0.1 to 24 keV)
Leuv1 - 10 to 36 nm
Leuv2 - 36 to 92 nm
Leuv - 10 to 92 nm
Lly - Lymann-alpha emission line

The user can also just get the X-ray luminosity as a float using Lx()

>>> Lx = mors.Lx(Mstar=1.0,Age=5000.0,Omega=10.0)

And similarly for the EUV and Lymann-alpha luminosities

>>> Leuv = mors.Leuv(Mstar=1.0,Age=5000.0,Omega=10.0)
>>> Lly = mors.Lly(Mstar=1.0,Age=5000.0,Omega=10.0)

The function Leuv() also takes the keyword argument 'band' which can be used to specify the band desired (=0 for 10-92 nm; =1 for 10-32 nm; =2 for 32-92 nm).

If the user wants a more detailed set of parameters for a given star, the ExtendedQuantities() function can be used and in this case, the star's mass, age, and envelope and core rotation rates must be specified.

>>> quantities = ExtendedQuantities(Mstar=1.0,Age=5000.0,OmegaEnv=10.0,OmegaCore=10.0)

This returns a larger dictionary with many other quantities, including some basic stellar properties such as radius and moment of inertia, the mass loss rate in the wind, the dipole field strength, and all of the torques acting on the star to influence its rotation. A list of the available parameters can be seen with

>>> list(quantities)

-------------------------------------------------------------------------------------------------------------

6. STELLAR EVOLUTION QUANTITIES

The rotation and activity evolution model requires that various basic stellar properties such as bolometric luminosity, radius, and convective turnover time, can be accessed at all ages for for all masses within the mass range considered (0.1 to 1.25 Msun). These are calculated using the evolutionary tracks fronm Spada et al. (2013) and the functions that do this are available to the user. First import the Mors package

>>> import Mors as mors

The stellar mass for a 1.0 Msun star with an age of 1000 Myr can be calculated using

>>> Rstar = mors.Rstar(1.0,1000.0)

The functions available are 

Rstar - radius (Rsun)
Lbol - bolometric luminosity (Lsun)
Teff - effective temperature (K)
Itotal - total moment of inertia in (g cm^2)
Icore - core moment of inertia in (g cm^2)
Ienv - envelope moment of inertia in (g cm^2)
Mcore - core mass in (Msun)
Menv - envelope mass in (Msun)
Rcore - core radius in (Rsun)
tau - convective turnover time (days)
dItotaldt - rate of change of total moment of inertia  (g cm^2 Myr^-1)
dIcoredt - rate of change of core moment of inertia  (g cm^2 Myr^-1)
dIenvdt - rate of change of envelope moment of inertia  (g cm^2 Myr^-1)
dMcoredt - rate of change of core mass (Msun Myr^-1)
dRcoredt - rate of change of core radius (Rsun Myr^-1)

Note that core and envelope have the definitions from the rotation model, so the 'core' is not the hydrogen burning core as it would typically be defined but everything interior to the outer convective zone, and the envelope is the outer convective zone. All of these functions take stellar mass in Msun and age in Myr. Alternatively, the user can call the function Value which takes the same two inputs and then the parameter name as a string. For example

>>> Rstar = mors.Value(1.0,1000.0,'Rstar')

which is equivalent to the above example. 

In all cases, multiple values for mass and age can be input using lists or numpy arrays and the functions will return numpy arrays holding the results for each input case. For example, if stellar mass is input as a 1D array or list of length 5, the result will be a length 5 numpy array. If both stellar mass and age are input as arrays, the result will be a 2-dimensional array of length len(Mstar)xlen(Age). Additionally, in the call to Value(), the string holding the parameter to calculate for can also be input as a list of strings, in which case the variable retuned will have an additional dimension with values for each of the input quantities. 

The first time one of these functions is called, the code loads all of the evolutionary tracks from the Spada et al. (2013) model. These will be written to the directors in the file SEmodels.pickle in the main directory of the evolutionary models and this file will be used to load the models in the future. This file can be safely deleted if the user wants since this will just cause the code to remake it next time it is needed. The code then does a linear interpolation between mass and age bins to get values at the input masses and ages. This can be time consuming, especially if an interpolation between mass bins is required, which will be the case if the input age is not an integer multiple of 0.05 Msun. The user can load a full evolutionary track for a specific stellar mass using LoadTrack(). For example, to load a track for a 1.02 Msun star, use

>>> mors.LoadTrack(1.02)

If a track for that specific mass is already loaded, this will do nothing. 

-------------------------------------------------------------------------------------------------------------

7. CLUSTER EVOLUTION CALCULATIONS

The code allows the user to calculate the evolution of stellar clusters in addition to single stars. This is done using the Cluster class. An instance of the Cluster class can be created in much the same way as an instance of the Star class using the same input keyword arguments. The only difference is that the masses and rotation rates of the stars should be given as arrays or lists

>>> Mstar = np.array( [ 1.0 , 0.5 , 0.75 ] )
>>> Omega = np.array( [ 10.0 , 10.0 , 10.0 ] )
>>> cluster = mors.Cluster(Mstar=Mstar,Omega=Omega)

This will create a cluster composed of three stars and immediately calculate evolutionary tracks for each. To see the progress of these calculations printed to the screen, use the verbose keyword argument

>>> cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)

As in the Star class, it is possible to specify the Age keyword argument. If this is not specified, the code will assume the Omega values are starting (1 Myr) rotation rates. If it is specified, the code will fit the rotation tracks such that they pass through the specified rotation rates at the speficied ages. If Age is specified as a single value (in Myr), it will be assumed that all stars have that age, and if it is specified as a numpy array then the array should have the same lengths as Mstar and Omega and the ages will be assumed individually for each star. As with the Star class, it is possible to input simulation parameters using the params keyword argument (see above).

As in the Star class, an instance of the Cluster class can be saved using the Save (or save) function

>>> cluster.Save()

By default, the data will be saved into a file called 'cluster.pickle' but using the user can specify the name o the file using the filename argument in the call to Save(). Saved clusters can be loaded using the Load() function.

>>> cluster = mors.Load("cluster.pickle")

This is advisable since the evolutionary tracks for clusters can take a lot of time to calculate if there are lots of stars.

The Star class instances for each star in the cluster is held in the stars list, which is an attribute of this class. To plot Lx tracks for each star for example, the user can use

>>> plt.plot(cluster.stars[0].AgeTrack,cluster.stars[0].LxTrack)
>>> plt.plot(cluster.stars[1].AgeTrack,cluster.stars[1].LxTrack)
>>> plt.plot(cluster.stars[2].AgeTrack,cluster.stars[2].LxTrack)

Alternatively, each star is an attribute of the class instance and has the name 'star' followed by its place in the list (starting at zero), so the above can be replaced with

>>> plt.plot(cluster.stars0.AgeTrack,cluster.stars0.LxTrack)
>>> plt.plot(cluster.stars1.AgeTrack,cluster.stars1.LxTrack)
>>> plt.plot(cluster.stars2.AgeTrack,cluster.stars2.LxTrack)

To get values for each star at a given age of a given parameter, the user can use the Values() function, specifying the age in Myr and a string for the quantity to retrieve. For example,

>>> Lx = cluster.Values(Age=100.0,Quantity='Lx')

This gives an array with Lx for each star in the cluster. Alternatively, as with the Star class, each quantity can be specified using its own function and the above is equivalent to.

>>> Lx = cluster.Lx(100.0)

So to plot the Lx distribution for a cluster at 100 Myr, the user can use

>>> plt.scatter(cluster.Mstar,cluster.Lx(100.0))

In Johnstone et al. (2020), many of the results are based on a composite cluster that is often referred to in the paper as the model cluster. This is composed of measured rotation distributions from several clusters with ages of ~150 Myr evolved back to 1 Myr. To get the masses and initial rotation rates of the model cluster, use the function ModelCluster().

>>> Mstar , Omega = ModelCluster()

If this model cluster is used in any research, please cite the papers listed in Table 1 of Johnstone et al. (2020) for the 150 Myr age bin as the original sources for the rotation measurements (note that the function ModelCluster does not return these measurements, but 1 Myr rotation rates for these stars derived from their measurments using the rotational evolution model of Johnstone et al. 2020). This cluster can then be evolved in the usual way

>>> Mstar , Omega = ModelCluster()
>>> cluster = mors.Cluster(Mstar=Mstar,Omega=Omega)



-------------------------------------------------------------------------------------------------------------

8. OTHER USEFUL FUNCTIONS






BELOW IS PART OF OLD DESCRIPTION AND MAYBE SHOULD BE REMOVED

-------------------------------------------------------------------------------------------------------------

RUNNING ROTATIONAL EVOLUTION MODEL

A rotational evolution model can be run using the function EvolveRotation. In the simplest case, a rotation model can be run simply by specifying the mass and initial rotation rate of the stars in the call to this function with no previous setup required. For example, to get the rotation track for a solar mass star with an initial rotation rate that is 10x the modern Sun's rotation (assumed throughout the code to be 2.67e-6 rad s^-1), use

>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0)

In this case, the evolutionary track starts are 1 Myr and ends at 5 Gyr. The user can specify starting and ending times in Myr using the AgeMin and AgeMax keyword arguments. 

>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0,AgeMin=1.0,AgeMax=100.0)

This does the same thing, but ends the time integration at 100 Myr. If core-envelope decoupling is being used, the initial core and envelope rotation rates can be specified separately using the OmegaEnv0 and OmegaCore0 arguments. 

The behavior of the code is determined by a large number of parameters that set which physical processes are included, how these processes are included, and other details about the run. The parameters that exist and their default values are held in the dictionary paramsDefault and can be seen using the function PrintParams.

>>> mors.PrintParams()

If the user wants the code to use a different set of parameters, a new parameter dictionary can be generated using the NewParams() function.

>>> myParams = mors.NewParams()

This will create a dictionary that the user can edit as they want. For example, to change the parameter 'param1' to 1.5

>>> myParams['param1'] = 1.5

Alternatively, this can also be done in the initial call to NewParams using keyword arguments

>>> myParams = mors.NewParams(param1=1.5)

These parameters can then be input into EvolveRotation using the params keyword argument.

>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0,params=myParams)

The user might also want to use a specific instance of the stellar evolution model, for example using different metallicities. This can be done by first creating and instance of the StarEvo class and then inputting it into EvolveRotation() using the StarEvo keyword argument.

>>> StarEvo = mors.StarEvo()
>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0,params=myParams,StarEvo=StarEvo)

-------------------------------------------------------------------------------------------------------------

LOOKING AT RESULTS FROM ROTATION MODEL

The function EvolveRotation returns a dictionary with evolutionary tracks for rotation. The main elements of this dictionary are Age in Myr, OmegaEnv in OmegaSun, and OmegaCore in OmegaSun. The following code will calculate and plot a rotation model for a solar mass star

>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0)
>>> plt.plot( tracks['Age'] , tracks['OmegaEnv'] , 'k' )
>>> plt.plot( tracks['Age'] , tracks['OmegaCore'] , 'k:' )
>>> plt.show()

By default, EvolveRotation() returns tracks for the envelope and core rotatation rates, so the dictionary that is returned will contain 'Age', 'dAge', 'OmegaEnv', and 'OmegaCore', where dAge is the length of each timestep (note that dAge[i]=Age[i]-Age[i-1] and dAge[0] is set to the starting value of dAge which for solvers that automatically determine step length will never be included). If core-envelope decoupling is not inlcuded (either because it was turned off as a parameter or because the mass of the star is below the fully convective boundary), then core and envelope rotation tracks will be identical.

The code is able to calculate a large number of rotation and activity related quantities. Many of these quantities are calculated by the code when running the rotation model, but their evolutionary tracks are not saved by default. The user can however can have these also returned using the ExtendedTracks parameter.

>>> myParams = mors.NewParams(ExtendedTracks=True)
>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0,params=myParams)



-------------------------------------------------------------------------------------------------------------


