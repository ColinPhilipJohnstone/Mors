
Mors - Model of rotational evolution

-------------------------------------------------------------------------------------------------------------

INSTALLATION

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
>>> myStar = mors.StellarEvoModel(starEvoDir=...)

where ... can be given as the path relative to the current directory. When this is done, no environmental variable needs to be set. 

-------------------------------------------------------------------------------------------------------------

BASIC STELLAR EVOLUTION

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

The function EvolveRotation returns a dictionary with evolutionary tracks for rotation. The main elements of this dictionary are Age in Myr, OmegaEnv in OmegaSun, and OmegaCore in OmegaSun. The following code will calculate and plot a rotation model for a solar mass star

>>> tracks = mors.EvolveRotation(Mstar=1.0,Omega0=10.0)
>>> plt.plot( tracks['Age'] , tracks['OmegaEnv'] )
>>> plt.show()




-------------------------------------------------------------------------------------------------------------






