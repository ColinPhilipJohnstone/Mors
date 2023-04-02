
"""
This code can be used to generate many of the figures appearing in Johnstone et al. (2020).

Note that many of the figures cannot be reproduced using this code, mostly since they contain data from
other studies (e.g. rotation or X-ray measurements for young clusters) which are not distributed with
this package. Also, this code does not attempt to format the plots in the same way as in the paper, mostly
for simplicity since many of these figures are very complex and the aim of this script is to demonstrate
how the package can be used and not to show how these figures were generated.

For Fig. 1 from the paper for example, use

$ python MakePlots.py Fig1

For all figures from the paper for example, use

$ python MakePlots.py All

"""

import Mors as mors
import Mors.constants as const

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

#==================================================================================================================

def Main(SystemArguments):
  
  """Just checks input parameters and then calls the right function."""
  
  # Make sure one argument was specified
  if not (len(SystemArguments) == 2 ):
    print("Give one and only one argument")
    print("e.g.")
    print(">>> python MakePlots.py Fig1")
    sys.exit()
  
  # Determine which figure to do
  if ((SystemArguments[1] == 'Fig1') or (SystemArguments[1] == 'All')):
    Fig1()
  if ((SystemArguments[1] == 'Fig2') or (SystemArguments[1] == 'All')):
    Fig2()
  if ((SystemArguments[1] == 'Fig3') or (SystemArguments[1] == 'All')):
    Fig3()
  if ((SystemArguments[1] == 'Fig4') or (SystemArguments[1] == 'All')):
    Fig4()
  if ((SystemArguments[1] == 'Fig5') or (SystemArguments[1] == 'All')):
    Fig5()
  if ((SystemArguments[1] == 'Fig6') or (SystemArguments[1] == 'All')):
    Fig6()
  if ((SystemArguments[1] == 'Fig7') or (SystemArguments[1] == 'All')):
    Fig7()
  if ((SystemArguments[1] == 'Fig8') or (SystemArguments[1] == 'All')):
    Fig8()
  if ((SystemArguments[1] == 'Fig9') or (SystemArguments[1] == 'All')):
    Fig9()
  if ((SystemArguments[1] == 'Fig10') or (SystemArguments[1] == 'All')):
    Fig10()
  if ((SystemArguments[1] == 'Fig11') or (SystemArguments[1] == 'All')):
    Fig11()
  if ((SystemArguments[1] == 'Fig12') or (SystemArguments[1] == 'All')):
    Fig12()
  if ((SystemArguments[1] == 'Fig13') or (SystemArguments[1] == 'All')):
    Fig13()
  if ((SystemArguments[1] == 'Fig14') or (SystemArguments[1] == 'All')):
    Fig14()
  if ((SystemArguments[1] == 'Fig15') or (SystemArguments[1] == 'All')):
    Fig15()
  if ((SystemArguments[1] == 'Fig16') or (SystemArguments[1] == 'All')):
    Fig16()
  if ((SystemArguments[1] == 'Fig17') or (SystemArguments[1] == 'All')):
    Fig17()
  if ((SystemArguments[1] == 'Fig18') or (SystemArguments[1] == 'All')):
    Fig18()
  if ((SystemArguments[1] == 'Fig19') or (SystemArguments[1] == 'All')):
    Fig19()
  if ((SystemArguments[1] == 'Fig20') or (SystemArguments[1] == 'All')):
    Fig20()
  
  return

#==================================================================================================================

def Fig1():
  
  """Plots the stuff from Fig 1."""
  
  #------------------------------------------------
  
  def UpperPanel():
    """Plots upper panel of Fig 1."""
    
    # Parameters
    nFrac = 100
    fracMin = 0.2
    fracMax = 1.0
    
    # Array for OmegaEnv/OmegaBreak
    frac = np.linspace(fracMin,fracMax,nFrac)
    
    # Just choose some stellar mass and age (should not really matter)
    Mstar = 1.0
    Age = 5000.0
    
    # Get radius 
    Rstar = mors.Rstar(Mstar,Age)
    
    # Get break-up rotation
    OmegaBreak = mors.OmegaBreak(Mstar,Rstar)
    
    # For each value, get mass loss factor
    MdotFactor = np.zeros(nFrac)
    for i in range(nFrac):
      MdotFactor[i] = mors.MdotFactor(Mstar,Rstar,frac[i]*OmegaBreak)
    
    # Make plot
    plt.figure()
    plt.plot( frac , MdotFactor )
    plt.xlabel("$\Omega_\mathrm{env} / \Omega_\mathrm{break}$")
    plt.ylabel("Mass loss factor")
    plt.yscale('log')
    plt.savefig("Fig1_upper.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def LowerPanel():
    """Plots lower panel of Fig 1."""
    
    # Ages to use
    AgeMin = 1.0
    AgeMax = 5000.0
    nAge = 100
    Age = np.logspace(np.log10(AgeMin),np.log10(AgeMax),nAge)
    
    # Masses to plot
    Mstar = np.array( [ 0.1 , 0.25 , 0.5 , 0.75 , 1.0 , 1.2 ] )
    
    # Array to hold results
    OmegaBreak = np.zeros((len(Mstar),nAge))
    
    # Loop over mass bins
    for i in range(len(Mstar)):
      Rstar = mors.Rstar(Mstar[i],Age)
      OmegaBreak[i,:] = mors.OmegaBreak(Mstar[i],Rstar[:])
    
    # Make plot
    plt.figure()
    for i in range(len(Mstar)):
      plt.plot( Age[:] , OmegaBreak[i,:] )
    plt.xlabel("Age (Myr)")
    plt.ylabel("Break-up rotation rate ($\Omega_\odot$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Fig1_lower.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Call panels
  UpperPanel()
  LowerPanel()
  
  return

#==================================================================================================================

def Fig2():
  print("Fig. 2 not included here.")
  return

#==================================================================================================================

def Fig3():
  
  #------------------------------------------------
  
  def Panel(Mstar):
    """Plots a panel from Fig. 3."""
    
    # Make stars
    starSlow = mors.Star( Mstar=Mstar , percentile='slow' )
    starMedium = mors.Star( Mstar=Mstar , percentile='medium' )
    starFast = mors.Star( Mstar=Mstar , percentile='fast' )
    
    # Make plot
    plt.figure()
    
    plt.plot( starSlow.AgeTrack , starSlow.OmegaEnvTrack , 'r' )
    plt.plot( starSlow.AgeTrack , starSlow.OmegaCoreTrack , 'r:' )
    
    plt.plot( starMedium.AgeTrack , starMedium.OmegaEnvTrack , 'g' )
    plt.plot( starMedium.AgeTrack , starMedium.OmegaCoreTrack , 'g:' )
    
    plt.plot( starFast.AgeTrack , starFast.OmegaEnvTrack , 'b' )
    plt.plot( starFast.AgeTrack , starFast.OmegaCoreTrack , 'b:' )
    
    plt.xlabel("Age (Myr)")
    plt.ylabel("Rotation rate ($\Omega_\odot$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.0,5000.0])
    plt.ylim([0.5,300.0])
    plt.savefig("Fig3_"+str(Mstar)+"Msun.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Panel for each mass bin
  Panel(1.0)
  Panel(0.75)
  Panel(0.5)
  Panel(0.25)
  
  return

#==================================================================================================================

def Fig4():
  print("Fig. 4 not included here.")
  return

#==================================================================================================================

def Fig5():
  """Plots stuff from Fig 5."""
  
  #------------------------------------------------
  
  def UpperPanel(Mstar,Age):
    """Plots upper panel of Fig 5."""
    
    # Get tauConv
    tauConv = mors.tauConv(Mstar,Age)
    
    # Make plot
    plt.figure()
    for i in range(len(Mstar)):
      plt.plot( Age[:] , tauConv[i,:] )
    plt.xlabel("Age (Myr)")
    plt.ylabel("Convective turnover time (days)")
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Fig5_upper.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def MiddlePanel(Mstar,Age):
    """Plots middle panel of Fig 5."""
    
    # Saturation rotation rate
    OmegaSat = mors.OmegaSat(Mstar=Mstar,Age=Age)
    
    # Make plot
    plt.figure()
    for i in range(len(Mstar)):
      plt.plot( Age[:] , OmegaSat[i,:] )
    plt.xlabel("Age (Myr)")
    plt.ylabel("Saturation rotation rate ($\Omega_\mathrm{sat}$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Fig5_middle.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def LowerPanel(Mstar,Age):
    """Plots lower panel of Fig 5."""
    
    # Saturation rotation rate
    OmegaSat = mors.OmegaSat(Mstar=Mstar,Age=Age)
    
    # Loop over masses and get Lx sat
    LxSat = np.zeros((len(Mstar),len(Age)))
    for iMstar in range(len(Mstar)):
      for iAge in range(len(Age)):
        LxSat[iMstar,iAge] = mors.Lx( Mstar=Mstar[iMstar] , Age=Age[iAge] , Omega=OmegaSat[iMstar,iAge] )
    
    # Make plot
    plt.figure()
    for i in range(len(Mstar)):
      plt.plot( Age[:] , LxSat[i,:] )
    plt.xlabel("Age (Myr)")
    plt.ylabel("Saturation $L_\mathrm{X}$ (erg s$^{-1}$ cm$^{-2}$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Fig5_lower.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Masses to plot
  Mstar = np.array( [ 0.1 , 0.3 , 0.4 , 0.6 , 0.8 , 1.0 , 1.1 , 1.2 ] )
  
  # Ages to use
  AgeMin = 1.0
  AgeMax = 5000.0
  nAge = 100
  Age = np.logspace(np.log10(AgeMin),np.log10(AgeMax),nAge)
  
  UpperPanel(Mstar,Age)
  MiddlePanel(Mstar,Age)
  LowerPanel(Mstar,Age)
  
  return

#==================================================================================================================

def Fig6():
  print("Fig. 6 not included here.")
  return

#==================================================================================================================

def Fig7():
  print("Fig. 7 not included here.")
  return

#==================================================================================================================

def Fig8():
  print("Fig. 8 not included here.")
  return

#==================================================================================================================

def Fig9():
  
  #------------------------------------------------
  
  def LeftColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for i in range(len(Ages)):
      
      # Start subplot
      plt.subplot(len(Ages),1,i+1)
      
      # Plot stars
      plt.scatter(cluster.Mstar,cluster.Values(Age=Ages[i],Quantity='OmegaEnv'),s=5)
      
      # Saturation rotation rate
      Mstar = np.linspace(0.1,1.2,100)
      OmegaSat = mors.OmegaSat(Mstar=Mstar,Age=Ages[i])
      plt.plot(Mstar,OmegaSat,'--k')
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_odot$)")
      plt.ylabel("Rotation rate ($\Omega_odot$)")
      plt.xlim([0.08,1.22])
      plt.ylim([0.1,150.0])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig9_leftcolumn.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def MiddleColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for i in range(len(Ages)):
      
      # Start subplot
      plt.subplot(len(Ages),1,i+1)
      
      # Get X-rays (including random scatter)
      Lx = cluster.Values(Age=Ages[i],Quantity='Lx')
      deltaLx = mors.XrayScatter(Lx)
      
      # Plot stars
      plt.scatter(cluster.Mstar,Lx+deltaLx,s=5)
      
      # Saturation Lx
      Mstar = np.linspace(0.1,1.2,100)
      OmegaSat = mors.OmegaSat(Mstar=Mstar,Age=Ages[i])
      LxSat = np.zeros(len(Mstar))
      for iMstar in range(len(Mstar)):
        LxSat[iMstar] = mors.Lx(Mstar=Mstar[iMstar],Omega=OmegaSat[iMstar],Age=Ages[i])
      plt.plot(Mstar,LxSat,'--k')
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_odot$)")
      plt.ylabel("X-ray luminosity (erg s$^{-1}$)")
      plt.xlim([0.08,1.22])
      plt.ylim([1.0e26,1.0e31])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig9_middlecolumn.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def RightColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for i in range(len(Ages)):
      
      # Start subplot
      plt.subplot(len(Ages),1,i+1)
      
      # Get X-rays (including random scatter)
      Lx = cluster.Values(Age=Ages[i],Quantity='Lx')
      deltaLx = mors.XrayScatter(Lx)
      
      # Get HZ distances for each star
      #aOrbHZ = mors.aOrbHZ(Mstar=cluster.Mstar)
      aOrbHZ = cluster.aOrbHZ
      
      # Get HZ flux
      FxHZ = ( Lx + deltaLx ) / ( 4.0 * const.Pi * (aOrbHZ['HZ']*const.AU)**2.0 )
      
      # Plot stars
      plt.scatter(cluster.Mstar,FxHZ,s=5)
      
      # Saturation FxHZ
      Mstar = np.linspace(0.1,1.2,100)
      OmegaSat = mors.OmegaSat(Mstar=Mstar,Age=Ages[i])
      LxSat = np.zeros(len(Mstar))
      for iMstar in range(len(Mstar)):
        LxSat[iMstar] = mors.Lx(Mstar=Mstar[iMstar],Omega=OmegaSat[iMstar],Age=Ages[i])
      aOrbHZ = mors.aOrbHZ(Mstar=Mstar)
      FxSat = LxSat / ( 4.0 * const.Pi * (aOrbHZ['HZ']*const.AU)**2.0 )
      plt.plot(Mstar,FxSat,'--k')
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_odot$)")
      plt.ylabel("X-ray flux in HZ (erg s$^{-1}$ cm$^{-2}$)")
      plt.xlim([0.08,1.22])
      plt.ylim([1.0e-1,1.0e5])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig9_lowercolumn.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Restore or calculate cluster
  if os.path.isfile('cluster.pickle'):
    cluster = mors.Load('cluster.pickle')
  else:
    Mstar , Omega = mors.ModelCluster()
    cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)
    cluster.save()
  
  # Ages to plot
  Ages = [ 2.0 , 10.0 , 100.0 , 500.0 , 1000.0 , 5000.0 ]
  
  # Make columns
  LeftColumn()
  MiddleColumn()
  RightColumn()
  
  return

#==================================================================================================================

def Fig10():
  
  # Restore or calculate cluster
  if os.path.isfile('cluster.pickle'):
    cluster = mors.Load('cluster.pickle')
  else:
    Mstar , Omega = mors.ModelCluster()
    cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)
    cluster.save()
  
  # Get saturation time for each star (assume Lx arbitrarily)
  AgeSat = cluster.ActivityLifetime(Quantity='Lx',Threshold='sat')
  
  # Make plot
  plt.figure()
  plt.scatter(cluster.Mstar,AgeSat,s=5)
  plt.xlabel("Stellar mass (M$_\odot$)")
  plt.ylabel("Saturation time (Myr)")
  plt.yscale('log')
  plt.savefig("Fig10.png")
  plt.close()
  
  return

#==================================================================================================================

def Fig11():
  """Makes panels for Fig. 11."""
  
  #------------------------------------------------
  
  def Panel(Mstar):
    """Plots a panel from Fig. 11."""
    
    # Make stars
    starSlow = mors.Star( Mstar=Mstar , percentile='slow' )
    starMedium = mors.Star( Mstar=Mstar , percentile='medium' )
    starFast = mors.Star( Mstar=Mstar , percentile='fast' )
    
    # Make plot
    plt.figure()
    
    plt.plot( starSlow.AgeTrack , starSlow.LxTrack , 'r' )    
    plt.plot( starMedium.AgeTrack , starMedium.LxTrack , 'g' )
    plt.plot( starFast.AgeTrack , starFast.LxTrack , 'b' )
    
    plt.xlabel("Age (Myr)")
    plt.ylabel("X-ray luminosity (erg s$^{-1}$)")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1.0,5000.0])
    plt.ylim([1.0e27,1.0e31])
    plt.savefig("Fig11_"+str(Mstar)+"Msun.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Panel for each mass bin
  Panel(1.0)
  Panel(0.75)
  Panel(0.5)
  Panel(0.25)
  
  return

#==================================================================================================================

def Fig12():
  
  #------------------------------------------------
  
  def LeftPanel():
    
    # Get saturation time for each star (assume Lx arbitrarily)
    AgeActive28 = cluster.ActivityLifetime(Quantity='Lx',Threshold=1.0e28)
    AgeActive29 = cluster.ActivityLifetime(Quantity='Lx',Threshold=1.0e29)
    AgeActive30 = cluster.ActivityLifetime(Quantity='Lx',Threshold=1.0e30)
    
    # Make plot
    plt.figure()
    plt.scatter(cluster.Mstar,AgeActive28,s=5,color='blue')
    plt.scatter(cluster.Mstar,AgeActive29,s=5,color='cyan')
    plt.scatter(cluster.Mstar,AgeActive30,s=5,color='orange')
    
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Active lifetime for $L_\mathrm{X}$ (Myr)")
    plt.yscale('log')
    plt.ylim([1.0,5000.0])
    plt.savefig("Fig12_leftpanel.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def RightPanel():
    
    # Get saturation time for each star (assume Lx arbitrarily)
    AgeActive10 = cluster.ActivityLifetime(Quantity='FxHZ',Threshold=10.0)
    AgeActive100 = cluster.ActivityLifetime(Quantity='FxHZ',Threshold=100.0)
    
    # Make plot
    plt.figure()
    plt.scatter(cluster.Mstar,AgeActive10,s=5,color='blue')
    plt.scatter(cluster.Mstar,AgeActive100,s=5,color='cyan')
    
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Active lifetime for $F_\mathrm{X}$ in HZ (Myr)")
    plt.yscale('log')
    plt.ylim([1.0,5000.0])
    plt.savefig("Fig12_rightpanel.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Restore or calculate cluster
  if os.path.isfile('cluster.pickle'):
    cluster = mors.Load('cluster.pickle')
  else:
    Mstar , Omega = mors.ModelCluster()
    cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)
    cluster.save()
  
  # Make panels
  LeftPanel()
  RightPanel()
  
  return

#==================================================================================================================

def Fig13():
  print("Fig. 13 not included here.")
  return

#==================================================================================================================

def Fig14():
  print("Fig. 14 not included here.")
  return

#==================================================================================================================

def Fig15():
  print("Fig. 15 not included here.")
  return

#==================================================================================================================

def Fig16():
  
  #------------------------------------------------
  
  def LeftColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for iAge in range(len(Ages)):
      
      # Loop over stars and sort into five groups
      color = []
      for iStar in range(cluster.nStars):
        if ( Fx[iAge,iStar] > Feuv[iAge,iStar] + Fly[iAge,iStar] ):
          color.append('blue')
        elif ( Fx[iAge,iStar] > Feuv[iAge,iStar] ):
          color.append('cyan')
        elif ( Feuv[iAge,iStar] > Fx[iAge,iStar] ) and ( Feuv[iAge,iStar] > Fly[iAge,iStar] ):
          color.append('green')
        elif ( Fly[iAge,iStar] > Fx[iAge,iStar] + Feuv[iAge,iStar] ):
          color.append('red')
        elif ( Fly[iAge,iStar] > Feuv[iAge,iStar] ):
          color.append('orange')
          
      # Start subplot
      plt.subplot(len(Ages),1,iAge+1)
      
      # Plot stars
      plt.scatter(cluster.Mstar,Fx[iAge,:],s=5,color=color)
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_\odot$)")
      plt.ylabel("$F_\mathrm{X}$ (erg s$^{-1}$ cm$^{-2}$)")
      plt.xlim([0.08,1.22])
      plt.ylim([1.0e4,2.0e8])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig16_leftcolumn.png")
    plt.close()
    
    
    return
  
  #------------------------------------------------
  
  def MiddleColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for iAge in range(len(Ages)):
      
      # Start subplot
      plt.subplot(len(Ages),1,iAge+1)
      
      # Plot stars
      plt.scatter(cluster.Mstar,Feuv[iAge,:]/Fx[iAge,:],s=5)
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_\odot$)")
      plt.ylabel("EUV to X-ray ratio")
      plt.xlim([0.08,1.22])
      plt.ylim([0.3,15.0])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig16_middlecolumn.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def RightColumn():
    
    # Start figure
    plt.figure(figsize=(7,5*len(Ages)))
    
    # Make each panel
    for iAge in range(len(Ages)):
      
      # Start subplot
      plt.subplot(len(Ages),1,iAge+1)
      
      # Plot stars
      plt.scatter(cluster.Mstar,Fly[iAge,:]/Fx[iAge,:],s=5)
      
      # Other stuff
      plt.yscale('log')
      plt.xlabel("Stellar Mass (M$_\odot$)")
      plt.ylabel(r"Ly-$\alpha$ to X-ray ratio")
      plt.xlim([0.08,1.22])
      plt.ylim([0.1,30.0])
    
    # Finalise
    plt.tight_layout()
    plt.savefig("Fig16_rightcolumn.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Set ages to plot
  Ages = np.array([2.0,100.0,1000.0,5000.0])
  
  # Restore or calculate cluster
  if os.path.isfile('cluster.pickle'):
    cluster = mors.Load('cluster.pickle')
  else:
    Mstar , Omega = mors.ModelCluster()
    cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)
    cluster.save()
  
  # Need to get at each age scattered values for each quantity
  # Start by looping over ages
  Fx = np.zeros((len(Ages),cluster.nStars))
  Feuv = np.zeros((len(Ages),cluster.nStars))
  Fly = np.zeros((len(Ages),cluster.nStars))
  for iAge in range(len(Ages)):
    
    # Get scatter values for all
    for iStar in range(cluster.nStars):
      
      Lxuv = mors.Lxuv(
        Mstar=cluster.stars[iStar].Mstar,
        Age=Ages[iAge],
        Omega=cluster.stars[iStar].Value(Age=Ages[iAge],Quantity='OmegaEnv')
        )
      deltaXUV = mors.XUVScatter(Lxuv)
      
      Fx[iAge,iStar] = Lxuv['Fx'] + deltaXUV['Fx']
      Feuv[iAge,iStar] = Lxuv['Feuv'] + deltaXUV['Feuv']
      Fly[iAge,iStar] = Lxuv['Fly'] + deltaXUV['Fly']
  
  # Plot
  LeftColumn()
  MiddleColumn()
  RightColumn()
  
  return

#==================================================================================================================

def Fig17():
  
  #------------------------------------------------
  
  def MainPanel():
    
    # Masses to include
    Mstar = np.array([0.1,0.3,0.6,0.8,1.0,1.2])
    
    # Ages to use
    AgeMin = 1.0
    AgeMax = 5000.0
    nAge = 100
    Age = np.logspace(np.log10(AgeMin),np.log10(AgeMax),nAge)
    
    # Make array to hold data
    LbolNorm = np.zeros((len(Mstar),nAge))
    for i in range(len(Mstar)):
      LbolNorm[i,:] = mors.Lbol(Mstar[i],Age) / mors.Lbol(Mstar[i],5000.0)
    
    # Make plot
    plt.figure()
    for i in range(len(Mstar)):
      plt.plot( Age[:] , LbolNorm[i,:] )
    plt.xlabel("Age (Myr)")
    plt.ylabel("Normalised bolometric flux")
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig("Fig17_main.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def Inset():
    
    # Masses and HZ boundaries
    Mstar = np.linspace(0.1,1.2,100)
    aOrbHZ = mors.aOrbHZ(Mstar=Mstar)
    
    # Make plot
    plt.figure()
    
    plt.plot( Mstar , aOrbHZ['HZ'] , 'k--' )
    
    plt.plot( Mstar , aOrbHZ['RecentVenus'] , 'c' )
    plt.plot( Mstar , aOrbHZ['RunawayGreenhouse'] , 'c' )
    plt.plot( Mstar , aOrbHZ['MoistGreenhouse'] , 'r' )
    plt.plot( Mstar , aOrbHZ['MaximumGreenhouse'] , 'b' )
    plt.plot( Mstar , aOrbHZ['EarlyMars'] , 'c' )

    
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("HZ distance (AU)")
    plt.yscale('log')
    plt.savefig("Fig17_inset.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  MainPanel()
  Inset()
  
  return

#==================================================================================================================

def Fig18():
  
  #------------------------------------------------
  
  def UpperLeft():
    
    # Get emitted energies
    energies = cluster.IntegrateEmission(AgeMin=1.0,AgeMax=1000.0,Band='XUV')
    
    # Make plot
    plt.figure()
    plt.scatter( cluster.Mstar , energies , s=5 )
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Emitted X-ray and EUV energy (erg)")
    plt.yscale('log')
    plt.savefig("Fig18_upperleft.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def UpperRight():
    
    # Get emitted energies
    energies = cluster.IntegrateEmission(AgeMin=1000.0,AgeMax=5000.0,Band='XUV')
    
    # Make plot
    plt.figure()
    plt.scatter( cluster.Mstar , energies , s=5 )
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Emitted X-ray and EUV energy (erg)")
    plt.yscale('log')
    plt.savefig("Fig18_upperright.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def LowerLeft():
    
    # Get emitted energies
    energies = cluster.IntegrateEmission(AgeMin=1.0,AgeMax=1000.0,Band='XUV',aOrb='HZ')
    
    # Make plot
    plt.figure()
    plt.scatter( cluster.Mstar , energies , s=5 )
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Integrated X-ray+EUV flux in HZ (erg cm$^{-2}$)")
    plt.yscale('log')
    plt.savefig("Fig18_lowerleft.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  def LowerRight():
    
    # Get emitted energies
    energies = cluster.IntegrateEmission(AgeMin=1000.0,AgeMax=5000.0,Band='XUV',aOrb='HZ')
    
    # Make plot
    plt.figure()
    plt.scatter( cluster.Mstar , energies , s=5 )
    plt.xlabel("Stellar mass (M$_\odot$)")
    plt.ylabel("Integrated X-ray+EUV flux in HZ (erg cm$^{-2}$)")
    plt.yscale('log')
    plt.savefig("Fig18_lowerright.png")
    plt.close()
    
    return
  
  #------------------------------------------------
  
  # Restore or calculate cluster
  if os.path.isfile('cluster.pickle'):
    cluster = mors.Load('cluster.pickle')
  else:
    Mstar , Omega = mors.ModelCluster()
    cluster = mors.Cluster(Mstar=Mstar,Omega=Omega,verbose=True)
    cluster.save()
  
  # Do plots
  UpperLeft()
  UpperRight()
  LowerLeft()
  LowerRight()
  
  return

#==================================================================================================================

def Fig19():
  print("Fig. 19 not included here.")
  return

#==================================================================================================================

def Fig20():
  print("Fig. 20 not included here.")
  return

#==================================================================================================================
#==================================================================================================================
#==================================================================================================================
#==================================================================================================================

Main(sys.argv)
