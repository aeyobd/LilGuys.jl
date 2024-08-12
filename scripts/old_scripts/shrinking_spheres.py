#!/usr/bin/env python3
# rerrani.github.io
# snippet of code that runs the shrinking spheres method
# on a list of gadget2 snapshots to trace the motion of
# the centre of density
import numpy as np
from glob import glob
import h5py
import re


def shrinksphere(xarray,yarray,zarray,vxarray, vyarray, vzarray, x0,y0,z0, r0, frac, Npart):
  # particle coordinate/velocity arrays x,y,z ...
  # initial center guess x0,y0,z0
  # initial radius to compute COM
  # fraction to decrease radius every iteration
  # min number of particles inside radius
  
  x = np.copy(xarray)
  y = np.copy(yarray)
  z = np.copy(zarray)
  vx = np.copy(vxarray)
  vy = np.copy(vyarray)
  vz = np.copy(vzarray)
  

  while True:
    r = np.sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )
    
    okIdx = np.where(r < r0)[0]
    if len(okIdx) < Npart: break

    x = x[okIdx]
    y = y[okIdx]
    z = z[okIdx]
    vx = vx[okIdx]
    vy = vy[okIdx]
    vz = vz[okIdx]
    
    x0 = np.average( x )
    y0 = np.average( y )
    z0 = np.average( z )
    
    r0 *= frac
  
  vx0 = np.average( vx )
  vy0 = np.average( vy )
  vz0 = np.average( vz )
  
  return x0,y0,z0, vx0,vy0,vz0
  
  
  

# run shrinking spheres on all snapshots in folder
if __name__ == "__main__": 

  G = 1.
  mU= 1.0e10   # Msol
  rU= 1.0      # kpc
  tU= 4.714e-3 # Gyrs
  vU= 207.4    # km/s


  # Track potential minimum
  orbitfile = open("rapha_ss_centres.csv", "w")

  print ("# time        x       y       z      ", file=orbitfile )
  filenames = np.array(glob("snapshot_*.hdf5"))
  print(filenames)
  index = [re.findall(r"\d+", f)[0] for f in filenames]
  filenames = filenames[np.argsort(index)]
  print(index)
  print(filenames)

  for filename in filenames:
    if filename[-4:] == "hdf5":
      f = h5py.File(filename, 'r')
      
      print (" ... working with file ", filename)
      
      header = f[u'Header']
      T = header.attrs["Time"]
      N = header.attrs["NumPart_ThisFile"][1]
      
      x = f['PartType1/Coordinates'][:,0]
      y = f['PartType1/Coordinates'][:,1]
      z = f['PartType1/Coordinates'][:,2]

      vx = f['PartType1/Velocities'][:,0]
      vy = f['PartType1/Velocities'][:,1]
      vz = f['PartType1/Velocities'][:,2]
      
      minIdx = np.argmin(f['PartType1/Potential'])
      xMin = x[minIdx]
      yMin = y[minIdx]
      zMin = z[minIdx]
      
      # new Center Of Density
      xCOD, yCOD, zCOD, vxCOD, vyCOD, vzCOD = shrinksphere(x,y,z,vx,vy,vz, xMin,yMin,zMin, 0.01/rU, 0.9, 50)

      #      1     2        3        4
      print (f"{T*tU},{xCOD*rU},{yCOD*rU},{zCOD*rU}", file=orbitfile )
    
      f.close()
      
  orbitfile.close()
