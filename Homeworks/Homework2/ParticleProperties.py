#!/usr/bin/env python
# coding: utf-8

# In[24]:


from ReadFile import Read
import numpy as np
import astropy.units as u

def ParticleInfo(filename, part_type, part_num):
    """
    ParticleInfo is a function that returns the distance, velocity, 
    and mass of a particle of any given type from a file. 

    Inputs: filename is the name of the file we want to read data from
            part_type is the type of particle we want to return the properties of 
            part_num is the number of the partcle we want to return the properties of

    Outputs: Mag of distance (astropy units kpc)
             Mag of velocity (astropy units km/s)
             Mass (astropy units of M_sun)
    """

    time,tot_particles,data = Read(filename) #calling Read(), a function we made previously to read the file
    
    index = np.where(data['type'] == part_type) #creating index to separate particle type specified by username

    #splitting data into relevant lists of masses, positions (x,y,z) and velocities (vx,vy,vz) and providing appropriate units
    #indexing the data so only the particles of the user-specified type are included in list

    #mass in 10**(10) M_sun
    m = data['m'][index]*u.M_sun*10**(10)

    #cartesian coordinates, in kpc
    x = data['x'][index]*u.kpc
    y = data['y'][index]*u.kpc
    z = data['z'][index]*u.kpc

    #velocity components, in km/s
    vx = data['vx'][index]*(u.km/u.s)
    vy = data['vy'][index]*(u.km/u.s)
    vz = data['vz'][index]*(u.km/u.s)

    #calculate magnitude of density (dist = sqrt(x**2 + y**2 + z**2)
    distance = np.sqrt(x[part_num-1]**2 + y[part_num-1]**2 + z[part_num-1]**2)
    
    #calculate magnitude of velocity (vel = vx**2 + vy**2 + vz**2)
    velocity = np.sqrt(vx[part_num-1]**2 + vy[part_num-1]**2 + vz[part_num-1]**2)

    distance = np.around(distance,3) #rounding distance to 3 sig figs
    velocity = np.around(velocity,3) #rounding velocity to 3 sig figs
    
    mass = m[part_num-1] #getting the particle mass
    
    return distance, velocity, mass #returning the distance, velocity, and the mass of the user-given particle
    
distance, velocity, mass = ParticleInfo("MW_000.txt",2,100) #calling function, giving it the file, the type of particle we want to calculate for (2 = disk in this case), and the number of the particle (100 in this case)

