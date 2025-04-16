#!/usr/bin/env python
# coding: utf-8

# In[234]:


from ReadFile import Read
import numpy as np
import astropy.units as u
import pandas as pd

def ComponentMass(mass_array):
    """
    ComponentMass is a function that returns the total mass of a user-specified galaxy component (disk, bulge, etc.)
    
    PARAMETERS
    ----------
       mass_array : 'np.ndarray' 
            an array that contains the masses of every particle in the 
            combined MW-M31 halo. 
    OUTPUTS
    -------
        M : 'float' the mass of the galaxy component in 10*10 M_Sun
    """
    time,tot_particles,data = Read(filename) #calling Read(), a function we made previously to read the file. the snapnumber, 
    #total number of particles, and the rest of the data array in the file are the outputs of Read().
    
    index = np.where(data['type'] == part_type) #creating index to separate the particles of the user-specified type out from the rest
    
    m = data['m'][index]*u.M_sun*10**(10) #getting the mass of each particle in the specified type, setting its units (10**10 M_sun) 
    #and storing as an array. 

    m = m/(10**(12)) #dividing by 10**12 to get m in 10**12 M_suns
    
    mtot = np.sum(m) #summing over all the masses in the mass array to get the total mass (mtot)

    mtot_round = np.round(mtot,3) #rounding the total mass to 3 decimal places

    return mtot_round #returning the rounded total mass



