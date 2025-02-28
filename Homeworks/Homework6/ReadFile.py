#!/usr/bin/env python
# coding: utf-8

# In[14]:


import numpy as np
import astropy.units as u

def Read(filename): #Read() is a function that takes the name of a file as an input, reads the first 2 lines of a file (the time and total number of particles) and retuns  the particle type, mass, x,y,z, vx,vy,vz as a data arrays
    file = open(filename,"r") #opening the file using read mode
    
    line1 = file.readline() #reading the first line of the file. This line contains the time.
    label,value = line1.split() #splitting the line into a list. The line contains both the label and the values of the time.
    time = float(value)*u.Myr #storing the time in terms of Myr 

    line2 = file.readline() #calling readline a second time reads the second line of the file. This is the number of particles
    label,value = line2.split() #splitting line into a list. This line contains a label and the total number of particles.
    tot_particles = float(value) #storing the particles

    file.close() #close the file

    data = np.genfromtxt(filename,dtype = None, names = True, skip_header = 4) #storing the remainder of the file.

    return time,tot_particles, data #returns the time, total number of particles, and the rest of the data from the file. 

#Read("MW_000.txt")
#print(data['type'][1])


# In[ ]:




