#!/usr/bin/env python
# coding: utf-8

# In[234]:


from ReadFile import Read
import numpy as np
import astropy.units as u
import pandas as pd

def ComponentMass(filename, part_type):
    """
    ComponentMass is a function that returns the total mass of a user-specified galaxy component (disk, bulge, etc.)
    
    Inputs: filename (string) is the name of the file we want to read data from
            part_type (integer) is the type of galaxy component we want to return the mass of 
            (1 = halo, 2 = disk, 3 = bulge)
            
    Outputs: M (astropy units Msun*10^12) the mass of the galaxy component
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

'''
# In[236]:


#calling ComponentMass() for every particle type in the Milky Way (1 = halo, 2 = disk, 3 = bulge)
MW_halo = ComponentMass("MW_000.txt",1)
MW_disk = ComponentMass("MW_000.txt",2)
MW_bulge = ComponentMass("MW_000.txt",3)

#calling ComponentMass() for every particle type in M31
M31_halo = ComponentMass("M31_000.txt",1)
M31_disk = ComponentMass("M31_000.txt",2)
M31_bulge = ComponentMass("M31_000.txt",3)

#calling ComponentMass() for every particle type in M33. M33 does not have a bulge -- only a disk and a halo
M33_halo = ComponentMass("M33_000.txt",1)
M33_disk = ComponentMass("M33_000.txt",2)

#totaling the component masses for each galaxy to get the total mass of each galaxy 
MWtot = MW_halo + MW_disk + MW_bulge 
M31tot = M31_halo + M31_disk + M31_bulge 
M33tot = M33_halo + M33_disk + M33_bulge 

#computing the total mass of the Local Group (MW + M31 + M33)
MLocal = MWtot + M31tot + M33tot 

#calculating f_bar, the baryon fraction for each galaxy and the Local Group. f_bar = stellar mass/total mass. 
#f_bar = (bulge + disk)/total mass. 
MWf_bar = (MW_bulge + MW_disk)/MWtot
M31f_bar = (M31_bulge + M31_disk)/M31tot
M33f_bar = M33_disk/M33tot #M33 does not have a bulge, so only the disk makes up the stellar mass
Localf_bar = (MW_bulge + MW_disk + M31_bulge + M31_disk + M33_disk)/MLocal

#rounding f_bar values to 3 significant figures
MWf_bar = np.round(MWf_bar,3)
M31f_bar = np.round(M31f_bar,3)
M33f_bar = np.round(M33f_bar,3)
Localf_bar = np.round(Localf_bar,3)


# In[237]:


#creating the table using the generated quantities for each galaxy and the Local Group
#"{:.2e}".format() formats the number as an exponent/in scientific notation. 

data = {'GALAXY NAME':["MW", "M31", "M33", "Local Group"], 
        "HALO MASS (10^12 Msun)":[(MW_halo.value), (M31_halo.value),(M33_halo.value), "-"],
        "DISK MASS (10^12 Msun)":[(MW_disk.value),(M31_disk.value),(M33_disk.value),"-"],
        "BULGE MASS (10^12 Msun)":[(MW_bulge.value),(M31_bulge.value),"-","-"],
        "TOTAL MASS (10^12 Msun)":[(MWtot.value),(M31tot.value),(M33tot.value),(MLocal.value)],
        "F_BAR":[(MWf_bar),(M31f_bar.value),(M33f_bar.value),(Localf_bar.value)]}
        #data is a dictionary, the column name is the key. 

mass_table = pd.DataFrame(data) #setting the data in the table using 
print(mass_table) #testing 
#table wraps and appears in 2 rows, but it's the same table. 


# In[240]:


mass_table.to_html('HW3table.html') #export table to html, where it can be saved as PDF 


# In[224]:


#comparing stellar masses of MW and M31, for part 4 question 2: 
MW_stellar = (MW_disk + MW_bulge)
M31_stellar = (M31_disk + M31_bulge)

print(MW_stellar)
print(M31_stellar)

#comparing MW and M31 dark matter (halo) mass ratio for part 4 question 3:
dark_ratio = M31_halo/MW_halo

print(dark_ratio)


# In[ ]:
'''



