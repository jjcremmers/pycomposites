#!/usr/bin/env python
# coding: utf-8

# This is the Python code is based on Example 3.6 of the book:
# 
# FIBER-REINFORCED COMPOSITES
# Materials, Manufacturing, and Design
# by: P.K. Mallick (2008) by Taylor & Francis Group, LLC
# 
# and is discussed during the lecture on classical laminate theory in the course:
# 
# Composite and Lightweight Materials (4MM00)
# 
# at Eindhoven University of Technology
# 
# This code:
# (C) Joris Remmers (2013-2023)

# ## Example 1 (Classical Laminate Theory)
# 
# Consider a fibre reinforced plastic that consists of uni-directional carbon fibres embedded in an epoxy matrix. The fibre volume fraction $V_f=0.6$. The properties of the transverse isoptric fibre are: $E_{\rm fL} = 220\,$GPa ; $E_{\rm fT} = 20\,$GPa ; $\nu_{\rm f} = 0.2$ ; $G_{\rm f} = 91.7\,$GPa.
# 
# The properties of the isotropic epoxy matrix are $E_{\rm m} = 3.6\,$GPa ; $\nu_{\rm f} = 0.35$ ; $G_{\rm f} = 1.33\,$GPa.
# 
# Determine the homogenised properties of the composite.
# 
# 

# ## Solution
# 
# Import the functions TransverseIstropic and mixMaterials from the composite module.

# In[1]:


from composite import TransverseIsotropic,mixMaterials


# Create two materials, carbon and epoxy with the correct properties. Note that carbon is transversely isotropic and epoxy is isotropic.

# In[8]:


carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)


# In[9]:


print("The properties of carbon are:\n",carbon)
print("The properties of epoxy are:\n",epoxy)


# The properties of the composites can be calculated using the mixMaterials function, which has the two materials and the volume faction as an input.

# In[10]:


udcomp = mixMaterials( carbon , epoxy , 0.6 )

print("Material properties of the composite material:\n\n",udcomp,"\n")

