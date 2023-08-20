#!/usr/bin/env python
# coding: utf-8

# This is the Python code is based on the examples in Chapter 3 of the book:
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

# ## Example 2 (Classical Laminate Theory)
# 
# Consider a fibre reinforced plastic that consists of uni-directional carbon fibres embedded in an epoxy matrix. 
# The fibre volume fraction $V_f=0.6$. The properties of the transverse isoptric fibre are: $E_{\rm fL} = 220\,$GPa ; 
# $E_{\rm fT} = 20\,$GPa ; $\nu_{\rm f} = 0.2$ ; $G_{\rm f} = 91.7\,$GPa.
# 
# The properties of the isotropic epoxy matrix are $E_{\rm m} = 3.6\,$GPa ; $\nu_{\rm f} = 0.35$ ; 
# $G_{\rm f} = 1.33\,$GPa.
# 
# Determine the ${\bf Q}$ matrix of this material and the matrices ${\bf \bar Q}(20)$ and  ${\bf \bar Q}(-20)$. Evaluate the results.
# 

# ## Solution

# In[13]:


from composite import TransverseIsotropic,mixMaterials


# For this carbon fiber composite material, the carbon fibers and the epoxy matrix are modeled as separate transversely isotropic materials. The T-300 carbon fibers have Young's Moduli $E_1=220\,$GPa,  $E_2=22\,$GPa, a Poisson's ratio $\nu_{12}=0.2$ and a Shear Modulus $G_{12}=91.7\,$GPa. The epoxy matrix may be considered isotropic with Young's modulus $E=3.6\,$GPa, Poisson's ratio $\nu=0.35$ and Shear Modulus $G_{12}=1.33\,$GPa. 

# In[24]:


carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)


# The properties of carbon and epoxy are:

# In[20]:


print(carbon)
print(epoxy)


# The composite consists of 60\% fibres and 40\% epoxy. The properties of the composite can be calculated by simple homogensiation.

# In[21]:


udcomp = mixMaterials( carbon , epoxy , 0.6 )
    


# The properties of the UD material are:  

# In[22]:


print(udcomp) 


# The stiffness matrix ${\bf Q}$ and the rotated matrices can be calculated using the equations in the book by Mallick.

# In[23]:


Q = udcomp.getQ()
print("The Q matrix of the composite is:\n\n",Q,"\n")


# In[10]:


Qbar = udcomp.getQbar( 20.0 )

print("The Q matrix of the composite under an angle of 20 degrees is:\n\n",
         Qbar,"\n")


# In[12]:


Qbar = udcomp.getQbar( -20.0 )

print("The Q matrix of the composite under an angle of minus 20 degrees is:\n\n",
         Qbar,"\n")


# Note that the 11, 12, 21, 22 and 66 terms of the ${\bf Q}_{20}$ matrix are identical to the terms in the ${\bf Q}_{-20}$ matrix. The 16 and 26 terms of these matrices have switched signs.

# In[ ]:




