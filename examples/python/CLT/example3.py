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

# ## Example 3 (Classical Laminate Theory)
# 
# Consider a fibre reinforced plastic that consists of uni-directional carbon fibres embedded in an epoxy matrix. 
# The fibre volume fraction $V_f=0.6$. The properties of the transverse isoptric fibre are: $E_{\rm fL} = 220\,$GPa ; 
# $E_{\rm fT} = 20\,$GPa ; $\nu_{\rm f} = 0.2$ ; $G_{\rm f} = 91.7\,$GPa.
# 
# The properties of the isotropic epoxy matrix are $E_{\rm m} = 3.6\,$GPa ; $\nu_{\rm f} = 0.35$ ; 
# $G_{\rm f} = 1.33\,$GPa.
# 
# The resulting UD composite material is used in a laminate consisting of 8 layers with a
# thickness $0.3\,$mm and stacking sequence $\lbrack 0,45,-45,90\rbrack_{\rm S}$. Calculate the ABD matrices of this laminate and evaluate the results.

# ## Solution

# In[29]:


from composite import TransverseIsotropic,mixMaterials


# For this carbon fiber composite material, the carbon fibers and the epoxy matrix are modeled as separate transversely isotropic materials. The T-300 carbon fibers have Young's Moduli $E_1=220\,$GPa,  $E_2=22\,$GPa, a Poisson's ratio $\nu_{12}=0.2$ and a Shear Modulus $G_{12}=91.7\,$GPa. The epoxy matrix may be considered isotropic with Young's modulus $E=3.6\,$GPa, Poisson's ratio $\nu=0.35$ and Shear Modulus $G_{12}=1.33\,$GPa. 

# In[30]:


carbon = TransverseIsotropic( [220e9,220e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)


# The properties of carbon and epoxy are:

# In[31]:


print(carbon)
print(epoxy)


# The composite consists of 60\% fibres and 40\% epoxy. The properties of the composite can be calculated by simple homogensiation.

# In[32]:


udcomp = mixMaterials( carbon , epoxy , 0.6 )
    


# Construct the laminate by adding the material and the different layers and calculate the $\bf{A},\bf{B}$ and ${\bf D}$ matrices.

# In[34]:


from composite import Laminate
lam = Laminate()

lam.addMaterial( 'UD' , udcomp )

orientations = [ 0. , 45. , -45. , 90 , 90 , -45. , 45. , 0. ]


# In[35]:


for angle in orientations:
  lam.addLayer( 'UD' , angle , 0.25e-3 )
    
print ("\nA matrix:\n",lam.getA())
print ("\nB matrix:\n",lam.getB())
print ("\nD matrix:\n",lam.getD())


# Note that for this stacking sequence, the B matrix is identical to zero (except some numerical noise).

# In[ ]:




