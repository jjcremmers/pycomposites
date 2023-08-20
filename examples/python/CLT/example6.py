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

# ## Example 6 (Classical Laminate Theory)
# 
# Consider a fibre reinforced plastic that consists of uni-directional carbon fibres embedded in an epoxy matrix. 
# The fibre volume fraction $V_f=0.6$. The properties of the transverse isoptric fibre are: $E_{\rm fL} = 220\,$GPa ; 
# $E_{\rm fT} = 20\,$GPa ; $\nu_{\rm f} = 0.2$ ; $G_{\rm f} = 91.7\,$GPa.
# 
# The properties of the isotropic epoxy matrix are $E_{\rm m} = 3.6\,$GPa ; $\nu_{\rm f} = 0.35$ ; 
# $G_{\rm f} = 1.33\,$GPa.
# 
# The resulting UD composite material is used in a two-layer unsymmetric $\lbrack 0 / 90\rbrack$ laminate where each layer has a thickness of $6\,$mm. Determine the curvatures of of this laminate when it is cooled down with a temperature drop $-100^{\circ}$\,C. Use $\alpha_{11}=-0.5\cdot 10^{-6}$ and $\alpha_{22}=12\cdot 10^{-6}$\,m/m per $^{\circ}$C as the thermal expansion coefficient of a single layer.

# ## Solution

# In[7]:


from composite import TransverseIsotropic,mixMaterials,Laminate


# For this carbon fiber composite material, the carbon fibers and the epoxy matrix are modeled as separate transversely isotropic materials. The T-300 carbon fibers have Young's Moduli $E_1=220\,$GPa,  $E_2=22\,$GPa, a Poisson's ratio $\nu_{12}=0.2$ and a Shear Modulus $G_{12}=91.7\,$GPa. The epoxy matrix may be considered isotropic with Young's modulus $E=3.6\,$GPa, Poisson's ratio $\nu=0.35$ and Shear Modulus $G_{12}=1.33\,$GPa. 

# In[8]:


carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)


# The properties of carbon and epoxy are:

# In[9]:


print(carbon)
print(epoxy)


# The composite consists of 60\% fibres and 40\% epoxy. The properties of the composite can be calculated by simple homogensiation.

# In[10]:


udcomp = mixMaterials( carbon , epoxy , 0.6 )
    


# Store the thermal expension coefficient of the material

# In[11]:


udcomp.setAlpha( [-5e-6 , 12e-6] )
print("The UD material properties are:",udcomp)


# In[12]:


lam = Laminate()

lam.addMaterial( 'UD' , udcomp )

lam.addLayer( 'UD' ,  0.0 , 6e-3 )
lam.addLayer( 'UD' , 90.0 , 6e-3 )


# The total system of equations is:
#   \begin{equation}
#   \begin{bmatrix}
#   {\bf N} \\
#   {\bf M}
#   \end{bmatrix}
#   =
#   \begin{bmatrix}
#   {\bf A} & {\bf B} \\
#   {\bf B} & {\bf D} \\  
#   \end{bmatrix}
#   \begin{bmatrix}
#   {\bf \epsilon} \\
#   {\bf \kappa}
#   \end{bmatrix}-
#     \begin{bmatrix}
#   {\bf T}^* \\
#   {\bf T}^{**}
#   \end{bmatrix}\Delta T
#   \end{equation}
#   By absence of external loads, this can be written as:
#   \begin{equation}
#     \begin{bmatrix}
#       {\bf A} & {\bf B} \\
#       {\bf B} & {\bf D} \\  
#     \end{bmatrix}
#     \begin{bmatrix}
#       {\bf \epsilon} \\
#       {\bf \kappa}
#     \end{bmatrix}
#     =
#     \begin{bmatrix}
#       {\bf T}^* \\
#       {\bf T}^{**}
#     \end{bmatrix}\Delta T
#   \end{equation}
#   Solving for ${\bf \epsilon}$ and ${\bf \kappa}$ gives:
#   \begin{equation}
#     \begin{bmatrix}
#       {\bf \epsilon} \\
#       {\bf \kappa}
#     \end{bmatrix}
#     =
#     \begin{bmatrix}
#       {\bf A}_1 & {\bf B}_1 \\
#       {\bf C}_1 & {\bf D}_1 \\  
#     \end{bmatrix}
#     \begin{bmatrix}
#       {\bf T}^* \\
#       {\bf T}^{**}
#     \end{bmatrix}\Delta T
#   \end{equation}
#   The curvature ${\bf \kappa}$ can be calculated as:
#   \begin{equation}
#     {\bf \kappa}= \left({\bf C}_1{\bf T}^{*} + {\bf D}_1 {\bf T}^{**}\right)  \Delta T
#   \end{equation}

# In[13]:


Ts  = lam.getTs()
Tss = lam.getTss()

print("Ts is  : ",Ts )
print("Tss is : ",Tss)


A1,B1,C1,D1 = lam.getInverseMatrices()

import numpy as np

deltaT = 100.0

kappa = ( np.dot( C1 , Ts ) + np.dot( D1 , Tss ) ) * deltaT

print("The curvature is: ",kappa)


# In[ ]:




