#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.6
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials,Laminate,stressTransformation

carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

print("The properties of carbon are:\n",carbon)
print("The properties of epoxy are:\n",epoxy)

udcomp = mixMaterials( carbon , epoxy , 0.6 )

lam = Laminate()

lam.addMaterial( 'UD' , udcomp )

orientations = [ 0. , 45. , -45. , 90 , 90 , -45. , 45. , 0. ]

for angle in orientations:
  lam.addLayer( 'UD' , angle , 0.25e-3 )

A1,B1,C1,D1 = lam.getInverseMatrices()

import numpy as np

N = np.array([1.0e4, 5.0e3, 6.0e3 ])
M = np.array([3.   , 10.  , 1.    ])

eps0  = np.dot( A1 , N ) + np.dot( B1 , M )
kappa = np.dot( C1 , N ) + np.dot( D1 , M )

sigmaplt = np.zeros(shape=(3,8))

for iLay,angle in enumerate(orientations):
  epsilon = eps0 + lam.getZcoord( iLay ) * kappa
  sigma   = np.dot( lam.getQbar(iLay) , epsilon )
  sigmaplt[:,iLay] = stressTransformation( sigma , angle )
  
import matplotlib.pyplot as plt

X = np.arange(8)
fig = plt.figure()

plt.bar(X - 0.26, sigmaplt[0], color = 'b', width = 0.25, label = "$\sigma_{11}$")
plt.bar(X       , sigmaplt[1], color = 'g', width = 0.25, label = "$\sigma_{22}$")
plt.bar(X + 0.26, sigmaplt[2], color = 'r', width = 0.25, label = "$\sigma_{12}$")

plt.xticks(X, orientations)

plt.xlabel("Layer orientations")
plt.ylabel("Stress [MPa]")
plt.title("Stress in each layer (12 coordinate system)")

plt.legend()
#plt.show()
plt.savefig('example4.png')
