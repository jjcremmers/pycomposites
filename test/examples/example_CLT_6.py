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

udcomp = mixMaterials( carbon , epoxy , 0.6 )

udcomp.setAlpha( [-5e-6 , 12e-6] )

lam = Laminate()

lam.addMaterial( 'UD' , udcomp )

lam.addLayer( 'UD' ,  0.0 , 6e-3 )
lam.addLayer( 'UD' , 90.0 , 6e-3 )

Ts  = lam.getTs()
Tss = lam.getTss()

print("Ts and Tss are : ",Ts,Tss)

A1,B1,C1,D1 = lam.getInverseMatrices()

import numpy as np

deltaT = 100.0

kappa = ( np.dot( C1 , Ts ) + np.dot( D1 , Tss ) ) * deltaT

print("The curvature is: ",kappa)
  