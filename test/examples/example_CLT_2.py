#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.6
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials

carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

udcomp = mixMaterials( carbon , epoxy , 0.6 )

Q = udcomp.getQ()

print("The Q matrix of the composite is:\n",Q,"\n")

Qbar = udcomp.getQbar( 20.0 )

print("The Q matrix of the composite under an angle of 20 degrees is:\n",
         Qbar,"\n")

Qbar = udcomp.getQbar( -20.0 )

print("The Q matrix of the composite under an angle of minus 20 degrees is:\n",
         Qbar,"\n")
