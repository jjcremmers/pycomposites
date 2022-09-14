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

print("The properties of carbon are:\n",carbon)
print("The properties of epoxy are:\n",epoxy)

udcomp = mixMaterials( carbon , epoxy , 0.6 )

print("Material properties of the composite material are:\n",udcomp)
