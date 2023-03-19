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
from numpy     import array,dot

carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

compmat = mixMaterials( carbon , epoxy , 0.6 )

print("Material properties of the composite material:\n\n",compmat,"\n")

print("Q matrix:\n\n",compmat.getQ(),"\n")
print("U[1..5]:\n\n",compmat.getU(),"\n")
print("Q matrix (45 degrees):\n\n",compmat.getQbar( 45. ),"\n")
print("Q matrix (-45 degrees):\n\n",compmat.getQbar( -45. ))

