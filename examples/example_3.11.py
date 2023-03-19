#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.11
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials,Laminate
from numpy     import array,dot

carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

compmat = mixMaterials( carbon , epoxy , 0.6 )

lam = Laminate()

lam.addMaterial( 'composite' , compmat )

total = []

for theta in range(91):
  lam.removeAllLayers()

  lam.addLayer( 'composite' , theta , 6e-3 )
  lam.addLayer( 'composite' , -theta , 6e-3 )
  lam.addLayer( 'composite' , -theta , 6e-3 )
  lam.addLayer( 'composite' , theta , 6e-3 )

  output = lam.getElastic()

  output.append(theta)

  total.append(output)

from pylab import plot, show, xlabel, ylabel

plot( [x[4] for x in total], [x[0] for x in total], 'r-' )
plot( [x[4] for x in total], [x[1] for x in total], 'b-' )
plot( [x[4] for x in total], [x[3] for x in total], 'g-' )
plot( [x[4] for x in total], [x[2]*5e10 for x in total], 'g-' )

show()
