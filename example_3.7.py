#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.7
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

lam.addLayer( 'composite' , 45. , 6e-3 )
lam.addLayer( 'composite' , -45. , 6e-3 )

print ("Laminate [45/-45]\n")
print ("A matrix:\n\n",lam.getA())
print ("\n")
print ("B matrix:\n\n",lam.getB())
print ("\n")
print ("D matrix:\n\n",lam.getD())
print ("\n")

lam.addLayer( 'composite' , -45. , 6e-3 )
lam.addLayer( 'composite' , 45. , 6e-3 )

print ("Laminate [45/-45/-45/45]\n")
print ("A matrix:\n\n",lam.getA())
print ("\n")
print ("B matrix:\n\n",lam.getB())
print ("\n")
print ("D matrix:\n\n",lam.getD())
print ("\n")

lam.removeAllLayers()

lam.addLayer( 'composite' , 0. , 6e-3 )
lam.addLayer( 'composite' , 90. , 6e-3 )
lam.addLayer( 'composite' , 0. , 6e-3 )
lam.addLayer( 'composite' , 0. , 6e-3 )
lam.addLayer( 'composite' , 90. , 6e-3 )
lam.addLayer( 'composite' , 0. , 6e-3 )

print ("Laminate [0/90/0/0/90/0]\n")
print ("A matrix:\n\n",lam.getA())
print ("\n")
print ("B matrix:\n\n",lam.getB())
print ("\n")
print ("D matrix:\n\n",lam.getD())
print ("\n")

A1,B1,C1,D1 = lam.getInverseMatrices()

print (1./A1[0,0]/lam.thick,12./D1[0,0]/(lam.thick**3),1./A1[0,0]/lam.thick/(12./D1[0,0]/(lam.thick**3))*100.)	
