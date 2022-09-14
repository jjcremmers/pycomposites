#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.6
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials,Laminate

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
    
print ("\nA matrix:\n",lam.getA())
print ("\nB matrix:\n",lam.getB())
print ("\nD matrix:\n",lam.getD())
 
