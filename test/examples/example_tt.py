#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.13
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite    import TransverseIsotropic,mixMaterials,Laminate
from numpy        import array,dot,zeros
from numpy.linalg import inv

carbon = TransverseIsotropic( [4.59e+11,1.44e+11],0.22,1.28e+11,[-3.6e-07,4.0e-06] )

lam = Laminate()

lam.addMaterial( 'carbon' , carbon )

lam.addLayer( 'carbon' , 0 , 6e-3 )
lam.addLayer( 'carbon' , 60 , 6e-3 )
lam.addLayer( 'carbon' , -60 , 6e-3 )
#lam.addLayer( 'carbon' , -45 , 6e-3 )

print(lam)

A = lam.getA()
B = lam.getB()
D = lam.getD()
Ts = lam.getTs()
print(A)
print("Ts",Ts)
print(lam.getElastic())
print("CTE",dot(inv(A),Ts))

