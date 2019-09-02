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

carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

compmat = mixMaterials( carbon , epoxy , 0.6 )

lam = Laminate()

lam.addMaterial( 'composite' , compmat )

lam.addLayer( 'composite' , -45 , 6e-3 )
lam.addLayer( 'composite' , 45 , 6e-3 )

A = lam.getA()
B = lam.getB()
D = lam.getD()

A1,B1,C1,D1 = lam.getInverseMatrices()

N = zeros(3)
N[0] = 1e5

eps0 = dot(A1,N)

print(eps0)

kappa = dot(C1,N)

print(kappa)

epsb = eps0-3e-3*kappa
epst = eps0+3e-3*kappa

print("eps bottom : ",epsb)
print("eps top    : ",epst)

stressb = dot(lam.getQbar(0),epsb)
stresst = dot(lam.getQbar(1),epst)

print(stressb,stresst)
