from composite import TransverseIsotropic
from composite import TransverseIsotropic
from math import cos
from numpy import zeros

as4 = TransverseIsotropic( [127600,11300],0.278,6000)

as4.setFailureProperties([1045.,844.,1044.,244.,244.],1.)

theta = 0.
sigma=zeros(3)
sigma[0]= -1.0
sigma[1] =0.5
print(sigma)

as4.setSLis( 95.1 )

fi = as4.getFILarc03( sigma )
