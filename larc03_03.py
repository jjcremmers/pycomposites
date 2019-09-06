from composite import TransverseIsotropic
from composite import stressTransformation
from math import cos,pi
from numpy import zeros

as4 = TransverseIsotropic( [127600,11300],0.278,6000)

as4.setFailureProperties([1045.,844.,1044.,244.,244.],1.)
as4.setSLis( 95.1 )

theta = zeros(100)
f     = zeros(100)

sigma=zeros(3)
sigma[0]= -10.0

for i in range(90):
  thdeg = i
  theta[i] = thdeg*pi/180

  fi = 0.
  j = 1
  while fi < 1.:
    sigapp = stressTransformation( j*sigma , theta[i] )
 
    failindx = as4.getFILarc03( sigapp )

    print(failindx)
    fi = max(failindx)

    j = j+1

  print(thdeg,j)
