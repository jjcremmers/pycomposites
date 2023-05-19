#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.11
#
#  (c) Joris Remmers, TU/e   2013 - 2017
#  
#  run with python 2.7
#------------------------------------------------------------------------------

from composite import TransverseIsotropic
from numpy import array,dot,linspace,zeros,pi,cos,sin
from pylab import *

glassepoxy = TransverseIsotropic( [39.0e9,8.6e9],0.28,3.254e9)

glassepoxy.setFailureProperties( [1080e6,620e6,39e6,128e6,89e6] )

theta = linspace(0,0.5*pi,10)
sf    = zeros(10)
sigmax = 1
sigma = zeros(3)


for i,t in enumerate(theta):
  sigma[0] = cos( t )**2
  sigma[1] = sin( t )**2
  sigma[2]  = -sin( t )*cos( t )

  sf[i] = glassepoxy.getSfTsaiWu( sigma )
  
  #print(sf[i])

figure()
plot(theta, sf, 'r')
xlabel('theta')
ylabel('sf')
title('title')
show()

