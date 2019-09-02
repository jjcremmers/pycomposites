#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.11
#
#  (c) Joris Remmers, TU/e   2013 - 2017
#  
#  run with python 2.7
#------------------------------------------------------------------------------

from material import TransverseIsotropic,mixMaterials,Laminate
from numpy import array,dot,linspace,zeros,pi,cos,sin
from pylab import *

glassepoxy = TransverseIsotropic( [39.0e9,8.6e9],0.28,3.254e9)

glassepoxy.setF1t( 1080e6 )
glassepoxy.setF1c( 620e6 )
glassepoxy.setF2t( 39e6 )
glassepoxy.setF2c( 128e6 )
glassepoxy.setF6( 89e6 )

theta = linspace(0,0.5*pi,1000)
sf    = zeros(1000)
sigmax = 1
sigma = zeros(3)

for i,t in enumerate(theta):
  sigma[0] = cos( t )**2
  sigma[1] = sin( t )**2
  sigma[2]  = -sin( t )*cos( t )

  sf[i] = glassepoxy.getSfTsaiWu( sigma )
  
  print sf[i]

figure()
plot(theta, sf, 'r')
xlabel('theta')
ylabel('sigma')
title('title')
show()

