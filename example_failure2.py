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

phi    = linspace(0,2*pi,10000);
sigma = zeros(3)

sig1plot = zeros(10000)
sig2plot = zeros(10000)

for i,p in enumerate(phi):
  sigma[0] = cos( p );
  sigma[1] = sin( p );
  
  sf,tt    = glassepoxy.getSfConservative( sigma )

  sig1plot[i] = sf*sigma[0];
  sig2plot[i] = sf*sigma[1];

figure()
plot(sig1plot, sig2plot, 'r')
xlabel('sigma1')
ylabel('sigma2')
title('title')
show()

