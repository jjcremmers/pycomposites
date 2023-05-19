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

glassepoxy.setFailureProperties( [1080e6,620e6,39e6,128e6,89e6] , 30e6 )
glassepoxy.setSLis(20e6)
print(glassepoxy)

n     = 1000;
phi   = linspace(0,2*pi,n);
sigma = zeros(3)

sig1plot = zeros(1000)
sig2plot = zeros(1000)

for i,p in enumerate(phi):
  sigma[0] = cos( p );
  sigma[1] = sin( p );
  
  tiny    = 1.0e-6
  fi      = 0.0
  sf_low  = 0.0
  sf_high = 1.0e15
  
  while fi > 1.+tiny or fi < 1.0-tiny:
    sf = 0.5*(sf_low+sf_high)
    
    fi = glassepoxy.getFILarc03( sf*sigma )
    
    if fi > 1.0:
      sf_high = sf
    else:
      sf_low  = sf
  
  sig1plot[i] = sf*sigma[0];
  sig2plot[i] = sf*sigma[1];

figure()
plot(sig1plot, sig2plot, 'r')
xlabel('sigma1')
ylabel('sigma2')
title('title')
show()

