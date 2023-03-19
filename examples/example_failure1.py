#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.11
#
#  (c) Joris Remmers, TU/e   2013 - 2017
#  
#  run with python 2.7
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials,Laminate
from numpy import array,dot,linspace,zeros,pi,cos,sin
from pylab import *

glassepoxy = TransverseIsotropic( [39.0e9,8.6e9],0.28,3.254e9)

glassepoxy.setFailureProperties( [1080e6,620e6,39e6,128e6,89e6] )

n = 1000

theta = linspace(0,0.5*pi,n)

splot = zeros(shape=(4,n))

sigma = zeros(3)
tiny  = 1.0e-6

for i,t in enumerate(theta):
  sigma[0] = cos( t )**2
  sigma[1] = sin( t )**2
  sigma[2]  = -sin( t )*cos( t )

  for j,ficrit in enumerate(["maxstress","maxstrain","tsaiwu"]):
    fi      = 0.0
    sf_low  = 0.0
    sf_high = 1.0e15
  
    while fi > 1.+tiny or fi < 1.0-tiny:
      sf = 0.5*(sf_low+sf_high)
    
      if ficrit == "maxstress":
        fi = glassepoxy.getFIMaximumStress( sf*sigma )        
      elif ficrit == "maxstrain":
        fi = glassepoxy.getFIMaximumStrain( sf*sigma )        
      elif ficrit == "tsaiwu":
        fi = glassepoxy.getFITsaiWu( sf*sigma )
           
      if fi > 1.0:
        sf_high = sf
      else:
        sf_low  = sf
  
    splot[j,i] = sf
  splot[3,i] = 180/pi*t
  

figure()
plot(splot[3,:] , splot[0,:] , 'r-.')
plot(splot[3,:] , splot[1,:] , 'b--')
plot(splot[3,:] , splot[2,:] , 'g')
xlabel('theta')
ylabel('sigma')
title('title')
show()

