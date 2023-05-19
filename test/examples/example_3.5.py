#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis (4MM00)
#
#  Example 3.5
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic
from numpy     import array,dot

stress = array([ 1.e9 , 0.5e9 , 0.0 ])

steel      = TransverseIsotropic( 207e9,0.33,1.0)
boronEpoxy = TransverseIsotropic( [207e9,19e9],0.21,6.4e9)

print("The strain in the steel is       : ",dot ( steel.getS()      , stress ),'\n')
print("The strain in the boron-epoxy is : ",dot ( boronEpoxy.getS() , stress ),'\n')
print("The strain in the boron-epoxy is : ",dot ( boronEpoxy.getSbar(45.) , stress ),'\n')

print(steel)

