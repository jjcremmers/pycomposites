#------------------------------------------------------------------------------
#  Composite and Lightweight Materials: Design and Analysis
#
#  Example 3.6
#
#  (c) Joris Remmers, TU/e   2013 - 2019
#  
#  run with python 3.x
#------------------------------------------------------------------------------

from composite import TransverseIsotropic,mixMaterials,Laminate,stressTransformation

carbon = TransverseIsotropic( [220e9,22e9],0.2,91.7e9)
epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)

print("The properties of carbon are:\n",carbon)
print("The properties of epoxy are:\n",epoxy)

udcomp = mixMaterials( carbon , epoxy , 0.6 )

lam = Laminate()

lam.addMaterial( 'UD' , udcomp )

total = []

for theta in range(91):
  lam.removeAllLayers()

  lam.addLayer( 'UD' ,  theta , 6e-3 )
  lam.addLayer( 'UD' , -theta , 6e-3 )
  lam.addLayer( 'UD' , -theta , 6e-3 )
  lam.addLayer( 'UD' ,  theta , 6e-3 )

  output = lam.getElastic()   #list of 4 values: Ex, Ey, Gxy and nuxy

  output.append(theta)        #the angle is added as 5th value

  total.append(output)    


import matplotlib.pyplot as plt

fig, ax1 = plt.subplots()

ax2 = ax1.twinx()

ax1.set_xlabel('layer angle $\\theta$')
ax1.set_ylabel('Stiffness [Pa]', color='b')
ax2.set_ylabel('Poisson ratio []', color='r')

ax1.plot( [x[4] for x in total], [x[0] for x in total], 'b-'  , label="$E_x$" )
ax1.plot( [x[4] for x in total], [x[1] for x in total], 'b--' , label="$E_y$" )
ax1.plot( [x[4] for x in total], [x[3] for x in total], 'b:'  , label="$G_{xy}$" )
ax2.plot( [x[4] for x in total], [x[2] for x in total], 'r-'  , label="$\\nu_{xy}$" )

ax1.legend(loc="center right")
ax2.legend(loc="upper right")

#plt.show()
plt.savefig('example5.png')
