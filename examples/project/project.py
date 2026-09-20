from pycomposites import TransverseIsotropic,Laminate

carbon = TransverseIsotropic( [62.45e9,62.452e9],0.037,3.71e9)
foam   = TransverseIsotropic( [1.,1.],0.037,1.)

print(carbon)


lam = Laminate()

lam.addMaterial( 'C' , carbon )
lam.addMaterial( 'F' , foam )


lam.addLayer( 'C' ,  0.0 , 2e-3 )
lam.addLayer( 'F' ,  0.0 , 0.005 )
lam.addLayer( 'C' ,  0.0 , 2e-3 )


help(Laminate)

D  = lam.getD()
Tss = lam.getTss()

F = 2500
L = 0.25
b = 0.2

delta = F*L*L*L/(48*b*D[0,0])

print("D is  : ",D )
print("de is : ",delta)