from material import TransverseIsotropic,mixMaterials

print "rt"

f = TransverseIsotropic( [5.],0.3,1.)
g = TransverseIsotropic( [5.,6.],0.3,2.)

h = mixMaterials(f,g,0.5)

print f.getQ()
print f.getU()

print f,g,h
