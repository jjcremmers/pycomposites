from composite    import TransverseIsotropic,Laminate
from numpy        import array,dot,zeros,linspace
from numpy.linalg import inv
from math         import sin,pi,sqrt
import numpy as np

from   mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from   matplotlib import cm


# Create a new material compUD with the correct input parameters

# In[ ]:


compUD = TransverseIsotropic( [130e9,7.2e9],0.337,4.2e9,[0.57e-6,35.1e-6],1514.)

print(compUD)

lam = Laminate()

lam.addMaterial( 'UCHSC200_SE84' , compUD )

lam.addLayer( 'UCHSC200_SE84' , 0  , 0.2e-3 )
lam.addLayer( 'UCHSC200_SE84' , 90 , 0.2e-3 )
lam.addLayer( 'UCHSC200_SE84' , 0  , 0.2e-3 )
lam.addLayer( 'UCHSC200_SE84' , 0  , 0.2e-3 )
lam.addLayer( 'UCHSC200_SE84' , 90 , 0.2e-3 )
lam.addLayer( 'UCHSC200_SE84' , 0  , 0.2e-3 )

print(lam)

D = lam.getD()

D1 = D[0,0]
D2 = D[1,1]
D3 = D[0,1]+2*D[2,2]

print("D1: ",D1," D2: ",D2," D3: ",D3)

F = 1200

a = 0.4
b = 0.6

x0 = 0.2
y0 = 0.55

nmax = 5

B = np.zeros((nmax+1,nmax+1))
A = np.zeros((nmax+1,nmax+1))

for m in range(1,nmax):
  for n in range(1,nmax):
    B[m,n] = 4*F*sin(m*pi*x0/a)*sin(n*pi*y0/b)/(a*b)
    A[m,n] = B[m,n]/(pi**4*(D1*(m/a)**4+2*D3*(m/a)**2*(n/b)**2+D2*(n/b)**4))

nplot = 5

xplot = linspace(0,a,nplot)
yplot = linspace(0,b,nplot)
    
q = np.zeros((nplot,nplot))
w = np.zeros((nplot,nplot))

x, y = np.meshgrid(xplot,yplot)
     
for m in range(1,nmax):
  for n in range(1,nmax):
    q += B[m,n]*sin(m*pi*x/a)*sin(n*pi*y/b)
    w += A[m,n]*sin(m*pi*x/a)*sin(n*pi*y/b)


fig = plt.figure()
ax = fig.gca(projection='3d')                


print (x)
print(xplot,yplot)

plt.title('Point load q(x,y)')
ax.plot_surface(x,y, q , alpha=0.8, cmap=cm.coolwarm)

ax.set_xlabel('x')
#ax.set_xlim(0, a)
ax.set_ylabel('y')
#ax.set_ylim(0, b)
ax.set_zlabel('q')

plt.savefig("Point_load.png")


fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')


plt.title('Vertical displacement w(x,y)')

ax1.plot_surface(x, y, w, alpha=0.8, cmap=cm.coolwarm)
cset = ax1.contour(x, y, w, zdir='z', offset=-0.05, cmap=cm.coolwarm)

ax1.set_xlabel('x')
#ax1.set_xlim(0, a)
ax1.set_ylabel('y')
#ax1.set_ylim(b, 0)
ax1.set_zlabel('w')


plt.savefig("Vertical_displacement.png")

plt.show()



