#==============================================================================
#  A set of libraries to perform Classical Laminate calculations for the Course
#
#  Composite and Lightweight Materials - Design and Analysis (4MM00)
#
# (c) Joris Remmers (2013-2019)
#==============================================================================

from numpy import zeros,ones,dot,transpose
from numpy.linalg import inv
from math import sin,cos,pi,sqrt,tan,atan

#------------------------------------------------------------------------------
#  class Transverse Isotropic.
#    constructor:
#    E     Young's modulus, can be a a list in case of an orthotropic material
#    nu12  Poisson ratio
#    G12   Shear modulus
#------------------------------------------------------------------------------

class TransverseIsotropic:

  def __init__( self , E , nu12 , G12 = 0. , alpha = 0. , rho = 0. ):

    if type(E) is list:
      if len(E) is 2:
        self.E1 = E[0]
        self.E2 = E[1]
      elif len(E) is 1:
        self.E1 = E[0]
        self.E2 = E[0]
      else:
        print('error')
    else:
      self.E1    = E
      self.E2    = E

    self.nu12  = nu12
    self.G12   = G12
    self.nu21  = self.E2/self.E1*self.nu12

    if type(alpha) is list:
      if len(alpha) is 2:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[1]
      elif len(alpha) is 1:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[0]
      else:
        print('error')
    else:
      self.alpha1    = alpha
      self.alpha2    = alpha

    self.rho = rho

#------------------------------------------------------------------------------
#  setFailureProperties
#    F is a vector with length 5 containing Xt,Xc,Yt,Yc,S
#------------------------------------------------------------------------------
  def setFailureProperties( self, F , Gfrac = 0 , alpha0deg = 53. ):

    if len(F) is 5:
      self.Xt = F[0]
      self.Xc = F[1]
      self.Yt = F[2]
      self.Yc = F[3]
      self.Sl = F[4]
    else:
      print('error')

    if type(Gfrac) is list:
      if len(Gfrac) is 2:
        self.GIc  = Gfrac[0]
        self.GIIc = Gfrac[1]
        self.g = self.GIc/self.GIIc
    else:
      self.GIc  = Gfrac
      self.GIIc = Gfrac
      self.g    = 1.0

    self.alpha0 = alpha0deg*pi/180

    self.cosa0  = cos( self.alpha0 )
    self.sina0  = sin( self.alpha0 )

    self.cos2a0 = cos( 2.0*self.alpha0 )
    self.tan2a0 = tan( 2.0*self.alpha0 )

    self.lam22 = 2.0*(1.0/self.E2-self.nu21*self.nu21/self.E1)
    self.lam44 = 1.0/self.G12

#------------------------------------------------------------------------------
#  setSLis
#    sets SLis parameter for LaRC03 model
#------------------------------------------------------------------------------

  def setSLis( self , SLis ):

    self.SLis = SLis

#------------------------------------------------------------------------------
#  print material properties
#------------------------------------------------------------------------------

  def __str__( self ):

    msg  = "  Elastic Properties:\n"
    msg += "  -----------------------------------------------------------\n"
    msg += "  E1     :  {:12.3e} , E2     :  {:12.3e} \n".format(self.E1,self.E2)
    msg += "  nu12   :  {:12.2f} , G12    :  {:12.3e} \n".format(self.nu12,self.G12)

    if self.rho > 0.:
      msg += "  rho    :  {:12.2f}\n".format(self.rho)

    if hasattr( self , "alpha1" ):
      msg += "\n  Thermal expansion coefficients:\n"
      msg += "  -----------------------------------------------------------\n"
      msg += "  alpha1 :  {:12.3e} , alpha2 :  {:12.3e} \n".format(self.alpha1,self.alpha2)

    if hasattr( self , "Xt" ):
      msg += "\n  Strengths and failure model parameters:\n"
      msg += "  -----------------------------------------------------------\n"
      msg += "  Xt     :  {:12.3e} , Xc     :  {:12.3e} \n".format(self.Xt,self.Xc)
      msg += "  Yt     :  {:12.3e} , Yc     :  {:12.3e} \n".format(self.Yt,self.Yc)
      msg += "  S      :  {:12.3e}\n".format(self.Sl)

    return msg

#------------------------------------------------------------------------------
#  getQ
#------------------------------------------------------------------------------

  def getQ( self ):

    if not hasattr( self , 'Q' ):
      self.Q = zeros( shape=(3,3) )

      self.Q[0,0] = self.E1/(1.-self.nu12*self.nu21)
      self.Q[0,1] = self.nu12*self.E2/(1.0-self.nu12*self.nu21)
      self.Q[1,1] = self.E2/(1.-self.nu12*self.nu21)
      self.Q[1,0] = self.Q[0,1]
      self.Q[2,2] = self.G12

    return self.Q

#------------------------------------------------------------------------------
#  getU
#------------------------------------------------------------------------------

  def getU( self ):

    if not hasattr( self , 'U' ):
      self.getQ()

      self.U = zeros(5)

      self.U[0] = 0.125*(3.*self.Q[0,0]+3.*self.Q[1,1]+2.*self.Q[0,1]+4.*self.Q[2,2])
      self.U[1] = 0.5*(self.Q[0,0]-self.Q[1,1])
      self.U[2] = 0.125*(self.Q[0,0]+self.Q[1,1]-2.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[3] = 0.125*(self.Q[0,0]+self.Q[1,1]+6.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[4] = 0.5*(self.U[0]-self.U[3])

    return self.U

#------------------------------------------------------------------------------
#  getS
#------------------------------------------------------------------------------

  def getS( self ):

    self.S = zeros( shape=(3,3) )

    self.S[0,0] = 1./self.E1
    self.S[0,1] = -self.nu12/self.E1
    self.S[1,1] = 1./self.E2
    self.S[1,0] = self.S[0,1]
    self.S[2,2] = 1./self.G12

    return self.S

#------------------------------------------------------------------------------
#  getV
#------------------------------------------------------------------------------

  def getV( self ):

    if not hasattr( self , 'V' ):
      self.getS()

      self.V = zeros(5)

      self.V[0] = 0.125*(3.*self.S[0,0]+3.*self.S[1,1]+2.*self.S[0,1]+self.S[2,2])
      self.V[1] = 0.5*(self.S[0,0]-self.S[1,1])
      self.V[2] = 0.125*(self.S[0,0]+self.S[1,1]-2.*self.S[0,1]-self.S[2,2])
      self.V[3] = 0.125*(self.S[0,0]+self.S[1,1]+6.*self.S[0,1]-self.S[2,2])
      self.V[4] = 2.*(self.V[0]-self.V[3])

    return self.V

#------------------------------------------------------------------------------
#  getQbar
#    Returns the Qbar matrix for a given theta (in degrees)
#------------------------------------------------------------------------------

  def getQbar( self , theta ):

    if not hasattr( self , 'U' ):
      self.getU()

    Qbar = zeros( shape=(3,3) )

    rad = theta*pi/180.

    s2 = sin(2.*rad)
    s4 = sin(4.*rad)

    c2 = cos(2.*rad)
    c4 = cos(4.*rad)

    Qbar[0,0] = self.U[0]+self.U[1]*c2+self.U[2]*c4
    Qbar[0,1] = self.U[3]-self.U[2]*c4
    Qbar[1,0] = Qbar[0,1]
    Qbar[1,1] = self.U[0]-self.U[1]*c2+self.U[2]*c4
    Qbar[0,2] = 0.5*self.U[1]*s2+self.U[2]*s4
    Qbar[1,2] = 0.5*self.U[1]*s2-self.U[2]*s4
    Qbar[2,0] = Qbar[0,2]
    Qbar[2,1] = Qbar[1,2]
    Qbar[2,2] = self.U[4]-self.U[2]*c4

    return Qbar

#------------------------------------------------------------------------------
#  getSbar
#    Returns the Sbar matrix for a given theta
#------------------------------------------------------------------------------

  def getSbar( self , theta ):

    if not hasattr( self , 'V' ):
      self.getV()

    Sbar = zeros( shape=(3,3) )

    rad = theta*pi/180.

    s2 = sin(2.*rad)
    s4 = sin(4.*rad)

    c2 = cos(2.*rad)
    c4 = cos(4.*rad)

    Sbar[0,0] = self.V[0]+self.V[1]*c2+self.V[2]*c4
    Sbar[0,1] = self.V[3]-self.V[2]*c4
    Sbar[1,0] = Sbar[0,1]
    Sbar[1,1] = self.V[0]-self.V[1]*c2+self.V[2]*c4
    Sbar[0,2] = self.V[1]*s2+2.*self.V[2]*s4
    Sbar[1,2] = self.V[1]*s2-2.*self.V[2]*s4
    Sbar[2,0] = Sbar[0,2]
    Sbar[2,1] = Sbar[1,2]
    Sbar[2,2] = self.V[4]-4.*self.V[2]*c4

    return Sbar

#------------------------------------------------------------------------------
#  getAlpha
#    Returns the thermal expansion vector alpha for a given theta
#------------------------------------------------------------------------------

  def getAlpha( self , theta ):

    alpha = zeros(3)

    rad = theta*pi/180.

    s = sin(rad)
    c = cos(rad)

    alpha[0] = self.alpha1*c*c+self.alpha2*s*s
    alpha[1] = self.alpha1*s*s+self.alpha2*c*c
    alpha[2] = 2.*s*c*(self.alpha1-self.alpha2)

    return alpha

#------------------------------------------------------------------------------
#  getFIMaximumStress
#------------------------------------------------------------------------------

  def getFIMaximumStress( self , sigma ):

    if sigma[0] > 0.:
      FI1 = sigma[0]/self.Xt
    elif sigma[0] < 0.:
      FI1 = abs(sigma[0]/self.Xc);
    else:
      FI1 = 0.;

    if sigma[1] > 0.:
      FI2 = sigma[1]/self.Yt
    elif sigma[1] < 0.:
      FI2 = abs(sigma[1]/self.Yc);
    else:
      FI2 = 0.

    if sigma[2] == 0.:
      FI6 = 0.
    else:
      FI6 = abs(sigma[2]/self.Sl)

    return max(FI1,FI2,FI6)

#------------------------------------------------------------------------------
#  getFIMaximumStrain
#------------------------------------------------------------------------------

  def getFIMaximumStrain( self , sigma ):

    eps1 = 1.0/self.E1*(sigma[0]-self.nu12*sigma[1])
    eps2 = 1.0/self.E2*(sigma[1]-self.nu21*sigma[0])

    if eps1 > 0.:
      FI1 = (sigma[0]-self.nu12*sigma[1])/self.Xt
    elif eps1 < 0.:
      FI1 = abs((sigma[0]-self.nu12*sigma[1])/self.Xc)
    else:
      FI1 = 0.

    if eps2 > 0.:
      FI2 = (sigma[1]-self.nu21*sigma[0])/self.Yt
    elif eps2 < 0.:
      FI2 = abs((sigma[1]-self.nu21*sigma[0])/self.Yc)
    else:
      FI2 = 0.

    if sigma[2] == 0.:
      FI6 = 0.
    else:
      FI6 = abs(sigma[2]/self.Sl)

    return max(FI1,FI2,FI6)

#------------------------------------------------------------------------------
#  getFITsaiWu
#    Returns the TsaiWu failure index
#------------------------------------------------------------------------------

  def getFITsaiWu( self , sigma ):

    f1  = 1.0/self.Xt - 1.0/self.Xc
    f2  = 1.0/self.Yt - 1.0/self.Yc
    f11 = 1.0/(self.Xt*self.Xc)
    f22 = 1.0/(self.Yt*self.Yc)
    f66 = 1.0/(self.Sl*self.Sl)

    f12 = -sqrt(f11*f22)/2

    a = f11*sigma[0]**2 + f22*sigma[1]**2 + f66*sigma[2]**2 + 2*f12*sigma[0]*sigma[1]
    b = f1*sigma[0] + f2*sigma[1]

    Discr = b*b + 4*a
    SfTsaiWu1 = (-b + sqrt(Discr)) / (2*a)
    SfTsaiWu2 = (-b - sqrt(Discr)) / (2*a)

    return 1.0/SfTsaiWu1

#------------------------------------------------------------------------------
#  getFIHashin73
#    Returns the Hashin73 failure index
#------------------------------------------------------------------------------

  def getFIHashin73( self , sigma ):

    FIf = 0.0
    FIm = 0.0

    if sigma[0] >= 0:
      FIf = ( sigma[0] / self.Xt )**2 + ( sigma[2] / self.Sl )**2
    else:
      FIf = -sigma[0] / self.Xc

    if sigma[1] >= 0:
      FIm = ( sigma[1] / self.Yt )**2 + ( sigma[2] / self.Sl )**2
    else:
      FIm = ( sigma[1] / self.Yc )**2 + ( sigma[2] / self.Sl )**2

    return max(FIf,FIm)

#------------------------------------------------------------------------------
#  getFIHashin80
#    Returns the Hashin80 failure index
#------------------------------------------------------------------------------

  def getFIHashin80( self , sigma ):

    FIf = 0.0
    FIm = 0.0

    if sigma[0] >= 0:
      FIf = ( sigma[0] / self.Xt )**2 + ( sigma[2] / self.Sl )**2
    else:
      FIf = -sigma[0] / self.Xc

    if sigma[1] >= 0:
      FIm = ( sigma[1] / self.Yt )**2 + ( sigma[2] / self.Sl )**2
    else:
      FIm = ( sigma[1] / (2*self.Sl) )**2 + (( self.Yc / (2*self.Sl) )**2-1.0)*sigma[1]/self.Yc+( sigma[2] / self.Sl )**2

    return max(FIf,FIm)

#------------------------------------------------------------------------------
#  getFILarc03
#    Returns the LaRC03 failure index (not finished yet)
#------------------------------------------------------------------------------

  def getFILarc03( self , sigma ):

    t = 0.1

    eps1 = 1.0/self.E1*(sigma[0]-self.nu12*sigma[1])
    epsFail1 = self.Xt / self.E1

    YTis = sqrt(8.0*self.GIc/(pi*t*self.lam22))

    if not hasattr( self , 'SLis' ):
      SLis = sqrt(8.0*self.GIIc/(pi*t*self.lam44))
    else:
      SLis = self.SLis

    etaT = -1./self.tan2a0
    etaL = -SLis*self.cos2a0/(self.Yc*self.cosa0*self.cosa0)

    tauTeff = 0.0
    tauLeff = 0.0

    for i in range(18):
      alp = i*5
      alpr = alp * pi/180

      aa = Macauley( -sigma[1]*cos(alpr)*(sin(alpr)-etaT*cos(alpr) ) )

      if (aa > tauTeff):
        tauTeff = aa

      aa = Macauley( cos(alpr)*abs(sigma[2]+etaL*sigma[1]*cos(alpr)) )

      if (aa > tauLeff):
        tauLeff = aa

    ST = self.Yc * self.cosa0 * ( self.sina0 + self.cosa0 / self.tan2a0 )

    c1 = SLis/self.Xc
    c2 = ( 1. - sqrt( 1. - 4.*(c1+etaL)*c1))/(2.0*(c1+etaL))
    phiC = atan( c2 )

    phi = (abs(sigma[2])+(self.G12-self.Xc)*phiC)/(self.G12+sigma[0]-sigma[1])

    sigmam = stressTransformation( sigma , phi )

    # Matrix cracking

    if sigma[1] >= 0:
      FIm = (1.0 - self.g)*(sigma[1]/YTis)+self.g*(sigma[1]/YTis)**2+(sigma[2]/SLis)**2
    else:
      if sigma[0] < self.Yc:
#        FIm = ( taumTeff / ST )**2 + ( taumLeff / SLis )**2
        FIm = ( tauTeff / ST )**2 + ( tauLeff / SLis )**2
      else:
        FIm = ( tauTeff / ST )**2 + ( tauLeff / SLis )**2

    # Fibre failure

    if sigma[0] >= 0:
      FIf = eps1 / epsFail1
    else:
      if sigmam[1] < 0:
        FIf = Macauley( ( abs( sigmam[2] ) + etaL*sigmam[1] ) / SLis )
      else:
        FIf = (1.0 - self.g)*(sigmam[1]/YTis)+self.g*(sigmam[1]/YTis)**2+(sigmam[2]/SLis)**2

    return max(FIf,FIm)

#==============================================================================
#  Class Layer
#==============================================================================

class Layer:

  def __init__( self , name , theta , thick ):

    self.name  = name
    self.theta = theta
    self.thick = thick

#==============================================================================
#  Class Laminate
#==============================================================================

class Laminate:

  def __init__( self ):

    self.materials = {}
    self.layers    = []

#------------------------------------------------------------------------------
#  prints the laminate (layers, thickness, material name)
#------------------------------------------------------------------------------

  def __str__( self ):

    msg  = "  Laminate properties\n"
    msg += "  -----------------------------------------------------------\n"
    msg += "  layer   thick orient.  material\n"
    msg += "  -----------------------------------------------------------\n"

    for i,lay in enumerate(self.layers):
      msg += "   {:4}   {:4}   {:4}   {:4}\n".format(i,lay.thick,lay.theta,lay.name)

    return msg

#------------------------------------------------------------------------------
#  addMaterial
#------------------------------------------------------------------------------

  def addMaterial( self , name , mat ):

    self.materials[name] = mat

#------------------------------------------------------------------------------
#  addLayer
#------------------------------------------------------------------------------

  def addLayer( self , name , theta , thick ):

    layer = Layer( name , theta , thick )

    self.layers.append( layer )

    self.h     = zeros( len(self.layers)+1 )
    self.thick = 0.

    for i,layer in enumerate(self.layers):
      self.h[i+1] = self.thick+layer.thick
      self.thick += layer.thick

    self.h += -0.5*self.thick*ones( len(self.h) )

#------------------------------------------------------------------------------
#  removeAllLayers
#    Erases all layers from the laminate
#------------------------------------------------------------------------------

  def removeAllLayers( self ):

    self.layers = []

#------------------------------------------------------------------------------
#  getA
#    Returns the A matrix
#------------------------------------------------------------------------------

  def getA( self ):

    self.A = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.A 	+= self.materials[name].getQbar( theta ) * (self.h[i+1]-self.h[i])

    return self.A

#------------------------------------------------------------------------------
#  getB
#    Returns the B matrix
#------------------------------------------------------------------------------

  def getB( self ):

    self.B = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.B += 0.5*self.materials[name].getQbar( theta ) * (self.h[i+1]**2-self.h[i]**2)

    return self.B

#------------------------------------------------------------------------------
#  getD
#    Returns the D matrix
#------------------------------------------------------------------------------

  def getD( self ):

    self.D = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.D += 1.0/3.0*self.materials[name].getQbar( theta ) * (self.h[i+1]**3-self.h[i]**3)

    return self.D

#------------------------------------------------------------------------------
#  getTs
#    Returns the thermal extension array Ts
#------------------------------------------------------------------------------

  def getTs( self ):

    self.Ts = zeros( 3 )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      Qalpha = dot( self.materials[name].getQbar( theta ) , self.materials[name].getAlpha( theta ) )

      self.Ts += Qalpha * (self.h[i+1]-self.h[i])

    return self.Ts

#------------------------------------------------------------------------------
#  getTss
#    Returns the thermal bending array Tss
#------------------------------------------------------------------------------

  def getTss( self ):

    self.Tss = zeros( 3 )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      Qalpha = dot( self.materials[name].getQbar( theta ) , self.materials[name].getAlpha( theta ) )

      self.Tss += 0.5 * Qalpha * (self.h[i+1]**2-self.h[i]**2)

    return self.Tss

#------------------------------------------------------------------------------
#  getRhoh
#    Returns the term rho*h for a lamminate
#------------------------------------------------------------------------------

  def getRhoh( self ):

    rhoh = 0.

    for i,layer in enumerate(self.layers):
      name  = layer.name

      rhoh += self.materials[name].rho * (self.h[i+1]-self.h[i])

    return rhoh

#------------------------------------------------------------------------------
#  getInverseMatrices
#    Returns the inverse matrices A1, B1, C1, D1 of a laminate
#------------------------------------------------------------------------------

  def getInverseMatrices( self ):

    self.getA()
    self.getB()
    self.getD()

    Ainv  = inv(self.A)
    Dstar = self.D-dot(self.B,dot(Ainv,self.B))
    Dsinv = inv(Dstar)

    self.A1 = Ainv + dot(Ainv,dot(self.B,dot(Dsinv,dot(self.B,Ainv))))
    self.B1 = -dot(Ainv,dot(self.B,Dsinv))
    self.C1 = transpose(self.B1)
    self.D1 = Dsinv

    return self.A1,self.B1,self.C1,self.D1

#------------------------------------------------------------------------------
#  getQbar
#    Returns Qbar for layer i in laminate
#------------------------------------------------------------------------------

  def getQbar( self , i ):

    name  = self.layers[i].name
    theta = self.layers[i].theta

    return self.materials[name].getQbar( theta )

#------------------------------------------------------------------------------
#  getElastic
#    Returns the apparent elastic properties of a laminate
#------------------------------------------------------------------------------

  def getElastic( self ):

    self.getA()

    Ex   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[1,1])
    Ey   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[0,0])
    nuxy = self.A[0,1]/self.A[1,1]
    Gxy  = self.A[2,2]/self.thick

    return [Ex,Ey,nuxy,Gxy]

#==============================================================================
#  Utility functions
#==============================================================================

#------------------------------------------------------------------------------
#  stressTransformation
#    Transforms stress from 12 coordinate system to xy coordinate system
#------------------------------------------------------------------------------

def stressTransformation( sigma , theta ):

  signew = zeros( 3 )

  rad = theta*pi/180.

  c = cos(rad)
  s = sin(rad)

  signew[0] = sigma[0]*c*c + sigma[1]*s*s + 2.*sigma[2]*c*s
  signew[1] = sigma[0]*s*s + sigma[1]*c*c - 2.*sigma[2]*c*s
  signew[2] = ( sigma[1] - sigma[0] )*c*s + sigma[2]*( c*c - s*s )

  return signew

#------------------------------------------------------------------------------
#  strainTransformation
#    Transforms strain from 12 coordinate system to xy coordinate system
#------------------------------------------------------------------------------

def strainTransformation( eps , theta ):

  epsnew = zeros( 3 )

  rad = theta*pi/180.

  c = cos(rad)
  s = sin(rad)

  epsnew[0] = eps[0]*c*c + eps[1]*s*s + eps[2]*c*s
  epsnew[1] = eps[0]*s*s + eps[1]*c*c - eps[2]*c*s
  epsnew[2] = 2.0*( eps[1] - eps[0] )*c*s + eps[2]*( c*c - s*s )

  return epsnew

#------------------------------------------------------------------------------
#  mixMaterials
#    Simple Voigt and Reuss volume averaging to determine the homogenised
#    elastic properties of two materials with volume fraction vf.
#------------------------------------------------------------------------------

def mixMaterials ( fibre , matrix , vf ):

  E1   = fibre.E1*vf+matrix.E1*(1.-vf)
  E2   = fibre.E1*matrix.E1/(fibre.E2*(1.0-vf)+matrix.E2*vf)
  nu12 = fibre.nu12*vf+matrix.nu12*(1.-vf)
  G12  = fibre.G12*matrix.G12/(fibre.G12*(1.0-vf)+matrix.G12*vf)

# Add alpha

  mat = TransverseIsotropic( [E1,E2] , nu12 , G12 )

  return mat

#------------------------------------------------------------------------------
#  Macauley
#------------------------------------------------------------------------------

def Macauley( x ):

  if x > 0:
    return x
  else:
    return 0.
