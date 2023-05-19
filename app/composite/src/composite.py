"""A set of classes and functions to perform Classical Laminate calculations

   This module contains a set of classes and functions to perform classical
   laminate calculations for the course:
   Composite and Lightweight Materials - Design and Analysis (4MM00)

  (c) Joris Remmers (2013-2023)"""

from numpy import zeros,ones,dot,transpose
from numpy.linalg import inv
from math import sin,cos,pi,sqrt,tan,atan

class TransverseIsotropic:

  '''A class to describe atransversely isotropic material model.'''

  def __init__( self , E , nu12 , G12 = 0. , alpha = 0. , rho = 0. ):
  
    '''Inits the class TransverseIsotropic.
    
       The material linear elastic properties of the transversely isotropic
       material model are set in this function.
       
       Args: 
         E:      Young's modulus in case of an orthtropic material this 
                 is a list containing the E1 and E2. When it is a float, E1 and
                 E2 are considered to be equal to E.
         nu12:   Poisson ratio. 
         G12:    Shear modulus. When not set, it is assumed that the material is
                 isotropic and G is calculated as G= E/(2(1+nu))
         alpha:  Thermal expansion coefficients - can be a list, initially 0
         rho:    Density.'''
   
    if type(E) == list:
      if len(E) == 2:
        self.E1 = E[0]
        self.E2 = E[1]
      elif len(E) == 1:
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

    if type(alpha) == list:
      if len(alpha) == 2:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[1]
      elif len(alpha) == 1:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[0]
      else:
        print('error')
    else:
      self.alpha1    = alpha
      self.alpha2    = alpha

    self.rho = rho
    
  def setAlpha( self , alpha ):
  
    if type(alpha) == list:
      if len(alpha) == 2:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[1]
      elif len(alpha) == 1:
        self.alpha1 = alpha[0]
        self.alpha2 = alpha[0]
      else:
        print('error')
    else:
      self.alpha1    = alpha
      self.alpha2    = alpha

  def setFailureProperties( self, F , Gfrac = 0 , alpha0deg = 53. ):

    '''Setting the failure properties of the material.
    
       Args:
         F     A list constaining the 5 main failure parameters:
               F=[Xt,Xc,Yt,Yc,Sl] where 
                 Xt (longitudinal tensile strength)
                 Xc (longitudinal compressive strength)
                 Yt (transverse tensile strength)
                 Yc (transverse compressive strength)
                 Sl (transverse shear strength)
                 
               For the Tsai-Wu failure criterion, the f12 parameter is
               entered as the 6th parameter in F. In that case:
               F=[Xt,Xc,Yt,Yc,Sl,f12]
         Gfrac:  Fracture toughness (needed for Larc03 model
         alpha0deg: alpgha0 in degrees (needed for Larc03)'''

    if len(F) == 5 or len(F) == 6:
      self.Xt = F[0]
      self.Xc = F[1]
      self.Yt = F[2]
      self.Yc = F[3]
      self.Sl = F[4]
      
      if len(F) == 6:
        self.f12 = F[5]
    else:
      print('error')

    if type(Gfrac) == list:
      if len(Gfrac) == 2:
        self.GIc  = Gfrac[0]
        self.GIIc = Gfrac[1]
    else:
      self.GIc  = Gfrac
      self.GIIc = Gfrac

    self.a0deg  = alpha0deg
    
    alpha0 = alpha0deg*pi/180

    self.cosa0  = cos( alpha0 )
    self.sina0  = sin( alpha0 )

    self.cos2a0 = cos( 2.0*alpha0 )
    self.tan2a0 = tan( 2.0*alpha0 )

  def setSLis( self , SLis ):

    '''Set the interlaminar shear strength as needed for the Larc03 model'''
    
    self.SLis = SLis

  def __str__( self ):
  
    '''Prints the properties in the TransvereIstropic material model.
    
       Usage:
       
         mat = TransverseIsotropic( args )
         
         print(mat)'''

    msg  = "\n  Elastic Properties:\n"
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
      msg += "  S      :  {:12.3e}".format(self.Sl)
      
      if hasattr( self, "f12" ):
        msg += " , f12    : {:12.3e} \n".format(self.f12)
      else:
        msg += "\n"
        
      if hasattr( self, "GIc" ):
        msg += "  GIc    :  {:12.3e} , GIIc   :  {:12.3e} \n".format(self.GIc,self.GIIc)
        msg += "  alpha0 :  {:12.3e}\n".format(self.a0deg)

    return msg

  def getQ( self ):
  
    '''Returns the Q matrix of the material model as a numpy matrix.'''
    
    if not hasattr( self , 'Q' ):
      self.Q = zeros( shape=(3,3) )

      self.Q[0,0] = self.E1/(1.-self.nu12*self.nu21)
      self.Q[0,1] = self.nu12*self.E2/(1.0-self.nu12*self.nu21)
      self.Q[1,1] = self.E2/(1.-self.nu12*self.nu21)
      self.Q[1,0] = self.Q[0,1]
      self.Q[2,2] = self.G12
  
    return self.Q

  def getU( self ):

    '''Returns the stiffness matrix invariant terms U[1..5] as a numpy array.'''
    
    if not hasattr( self , 'U' ):
      self.getQ()

      self.U = zeros(5)

      self.U[0] = 0.125*(3.*self.Q[0,0]+3.*self.Q[1,1]+2.*self.Q[0,1]+4.*self.Q[2,2])
      self.U[1] = 0.5*(self.Q[0,0]-self.Q[1,1])
      self.U[2] = 0.125*(self.Q[0,0]+self.Q[1,1]-2.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[3] = 0.125*(self.Q[0,0]+self.Q[1,1]+6.*self.Q[0,1]-4.*self.Q[2,2])
      self.U[4] = 0.5*(self.U[0]-self.U[3])

    return self.U

  def getS( self ):

    '''Returns the compliance matrix of the material as a numpy array.'''
    
    self.S = zeros( shape=(3,3) )

    self.S[0,0] = 1./self.E1
    self.S[0,1] = -self.nu12/self.E1
    self.S[1,1] = 1./self.E2
    self.S[1,0] = self.S[0,1]
    self.S[2,2] = 1./self.G12

    return self.S

  def getV( self ):
 
    '''Returns the compliance matrix invariant terms V[1..5] as a numpy array.'''
    
    if not hasattr( self , 'V' ):
      self.getS()

      self.V = zeros(5)

      self.V[0] = 0.125*(3.*self.S[0,0]+3.*self.S[1,1]+2.*self.S[0,1]+self.S[2,2])
      self.V[1] = 0.5*(self.S[0,0]-self.S[1,1])
      self.V[2] = 0.125*(self.S[0,0]+self.S[1,1]-2.*self.S[0,1]-self.S[2,2])
      self.V[3] = 0.125*(self.S[0,0]+self.S[1,1]+6.*self.S[0,1]-self.S[2,2])
      self.V[4] = 2.*(self.V[0]-self.V[3])

    return self.V

  def getQbar( self , theta ):
  
    '''Calculates and returns the global stiffness matrix Q for a given theta angle.
    
       Args:
         theta:    Angle of the fibre direction with respect to the global x-axis in degrees.'''
         
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

  def getSbar( self , theta ):

    '''Calculates and returns the global compliance matrix S for a given theta angle.
    
       Args:
         theta:    Angle of the fibre direction with respect to the global x axis in degrees.'''
         
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

  def getAlpha( self , theta ):
  
    '''Calculates the thermal expansion vector alpha for a given angle.
    
       Args:
         theta:  Angle of the fibre direction with respect to the global x axis in degrees.'''

    alpha = zeros(3)

    rad = theta*pi/180.

    s = sin(rad)
    c = cos(rad)
 
    alpha[0] = self.alpha1*c*c+self.alpha2*s*s
    alpha[1] = self.alpha1*s*s+self.alpha2*c*c
    alpha[2] = 2.*s*c*(self.alpha1-self.alpha2)

    return alpha

  def getFIMaximumStress( self , sigma ):
  
    '''Calculates and return the Failure index for a given stress state according to
       the maximum stress criterion
       
       Args:
         sigma:    The current stress state.'''
    
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

  def getFIMaximumStrain( self , sigma ):
  
    '''Calculates and return the Failure index for a given stress state according to
       the maximum strain criterion
       
       Args:
         sigma:    The current stress state.'''
    
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

  def getFITsaiWu( self , sigma ):

    '''Calculates and return the Failure index for a given stress state according to
       the maximum Tsai-Wu criterion
       
       Args:
         sigma:    The current stress state.'''
             
    f1  = 1.0/self.Xt - 1.0/self.Xc
    f2  = 1.0/self.Yt - 1.0/self.Yc
    f11 = 1.0/(self.Xt*self.Xc)
    f22 = 1.0/(self.Yt*self.Yc)
    f66 = 1.0/(self.Sl*self.Sl)
    
    if hasattr( self , 'f12'):
      f12 = self.f12
    else:
      f12 = -sqrt(f11*f22)/2
        
    a = f11*sigma[0]**2 + f22*sigma[1]**2 + f66*sigma[2]**2 + 2*f12*sigma[0]*sigma[1]
    b = f1*sigma[0] + f2*sigma[1]
    
    Discr = b*b + 4*a
    SfTsaiWu1 = (-b + sqrt(Discr)) / (2*a)
    SfTsaiWu2 = (-b - sqrt(Discr)) / (2*a)
    
    return 1.0/SfTsaiWu1

  def getFIHashin73( self , sigma ):
  
    '''Calculates and return the Failure index for a given stress state according to
       the Hashin 73 criterion
       
       Args:
         sigma:    The current stress state.'''

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

  def getFIHashin80( self , sigma ):
  
    '''Calculates and return the Failure index for a given stress state according to
       the maximum Hashin80 criterion
       
       Args:
         sigma:    The current stress state.'''

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

  def getFILarc03( self , sigma ):
  
    '''Calculates and return the Failure index for a given stress state according to
       the maximum Larc03 criterion.
       
       Args:
         sigma:    The current stress state.'''

    t = 0.1
    
    lam22 = 2.0*(1.0/self.E2-self.nu21*self.nu21/self.E1)
    lam44 = 1.0/self.G12  
    
    g     = self.GIc/self.GIIc  

    eps1     = 1.0/self.E1*(sigma[0]-self.nu12*sigma[1])
    epsFail1 = self.Xt / self.E1
   
    YTis     = sqrt(8.0*self.GIc/(pi*t*lam22))

    if not hasattr( self , 'SLis' ):
      SLis = sqrt(8.0*self.GIIc/(pi*t*lam44))
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
   
      aa = Macauley( cos(alpr)*(abs(sigma[2])+etaL*sigma[1]*cos(alpr)))

      if (aa > tauLeff):
        tauLeff = aa

    ST = self.Yc * self.cosa0 * ( self.sina0 + self.cosa0 / self.tan2a0 )

    c1 = self.SLis/self.Xc
    c2 = ( 1. - sqrt( 1. - 4.*(c1+etaL)*c1))/(2.0*(c1+etaL))
    
    phiC = atan( c2 )

    phi = (abs(sigma[2])+(self.G12-self.Xc)*phiC)/(self.G12+sigma[0]-sigma[1])

    sigmam = stressTransformation( sigma , 180*phi/pi )

    # --- Matrix cracking ---

    if sigma[1] >= 0:
      FIm = (1.0 - g)*(sigma[1]/YTis)+g*(sigma[1]/YTis)**2+(sigma[2]/SLis)**2
    else:
      if sigma[0] < self.Yc:
        FIm = ( tauTeff / ST )**2 + ( tauLeff / SLis )**2
      else:
        FIm = ( tauTeff / ST )**2 + ( tauLeff / SLis )**2
     
    # --- Fibre failure ---
    
    if sigma[0] >= 0:
      FIf = eps1 / epsFail1
    else:
      if sigmam[1] < 0:
        FIf = Macauley( ( abs( sigmam[2] ) + etaL*sigmam[1] ) / SLis )
      else:
        FIf = (1.0 - g)*(sigmam[1]/YTis)+g*(sigmam[1]/YTis)**2+(sigmam[2]/SLis)**2

    return max(FIf,FIm)

#===============================================================================
#  Class Layer
#===============================================================================

class Layer:

  '''The layer class specifies the properties of a single layer.
  
     These properties are the name of the transversely isotropitc
     material, the orientation of the fibre direction with respect to the global
     x-axis and the thickness of the layer.'''

  def __init__( self , name , theta , thick ):
   
    '''Inits the class layer.
 
       Args:
          name:     name of the material model used.
          theta:    theorientation of the fibre direction in degrees.
          thick:    the thickness of the layer.'''
      
    self.name  = name
    self.theta = theta
    self.thick = thick

#===============================================================================
#  Class Laminate
#===============================================================================

class Laminate:
  
  '''A class to describe the stacking sequency of a laminate.'''
 
  def __init__( self ): 

    self.materials = {}
    self.layers    = []

  def __str__( self ):
  
    '''Prints information about the created laminate structure.'''

    msg  = "  Laminate properties\n"
    msg += "  -----------------------------------------------------------\n"
    msg += "  layer   thick orient.  material\n"
    msg += "  -----------------------------------------------------------\n"

    for i,lay in enumerate(self.layers):
      msg += "   {:4}   {:4}   {:4}   {:4}\n".format(i,lay.thick,lay.theta,lay.name)

    return msg

  def addMaterial( self , name , mat ):

    '''Adds a material to the laminate class.
    
       Adding a material model with a name to the created laminate structure.
    
       Args:
         name:      Name of the material model. This will be used as an identifier.
         mat:       The instance of the transverse material model.
         
       Usage:
         
         mat = TransverseIsotropic(..args..)
         
         lam = Laminate()
         
         lam.addMaterial( "glassfibre", mat )'''  
         
    self.materials[name] = mat

  def addLayer( self , name , theta , thick ):

    '''Adding a layer to the created laminate structure
    
       Args:
         name:    name of the material type
         theta:   angle with respect to x-axis in degrees.
         thick:   thickness of the layer.'''
                
    layer = Layer( name , theta , thick )

    self.layers.append( layer )

    self.h     = zeros( len(self.layers)+1 )
    self.thick = 0.
    
    for i,layer in enumerate(self.layers):
      self.h[i+1] = self.thick+layer.thick
      self.thick += layer.thick

    self.h += -0.5*self.thick*ones( len(self.h) )

  def removeAllLayers( self ):

    '''Erases all layers from the current laminate.'''
    
    self.layers = []

  def getA( self ):
  
    '''Return the extensional stiffness matrix A for the given laminate as
       a (3x3) numpy matrix.'''

    self.A = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.A 	+= self.materials[name].getQbar( theta ) * (self.h[i+1]-self.h[i])

    return self.A

  def getB( self ):
  
    '''Return the coupling stiffness matrix B for the given laminate as
       a (3x3) numpy matrix.'''  

    self.B = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.B += 0.5*self.materials[name].getQbar( theta ) * (self.h[i+1]**2-self.h[i]**2)

    return self.B

  def getD( self ):
  
    '''Return the bending stiffness matrix D for the given laminate as
       a (3x3) numpy matrix.'''

    self.D = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.D += 1.0/3.0*self.materials[name].getQbar( theta ) * (self.h[i+1]**3-self.h[i]**3)

    return self.D

  def getTs( self ):

    '''Calculates and returns the thermal exansion array T* pf the laminate as
       a numpy array of length 3x1.'''
       
    self.Ts = zeros( 3 )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      Qalpha = dot( self.materials[name].getQbar( theta ) , self.materials[name].getAlpha( theta ) )

      self.Ts += Qalpha * (self.h[i+1]-self.h[i])

    return self.Ts

  def getTss( self ):
  
    '''Calculates and returns the thermal exansion array T** of the laminate as
       a numpy array of length 3x1.'''

    self.Tss = zeros( 3 )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      Qalpha = dot( self.materials[name].getQbar( theta ) , self.materials[name].getAlpha( theta ) )

      self.Tss += 0.5 * Qalpha * (self.h[i+1]**2-self.h[i]**2)

    return self.Tss

  def getRhoh( self ):

    '''Calculates and returns the mass per unit area of the laminate (rho*h) fo the laminate.'''
    
    rhoh = 0.

    for i,layer in enumerate(self.layers):
      name  = layer.name

      rhoh += self.materials[name].rho * (self.h[i+1]-self.h[i])

    return rhoh
    
  def getZcoord( self , j ):

    '''Calculates and returns the z coordinate of layer j. The z coordinate is
       defined as the centre of the layer.
       
       Args:
         j     layer number.'''  
         
    return 0.5*(self.h[j] + self.h[j+1])  

  def getLayerBounds( self , j ):

    '''Returns the z coordinates of the top and the bottom of layer j.      
       
       Args:
         j     layer number.'''  
         
    return (self.h[j],self.h[j+1]) 
    
  def getInverseMatrices( self ):
  
    '''Calculate and returns the four inverse matrices A1, B1, C1, D1 as 4
       numpy arrays of shape (3x3).'''

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

  def getQbar( self , j ):

    '''Returns the stiffness matrix Qbar for a given layer.
    
       Args:
         j    layer number.'''
         
    name  = self.layers[j].name
    theta = self.layers[j].theta

    return self.materials[name].getQbar( theta )

  def getElastic( self ):
  
    '''Calculates and returns the 4 apparent elastic properties of the laminate as
       a list: [Ex,Ey,nuxy,Gxy].'''
    
    self.getA()

    Ex   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[1,1])
    Ey   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[0,0])
    nuxy = self.A[0,1]/self.A[1,1]
    Gxy  = self.A[2,2]/self.thick

    return [Ex,Ey,nuxy,Gxy]

#==============================================================================
#  Utility functions
#==============================================================================

def stressTransformation( sigma , theta ):

  '''Function to transform stress from the 12 coordinate system 
     (aligned with the fibre direction) to xy coordinate system (global axis)
     
     Args:
       sigma:     The stress state in the 12 frame of reference
       theta:     the angle between 122 and xy in degrees.'''
       
  signew = zeros( 3 )

  rad = theta*pi/180.

  c = cos(rad)
  s = sin(rad)

  signew[0] = sigma[0]*c*c + sigma[1]*s*s + 2.*sigma[2]*c*s
  signew[1] = sigma[0]*s*s + sigma[1]*c*c - 2.*sigma[2]*c*s
  signew[2] = ( sigma[1] - sigma[0] )*c*s + sigma[2]*( c*c - s*s )
   
  return signew

def strainTransformation( eps , theta ):

  '''Function to transform strain from the 12 coordinate system 
     (aligned with the fibre direction) to xy coordinate system (global axis)
     
     Args:
       eps:       The stress state in the 12 frame of reference
       theta:     the angle between the 12 and xy coordinate frame in degrees.'''

  epsnew = zeros( 3 )

  rad = theta*pi/180.

  c = cos(rad)
  s = sin(rad)

  epsnew[0] = eps[0]*c*c + eps[1]*s*s + eps[2]*c*s
  epsnew[1] = eps[0]*s*s + eps[1]*c*c - eps[2]*c*s
  epsnew[2] = 2.0*( eps[1] - eps[0] )*c*s + eps[2]*( c*c - s*s )

  return epsnew
  
def mixMaterials ( fibre , matrix , vf ):

  '''Simple Voigt and Reuss volume averaging to determine the 
     homogenised elastic properties of two materials with fibre 
     volume fraction vf.
     
     Args:
       fibre:   the material model of the fibre material
       matrix:  the material model of the matrix material
       vf:      the fibre volume fraction.'''

  E1   = fibre.E1*vf+matrix.E1*(1.-vf)
  E2   = fibre.E1*matrix.E1/(fibre.E2*(1.0-vf)+matrix.E2*vf)
  nu12 = fibre.nu12*vf+matrix.nu12*(1.-vf)
  G12  = fibre.G12*matrix.G12/(fibre.G12*(1.0-vf)+matrix.G12*vf)

  mat = TransverseIsotropic( [E1,E2] , nu12 , G12 )

  return mat

def Macauley( x ):

  '''The Macauley operater. Returns argument x when x > 0. Otherwise 0.'''

  if x > 0:
    return x
  else:
    return 0.
