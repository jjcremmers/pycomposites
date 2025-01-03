import numpy as np
from numpy import zeros,ones,dot,transpose
from numpy.linalg import inv
from math import sin,cos,pi,sqrt,tan,atan

'''
A set of classes and functions to perform Classical Laminate calculations

This module contains a set of classes and functions to perform classical
laminate calculations for the course:
Composite and Lightweight Materials - Design and Analysis (4MM00)

(c) Joris Remmers (2013-2023)
'''

class TransverseIsotropic:

    """
    A class to represent a transversely isotropic material model.

    This class models materials that exhibit isotropic properties in a plane
    perpendicular to a single axis of symmetry and anisotropic behavior along
    the axis of symmetry.
    
    Args:
        E (float or list[float]): Young's modulus. 
            If `E` is a list, 
            it should contain `E1` and `E2` for an orthotropic material. 
            If `E` is a float, `E1` and `E2` are considered equal to `E`.    
        nu12 (float): Poisson's ratio.
        G12 (float, optional): Shear modulus.  Defaults to 0.
            If not specified, it is assumed 
            that the material is isotropic, and `G` is 
            calculated as 
            $$
            G = E / (2 * (1 + nu))
            $$
        alpha (float or list[float], optional): Thermal expansion coefficients. Can be a single value or a list. Defaults to 0.
        rho (float, optional): Density of the material. Defaults to 0.                       
    """

    def __init__(self, E: float | list[float], nu12: float, 
                 G12: float = 0.0, alpha: float | list[float] = 0.0, 
                 rho: float = 0.0) -> None:
                 
        """
        Initialize the TransverseIsotropic material model.
        """
   
        if type(E) == list:
            if len(E) == 2:
                self.E1 = E[0]
                self.E2 = E[1]
            elif len(E) == 1:
                self.E1 = E[0]
                self.E2 = E[0]
            else:
                raise ValueError("E list must have one or two elements.")
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
                raise ValueError("alpha list must have one or two elements.")
        else:
            self.alpha1    = alpha
            self.alpha2    = alpha

        self.rho = rho

#-
#
#-
   
    def setAlpha(self, alpha: float | list[float]) -> None:
    
        """
        Set the thermal expansion coefficients (alpha) for the material.

        :param alpha: Thermal expansion coefficients. Can be a single float value 
                      or a list of one or two float values. If a single value is 
                      provided, it is assigned to both `alpha1` and `alpha2`. 
                      If a list of two values is provided, the first value is 
                      assigned to `alpha1` and the second to `alpha2`.
        :type alpha:  float or list[float]
        :raises ValueError: If the list has more than two elements.
        """
        
        if isinstance(alpha, list):
            if len(alpha) == 2:
                self.alpha1 = alpha[0]
                self.alpha2 = alpha[1]
            elif len(alpha) == 1:
                self.alpha1 = alpha[0]
                self.alpha2 = alpha[0]
            else:
                raise ValueError("alpha list must have one or two elements.")
        else:
            self.alpha1 = alpha
            self.alpha2 = alpha

#
#
#

    def setFailureProperties(self, F: list[float], Gfrac: list[float] = 0.0, alpha0deg: float = 53.0) -> None:

        """
        Set the failure properties of the material.

        Args:
            F (list[float]): A list containing the main failure parameters:
                             F = [Xt, Xc, Yt, Yc, Sl] where:
                                  - Xt: Longitudinal tensile strength
                                  - Xc: Longitudinal compressive strength
                                  - Yt: Transverse tensile strength
                                  - Yc: Transverse compressive strength
                                  - Sl: Transverse shear strength
                              For the Tsai-Wu failure criterion, f12 can be provided as the 6th parameter:
                              F = [Xt, Xc, Yt, Yc, Sl, f12].                              
            Gfrac (list[float], optional): Fracture toughness. Can be a single value or a list of two values 
                [GIc, GIIc]. Defaults to 0.                                                   
            alpha0deg (float, optional): Alpha0 in degrees, needed for the Larc03 model. Defaults to 53.0.
        Raises:
            ValueError: If `F` has an invalid length.            
        """
        
        if len(F) == 5 or len(F) == 6:
            self.Xt, self.Xc, self.Yt, self.Yc, self.Sl = F[:5]
            if len(F) == 6:
                self.f12 = F[5]
        else:
            raise ValueError("F must contain 5 or 6 elements.")

        if isinstance(Gfrac, list):
            if len(Gfrac) == 2:
                self.GIc, self.GIIc = Gfrac
            else:
                raise ValueError("Gfrac list must contain exactly 2 elements.")
        else:
            self.GIc = GIIc = Gfrac

        self.a0deg = alpha0deg
        alpha0 = alpha0deg * pi / 180

        self.cosa0 = cos(alpha0)
        self.sina0 = sin(alpha0)
        self.cos2a0 = cos(2.0 * alpha0)
        self.tan2a0 = tan(2.0 * alpha0)
        
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setSLis(self, SLis: float) -> None:
        
        """
        Set the interlaminar shear strength (SLis) for the Larc03 model.

        :param SLis: Interlaminar shear strength.
        :type SLis: float
        """
        
        self.SLis = SLis

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __str__(self) -> str:
        
        """
        Print the properties of the TransverseIsotropic material model.

        :return: Formatted string containing the material properties.
        :rtype: str
        """
        
        msg = "\n  Elastic Properties:\n"
        msg += "  -----------------------------------------------------------\n"
        msg += f"  E1     :  {self.E1:12.3e} , E2     :  {self.E2:12.3e} \n"
        msg += f"  nu12   :  {self.nu12:12.2f} , G12    :  {self.G12:12.3e} \n"

        if self.rho > 0.0:
            msg += f"  rho    :  {self.rho:12.2f}\n"

        if hasattr(self, "alpha1"):
            msg += "\n  Thermal expansion coefficients:\n"
            msg += "  -----------------------------------------------------------\n"
            msg += f"  alpha1 :  {self.alpha1:12.3e} , alpha2 :  {self.alpha2:12.3e} \n"

        if hasattr(self, "Xt"):
            msg += "\n  Strengths and failure model parameters:\n"
            msg += "  -----------------------------------------------------------\n"
            msg += f"  Xt     :  {self.Xt:12.3e} , Xc     :  {self.Xc:12.3e} \n"
            msg += f"  Yt     :  {self.Yt:12.3e} , Yc     :  {self.Yc:12.3e} \n"
            msg += f"  S      :  {self.Sl:12.3e}"
            if hasattr(self, "f12"):
                msg += f" , f12    : {self.f12:12.3e} \n"
            else:
                msg += "\n"

            if hasattr(self, "GIc"):
                msg += f"  GIc    :  {self.GIc:12.3e} , GIIc   :  {self.GIIc:12.3e} \n"
                msg += f"  alpha0 :  {self.a0deg:12.3e}\n"

        return msg

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getQ(self) -> np.ndarray:
        
        """
        Calculates and returns the stiffness matrix Q of the material model.

        :return: A 3x3 numpy array representing the stiffness matrix Q.
        """
        
        if not hasattr(self, 'Q'):
            self.Q = np.zeros((3, 3))

            factor = 1.0 / (1.0 - self.nu12 * self.nu21)
            self.Q[0, 0] = self.E1 * factor
            self.Q[0, 1] = self.nu12 * self.E2 * factor
            self.Q[1, 1] = self.E2 * factor
            self.Q[1, 0] = self.Q[0, 1]
            self.Q[2, 2] = self.G12

        return self.Q

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getU(self) -> np.ndarray:

        """
        Returns the stiffness matrix invariant terms U[1..5].

        :return: A numpy array of size 5 containing the invariant terms.
        """

        if not hasattr(self, 'U'):
            self.getQ()

            self.U = np.zeros(5)
            self.U[0] = 0.125 * (3.0 * self.Q[0, 0] + 3.0 * self.Q[1, 1] + 2.0 * self.Q[0, 1] + 4.0 * self.Q[2, 2])
            self.U[1] = 0.5 * (self.Q[0, 0] - self.Q[1, 1])
            self.U[2] = 0.125 * (self.Q[0, 0] + self.Q[1, 1] - 2.0 * self.Q[0, 1] - 4.0 * self.Q[2, 2])
            self.U[3] = 0.125 * (self.Q[0, 0] + self.Q[1, 1] + 6.0 * self.Q[0, 1] - 4.0 * self.Q[2, 2])
            self.U[4] = 0.5 * (self.U[0] - self.U[3])

        return self.U

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getS(self) -> np.ndarray:

        """
        Returns the compliance matrix S of the material.

        :return: A 3x3 numpy array representing the compliance matrix S.
        """

        S = np.zeros((3, 3))

        S[0, 0] = 1.0 / self.E1
        S[0, 1] = -self.nu12 / self.E1
        S[1, 1] = 1.0 / self.E2
        S[1, 0] = S[0, 1]
        S[2, 2] = 1.0 / self.G12

        return S

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getV(self) -> np.ndarray:

        """
        Returns the compliance matrix invariant terms V[1..5].

        :return: A numpy array of size 5 containing the invariant terms.
        """

        if not hasattr(self, 'V'):
            S = self.getS()

            self.V = np.zeros(5)
            self.V[0] = 0.125 * (3.0 * S[0, 0] + 3.0 * S[1, 1] + 2.0 * S[0, 1] + S[2, 2])
            self.V[1] = 0.5 * (S[0, 0] - S[1, 1])
            self.V[2] = 0.125 * (S[0, 0] + S[1, 1] - 2.0 * S[0, 1] - S[2, 2])
            self.V[3] = 0.125 * (S[0, 0] + S[1, 1] + 6.0 * S[0, 1] - S[2, 2])
            self.V[4] = 2.0 * (self.V[0] - self.V[3])

        return self.V

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getQbar(self, theta: float) -> np.ndarray:

        """
        Calculates the global stiffness matrix Q for a given fiber angle.

        :param theta: Angle of the fiber direction with respect to the global x-axis, in degrees.
        :return: A 3x3 numpy array representing the global stiffness matrix Q.
        """

        if not hasattr(self, 'U'):
            self.getU()

        Qbar = np.zeros((3, 3))

        rad = theta * pi / 180.0
        s2 = sin(2.0 * rad)
        s4 = sin(4.0 * rad)
        c2 = cos(2.0 * rad)
        c4 = cos(4.0 * rad)

        Qbar[0, 0] = self.U[0] + self.U[1] * c2 + self.U[2] * c4
        Qbar[0, 1] = self.U[3] - self.U[2] * c4
        Qbar[1, 0] = Qbar[0, 1]
        Qbar[1, 1] = self.U[0] - self.U[1] * c2 + self.U[2] * c4
        Qbar[0, 2] = 0.5 * self.U[1] * s2 + self.U[2] * s4
        Qbar[1, 2] = 0.5 * self.U[1] * s2 - self.U[2] * s4
        Qbar[2, 0] = Qbar[0, 2]
        Qbar[2, 1] = Qbar[1, 2]
        Qbar[2, 2] = self.U[4] - self.U[2] * c4

        return Qbar

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getSbar(self, theta: float) -> np.ndarray:

        """
        Calculates the global compliance matrix S for a given fiber angle.

        :param theta: Angle of the fiber direction with respect to the global x-axis, in degrees.
        :return: A 3x3 numpy array representing the global compliance matrix S.
        """

        if not hasattr(self, 'V'):
            self.getV()

        Sbar = np.zeros((3, 3))

        rad = theta * pi / 180.0
        s2 = sin(2.0 * rad)
        s4 = sin(4.0 * rad)
        c2 = cos(2.0 * rad)
        c4 = cos(4.0 * rad)

        Sbar[0, 0] = self.V[0] + self.V[1] * c2 + self.V[2] * c4
        Sbar[0, 1] = self.V[3] - self.V[2] * c4
        Sbar[1, 0] = Sbar[0, 1]
        Sbar[1, 1] = self.V[0] - self.V[1] * c2 + self.V[2] * c4
        Sbar[0, 2] = self.V[1] * s2 + 2.0 * self.V[2] * s4
        Sbar[1, 2] = self.V[1] * s2 - 2.0 * self.V[2] * s4
        Sbar[2, 0] = Sbar[0, 2]
        Sbar[2, 1] = Sbar[1, 2]
        #Sbar[2, 2] TO HERE

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getAlpha(self, theta: float) -> np.ndarray:
        
        """
        Calculate the thermal expansion vector alpha for a given angle.

        Args:
            theta (float): Angle of the fibre direction with respect to the global x-axis in degrees.

        Returns:
            np.ndarray: A numpy array representing the thermal expansion vector.
        """
        
        alpha = np.zeros(3)
        rad = theta * pi / 180

        s = sin(rad)
        c = cos(rad)

        alpha[0] = self.alpha1 * c * c + self.alpha2 * s * s
        alpha[1] = self.alpha1 * s * s + self.alpha2 * c * c
        alpha[2] = 2.0 * s * c * (self.alpha1 - self.alpha2)

        return alpha

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getFIMaximumStress(self, sigma: np.ndarray ) -> float:
        
        """Calculates and returns the Failure index for a given stress state according to the maximum stress criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the maximum stress criterion.
        """
        
        if sigma[0] > 0.:
            FI1 = sigma[0] / self.Xt
        elif sigma[0] < 0.:
            FI1 = abs(sigma[0] / self.Xc)
        else:
            FI1 = 0.

        if sigma[1] > 0.:
            FI2 = sigma[1] / self.Yt
        elif sigma[1] < 0.:
            FI2 = abs(sigma[1] / self.Yc)
        else:
            FI2 = 0.

        if sigma[2] == 0.:
            FI6 = 0.
        else:
            FI6 = abs(sigma[2] / self.Sl)

        return max(FI1, FI2, FI6)

    def getFIMaximumStrain(self, sigma: np.ndarray ) -> float:
        """Calculates and returns the Failure index for a given stress state according to the maximum strain criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the maximum strain criterion.
        """
        eps1 = 1.0 / self.E1 * (sigma[0] - self.nu12 * sigma[1])
        eps2 = 1.0 / self.E2 * (sigma[1] - self.nu21 * sigma[0])

        if eps1 > 0.:
            FI1 = (sigma[0] - self.nu12 * sigma[1]) / self.Xt
        elif eps1 < 0.:
            FI1 = abs((sigma[0] - self.nu12 * sigma[1]) / self.Xc)
        else:
            FI1 = 0.

        if eps2 > 0.:
            FI2 = (sigma[1] - self.nu21 * sigma[0]) / self.Yt
        elif eps2 < 0.:
            FI2 = abs((sigma[1] - self.nu21 * sigma[0]) / self.Yc)
        else:
            FI2 = 0.

        if sigma[2] == 0.:
            FI6 = 0.
        else:
            FI6 = abs(sigma[2] / self.Sl)

        return max(FI1, FI2, FI6)

    def getFITsaiWu(self, sigma: np.ndarray ) -> float:
        """Calculates and returns the Failure index for a given stress state according to the maximum Tsai-Wu criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the Tsai-Wu criterion.
        """
        f1 = 1.0 / self.Xt - 1.0 / self.Xc
        f2 = 1.0 / self.Yt - 1.0 / self.Yc
        f11 = 1.0 / (self.Xt * self.Xc)
        f22 = 1.0 / (self.Yt * self.Yc)
        f66 = 1.0 / (self.Sl * self.Sl)

        if hasattr(self, 'f12'):
            f12 = self.f12
        else:
            f12 = -sqrt(f11 * f22) / 2

        a = f11 * sigma[0]**2 + f22 * sigma[1]**2 + f66 * sigma[2]**2 + 2 * f12 * sigma[0] * sigma[1]
        b = f1 * sigma[0] + f2 * sigma[1]

        Discr = b * b + 4 * a
        SfTsaiWu1 = (-b + sqrt(Discr)) / (2 * a)
        SfTsaiWu2 = (-b - sqrt(Discr)) / (2 * a)

        return 1.0 / SfTsaiWu1

    def getFIHashin73(self, sigma: np.ndarray ) -> float:
        """Calculates and returns the Failure index for a given stress state according to the Hashin 73 criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the Hashin 73 criterion.
        """
        FIf = 0.0
        FIm = 0.0

        if sigma[0] >= 0:
            FIf = (sigma[0] / self.Xt)**2 + (sigma[2] / self.Sl)**2
        else:
            FIf = -sigma[0] / self.Xc

        if sigma[1] >= 0:
            FIm = (sigma[1] / self.Yt)**2 + (sigma[2] / self.Sl)**2
        else:
            FIm = (sigma[1] / self.Yc)**2 + (sigma[2] / self.Sl)**2

        return max(FIf, FIm)

    def getFIHashin80(self, sigma: np.ndarray ) -> float:
        """Calculates and returns the Failure index for a given stress state according to the Hashin 80 criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the Hashin 80 criterion.
        """
        FIf = 0.0
        FIm = 0.0

        if sigma[0] >= 0:
            FIf = (sigma[0] / self.Xt)**2 + (sigma[2] / self.Sl)**2
        else:
            FIf = -sigma[0] / self.Xc

        if sigma[1] >= 0:
            FIm = (sigma[1] / self.Yt)**2 + (sigma[2] / self.Sl)**2
        else:
            FIm = (sigma[1] / (2 * self.Sl))**2 + ((self.Yc / (2 * self.Sl))**2 - 1.0) * sigma[1] / self.Yc + (sigma[2] / self.Sl)**2

        return max(FIf, FIm)

    def getFILarc03(self, sigma: np.ndarray ) -> float:
        """Calculates and returns the Failure index for a given stress state according to the Larc03 criterion.

        Args:
            sigma: The current stress state.

        Returns:
            The Failure index based on the Larc03 criterion.
        """
        t = 0.1

        lam22 = 2.0 * (1.0 / self.E2 - self.nu21 * self.nu21 / self.E1)
        lam44 = 1.0 / self.G12

        g = self.GIc / self.GIIc

        eps1 = 1.0 / self.E1 * (sigma[0] - self.nu12 * sigma[1])
        epsFail1 = self.Xt / self.E1

        YTis = sqrt(8.0 * self.GIc / (pi * t * lam22))

        if not hasattr(self, 'SLis'):
            SLis = sqrt(8.0 * self.GIIc / (pi * t * lam44))
        else:
            SLis = self.SLis

        etaT = -1. / self.tan2a0
        etaL = -SLis * self.cos2a0 / (self.Yc * self.cosa0 * self.cosa0)

        tauTeff = 0.0
        tauLeff = 0.0

        for i in range(18):
            alp = i * 5
            alpr = alp * pi / 180

            aa = Macauley(-sigma[1] * cos(alpr) * (sin(alpr) - etaT * cos(alpr)))

            if aa > tauTeff:
                tauTeff = aa

            aa = Macauley(cos(alpr) * (abs(sigma[2]) + etaL * sigma[1] * cos(alpr)))

            if aa > tauLeff:
                tauLeff = aa

        ST = self.Yc * self.cosa0 * (self.sina0 + self.cosa0 / self.tan2a0)

        c1 = self.SLis / self.Xc
        c2 = (1. - sqrt(1. - 4. * (c1 + etaL) * c1)) / (2.0 * (c1 + etaL))

        phiC = atan(c2)

        phi = (abs(sigma[2]) + (self.G12 - self.Xc) * phiC) / (self.G12 + sigma[0] - sigma[1])

        sigmam = stressTransformation(sigma, 180 * phi / pi)

        # --- Matrix cracking ---
        if sigma[1] >= 0:
            FIm = (1.0 - g) * (sigma[1] / YTis) + g * (sigma[1] / YTis)**2 + (sigma[2] / SLis)**2
        else:
            if sigma[0] < self.Yc:
                FIm = (tauTeff / ST)**2 + (tauLeff / SLis)**2
            else:
                FIm = (tauTeff / ST)**2 + (tauLeff / SLis)**2

        # --- Fibre failure ---
        if sigma[0] >= 0:
            FIf = eps1 / epsFail1
        else:
            if sigmam[1] < 0:
                FIf = Macauley((abs(sigmam[2]) + etaL * sigmam[1]) / SLis)
            else:
                FIf = (1.0 - g) * (sigmam[1] / YTis) + g * (sigmam[1] / YTis)**2 + (sigmam[2] / SLis)**2

        return max(FIf, FIm)

#===============================================================================
#  Class Layer
#===============================================================================

class Layer:

  '''
  The layer class specifies the properties of a single layer.
  
  These properties are the name of the transversely isotropitc
  material, the orientation of the fibre direction with respect to the global
  x-axis and the thickness of the layer.
     
  Args:
    name:     name of the material model used.
    theta:    theorientation of the fibre direction in degrees.
    thick:    the thickness of the layer.         
  '''

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def __init__( self , name , theta , thick ):
   
    '''
    Inits the class layer.
    '''
      
    self.name  = name
    self.theta = theta
    self.thick = thick

#===============================================================================
#  Class Laminate
#===============================================================================

class Laminate:
  
  '''A class to describe the stacking sequency of a laminate.'''

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
 
  def __init__( self ): 

    self.materials = {}
    self.layers    = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def __str__( self ):
  
    '''Prints information about the created laminate structure.'''

    msg  = "  Laminate properties\n"
    msg += "  -----------------------------------------------------------\n"
    msg += "  layer   thick orient.  material\n"
    msg += "  -----------------------------------------------------------\n"

    for i,lay in enumerate(self.layers):
      msg += "   {:4}   {:4}   {:4}   {:4}\n".format(i,lay.thick,lay.theta,lay.name)

    return msg

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def addMaterial( self , name , mat ):

    """ Adds a material to the laminate class.
    
        Adding a material model with a name to the created laminate structure.
    
        Args:
          name:      Name of the material model. This will be used as an identifier.
          mat:       The instance of the transverse material model.
         
        Example:
         
          mat = TransverseIsotropic(..args..)      
          lam = Laminate()       
          lam.addMaterial( "glassfibre", mat )
         
          Here's an example of how to use this function:

          >>> result = function_name(5, 'example')
          >>> print(result)
         True
    """  
         
    self.materials[name] = mat

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  
  def removeAllLayers( self ):

    '''Erases all layers from the current laminate.'''
    
    self.layers = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getA( self ):
  
    '''Return the extensional stiffness matrix A for the given laminate as
       a (3x3) numpy matrix.'''

    self.A = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.A 	+= self.materials[name].getQbar( theta ) * (self.h[i+1]-self.h[i])

    return self.A

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getB( self ):
  
    '''Return the coupling stiffness matrix B for the given laminate as
       a (3x3) numpy matrix.'''  

    self.B = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.B += 0.5*self.materials[name].getQbar( theta ) * (self.h[i+1]**2-self.h[i]**2)

    return self.B

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getD( self ):
  
    '''Return the bending stiffness matrix D for the given laminate as
       a (3x3) numpy matrix.'''

    self.D = zeros( shape = ( 3,3) )

    for i,layer in enumerate(self.layers):
      name  = layer.name
      theta = layer.theta

      self.D += 1.0/3.0*self.materials[name].getQbar( theta ) * (self.h[i+1]**3-self.h[i]**3)

    return self.D

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getRhoh( self ):

    '''Calculates and returns the mass per unit area of the laminate (rho*h) fo the laminate.'''
    
    rhoh = 0.

    for i,layer in enumerate(self.layers):
      name  = layer.name

      rhoh += self.materials[name].rho * (self.h[i+1]-self.h[i])

    return rhoh

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
  def getZcoord( self , j : int ) -> float:

    '''
    Calculates and returns the z coordinate of layer j. The z coordinate is
    defined as the centre of the layer.
       
    Args:
      j     layer ID.
         
    Returns:
      z (float): z coordinate of the layer. 
    '''  
         
    return 0.5*(self.h[j] + self.h[j+1])  

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getLayerBounds( self , j ):

    '''Returns the z coordinates of the top and the bottom of layer j.      
       
       Args:
         j     layer number.'''  
         
    return (self.h[j],self.h[j+1]) 

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getQbar( self , j : int ) -> np.ndarray:

    '''Returns the stiffness matrix Qbar for a given layer.
    
       Args:
         j    layer number.'''
         
    name  = self.layers[j].name
    theta = self.layers[j].theta

    return self.materials[name].getQbar( theta )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

  def getElastic( self ) -> list[float]:
  
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

def stressTransformation( sigma : np.ndarray , theta : float ) -> np.ndarray:

    '''
    Function to transform stress from the 12 coordinate system 
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def strainTransformation( eps : np.ndarray , theta : float ) -> np.ndarray:

    '''
    Function to transform strain from the 12 coordinate system 
    (aligned with the fibre direction) to xy coordinate system (global axis)
     
    Args:
        eps:       The stress state in the 12 frame of reference
        theta:     the angle between the 12 and xy coordinate frame in degrees.
    '''

    epsnew = zeros( 3 )

    rad = theta*pi/180.

    c = cos(rad)
    s = sin(rad)

    epsnew[0] = eps[0]*c*c + eps[1]*s*s + eps[2]*c*s
    epsnew[1] = eps[0]*s*s + eps[1]*c*c - eps[2]*c*s
    epsnew[2] = 2.0*( eps[1] - eps[0] )*c*s + eps[2]*( c*c - s*s )

    return epsnew

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
  

def mixMaterials ( fibre : TransverseIsotropic , matrix  : TransverseIsotropic , vf : float ) ->  TransverseIsotropic:

    '''
    Simple Voigt and Reuss volume averaging to determine the 
    homogenised elastic properties of two materials with fibre 
    volume fraction vf.
     
    Args:
        fibre(TransverseIsotropic):   the material model of the fibre material
        matrix(TransverseIsotropic):  the material model of the matrix material
        vf(float):      the fibre volume fraction.
        
    Returns:
        TransverseIsotropic:   The matrials
        
    Raises:      
        ValueError: when vf is not between 0 and 1.
    '''
       
    if vf < 0.0 or vf > 1.0:
        raise RuntimeError('Volumefraction must be between 0.0 and 1.0.')       

    vm   = 1.0-vf
  
    E1   = fibre.E1*vf+matrix.E1*vm
    E2   = fibre.E2*matrix.E2/(fibre.E2*vm+matrix.E2*vf)
    nu12 = fibre.nu12*vf+matrix.nu12*vm
    G12  = fibre.G12*matrix.G12/(fibre.G12*vm+matrix.G12*vf)

    mat = TransverseIsotropic( [E1,E2] , nu12 , G12 )

    return mat

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

def Macauley( x ):

  '''The Macauley operater. Returns argument x when x > 0. Otherwise 0.'''

  if x > 0:
    return x
  else:
    return 0.
