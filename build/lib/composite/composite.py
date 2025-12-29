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
    Material model representing a transversely isotropic (or orthotropic) solid
    in plane stress conditions.

    Transversely isotropic materials are isotropic in one plane and have
    distinct properties along a single axis of symmetry (e.g. fiber-reinforced
    composites).

    Parameters
    ----------
    E : float or list of float
        Young's modulus. If a list is given:
          - [E1, E2] for orthotropic materials
          - [E] for isotropic materials
        If a single float is given, both E1 and E2 are set to that value.
    nu12 : float
        Poisson's ratio ν₁₂.
    G12 : float, optional
        In-plane shear modulus. Default is 0. If not specified, it is assumed
        isotropic and computed as G = E / [2(1 + ν)].
    alpha : float or list of float, optional
        Thermal expansion coefficient(s). Either:
          - single float, applied in both directions, or
          - [α1, α2]. Default is 0.
    rho : float, optional
        Density. Default is 0.

    Notes
    -----
    Provides methods to:
      - compute stiffness and compliance matrices (Q, S)
      - compute rotated properties (Q̅, S̅)
      - evaluate failure indices (Maximum Stress, Maximum Strain, Tsai-Wu,
        Hashin, Larc03).
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
   
    def setAlpha(self, alpha: float | list[float]) -> None:
    
        """
        Set thermal expansion coefficients.

        Parameters
        ----------
        alpha : float or list of float
            Thermal expansion coefficients. Can be:
              - single float (applied to both directions)
              - [α1, α2]

        Raises
        ------
        ValueError
            If a list with more than two elements is provided.
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def setFailureProperties(self, F: list[float], Gfrac: list[float] = 0.0, alpha0deg: float = 53.0) -> None:

        """
        Set material failure properties.

        Parameters
        ----------
        F : list of float
            Failure parameters. Expected length:
              - 5 values: [Xt, Xc, Yt, Yc, Sl]
              - 6 values: [Xt, Xc, Yt, Yc, Sl, f12] (for Tsai-Wu)
        Gfrac : float or list of float, optional
            Fracture toughness values [GIc, GIIc]. Default is 0.
        alpha0deg : float, optional
            Fiber orientation angle (degrees) for Larc03 model. Default is 53°.

        Raises
        ------
        ValueError
            If F does not have length 5 or 6, or if Gfrac is invalid.
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
        Set interlaminar shear strength (for Larc03).

        Parameters
        ----------
        SLis : float
            Interlaminar shear strength.
        """
        
        self.SLis = SLis

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __str__(self) -> str:
        
        """
        Return formatted string with material properties.

        Returns
        -------
        str
            Elastic, thermal, and strength properties (if defined).
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

            if hasattr(self, "GIc") and hasattr(self, "GIIc"):
                msg += f"  GIc    :  {self.GIc:12.3e} , GIIc   :  {self.GIIc:12.3e} \n"
                msg += f"  alpha0 :  {self.a0deg:12.3e}\n"

        return msg

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getQ(self) -> np.ndarray:
        
        """
        Compute reduced stiffness matrix Q.

        Returns
        -------
        np.ndarray
            3x3 reduced stiffness matrix Q.
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
        Compute stiffness invariants U.

        Returns
        -------
        np.ndarray
            Array [U1..U5].
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
        Compute reduced compliance matrix S.

        Returns
        -------
        np.ndarray
            3x3 compliance matrix S.
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
        Compute compliance invariants V.

        Returns
        -------
        np.ndarray
            Array [V1..V5].
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
        Compute global stiffness matrix Q̅ for a ply angle.

        Parameters
        ----------
        theta : float
            Ply angle (degrees).

        Returns
        -------
        np.ndarray
            3x3 rotated stiffness matrix Q̅.
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
        Compute global compliance matrix S̅ for a ply angle.

        Parameters
        ----------
        theta : float
            Ply angle (degrees).

        Returns
        -------
        np.ndarray
            3x3 rotated compliance matrix S̅.
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
        Sbar[2, 2] = self.V[4] - self.V[2] * c4

        return Sbar

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getAlpha(self, theta: float) -> np.ndarray:
        
        """
        Compute thermal expansion vector in global axes.

        Parameters
        ----------
        theta : float
            Ply angle (degrees).

        Returns
        -------
        np.ndarray
            Vector [αx, αy, αxy].
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
        
        """
        Compute failure index (Maximum Stress criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getFIMaximumStrain(self, sigma: np.ndarray ) -> float:
        
        """
        Compute failure index (Maximum Strain criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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
        
        """
        Compute failure index (Tsai-Wu criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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
        
        """
        Compute failure index (Hashin 1973 criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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
        
        """
        Compute failure index (Hashin 1980 criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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
        
        """
        Compute failure index (Larc03 criterion).

        Parameters
        ----------
        sigma : np.ndarray
            Stress vector [σx, σy, τxy].

        Returns
        -------
        float
            Failure index.
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

    """
    Composite layer definition.

    Parameters
    ----------
    name : str
        Identifier for the material.
    theta : float
        Ply orientation (degrees).
    thick : float
        Ply thickness.
    """

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __init__(self, name: str, theta: float, thick: float) -> None:

   
        """
        Initialize a composite layer.

        Parameters
        ----------
        name : str
            Identifier for the material model.
        theta : float
            Ply angle (degrees).
        thick : float
            Ply thickness.
        """
      
        self.name  = name
        self.theta = theta
        self.thick = thick

#===============================================================================
#  Class Laminate
#===============================================================================

class Laminate:
  
    """"
    Composite laminate consisting of multiple layers.

    Provides methods to:
      - add materials and layers
      - compute A, B, D stiffness matrices
      - compute effective laminate properties
      - extract geometry and z-coordinates
    """

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
 
    def __init__( self ) -> None: 

        """
        Initialize an empty laminate.
        """
        
        self.materials = {}
        self.layers    = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def __str__( self ) -> str:
  
        """
        Return formatted string with laminate stacking information.

        Returns
        -------
        str
            Table of layer indices, thicknesses, orientations, and materials.
        """

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

    def addMaterial(self, name: str, mat: TransverseIsotropic) -> None:


        """
        Add a material to the laminate.

        Parameters
        ----------
        name : str
            Identifier for the material model.
        mat : TransverseIsotropic
            Instance of the material model.

        Examples
        --------
        >>> mat = TransverseIsotropic(70e9, 0.2, 5e9)
        >>> lam = Laminate()
        >>> lam.addMaterial("glassfibre", mat)
        """
         
        self.materials[name] = mat

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def addLayer(self, name: str, theta: float, thick: float) -> None:
        
        """
        Add a ply (layer) to the laminate.

        Parameters
        ----------
        name : str
            Material identifier.
        theta : float
            Ply angle (degrees).
        thick : float
            Ply thickness.
        """
                
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
  
    def removeAllLayers(self) -> None:
        
        """
        Remove all layers from the laminate.
        """
    
        self.layers = []

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getA(self) -> np.ndarray:
        
        """
        Compute extensional stiffness matrix A.

        Returns
        -------
        np.ndarray
            3x3 matrix A.
        """

        self.A = zeros( shape = ( 3,3) )

        for i,layer in enumerate(self.layers):
            name  = layer.name
            theta = layer.theta

            self.A += self.materials[name].getQbar( theta ) * (self.h[i+1]-self.h[i])

        return self.A

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getB(self) -> np.ndarray:
        
        """
        Compute coupling stiffness matrix B.

        Returns
        -------
        np.ndarray
            3x3 matrix B.
        """
        self.B = zeros( shape = ( 3,3) )

        for i,layer in enumerate(self.layers):
            name  = layer.name
            theta = layer.theta

            self.B += 0.5*self.materials[name].getQbar( theta ) * (self.h[i+1]**2-self.h[i]**2)

        return self.B

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getD(self) -> np.ndarray:
    
        """
        Compute bending stiffness matrix D.

        Returns
        -------
        np.ndarray
            3x3 matrix D.
        """
        self.D = zeros( shape = ( 3,3) )

        for i,layer in enumerate(self.layers):
            name  = layer.name
            theta = layer.theta

            self.D += 1.0/3.0*self.materials[name].getQbar( theta ) * (self.h[i+1]**3-self.h[i]**3)

        return self.D

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getTs(self) -> np.ndarray:
        
        """
        Compute thermal resultants T*.

        Returns
        -------
        np.ndarray
            3x1 vector of thermal forces.
        """
       
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

    def getTss(self) -> np.ndarray:
        
        """
        Compute thermal resultants T**.

        Returns
        -------
        np.ndarray
            3x1 vector of thermal moments.
        """

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

    def getRhoh(self) -> float:
        
        """
        Compute mass per unit area of the laminate.

        Returns
        -------
        float
            Areal mass density (ρh).
        """
    
        rhoh = 0.

        for i,layer in enumerate(self.layers):
            name  = layer.name

            rhoh += self.materials[name].rho * (self.h[i+1]-self.h[i])

        return rhoh

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getZcoord(self, j: int) -> float:
        """
        Get z-coordinate of the centroid of a given layer.

        Parameters
        ----------
        j : int
            Layer index.

        Returns
        -------
        float
            z-coordinate of the layer centroid.
        """
         
        return 0.5*(self.h[j] + self.h[j+1])  

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getLayerBounds(self, j: int) -> tuple[float, float]:
        """
        Get top and bottom z-coordinates of a given layer.

        Parameters
        ----------
        j : int
            Layer index.

        Returns
        -------
        tuple[float, float]
            (z_top, z_bottom).
        """

        return (self.h[j], self.h[j+1])

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
    
    def getInverseMatrices(self) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        
        """
        Compute inverse stiffness matrices A1, B1, C1, D1.

        Returns
        -------
        tuple of np.ndarray
            A1, B1, C1, D1 (all 3x3).
        """
        
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

    def getQbar(self, j: int) -> np.ndarray:
        
        """
        Get Q̅ matrix for a specific layer.

        Parameters
        ----------
        j : int
            Layer index.

        Returns
        -------
        np.ndarray
            3x3 Q̅ matrix for the layer.
        """
         
        name  = self.layers[j].name
        theta = self.layers[j].theta

        return self.materials[name].getQbar( theta )

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------

    def getElastic(self) -> list[float]:
    
        """
        Compute effective laminate properties.

        Returns
        -------
        list of float
            [Ex, Ey, νxy, Gxy].
        """
    
        self.getA()

        Ex   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[1,1])
        Ey   = (self.A[0,0]*self.A[1,1]-self.A[0,1]*self.A[0,1])/(self.thick*self.A[0,0])
        nuxy = self.A[0,1]/self.A[1,1]
        Gxy  = self.A[2,2]/self.thick

        return [Ex,Ey,nuxy,Gxy]

#==============================================================================
#  Utility functions
#==============================================================================

def stressTransformation(sigma: np.ndarray, theta: float) -> np.ndarray:

    """
    Transform stresses from global (x–y) to local (1–2) axes.

    Parameters
    ----------
    sigma : np.ndarray
        Stress vector [σ1, σ2, τ12] in local coordinates.
    theta : float
        Ply angle (degrees).

    Returns
    -------
    np.ndarray
        Stress vector [σx, σy, τxy] in global coordinates.
    """
       
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

def strainTransformation(eps: np.ndarray, theta: float) -> np.ndarray:

    """
    Transform strains from global (x–y) to local (1–2) axes.

    Parameters
    ----------
    eps : np.ndarray
        Strain vector [ε1, ε2, γ12] in local coordinates.
    theta : float
        Ply angle (degrees).

    Returns
    -------
    np.ndarray
        Strain vector [εx, εy, γxy] in global coordinates.
    """

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
  

def mixMaterials(fibre: TransverseIsotropic,
                 matrix: TransverseIsotropic,
                 vf: float) -> TransverseIsotropic:
                 
    """
    Homogenize fiber and matrix into an equivalent ply (rule of mixtures).

    Parameters
    ----------
    fibre : TransverseIsotropic
        Fiber material.
    matrix : TransverseIsotropic
        Matrix material.
    vf : float
        Fiber volume fraction (0 ≤ vf ≤ 1).

    Returns
    -------
    TransverseIsotropic
        Homogenized material model.

    Raises
    ------
    RuntimeError
        If vf is outside [0, 1].
    """
       
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

def Macauley(x: float) -> float:

    """
    Macauley operator ⟨x⟩.

    Parameters
    ----------
    x : float
        Input value.

    Returns
    -------
    float
        x if x > 0, else 0.
    """

    if x > 0:
        return x
    else:
        return 0.
