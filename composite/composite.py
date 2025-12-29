import numpy as np
from numpy import zeros,ones
from numpy.linalg import inv
from math import sin,cos,pi,sqrt,tan,atan

"""
PyComposites - Classical Laminate Theory Analysis Module
=========================================================

A comprehensive Python package for performing Classical Laminate Theory (CLT) 
calculations on composite materials and laminates.

This module provides tools for analyzing transversely isotropic materials and 
composite laminates under plane stress conditions, including:

- Material property definitions and transformations
- Stiffness and compliance matrix calculations  
- Laminate stacking sequence analysis
- Thermal effects and expansion coefficients
- Multiple failure criteria (Maximum Stress, Maximum Strain, Tsai-Wu, Hashin, Larc03)
- Effective laminate properties

Classes
-------
TransverseIsotropic
    Material model for orthotropic/transversely isotropic materials
Layer
    Individual ply definition within a laminate
Laminate
    Multi-layer composite laminate with full CLT analysis

Functions
---------
stressTransformation(sigma, theta)
    Transform stress vectors between coordinate systems
strainTransformation(eps, theta)
    Transform strain vectors between coordinate systems
mixMaterials(fibre, matrix, vf)
    Homogenize fiber and matrix using rule of mixtures
Macauley(x)
    Macauley bracket operator <x>

Examples
--------
>>> # Create a carbon fiber composite material
>>> carbon = TransverseIsotropic([140e9, 10e9], 0.3, 5e9, alpha=[1e-6, 25e-6])
>>> 
>>> # Define a symmetric laminate [0/90]s
>>> lam = Laminate()
>>> lam.addMaterial("carbon", carbon)
>>> lam.addLayer("carbon", 0, 0.125)
>>> lam.addLayer("carbon", 90, 0.125)
>>> lam.addLayer("carbon", 90, 0.125)
>>> lam.addLayer("carbon", 0, 0.125)
>>> 
>>> # Calculate effective properties
>>> Ex, Ey, nuxy, Gxy = lam.getElastic()
>>> print(f"Ex = {Ex/1e9:.2f} GPa, Ey = {Ey/1e9:.2f} GPa")

Notes
-----
This module was developed for the course:
Composite and Lightweight Materials - Design and Analysis (4MM00)

Author: Joris Remmers
Copyright (c) 2013-2025
"""

class TransverseIsotropic:

    """
    Material model representing a transversely isotropic (or orthotropic) solid
    in plane stress conditions.

    Transversely isotropic materials are isotropic in one plane and have
    distinct properties along a single axis of symmetry (e.g. fiber-reinforced
    composites). This class provides a complete framework for analyzing such
    materials including elastic properties, thermal effects, and failure criteria.

    Parameters
    ----------
    E : float or list of float
        Young's modulus (Pa). If a list is given:
          - [E1, E2] for orthotropic materials (fiber and transverse directions)
          - [E] for isotropic materials
        If a single float is given, both E1 and E2 are set to that value.
    nu12 : float
        Poisson's ratio ν₁₂ (dimensionless, typically 0.2-0.35).
    G12 : float, optional
        In-plane shear modulus G₁₂ (Pa). Default is 0.
    alpha : float or list of float, optional
        Thermal expansion coefficient(s) (1/K or 1/°C). Can be:
          - single float (applied to both α₁ and α₂)
          - [α1, α2] for different longitudinal and transverse values
        Default is 0.
    rho : float, optional
        Material density (kg/m³). Default is 0.

    Attributes
    ----------
    E1, E2 : float
        Young's moduli in fiber and transverse directions
    nu12, nu21 : float
        Major and minor Poisson's ratios
    G12 : float
        In-plane shear modulus
    alpha1, alpha2 : float
        Thermal expansion coefficients
    rho : float
        Material density
    Q : ndarray
        Reduced stiffness matrix (computed on first access)
    U, V : ndarray
        Stiffness and compliance invariants (computed on first access)

    Methods
    -------
    getQ()
        Compute reduced stiffness matrix Q
    getS()
        Compute compliance matrix S
    getQbar(theta)
        Compute rotated stiffness matrix Q̅
    getSbar(theta)
        Compute rotated compliance matrix S̅
    getAlpha(theta)
        Compute rotated thermal expansion vector
    setFailureProperties(F, Gfrac, alpha0deg)
        Define failure properties for failure criteria
    getFIMaximumStress(sigma)
        Evaluate Maximum Stress failure criterion
    getFIMaximumStrain(sigma)
        Evaluate Maximum Strain failure criterion
    getFITsaiWu(sigma)
        Evaluate Tsai-Wu failure criterion
    getFIHashin73(sigma)
        Evaluate Hashin 1973 failure criterion
    getFIHashin80(sigma)
        Evaluate Hashin 1980 failure criterion
    getFILarc03(sigma)
        Evaluate Larc03 failure criterion

    Examples
    --------
    Create a carbon fiber/epoxy composite material:
    
    >>> mat = TransverseIsotropic(
    ...     E=[140e9, 10e9],      # E1=140 GPa, E2=10 GPa
    ...     nu12=0.3,              # Poisson's ratio
    ...     G12=5e9,               # Shear modulus 5 GPa
    ...     alpha=[1e-6, 25e-6],   # Thermal expansion coefficients
    ...     rho=1600               # Density in kg/m³
    ... )
    
    Define failure properties:
    
    >>> mat.setFailureProperties(
    ...     F=[2000e6, 1200e6, 50e6, 200e6, 80e6],  # Xt, Xc, Yt, Yc, S
    ...     Gfrac=[500, 1000],                       # GIc, GIIc
    ...     alpha0deg=53                             # Fracture angle
    ... )
    
    Calculate stiffness matrix at 45°:
    
    >>> Qbar = mat.getQbar(45.0)
    >>> print(Qbar)
    
    Check failure at a stress state:
    
    >>> import numpy as np
    >>> sigma = np.array([500e6, 20e6, 10e6])  # [σ1, σ2, τ12]
    >>> FI = mat.getFITsaiWu(sigma)
    >>> if FI >= 1.0:
    ...     print("Material has failed!")

    Notes
    -----
    - All stress/strain calculations assume plane stress conditions
    - Failure indices > 1.0 indicate failure
    - Material axes: 1 = fiber direction, 2 = transverse direction
    - Angles are measured in degrees, positive counterclockwise
    - The relationship ν₂₁ = (E₂/E₁)ν₁₂ is automatically enforced
    
    See Also
    --------
    Laminate : Multi-layer laminate analysis
    mixMaterials : Homogenization using rule of mixtures
    """

    def __init__(self, E: float | list[float], nu12: float, 
                 G12: float = 0.0, alpha: float | list[float] = 0.0, 
                 rho: float = 0.0) -> None:
                 
        """
        Initialize the TransverseIsotropic material model.
        
        Parameters
        ----------
        E : float or list of float
            Young's modulus (Pa). Either scalar or [E1, E2].
        nu12 : float
            Major Poisson's ratio (dimensionless).
        G12 : float, optional
            Shear modulus (Pa). Default is 0.
        alpha : float or list of float, optional
            Thermal expansion coefficient(s) (1/K). Default is 0.
        rho : float, optional
            Density (kg/m³). Default is 0.
            
        Raises
        ------
        ValueError
            If E or alpha lists have invalid length (not 1 or 2 elements).
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
            Thermal expansion coefficients (1/K or 1/°C). Can be:
              - single float (applied to both α₁ and α₂)
              - [α1, α2] for different longitudinal and transverse values

        Raises
        ------
        ValueError
            If a list with more than two elements is provided.
            
        Examples
        --------
        >>> mat = TransverseIsotropic(140e9, 0.3, 5e9)
        >>> mat.setAlpha(1e-6)              # Set both to 1e-6
        >>> mat.setAlpha([1e-6, 25e-6])     # Set α₁=1e-6, α₂=25e-6
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
        Set material failure properties for failure criteria evaluation.

        Parameters
        ----------
        F : list of float
            Failure strength parameters (Pa). Expected length:
              - 5 values: [Xt, Xc, Yt, Yc, Sl]
                * Xt: Longitudinal tensile strength
                * Xc: Longitudinal compressive strength
                * Yt: Transverse tensile strength
                * Yc: Transverse compressive strength
                * Sl: In-plane shear strength
              - 6 values: [Xt, Xc, Yt, Yc, Sl, f12]
                * f12: Tsai-Wu interaction coefficient (dimensionless)
        Gfrac : float or list of float, optional
            Fracture toughness values (J/m² or N/m):
              - scalar: same value for GIc and GIIc
              - [GIc, GIIc]: Mode I and Mode II fracture toughness
            Default is 0.
        alpha0deg : float, optional
            Fracture plane angle (degrees) for Larc03 model. Default is 53°.

        Raises
        ------
        ValueError
            If F does not have 5 or 6 elements, or Gfrac list is invalid.
            
        Examples
        --------
        >>> mat = TransverseIsotropic(140e9, 0.3, 5e9)
        >>> # Set basic failure properties
        >>> mat.setFailureProperties([2000e6, 1200e6, 50e6, 200e6, 80e6])
        >>> 
        >>> # Include Tsai-Wu interaction coefficient
        >>> mat.setFailureProperties(
        ...     [2000e6, 1200e6, 50e6, 200e6, 80e6, -0.5],
        ...     Gfrac=[500, 1000],
        ...     alpha0deg=53
        ... )

        Notes
        -----
        - All strength values should be positive (absolute values)
        - GIc and GIIc are used in Larc03 failure criterion
        - f12 is the interaction term in Tsai-Wu criterion, typically negative
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
        Set interlaminar shear strength for Larc03 failure criterion.

        Parameters
        ----------
        SLis : float
            Interlaminar shear strength (Pa).
            
        Examples
        --------
        >>> mat = TransverseIsotropic(140e9, 0.3, 5e9)
        >>> mat.setFailureProperties([2000e6, 1200e6, 50e6, 200e6, 80e6])
        >>> mat.setSLis(60e6)  # 60 MPa interlaminar shear strength
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
        Compute reduced stiffness matrix Q in material coordinates (plane stress).

        The reduced stiffness matrix relates stresses to strains in the material
        principal axes: {σ} = [Q]{ε}
        
        Returns
        -------
        np.ndarray
            3x3 reduced stiffness matrix Q with structure:
            [[Q11, Q12,   0],
             [Q12, Q22,   0],
             [  0,   0, Q66]]
            Units: Pa
            
        Notes
        -----
        The matrix is cached after first computation for efficiency.
        Matrix is symmetric with Q12 = Q21.
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
        Compute reduced compliance matrix S in material coordinates (plane stress).

        The compliance matrix relates strains to stresses in the material
        principal axes: {ε} = [S]{σ}
        
        Returns
        -------
        np.ndarray
            3x3 compliance matrix S = Q⁻¹. Units: 1/Pa
            
        Notes
        -----
        S is the inverse of the stiffness matrix Q.
        Matrix is symmetric with S12 = S21.
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
        Compute rotated (global) stiffness matrix Q̅ for a given ply angle.

        Transforms the reduced stiffness matrix from material coordinates
        (1-2) to global coordinates (x-y) using the transformation angle θ.

        Parameters
        ----------
        theta : float
            Ply orientation angle (degrees). Measured counterclockwise from
            global x-axis to material 1-axis.

        Returns
        -------
        np.ndarray
            3x3 rotated stiffness matrix Q̅(θ) in global coordinates. Units: Pa
            
        Examples
        --------
        >>> mat = TransverseIsotropic([140e9, 10e9], 0.3, 5e9)
        >>> Qbar_0 = mat.getQbar(0)      # 0° ply
        >>> Qbar_90 = mat.getQbar(90)    # 90° ply
        >>> Qbar_45 = mat.getQbar(45)    # 45° ply
        
        Notes
        -----
        Uses stiffness invariants (U) for efficient computation.
        At θ=0°, Q̅ = Q (material and global axes aligned).
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
    Individual composite ply (layer) within a laminate.

    A Layer represents a single ply with specific material, orientation,
    and thickness. Multiple layers are combined to form a Laminate.

    Parameters
    ----------
    name : str
        Identifier for the material model (must match a material added
        to the parent Laminate).
    theta : float
        Ply orientation angle (degrees). Measured counterclockwise from
        global x-axis to material fiber direction.
    thick : float
        Ply thickness (m or consistent units).

    Attributes
    ----------
    name : str
        Material identifier
    theta : float
        Ply angle in degrees
    thick : float
        Ply thickness
        
    Examples
    --------
    >>> layer = Layer("carbon", 45.0, 0.000125)  # 45° ply, 0.125mm thick
    
    Notes
    -----
    Layers are typically created through Laminate.addLayer() rather than
    directly instantiating this class.
    
    See Also
    --------
    Laminate.addLayer : Add a layer to a laminate
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
  
    """
    Multi-layer composite laminate for Classical Laminate Theory analysis.

    A Laminate consists of multiple layers (plies) stacked together, each
    with potentially different materials and orientations. This class provides
    comprehensive CLT analysis including stiffness matrices, effective properties,
    thermal effects, and layer-by-layer stress/strain recovery.

    Attributes
    ----------
    materials : dict
        Dictionary of material models, keyed by name
    layers : list of Layer
        Ordered list of layers from bottom to top
    thick : float
        Total laminate thickness (m)
    h : ndarray
        z-coordinates of layer interfaces from laminate midplane (m)
    A, B, D : ndarray (computed on demand)
        Extensional, coupling, and bending stiffness matrices

    Methods
    -------
    addMaterial(name, mat)
        Add a material model to the laminate material library
    addLayer(name, theta, thick)
        Add a ply to the laminate stacking sequence
    removeAllLayers()
        Clear all layers from the laminate
    getA(), getB(), getD()
        Compute stiffness matrices [A], [B], [D]
    getTs(), getTss()
        Compute thermal force and moment resultants
    getInverseMatrices()
        Compute inverse stiffness matrices [A*], [B*], [C*], [D*]
    getElastic()
        Compute effective laminate elastic properties
    getRhoh()
        Compute areal mass density
    getZcoord(j)
        Get z-coordinate of layer j centroid
    getLayerBounds(j)
        Get top and bottom z-coordinates of layer j
    getQbar(j)
        Get rotated stiffness matrix for layer j

    Examples
    --------
    Create a quasi-isotropic [0/±60]s laminate:
    
    >>> # Define material
    >>> mat = TransverseIsotropic([140e9, 10e9], 0.3, 5e9)
    >>> 
    >>> # Build laminate
    >>> lam = Laminate()
    >>> lam.addMaterial("carbon", mat)
    >>> lam.addLayer("carbon", 0, 0.125e-3)
    >>> lam.addLayer("carbon", 60, 0.125e-3)
    >>> lam.addLayer("carbon", -60, 0.125e-3)
    >>> lam.addLayer("carbon", -60, 0.125e-3)
    >>> lam.addLayer("carbon", 60, 0.125e-3)
    >>> lam.addLayer("carbon", 0, 0.125e-3)
    >>> 
    >>> # Get stiffness matrices
    >>> A = lam.getA()  # Extensional stiffness
    >>> B = lam.getB()  # Coupling stiffness (should be ~0 for symmetric)
    >>> D = lam.getD()  # Bending stiffness
    >>> 
    >>> # Get effective properties
    >>> Ex, Ey, nuxy, Gxy = lam.getElastic()
    >>> print(f"Ex = {Ex/1e9:.1f} GPa")
    >>> print(f"Ey = {Ey/1e9:.1f} GPa")

    Notes
    -----
    - Layers are added from bottom (-z) to top (+z)
    - The laminate midplane is at z=0
    - For symmetric laminates, B ≈ 0 (no extension-bending coupling)
    - Units must be consistent throughout (typically SI: Pa, m, kg)
    
    See Also
    --------
    TransverseIsotropic : Material model for individual plies
    Layer : Individual ply definition
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

        self.h     = np.zeros( len(self.layers)+1 )
        self.thick = 0.
    
        for i,layer in enumerate(self.layers):
            self.h[i+1] = self.thick+layer.thick
            self.thick += layer.thick

        self.h += -0.5*self.thick*np.ones( len(self.h) )

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

        self.A = np.zeros( shape = ( 3,3) )

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
        self.B = np.zeros( shape = ( 3,3) )

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
        self.D = np.zeros( shape = ( 3,3) )

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
       
        self.Ts = np.zeros( 3 )

        for i,layer in enumerate(self.layers):
            name  = layer.name
            theta = layer.theta

            Qalpha = self.materials[name].getQbar( theta ) @ self.materials[name].getAlpha( theta )

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

        self.Tss = np.zeros( 3 )

        for i,layer in enumerate(self.layers):
            name  = layer.name
            theta = layer.theta

            Qalpha = self.materials[name].getQbar( theta ) @ self.materials[name].getAlpha( theta )

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
        Dstar = self.D - self.B @ Ainv @ self.B
        Dsinv = inv(Dstar)

        self.A1 = Ainv + Ainv @ self.B @ Dsinv @ self.B @ Ainv
        self.B1 = -Ainv @ self.B @ Dsinv
        self.C1 = self.B1.T
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
    Transform stress components between coordinate systems.

    Rotates stress tensor from one coordinate system to another using
    the transformation angle θ.

    Parameters
    ----------
    sigma : np.ndarray
        Stress vector [σ₁, σ₂, τ₁₂] in the original coordinate system (Pa).
    theta : float
        Rotation angle (degrees). Positive counterclockwise.

    Returns
    -------
    np.ndarray
        Transformed stress vector [σ'₁, σ'₂, τ'₁₂] in rotated coordinates (Pa).
        
    Examples
    --------
    >>> import numpy as np
    >>> # Pure σx stress
    >>> sigma_x = np.array([100e6, 0, 0])
    >>> # Rotate by 45°
    >>> sigma_45 = stressTransformation(sigma_x, 45)
    >>> print(sigma_45)  # [50e6, 50e6, 25e6]
    
    Notes
    -----
    Uses the standard stress transformation equations:
    σ'₁ = σ₁cos²θ + σ₂sin²θ + 2τ₁₂sinθcosθ
    σ'₂ = σ₁sin²θ + σ₂cos²θ - 2τ₁₂sinθcosθ  
    τ'₁₂ = (σ₂ - σ₁)sinθcosθ + τ₁₂(cos²θ - sin²θ)
    
    See Also
    --------
    strainTransformation : Transform strain components
    """
       
    signew = np.zeros( 3 )

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
    Transform strain components between coordinate systems.

    Rotates strain tensor from one coordinate system to another using
    the transformation angle θ. Note the factor of 2 difference in
    shear strain transformation compared to stress.

    Parameters
    ----------
    eps : np.ndarray
        Strain vector [ε₁, ε₂, γ₁₂] in the original coordinate system.
        Engineering shear strain γ = 2ε₁₂ (dimensionless).
    theta : float
        Rotation angle (degrees). Positive counterclockwise.

    Returns
    -------
    np.ndarray
        Transformed strain vector [ε'₁, ε'₂, γ'₁₂] in rotated coordinates.
        
    Examples
    --------
    >>> import numpy as np
    >>> eps = np.array([0.001, 0.0, 0.0])  # Pure ε₁ strain
    >>> eps_90 = strainTransformation(eps, 90)
    >>> print(eps_90)  # [0.0, 0.001, 0.0]
    
    Notes
    -----
    Uses standard strain transformation with engineering shear strain:
    ε'₁ = ε₁cos²θ + ε₂sin²θ + γ₁₂sinθcosθ
    ε'₂ = ε₁sin²θ + ε₂cos²θ - γ₁₂sinθcosθ
    γ'₁₂ = 2(ε₂ - ε₁)sinθcosθ + γ₁₂(cos²θ - sin²θ)
    
    See Also
    --------
    stressTransformation : Transform stress components
    """

    epsnew = np.zeros( 3 )

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
    Homogenize fiber and matrix into equivalent ply using micromechanics.

    Applies classical rule of mixtures and inverse rule of mixtures to
    compute effective ply properties from constituent fiber and matrix
    properties.

    Parameters
    ----------
    fibre : TransverseIsotropic
        Fiber material properties (typically stiffer, e.g., carbon, glass).
    matrix : TransverseIsotropic
        Matrix material properties (typically softer, e.g., epoxy, PEEK).
    vf : float
        Fiber volume fraction, 0 ≤ vf ≤ 1 (dimensionless).
        For example, vf=0.6 means 60% fibers, 40% matrix.

    Returns
    -------
    TransverseIsotropic
        Homogenized composite material with effective properties:
        - E₁: longitudinal modulus (rule of mixtures)
        - E₂: transverse modulus (inverse rule of mixtures)
        - ν₁₂: Poisson's ratio (rule of mixtures)
        - G₁₂: shear modulus (inverse rule of mixtures)

    Raises
    ------
    RuntimeError
        If vf is outside the valid range [0, 1].
        
    Examples
    --------
    >>> # Carbon fiber properties
    >>> carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
    >>> 
    >>> # Epoxy matrix properties  
    >>> epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
    >>> 
    >>> # Create composite with 60% fiber volume fraction
    >>> composite = mixMaterials(carbon, epoxy, 0.6)
    >>> print(f"E1 = {composite.E1/1e9:.1f} GPa")
    >>> print(f"E2 = {composite.E2/1e9:.1f} GPa")

    Notes
    -----
    Micromechanics equations used:
    - E₁ = Ef·vf + Em·vm (rule of mixtures)
    - 1/E₂ = vf/Ef + vm/Em (inverse rule of mixtures)
    - ν₁₂ = νf·vf + νm·vm
    - 1/G₁₂ = vf/Gf + vm/Gm
    
    where vm = 1 - vf is the matrix volume fraction.
    
    These are simplified analytical models. For more accurate predictions,
    consider finite element micromechanics or semi-empirical models.
    
    See Also
    --------
    TransverseIsotropic : Material model class
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
    Macauley bracket operator ⟨x⟩ = max(x, 0).

    Returns the positive part of a value, commonly used in damage mechanics
    and failure criteria to distinguish tension from compression.

    Parameters
    ----------
    x : float
        Input value (any real number).

    Returns
    -------
    float
        x if x > 0, else 0.
        
    Examples
    --------
    >>> Macauley(5.0)
    5.0
    >>> Macauley(-3.0)
    0.0
    >>> Macauley(0.0)
    0.0
    
    Notes
    -----
    Also known as the ramp function or positive part function.
    Mathematical definition: ⟨x⟩ = (x + |x|)/2 = max(x, 0)
    
    Used in Larc03 and other failure criteria to handle compression
    vs. tension differently.
    """

    if x > 0:
        return x
    else:
        return 0.
