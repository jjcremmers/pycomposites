"""
Unit tests for example scripts from examples/python directory.

These tests verify that the example scripts execute successfully and
produce expected results.

(c) Joris Remmers (2025)
"""

import unittest
import numpy as np

from composite import (
    TransverseIsotropic, Laminate, mixMaterials,
    stressTransformation, strainTransformation
)


class TestCLTExample1(unittest.TestCase):
    """Test Example 1 - Material Homogenization."""

    def test_homogenization_carbon_epoxy(self) -> None:
        """Test homogenization of carbon fibers in epoxy matrix."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials with 60% fiber volume fraction
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Verify properties are within expected range
        self.assertGreater(udcomp.E1, 100e9)
        self.assertGreater(udcomp.E2, 5e9)
        self.assertGreater(udcomp.G12, 2e9)
        self.assertGreater(udcomp.nu12, 0.1)
        self.assertLess(udcomp.nu12, 0.4)


class TestCLTExample2(unittest.TestCase):
    """Test Example 2 - Q and Qbar Matrices."""

    def test_Q_and_Qbar_matrices(self) -> None:
        """Test computation of Q and Qbar matrices at different angles."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Get Q matrix
        Q = udcomp.getQ()
        self.assertEqual(Q.shape, (3, 3))
        self.assertGreater(Q[0, 0], 0)
        
        # Get Qbar at different angles
        Qbar_20 = udcomp.getQbar(20.0)
        Qbar_minus_20 = udcomp.getQbar(-20.0)
        
        # Verify symmetry properties
        # Q11, Q12, Q22, Q66 should be identical for +/- angles
        self.assertAlmostEqual(Qbar_20[0, 0], Qbar_minus_20[0, 0], delta=1e6)
        self.assertAlmostEqual(Qbar_20[0, 1], Qbar_minus_20[0, 1], delta=1e6)
        self.assertAlmostEqual(Qbar_20[1, 1], Qbar_minus_20[1, 1], delta=1e6)
        self.assertAlmostEqual(Qbar_20[2, 2], Qbar_minus_20[2, 2], delta=1e6)
        
        # Q16 and Q26 should have opposite signs
        self.assertAlmostEqual(Qbar_20[0, 2], -Qbar_minus_20[0, 2], delta=1e6)
        self.assertAlmostEqual(Qbar_20[1, 2], -Qbar_minus_20[1, 2], delta=1e6)


class TestCLTExample3(unittest.TestCase):
    """Test Example 3 - ABD Matrices for Symmetric Laminate."""

    def test_ABD_symmetric_laminate(self) -> None:
        """Test ABD matrices for [0/45/-45/90]s laminate."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 220e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Create laminate
        lam = Laminate()
        lam.addMaterial('UD', udcomp)
        
        # Add layers with symmetric stacking sequence
        orientations = [0., 45., -45., 90., 90., -45., 45., 0.]
        for angle in orientations:
            lam.addLayer('UD', angle, 0.25e-3)
        
        # Get ABD matrices
        A = lam.getA()
        B = lam.getB()
        D = lam.getD()
        
        # Verify shapes
        self.assertEqual(A.shape, (3, 3))
        self.assertEqual(B.shape, (3, 3))
        self.assertEqual(D.shape, (3, 3))
        
        # B matrix should be nearly zero for symmetric laminate
        self.assertTrue(np.allclose(B, np.zeros((3, 3)), atol=1e-6))
        
        # A and D should be positive definite (diagonal elements positive)
        self.assertGreater(A[0, 0], 0)
        self.assertGreater(A[1, 1], 0)
        self.assertGreater(A[2, 2], 0)
        self.assertGreater(D[0, 0], 0)
        self.assertGreater(D[1, 1], 0)
        self.assertGreater(D[2, 2], 0)


class TestCLTExample4(unittest.TestCase):
    """Test Example 4 - Stress Analysis in Layers."""

    def test_stress_recovery_in_layers(self) -> None:
        """Test stress recovery in each layer of a loaded laminate."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Create laminate
        lam = Laminate()
        lam.addMaterial('UD', udcomp)
        
        # Add layers
        orientations = [0., 45., -45., 90., 90., -45., 45., 0.]
        for angle in orientations:
            lam.addLayer('UD', angle, 0.3e-3)
        
        # Define loads: N = [10, 5, 6] kN/m, M = [3, 10, 1] Nm/m
        N = np.array([10e3, 5e3, 6e3])
        M = np.array([3.0, 10.0, 1.0])
        
        # Get inverse matrices
        A1, B1, C1, D1 = lam.getInverseMatrices()
        
        # Calculate midplane strains and curvatures
        eps0 = A1 @ N + B1 @ M
        kappa = C1 @ N + D1 @ M
        
        # Verify results are reasonable (not NaN or infinite)
        self.assertTrue(np.all(np.isfinite(eps0)))
        self.assertTrue(np.all(np.isfinite(kappa)))
        
        # Calculate stresses in first layer (layer 0)
        z = lam.getZcoord(0)
        eps_global = eps0 + z * kappa
        Qbar = lam.getQbar(0)
        sigma_global = Qbar @ eps_global
        
        # Verify stresses are finite
        self.assertTrue(np.all(np.isfinite(sigma_global)))


class TestCLTExample5(unittest.TestCase):
    """Test Example 5 - Effective Properties vs Angle."""

    def test_effective_properties_angle_ply(self) -> None:
        """Test effective properties for [θ/-θ]s laminate at various angles."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Test at a few representative angles
        test_angles = [0, 15, 30, 45, 60, 75, 90]
        
        for theta in test_angles:
            lam = Laminate()
            lam.addMaterial('UD', udcomp)
            
            # Create [θ/-θ]s laminate
            lam.addLayer('UD', theta, 1.5e-3)
            lam.addLayer('UD', -theta, 1.5e-3)
            lam.addLayer('UD', -theta, 1.5e-3)
            lam.addLayer('UD', theta, 1.5e-3)
            
            # Get effective properties
            Ex, Ey, nuxy, Gxy = lam.getElastic()
            
            # Verify all properties are positive and finite
            self.assertGreater(Ex, 0)
            self.assertGreater(Ey, 0)
            self.assertGreater(Gxy, 0)
            self.assertTrue(np.isfinite(Ex))
            self.assertTrue(np.isfinite(Ey))
            self.assertTrue(np.isfinite(nuxy))
            self.assertTrue(np.isfinite(Gxy))
            
            # Verify Poisson's ratio is reasonable
            self.assertGreater(nuxy, 0)

class TestCLTExample6(unittest.TestCase):
    """Test Example 6 - Thermal Curvature."""

    def test_thermal_curvature_unsymmetric_laminate(self) -> None:
        """Test thermal curvature of [0/90] unsymmetric laminate."""
        # Create materials
        carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        # Mix materials
        udcomp = mixMaterials(carbon, epoxy, 0.6)
        
        # Set thermal expansion coefficients
        udcomp.setAlpha([-5e-7, 12e-6])
        
        # Create laminate
        lam = Laminate()
        lam.addMaterial('UD', udcomp)
        lam.addLayer('UD', 0.0, 6e-3)
        lam.addLayer('UD', 90.0, 6e-3)
        
        # Get ABD matrices
        A = lam.getA()
        B = lam.getB()
        D = lam.getD()
        
        # Get thermal forces and moments
        Ts = lam.getTs()
        Tss = lam.getTss()
        
        # Temperature change
        dT = -100.0
        
        # Thermal loads
        NT = Ts * dT
        MT = Tss * dT
        
        # Get inverse matrices
        A1, B1, C1, D1 = lam.getInverseMatrices()
        
        # Calculate thermal strains and curvatures
        eps0_thermal = A1 @ NT + B1 @ MT
        kappa_thermal = C1 @ NT + D1 @ MT
        
        # Verify results are finite
        self.assertTrue(np.all(np.isfinite(eps0_thermal)))
        self.assertTrue(np.all(np.isfinite(kappa_thermal)))
        
        # For unsymmetric laminate, B is not zero
        self.assertFalse(np.allclose(B, np.zeros((3, 3)), atol=1e-6))
        
        # Curvature should be non-zero due to thermal mismatch
        self.assertGreater(np.linalg.norm(kappa_thermal), 1e-10)


class TestStressStrainTransformations(unittest.TestCase):
    """Test stress and strain transformation functions."""

    def test_stress_transformation_identity(self) -> None:
        """Test stress transformation at 0 degrees (identity)."""
        sigma = np.array([100e6, 50e6, 25e6])
        sigma_0 = stressTransformation(sigma, 0.0)
        self.assertTrue(np.allclose(sigma, sigma_0, atol=1e-6))

    def test_stress_transformation_90_degrees(self) -> None:
        """Test stress transformation at 90 degrees."""
        sigma = np.array([100e6, 0.0, 0.0])
        sigma_90 = stressTransformation(sigma, 90.0)
        # σx at 0° becomes σy at 90°
        self.assertAlmostEqual(sigma_90[1], 100e6, delta=1e3)
        self.assertAlmostEqual(sigma_90[0], 0.0, delta=1e3)

    def test_strain_transformation_identity(self) -> None:
        """Test strain transformation at 0 degrees (identity)."""
        eps = np.array([0.001, 0.0005, 0.0002])
        eps_0 = strainTransformation(eps, 0.0)
        self.assertTrue(np.allclose(eps, eps_0, atol=1e-9))

    def test_strain_transformation_90_degrees(self) -> None:
        """Test strain transformation at 90 degrees."""
        eps = np.array([0.001, 0.0, 0.0])
        eps_90 = strainTransformation(eps, 90.0)
        # εx at 0° becomes εy at 90°
        self.assertAlmostEqual(eps_90[1], 0.001, delta=1e-9)
        self.assertAlmostEqual(eps_90[0], 0.0, delta=1e-9)


class TestPlateEquationsExamples(unittest.TestCase):
    """Test examples from PlateEquations directory."""

    def test_basic_plate_stiffness(self) -> None:
        """Test basic plate stiffness calculations."""
        # Create a simple isotropic material
        aluminum = TransverseIsotropic(70e9, 0.33, 26e9)
        
        # Create a single-layer "laminate"
        lam = Laminate()
        lam.addMaterial('Al', aluminum)
        lam.addLayer('Al', 0, 1e-3)
        
        # Get stiffness matrices
        D = lam.getD()
        
        # Verify D matrix is positive definite
        self.assertGreater(D[0, 0], 0)
        self.assertGreater(D[1, 1], 0)
        self.assertGreater(D[2, 2], 0)
        
        # For isotropic material at 0°, D11 should equal D22
        self.assertAlmostEqual(D[0, 0], D[1, 1], delta=D[0, 0] * 0.01)


if __name__ == '__main__':
    unittest.main()
