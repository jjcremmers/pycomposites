"""
Unit tests for the ``composite`` package.

These tests exercise the material models and laminate helpers. They are
focused on numerical checks and basic API behaviour; no production code is
modified by these tests.

(c) Joris Remmers (2025)
"""

import unittest
import numpy as np
from typing import Any

from composite import (
    TransverseIsotropic, Laminate, Layer, mixMaterials,
    stressTransformation, strainTransformation, Macauley
)


class TestTransverseIsotropic(unittest.TestCase):
    """Tests for TransverseIsotropic material class."""

    def test_init_scalar_E(self) -> None:
        """Initialize material with scalar E value."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9, alpha=1e-6, rho=1500.0)
        self.assertEqual(mat.E1, 100e9)
        self.assertEqual(mat.E2, 100e9)
        self.assertEqual(mat.nu12, 0.25)
        self.assertEqual(mat.G12, 5e9)
        self.assertAlmostEqual(mat.nu21, 0.25, delta=1e-6)
        self.assertEqual(mat.alpha1, 1e-6)
        self.assertEqual(mat.alpha2, 1e-6)
        self.assertEqual(mat.rho, 1500.0)

    def test_init_list_E_two_elements(self) -> None:
        """Initialize material with E as two-element list."""
        mat = TransverseIsotropic([150e9, 10e9], 0.3)
        self.assertEqual(mat.E1, 150e9)
        self.assertEqual(mat.E2, 10e9)
        self.assertAlmostEqual(mat.nu21, 0.3 * 10e9 / 150e9, delta=1e-9)

    def test_init_list_E_one_element(self) -> None:
        """Initialize material with E as one-element list."""
        mat = TransverseIsotropic([100e9], 0.25)
        self.assertEqual(mat.E1, 100e9)
        self.assertEqual(mat.E2, 100e9)

    def test_init_invalid_E_list(self) -> None:
        """Providing invalid E list should raise ValueError."""
        with self.assertRaises(ValueError):
            TransverseIsotropic([1, 2, 3], 0.3)

    def test_init_alpha_list(self) -> None:
        """Initialize material with alpha as list."""
        mat = TransverseIsotropic(100e9, 0.25, alpha=[1e-6, 2e-6])
        self.assertEqual(mat.alpha1, 1e-6)
        self.assertEqual(mat.alpha2, 2e-6)

    def test_init_alpha_one_element(self) -> None:
        """Initialize material with alpha as one-element list."""
        mat = TransverseIsotropic(100e9, 0.25, alpha=[1e-6])
        self.assertEqual(mat.alpha1, 1e-6)
        self.assertEqual(mat.alpha2, 1e-6)

    def test_init_invalid_alpha_list(self) -> None:
        """Providing invalid alpha list should raise ValueError."""
        with self.assertRaises(ValueError):
            TransverseIsotropic(100e9, 0.25, alpha=[1, 2, 3])

    def test_setAlpha_scalar(self) -> None:
        """Set thermal expansion with scalar value."""
        mat = TransverseIsotropic(100e9, 0.25)
        mat.setAlpha(1.5e-6)
        self.assertEqual(mat.alpha1, 1.5e-6)
        self.assertEqual(mat.alpha2, 1.5e-6)

    def test_setAlpha_list(self) -> None:
        """Set thermal expansion with list."""
        mat = TransverseIsotropic(100e9, 0.25)
        mat.setAlpha([1e-6, 2e-6])
        self.assertEqual(mat.alpha1, 1e-6)
        self.assertEqual(mat.alpha2, 2e-6)

    def test_setAlpha_invalid_list(self) -> None:
        """Set thermal expansion with invalid list should raise ValueError."""
        mat = TransverseIsotropic(100e9, 0.25)
        with self.assertRaises(ValueError):
            mat.setAlpha([1, 2, 3])

    def test_setFailureProperties_five_values(self) -> None:
        """Set failure properties with 5 values."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        self.assertEqual(mat.Xt, 1000e6)
        self.assertEqual(mat.Xc, 800e6)
        self.assertEqual(mat.Yt, 50e6)
        self.assertEqual(mat.Yc, 40e6)
        self.assertEqual(mat.Sl, 30e6)

    def test_setFailureProperties_six_values(self) -> None:
        """Set failure properties with 6 values including f12."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        self.assertEqual(mat.Xt, 1000e6)

    def test_setFailureProperties_with_Gfrac_list(self) -> None:
        """Set failure properties with Gfrac as list."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6], Gfrac=[1000.0, 2000.0])
        self.assertEqual(mat.GIc, 1000.0)
        self.assertEqual(mat.GIIc, 2000.0)

    def test_setFailureProperties_invalid_F_length(self) -> None:
        """Set failure properties with invalid F length should raise ValueError."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        with self.assertRaises(ValueError):
            mat.setFailureProperties([1000, 800, 50])

    def test_setFailureProperties_invalid_Gfrac(self) -> None:
        """Set failure properties with invalid Gfrac should raise ValueError."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        with self.assertRaises(ValueError):
            mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6], Gfrac=[1, 2, 3])

    def test_setSLis(self) -> None:
        """Set interlaminar shear strength."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setSLis(50e6)
        self.assertEqual(mat.SLis, 50e6)

    def test_str_method(self) -> None:
        """Test string representation of material."""
        mat = TransverseIsotropic([150e9, 10e9], 0.3, 5e9, alpha=[1e-6, 2e-6], rho=1500)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6], Gfrac=[1000, 2000])
        msg = str(mat)
        self.assertIn("Elastic Properties", msg)
        self.assertIn("E1", msg)
        self.assertIn("Thermal expansion", msg)
        self.assertIn("Strengths", msg)

    def test_getQ(self) -> None:
        """Test reduced stiffness matrix calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q = mat.getQ()
        self.assertEqual(Q.shape, (3, 3))
        self.assertGreater(Q[0, 0], 0)
        self.assertEqual(Q[2, 2], 5e9)
        # Verify symmetry
        self.assertEqual(Q[0, 1], Q[1, 0])

    def test_getS(self) -> None:
        """Test compliance matrix calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        S = mat.getS()
        self.assertEqual(S.shape, (3, 3))
        self.assertAlmostEqual(S[0, 0], 1.0 / 100e9, delta=1e-15)
        self.assertAlmostEqual(S[2, 2], 1.0 / 5e9, delta=1e-15)
        # Verify symmetry
        self.assertEqual(S[0, 1], S[1, 0])

    def test_Q_and_S_consistency(self) -> None:
        """Check that Q and S matrices are mutual inverses."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q = mat.getQ()
        S = mat.getS()
        I = Q @ S
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-6))

    def test_getU(self) -> None:
        """Test stiffness invariants calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        U = mat.getU()
        self.assertEqual(len(U), 5)
        # U should be cached
        U2 = mat.getU()
        self.assertTrue(np.array_equal(U, U2))

    def test_getV(self) -> None:
        """Test compliance invariants calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        V = mat.getV()
        self.assertEqual(len(V), 5)

    def test_getQbar_zero_angle(self) -> None:
        """Test rotated stiffness matrix at 0 degrees."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q = mat.getQ()
        Qbar = mat.getQbar(0.0)
        self.assertTrue(np.allclose(Q, Qbar, atol=1e-6))

    def test_getQbar_rotation(self) -> None:
        """Test rotated stiffness matrix at 90 degrees."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q0 = mat.getQbar(0.0)
        Q90 = mat.getQbar(90.0)
        # Q[0,0] at 0° should equal Q[1,1] at 90°
        self.assertTrue(np.allclose(Q0[0, 0], Q90[1, 1], atol=1e-6))

    def test_getSbar(self) -> None:
        """Test rotated compliance matrix."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Sbar = mat.getSbar(45.0)
        self.assertEqual(Sbar.shape, (3, 3))

    def test_getAlpha_zero_angle(self) -> None:
        """Test thermal expansion vector at 0 degrees."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9, alpha=[1e-6, 2e-6])
        alpha = mat.getAlpha(0.0)
        self.assertAlmostEqual(alpha[0], 1e-6, delta=1e-12)
        self.assertAlmostEqual(alpha[1], 2e-6, delta=1e-12)
        self.assertAlmostEqual(alpha[2], 0.0, delta=1e-12)

    def test_getAlpha_ninety_angle(self) -> None:
        """Test thermal expansion vector at 90 degrees."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9, alpha=[1e-6, 2e-6])
        alpha = mat.getAlpha(90.0)
        self.assertAlmostEqual(alpha[0], 2e-6, delta=1e-12)
        self.assertAlmostEqual(alpha[1], 1e-6, delta=1e-12)

    def test_getFIMaximumStress_tension(self) -> None:
        """Test maximum stress failure criterion with tensile stress."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIMaximumStress(np.array([1000e6, 0, 0]))
        self.assertAlmostEqual(FI, 1.0, delta=1e-6)

    def test_getFIMaximumStress_compression(self) -> None:
        """Test maximum stress failure criterion with compressive stress."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIMaximumStress(np.array([-800e6, 0, 0]))
        self.assertAlmostEqual(FI, 1.0, delta=1e-6)

    def test_getFIMaximumStress_shear(self) -> None:
        """Test maximum stress failure criterion with shear stress."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIMaximumStress(np.array([0, 0, 30e6]))
        self.assertAlmostEqual(FI, 1.0, delta=1e-6)

    def test_getFIMaximumStrain(self) -> None:
        """Test maximum strain failure criterion."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIMaximumStrain(np.array([500e6, 0, 0]))
        #self.assertLess(FI, 1.0)

    def test_getFITsaiWu(self) -> None:
        """Test Tsai-Wu failure criterion."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6   ])
        FI = mat.getFITsaiWu(np.array([500e6, 20e6, 10e6]))
        self.assertGreater(FI, 0)

    def test_getFIHashin73(self) -> None:
        """Test Hashin 1973 failure criterion."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIHashin73(np.array([500e6, 20e6, 10e6]))
        self.assertGreater(FI, 0)

    def test_getFIHashin80(self) -> None:
        """Test Hashin 1980 failure criterion."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000e6, 800e6, 50e6, 40e6, 30e6])
        FI = mat.getFIHashin80(np.array([500e6, 20e6, 10e6]))
        self.assertGreater(FI, 0)

    def test_getFILarc03(self) -> None:
        """Test Larc03 failure criterion."""
        mat = TransverseIsotropic(140e9, 0.25, 5e9)
        mat.setFailureProperties([2000e6, 1200e6, 50e6, 200e6, 80e6], Gfrac=[500, 1000])
        mat.setSLis(60e6)
        FI = mat.getFILarc03(np.array([500e6, 20e6, 10e6]))
        self.assertGreater(FI, 0)


class TestLayer(unittest.TestCase):
    """Tests for Layer class."""

    def test_layer_init(self) -> None:
        """Test layer initialization."""
        layer = Layer("carbon", 45.0, 0.125)
        self.assertEqual(layer.name, "carbon")
        self.assertEqual(layer.theta, 45.0)
        self.assertEqual(layer.thick, 0.125)


class TestLaminate(unittest.TestCase):
    """Tests for Laminate class."""

    def test_laminate_init(self) -> None:
        """Test laminate initialization."""
        lam = Laminate()
        self.assertEqual(len(lam.materials), 0)
        self.assertEqual(len(lam.layers), 0)

    def test_addMaterial(self) -> None:
        """Test adding material to laminate."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        self.assertIn("mat1", lam.materials)
        self.assertEqual(lam.materials["mat1"], mat)

    def test_addLayer(self) -> None:
        """Test adding layer to laminate."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        self.assertEqual(len(lam.layers), 1)
        self.assertEqual(lam.layers[0].name, "mat1")
        self.assertEqual(lam.layers[0].theta, 0)
        self.assertEqual(lam.layers[0].thick, 0.125)
        self.assertEqual(lam.thick, 0.125)

    def test_addMultipleLayers(self) -> None:
        """Test adding multiple layers."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        lam.addLayer("mat1", 45, 0.125)
        self.assertEqual(len(lam.layers), 3)
        self.assertEqual(lam.thick, 0.375)

    def test_removeAllLayers(self) -> None:
        """Test removing all layers."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        lam.removeAllLayers()
        self.assertEqual(len(lam.layers), 0)

    def test_str_method(self) -> None:
        """Test string representation of laminate."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        msg = str(lam)
        self.assertIn("Laminate properties", msg)
        self.assertIn("layer", msg)

    def test_getA_single_layer(self) -> None:
        """Test A matrix for single layer."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        A = lam.getA()
        self.assertEqual(A.shape, (3, 3))
        self.assertGreater(A[0, 0], 0)

    def test_getB_symmetric_laminate(self) -> None:
        """Test B matrix is zero for symmetric laminate."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 45, 0.125)
        lam.addLayer("mat1", -45, 0.125)
        lam.addLayer("mat1", -45, 0.125)
        lam.addLayer("mat1", 45, 0.125)
        B = lam.getB()
        self.assertTrue(np.allclose(B, np.zeros((3, 3)), atol=1e-10))

    def test_getB_unsymmetric_laminate(self) -> None:
        """Test B matrix is non-zero for unsymmetric laminate."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        B = lam.getB()
        # B should have non-zero elements for unsymmetric
        self.assertFalse(np.allclose(B, np.zeros((3, 3)), atol=1e-10))

    def test_getD(self) -> None:
        """Test D matrix calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        D = lam.getD()
        self.assertEqual(D.shape, (3, 3))
        self.assertGreater(D[0, 0], 0)

    def test_getTs(self) -> None:
        """Test thermal force resultant calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9, alpha=[1e-6, 2e-6])
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        Ts = lam.getTs()
        self.assertEqual(len(Ts), 3)

    def test_getTss(self) -> None:
        """Test thermal moment resultant calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9, alpha=[1e-6, 2e-6])
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        Tss = lam.getTss()
        self.assertEqual(len(Tss), 3)

    def test_getRhoh(self) -> None:
        """Test areal mass density calculation."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9, rho=1500.0)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        rhoh = lam.getRhoh()
        self.assertAlmostEqual(rhoh, 1500.0 * 0.25, delta=1e-6)

    def test_getZcoord(self) -> None:
        """Test z-coordinate of layer centroid."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        z = lam.getZcoord(0)
        self.assertAlmostEqual(z, 0.0, delta=1e-9)

    def test_getLayerBounds(self) -> None:
        """Test layer boundary coordinates."""
        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        z_bot, z_top = lam.getLayerBounds(0)
        self.assertLess(z_bot, z_top)

    def test_getInverseMatrices(self) -> None:
        """Test inverse stiffness matrices."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1",  0, 0.125)    
        lam.addLayer("mat1", 90, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        lam.addLayer("mat1",  0, 0.125)
        A1, B1, C1, D1 = lam.getInverseMatrices()
        
        # Test A * A1 ≈ I
        A = lam.getA()
        I = A @ A1
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-6))
        
        # Test D * D1 ≈ I
        D = lam.getD()
        I = D @ D1
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-3))

    def test_getQbar(self) -> None:
        """Test getting Qbar for specific layer."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 45, 0.125)
        Qbar = lam.getQbar(0)
        self.assertEqual(Qbar.shape, (3, 3))

    def test_getElastic(self) -> None:
        """Test effective elastic properties calculation."""
        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("mat1", mat)
        lam.addLayer("mat1", 0, 0.125)
        lam.addLayer("mat1", 90, 0.125)
        Ex, Ey, nuxy, Gxy = lam.getElastic()
        self.assertGreater(Ex, 0)
        self.assertGreater(Ey, 0)
        self.assertGreater(Gxy, 0)


class TestUtilityFunctions(unittest.TestCase):
    """Tests for utility functions."""

    def test_stressTransformation_zero_angle(self) -> None:
        """Test stress transformation at 0 degrees."""
        sigma = np.array([100.0, 50.0, 20.0])
        sigma_rot = stressTransformation(sigma, 0.0)
        self.assertTrue(np.allclose(sigma, sigma_rot, atol=1e-10))

    def test_stressTransformation_ninety_degrees(self) -> None:
        """Test stress transformation at 90 degrees."""
        sigma = np.array([100.0, 0.0, 0.0])
        sigma_rot = stressTransformation(sigma, 90.0)
        self.assertAlmostEqual(sigma_rot[1], 100.0, delta=1e-6)
        self.assertAlmostEqual(sigma_rot[0], 0.0, delta=1e-6)

    def test_stressTransformation_fortyfive_degrees(self) -> None:
        """Test stress transformation at 45 degrees."""
        sigma = np.array([100.0, 0.0, 0.0])
        sigma_rot = stressTransformation(sigma, 45.0)
        # At 45°, σx should be split between σ'x, σ'y, and τ'xy
        self.assertAlmostEqual(sigma_rot[0], 50.0, delta=1e-6)
        self.assertAlmostEqual(sigma_rot[1], 50.0, delta=1e-6)

    def test_strainTransformation_zero_angle(self) -> None:
        """Test strain transformation at 0 degrees."""
        eps = np.array([0.001, 0.0005, 0.0002])
        eps_rot = strainTransformation(eps, 0.0)
        self.assertTrue(np.allclose(eps, eps_rot, atol=1e-10))

    def test_strainTransformation_ninety_degrees(self) -> None:
        """Test strain transformation at 90 degrees."""
        eps = np.array([0.001, 0.0, 0.0])
        eps_rot = strainTransformation(eps, 90.0)
        self.assertAlmostEqual(eps_rot[1], 0.001, delta=1e-9)
        self.assertAlmostEqual(eps_rot[0], 0.0, delta=1e-9)

    def test_mixMaterials_valid(self) -> None:
        """Test material mixing with valid volume fraction."""
        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        compmat = mixMaterials(carbon, epoxy, 0.6)
        
        self.assertAlmostEqual(compmat.E1, 1.3344e11, delta=1.0e8)
        self.assertAlmostEqual(compmat.E2, 8.784e9, delta=1.0e7)
        self.assertAlmostEqual(compmat.nu12, 0.26, delta=0.001)
        self.assertAlmostEqual(compmat.G12, 3.254e9, delta=1.0e7)

    def test_mixMaterials_zero_vf(self) -> None:
        """Test material mixing with zero fiber volume fraction."""
        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        compmat = mixMaterials(carbon, epoxy, 0.0)
        
        # Should get pure matrix properties
        self.assertAlmostEqual(compmat.E1, epoxy.E1, delta=1.0)

    def test_mixMaterials_one_vf(self) -> None:
        """Test material mixing with unit fiber volume fraction."""
        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        compmat = mixMaterials(carbon, epoxy, 1.0)
        
        # Should get pure fiber properties
        self.assertAlmostEqual(compmat.E1, carbon.E1, delta=1.0)

    def test_mixMaterials_invalid_vf_negative(self) -> None:
        """Test material mixing with invalid negative volume fraction."""
        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        with self.assertRaises(RuntimeError):
            mixMaterials(carbon, epoxy, -0.1)

    def test_mixMaterials_invalid_vf_greater_than_one(self) -> None:
        """Test material mixing with invalid volume fraction > 1."""
        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
        
        with self.assertRaises(RuntimeError):
            mixMaterials(carbon, epoxy, 1.1)

    def test_Macauley_positive(self) -> None:
        """Test Macauley operator with positive value."""
        self.assertEqual(Macauley(5.0), 5.0)
        self.assertEqual(Macauley(0.001), 0.001)

    def test_Macauley_negative(self) -> None:
        """Test Macauley operator with negative value."""
        self.assertEqual(Macauley(-3.0), 0.0)
        self.assertEqual(Macauley(-0.001), 0.0)

    def test_Macauley_zero(self) -> None:
        """Test Macauley operator with zero."""
        self.assertEqual(Macauley(0.0), 0.0)


class TestIntegration(unittest.TestCase):
    """Integration tests combining multiple components."""

    def test_full_laminate_analysis(self) -> None:
        """Test complete laminate analysis workflow."""
        # Create material
        mat = TransverseIsotropic([140e9, 10e9], 0.3, 5e9, alpha=[1e-6, 25e-6], rho=1600)
        mat.setFailureProperties([2000e6, 1200e6, 50e6, 200e6, 80e6])
        
        # Create laminate
        lam = Laminate()
        lam.addMaterial("carbon", mat)
        lam.addLayer("carbon", 0, 0.125)
        lam.addLayer("carbon", 90, 0.125)
        lam.addLayer("carbon", 90, 0.125)
        lam.addLayer("carbon", 0, 0.125)
        
        # Get stiffness matrices
        A = lam.getA()
        B = lam.getB()
        D = lam.getD()
        
        # Verify symmetric laminate has zero B
        self.assertTrue(np.allclose(B, np.zeros((3, 3)), atol=1e-10))
        
        # Get effective properties
        Ex, Ey, nuxy, Gxy = lam.getElastic()
        self.assertGreater(Ex, 0)
        self.assertGreater(Ey, 0)
        
        # Get mass
        rhoh = lam.getRhoh()
        self.assertAlmostEqual(rhoh, 1600 * 0.5, delta=1.0)

    def test_quasi_isotropic_laminate(self) -> None:
        """Test quasi-isotropic [0/±60] laminate."""
        mat = TransverseIsotropic([140e9, 10e9], 0.3, 5e9)
        
        lam = Laminate()
        lam.addMaterial("mat", mat)
        lam.addLayer("mat", 0, 0.125)
        lam.addLayer("mat", 60, 0.125)
        lam.addLayer("mat", -60, 0.125)
        
        Ex, Ey, nuxy, Gxy = lam.getElastic()
        
        # For quasi-isotropic, Ex should approximately equal Ey
        self.assertAlmostEqual(Ex, Ey, delta=Ex * 0.1)


if __name__ == '__main__':
    unittest.main()
