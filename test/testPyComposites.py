"""
Unit tests for the ``composite`` package.

These tests exercise the material models and laminate helpers. They are
focused on numerical checks and basic API behaviour; no production code is
modified by these tests.

(c) Joris Remmers (2025)
"""

import unittest, os
import numpy as np
from typing import Any

from composite import (
    TransverseIsotropic, Laminate, mixMaterials,
    stressTransformation, strainTransformation, Macauley
)


class PyCompositeTesting(unittest.TestCase):
    """Collection of unit tests for core composite functionality.

    Each test verifies a small, self-contained behaviour (initialisation,
    matrix generation, transformations, failure checks, etc.). Tests only
    read and assert on outputs and do not modify library behaviour.
    """

    def testMixMaterials(self) -> None:
        """Mix two materials and check homogenised properties."""

        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)
        epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)

        compmat = mixMaterials(carbon, epoxy, 0.6)

        self.assertAlmostEqual(compmat.E1, 1.3344e+11, delta=1.0e8)
        self.assertAlmostEqual(compmat.E2, 8.784e+09, delta=1.0e7)
        self.assertAlmostEqual(compmat.nu12, 0.26, delta=0.001)
        self.assertAlmostEqual(compmat.G12, 3.254e+09, delta=1.0e7)

    def testSetAlpha(self) -> None:
        """Ensure ``setAlpha`` accepts scalar and list inputs correctly."""

        carbon = TransverseIsotropic(220e9, 0.2, 91.7e9)

        carbon.setAlpha(1.0e-6)

        self.assertAlmostEqual(carbon.alpha1, 1.0e-6, delta=1.0e-12)
        self.assertAlmostEqual(carbon.alpha2, 1.0e-6, delta=1.0e-12)

        carbon.setAlpha([2.0e-6, 3.0e-6])

        self.assertAlmostEqual(carbon.alpha1, 2.0e-6, delta=1.0e-12)
        self.assertAlmostEqual(carbon.alpha2, 3.0e-6, delta=1.0e-12)

    def test_material_init_scalar(self) -> None:
        """Initialise material with scalar transverse property values."""

        mat = TransverseIsotropic(100e9, 0.25)
        self.assertEqual(mat.E1, 100e9)
        self.assertEqual(mat.E2, 100e9)
        self.assertAlmostEqual(mat.nu21, 0.25, delta=1e-6)

    def test_material_init_list(self) -> None:
        """Initialise material with E provided as a two-element list."""

        mat = TransverseIsotropic([150e9, 10e9], 0.3)
        self.assertEqual(mat.E1, 150e9)
        self.assertEqual(mat.E2, 10e9)

    def test_invalid_E_list(self) -> None:
        """Providing an invalid E list should raise ValueError."""

        with self.assertRaises(ValueError):
            TransverseIsotropic([1, 2, 3], 0.3)

    def test_setAlpha(self) -> None:
        """Set thermal expansion coefficients using a list input."""

        mat = TransverseIsotropic(70e9, 0.2)
        mat.setAlpha([1e-6, 2e-6])
        self.assertEqual(mat.alpha1, 1e-6)
        self.assertEqual(mat.alpha2, 2e-6)

    def test_getQ_and_S_consistency(self) -> None:
        """Check that Q and S matrices are mutual inverses (plane stress)."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q = mat.getQ()
        S = mat.getS()
        # Q * S ≈ identity (in plane stress)
        I = Q @ S
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-6))

    def test_getQbar_rotation(self) -> None:
        """Verify reduced stiffness matrix rotates consistently for 0/90°."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        Q0 = mat.getQbar(0.0)
        Q90 = mat.getQbar(90.0)
        self.assertTrue(np.allclose(Q0[0, 0], Q90[1, 1], atol=1e-6))

    def test_getAlpha_rotation(self) -> None:
        """Check rotated thermal expansion vector components."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9, alpha=[1e-6, 2e-6])
        a0 = mat.getAlpha(0.0)
        a90 = mat.getAlpha(90.0)
        self.assertAlmostEqual(a0[0], 1e-6)
        self.assertAlmostEqual(a90[1], 1e-6)

    def test_laminate_ABD(self) -> None:
        """Assemble A, B, D matrices and verify B is zero for symmetric layup."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("m", mat)
        lam.addLayer("m", 0, 0.1)
        lam.addLayer("m", 0, 0.1)
        A = lam.getA()
        B = lam.getB()
        D = lam.getD()
        self.assertTrue(np.allclose(B, np.zeros((3, 3)), atol=1e-12))

    def test_laminate_inverseABD(self) -> None:
        """Validate the inverse matrices returned by laminate routines."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("m", mat)
        lam.addLayer("m", 0, 0.1)
        lam.addLayer("m", 0, 0.1)
        A = lam.getA()
        D = lam.getD()

        (A1, B1, C1, D1) = lam.getInverseMatrices()

        I = A1 @ A
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-6))

        I = D1 @ D
        self.assertTrue(np.allclose(I, np.eye(3), atol=1e-6))

    def test_elastic_effective(self) -> None:
        """Compute effective elastic moduli from a small laminate."""

        mat = TransverseIsotropic([100e9, 10e9], 0.25, 5e9)
        lam = Laminate()
        lam.addMaterial("m", mat)
        lam.addLayer("m", 0, 0.1)
        lam.addLayer("m", 90, 0.1)
        Ex, Ey, nuxy, Gxy = lam.getElastic()
        self.assertGreater(Ex, 1e9)
        self.assertGreater(Ey, 1e9)

    def test_failure_max_stress(self) -> None:
        """Check maximum-stress failure index equals 1 at the material strength."""

        mat = TransverseIsotropic(100e9, 0.25, 5e9)
        mat.setFailureProperties([1000, 800, 50, 40, 30])
        FI = mat.getFIMaximumStress(np.array([1000, 0, 0]))
        self.assertAlmostEqual(FI, 1.0)

    def test_stressTransformation(self) -> None:
        """Stress transformation should rotate a stress vector by given angle."""

        sigma = np.array([1.0, 0.0, 0.0])
        sig90 = stressTransformation(sigma, 90)
        self.assertAlmostEqual(sig90[1], 1.0, places=6)

    def test_Macauley(self) -> None:
        """Macauley returns the positive part of a scalar (max(x, 0))."""

        self.assertEqual(Macauley(5.0), 5.0)
        self.assertEqual(Macauley(-3.0), 0.0)


if __name__ == '__main__':
    unittest.main()
