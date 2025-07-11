"""
  Unittests
  
  (c) Joris Remmers (2025)
  
"""
import unittest,os
import numpy as np
from pycomposites.composite import TransverseIsotropic,mixMaterials, Laminate, stressTransformation

  
class PyCompositeTesting(unittest.TestCase):
        
    def testMixMaterials(self): 
    
        carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
        epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)
        
        compmat = mixMaterials( carbon , epoxy , 0.6 )
                
        self.assertAlmostEqual( compmat.E1   , 1.3344e+11 , delta=1.0e8 )
        self.assertAlmostEqual( compmat.E2   , 8.784e+09  , delta=1.0e7 )   
        self.assertAlmostEqual( compmat.nu12 , 0.26       , delta=0.001 )
        self.assertAlmostEqual( compmat.G12  , 3.254e+09  , delta=1.0e7 ) 
        
    def testSetAlpha(self): 
    
        carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
        
        carbon.setAlpha( 1.0e-6 )
        
        self.assertAlmostEqual( carbon.alpha1 , 1.0e-6 , delta=1.0e-12 )
        self.assertAlmostEqual( carbon.alpha2 , 1.0e-6 , delta=1.0e-12 )        

        carbon.setAlpha( [ 2.0e-6 , 3.0e-6 ] )
        
        self.assertAlmostEqual( carbon.alpha1 , 2.0e-6 , delta=1.0e-12 )

        self.assertAlmostEqual( carbon.alpha2 , 3.0e-6 , delta=1.0e-12 )         

    def testAddlayer(self): 
    
        carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
        epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)
        
        compmat = mixMaterials( carbon , epoxy , 0.6 )

        lam = Laminate()
        
        lam.addMaterial('comp', compmat)
        lam.addLayer('comp', 35.0, 0.2e-3)
        lam.addLayer( 'comp' , 45 , 6e-3 )

        self.assertEqual(len(lam.layers), 2)

        first_layer = lam.layers[0]
        self.assertEqual(first_layer.name, 'comp')
        self.assertAlmostEqual(first_layer.theta, 35.0)
        self.assertAlmostEqual(first_layer.thick, 0.2e-3)

        second_layer = lam.layers[1]
        self.assertEqual(second_layer.name, 'comp')
        self.assertAlmostEqual(second_layer.theta, 45.0)
        self.assertAlmostEqual(second_layer.thick, 6e-3)

    def testStressTransformation(self):

        sigma = np.array([100, 50, 25])
        result = stressTransformation(sigma, 0)
        np.testing.assert_allclose(result, sigma, atol=1e-12)

        sigma = np.array([120, 30, 0])
        result = stressTransformation(sigma, 90)
        expected = np.array([30, 120, 0])
        print(result)
    
        np.testing.assert_allclose(result, expected, atol=1e-6)


    def testElasticCalculation(self):
        carbon = TransverseIsotropic( 220e9,0.2,91.7e9)
        epoxy  = TransverseIsotropic( 3.6e9,0.35,1.33e9)
        
        compmat = mixMaterials( carbon , epoxy , 0.6 )

        lam = Laminate()
        
        lam.addMaterial('comp', compmat)
        lam.addLayer('comp', 0, 0.2e-3)
                
        elastic = lam.getElastic()
        self.assertEqual(len(elastic), 4)
        self.assertTrue(all(e > 0 for e in elastic))

if __name__ == '__main__':
    unittest.main()
