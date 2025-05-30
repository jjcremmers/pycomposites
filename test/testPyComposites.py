"""A set of classes and functions to read and transform hdf5 files in the NDF 
    file format.

  (c) Joris Remmers (2021-2023)
  
"""
import unittest,os
from pycomposites.composite import TransverseIsotropic,mixMaterials

  
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
                
if __name__ == '__main__':
    unittest.main()
