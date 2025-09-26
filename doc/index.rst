.. PyComposites documentation master file, created by
   sphinx-quickstart on Sun Dec 29 09:56:23 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PyComposites
============

The ``composite`` module provides a compact framework for performing 
**Classical Laminate Theory (CLT)** calculations. It was originally 
developed for the TU/e course *Composite and Lightweight Materials – 
Design and Analysis (4MM00)*, but is broadly applicable to engineering 
and research involving laminated composite structures.

Core Components
---------------

- **``TransverseIsotropic``**
  
  Represents a transversely isotropic unidirectional ply. It stores 
  elastic constants (E, ν, G), density, thermal expansion coefficients, 
  and failure properties. The class provides methods to compute stiffness 
  (Q) and compliance (S) matrices, invariant forms (U, V), rotated 
  stiffness/compliance (Q̅, S̅), and thermal expansion vectors. Several 
  classical failure criteria are implemented, including Maximum Stress, 
  Maximum Strain, Tsai–Wu, Hashin, and Larc03.

- **``Layer``**
  
  A lightweight container describing a single ply, defined by its 
  material, orientation, and thickness.

- **``Laminate``**
  
  Defines a stacking sequence of multiple layers. It computes geometric 
  layer boundaries, laminate thickness, stiffness matrices (A, B, D), 
  thermal resultants, effective elastic constants, and inverse matrices 
  for solving laminate constitutive equations.

- **Utility functions**
  
  The module provides helper routines for stress and strain transformations 
  between local and global coordinates, homogenization of fiber–matrix 
  systems via rule-of-mixtures (``mixMaterials``), and the Macauley operator 
  ⟨x⟩ used in failure models.

Purpose and Scope
-----------------

The ``composite`` module is designed with **clarity, educational value, 
and extendibility** in mind. It enables students and engineers to analyze 
composite laminates from the ply level to the full laminate response. 
Its simple structure makes it ideal for classroom demonstrations, 
student assignments, and rapid prototyping of research ideas. 
While lightweight by design, it captures the essential mechanics of 
laminated composites and serves as a solid foundation for more advanced 
analytical or numerical models.


.. toctree::
   :maxdepth: 2
   
   install
   examples
   
   :caption: API Reference
   
   composite



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
