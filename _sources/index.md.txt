# PyComposites

```{image} img/pycomposites_logo.png
:alt: PyComposites logo
:class: pycomposites-logo
:width: 360px
```

The `pycomposites` package provides a compact framework for performing
**Classical Laminate Theory (CLT)** calculations. It was originally developed
for the TU/e course *Composite and Lightweight Materials - Design and Analysis
(4MM00)*, but is broadly applicable to engineering and research involving
laminated composite structures.

## Core Components

- **`TransverseIsotropic`**

  Represents a transversely isotropic unidirectional ply. It stores elastic
  constants, density, thermal expansion coefficients, and failure properties.
  The class provides methods to compute stiffness and compliance matrices,
  invariant forms, rotated stiffness/compliance, and thermal expansion vectors.
  Several classical failure criteria are implemented, including Maximum Stress,
  Maximum Strain, Tsai-Wu, Hashin, and Larc03.

- **`Layer`**

  A lightweight container describing a single ply, defined by its material,
  orientation, and thickness.

- **`Laminate`**

  Defines a stacking sequence of multiple layers. It computes geometric layer
  boundaries, laminate thickness, stiffness matrices, thermal resultants,
  effective elastic constants, and inverse matrices for solving laminate
  constitutive equations.

- **Utility functions**

  Helper routines for stress and strain transformations between local and
  global coordinates, homogenization of fiber-matrix systems via
  rule-of-mixtures (`mixMaterials`), and the Macauley operator used in failure
  models.

## Purpose and Scope

The `pycomposites` package is designed with clarity, educational value, and
extendibility in mind. It enables students and engineers to analyze composite
laminates from the ply level to the full laminate response. Its simple
structure makes it useful for classroom demonstrations, student assignments,
and rapid prototyping of research ideas.

```{toctree}
:maxdepth: 2
:caption: Contents

install
examples
```

```{toctree}
:maxdepth: 2
:caption: API Reference

composite
```
