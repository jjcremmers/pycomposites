# API Reference

This page documents the public Python API exposed by `pycomposites`.

## Classes

### TransverseIsotropic

Material model for a transversely isotropic or orthotropic ply in plane stress.

```python
TransverseIsotropic(E, nu12, G12=0.0, alpha=0.0, rho=0.0)
```

#### TransverseIsotropic methods

```{eval-rst}
.. automethod:: pycomposites.composite.TransverseIsotropic.setAlpha

.. automethod:: pycomposites.composite.TransverseIsotropic.setFailureProperties

.. automethod:: pycomposites.composite.TransverseIsotropic.setSLis

.. automethod:: pycomposites.composite.TransverseIsotropic.getQ

.. automethod:: pycomposites.composite.TransverseIsotropic.getS

.. automethod:: pycomposites.composite.TransverseIsotropic.getQbar

.. automethod:: pycomposites.composite.TransverseIsotropic.getSbar

.. automethod:: pycomposites.composite.TransverseIsotropic.getAlpha

.. automethod:: pycomposites.composite.TransverseIsotropic.getFIMaximumStress

.. automethod:: pycomposites.composite.TransverseIsotropic.getFIMaximumStrain

.. automethod:: pycomposites.composite.TransverseIsotropic.getFITsaiWu

.. automethod:: pycomposites.composite.TransverseIsotropic.getFIHashin73

.. automethod:: pycomposites.composite.TransverseIsotropic.getFIHashin80

.. automethod:: pycomposites.composite.TransverseIsotropic.getFILarc03
```

### Layer

Container describing one laminate ply.

```python
Layer(name, theta, thick)
```

### Laminate

Stacking sequence of multiple plies with Classical Laminate Theory helpers.

```python
Laminate()
```

#### Laminate methods

```{eval-rst}
.. automethod:: pycomposites.composite.Laminate.addMaterial

.. automethod:: pycomposites.composite.Laminate.addLayer

.. automethod:: pycomposites.composite.Laminate.removeAllLayers

.. automethod:: pycomposites.composite.Laminate.getA

.. automethod:: pycomposites.composite.Laminate.getB

.. automethod:: pycomposites.composite.Laminate.getD

.. automethod:: pycomposites.composite.Laminate.getTs

.. automethod:: pycomposites.composite.Laminate.getTss

.. automethod:: pycomposites.composite.Laminate.getRhoh

.. automethod:: pycomposites.composite.Laminate.getZcoord

.. automethod:: pycomposites.composite.Laminate.getLayerBounds

.. automethod:: pycomposites.composite.Laminate.getInverseMatrices

.. automethod:: pycomposites.composite.Laminate.getQbar

.. automethod:: pycomposites.composite.Laminate.getElastic
```

## Utility Functions

```{eval-rst}
.. autofunction:: pycomposites.composite.stressTransformation

.. autofunction:: pycomposites.composite.strainTransformation

.. autofunction:: pycomposites.composite.mixMaterials

.. autofunction:: pycomposites.composite.Macauley
```
