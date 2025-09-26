.. toctree::
   :maxdepth: 2

Classical Laminate Theory
=========================

.. raw:: html

   <div class="blue-block">
       The examples in are obtained from abcd
   </div>
   
   <div class="green-block">
       The examples in are obtained from.
   </div>
   
   <div class="red-block">
       The examples in are obtained from.
   </div>



Example: Homogenisation
-----------------------

**Problem Statement**

Consider a fibre-reinforced plastic that consists of uni-directional 
carbon fibres embedded in an epoxy matrix. The fibre volume fraction 
:math:`V_f = 0.6`. The properties of the transversely isotropic fibre are: 
:math:`E_{fL} = 220 \, \text{GPa}`, :math:`E_{fT} = 20 \, \text{GPa}`, 
:math:`\nu_f = 0.2`, :math:`G_{\text{f}} = 91.7 \, \text{GPa}`.

The properties of the isotropic epoxy matrix are: 
:math:`E_{\text{m}} = 3.6 \, \text{GPa}`, :math:`\nu_{\text{m}} = 0.35`, 
:math:`G_{\text{m}} = 1.33 \, \text{GPa}`.

Determine the homogenised properties of the composite.

**Solution**

.. currentmodule:: pycomposites.composite
   
Import the functions :py:func:`TransverseIsotropic` and :py:func:`mixMaterials` from the 
`composite` module:

.. code-block:: python

   from composite import TransverseIsotropic, mixMaterials

Create two materials, `carbon` and `epoxy`, with the correct properties. Note 
that `carbon` is transversely isotropic and `epoxy` is isotropic.

.. code-block:: python

   carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
   epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)

   print("The properties of carbon are:\n", carbon)
   print("The properties of epoxy are:\n", epoxy)

The properties of carbon are:

.. code-block:: text

   Elastic Properties:
   -----------------------------------------------------------
   E1     :     2.200e+11 , E2     :     2.200e+10 
   nu12   :          0.20 , G12    :     9.170e+10 

   Thermal expansion coefficients:
   -----------------------------------------------------------
   alpha1 :     0.000e+00 , alpha2 :     0.000e+00 

The properties of epoxy are:

.. code-block:: text

   Elastic Properties:
   -----------------------------------------------------------
   E1     :     3.600e+09 , E2     :     3.600e+09 
   nu12   :          0.35 , G12    :     1.330e+09 

   Thermal expansion coefficients:
   -----------------------------------------------------------
   alpha1 :     0.000e+00 , alpha2 :     0.000e+00 

The properties of the composite can be calculated using the 
:py:func:`mixMaterials` 
function, which takes the two materials and the volume fraction as input.

.. code-block:: python

   udcomp = mixMaterials(carbon, epoxy, 0.6)

   print("Material properties of the composite material:\n\n", udcomp, "\n")

Material properties of the composite material:

.. code-block:: text

   Elastic Properties:
   -----------------------------------------------------------
   E1     :     1.334e+11 , E2     :     7.226e+09 
   nu12   :          0.26 , G12    :     3.254e+09 

   Thermal expansion coefficients:
   -----------------------------------------------------------
   alpha1 :     0.000e+00 , alpha2 :     0.000e+00 
   
   
Example 2: Classical Laminate Theory
------------------------------------

**Problem Statement**

Consider a fibre-reinforced plastic consisting of uni-directional carbon fibres embedded in an epoxy matrix. The fibre volume fraction is :math:`V_f = 0.6`. The properties of the transversely isotropic fibre are:

- :math:`E_f^L = 220 \, \text{GPa}`
- :math:`E_f^T = 20 \, \text{GPa}`
- :math:`\nu_f = 0.2`
- :math:`G_f = 91.7 \, \text{GPa}`

The properties of the isotropic epoxy matrix are:

- :math:`E_m = 3.6 \, \text{GPa}`
- :math:`\nu_m = 0.35`
- :math:`G_m = 1.33 \, \text{GPa}`

Determine the :math:`\mathbf{Q}` matrix of this material, as well as the rotated matrices :math:`\mathbf{Q}^{20}` and :math:`\mathbf{Q}^{-20}`. Evaluate the results.

**Solution**

Import the necessary classes and functions:

.. code-block:: python

    from composite import TransverseIsotropic, mixMaterials

For this carbon fibre composite material, the carbon fibres and the epoxy matrix are 
modeled as separate transversely isotropic materials. The T-300 carbon fibres have 
the following properties: :math:`E_1 = 220 \, \text{GPa}`, :math:`E_2 = 22 \, \text{GPa}`, 
:math:`\nu_{12} = 0.2` and :math:`G_{12} = 91.7 \, \text{GPa}`. 
The epoxy matrix is isotropic with the following properties: :math:`E = 3.6 \, \text{GPa}`,
:math:`\nu = 0.35` and :math:`G = 1.33 \, \text{GPa}`. 

Define the materials in Python using the function
:py:func:`TransverseIsotropic`:

.. code-block:: python

    carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
    epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)

The properties of the carbon fibre and epoxy matrix are:

.. code-block:: python

    print(carbon)
    print(epoxy)

Output:

::

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     2.200e+11 , E2     :     2.200e+11 
    nu12   :          0.20 , G12    :     9.170e+10 

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     3.600e+09 , E2     :     3.600e+09 
    nu12   :          0.35 , G12    :     1.330e+09 

---

**Composite Material Properties**

The composite consists of 60% fibres and 40% epoxy. Use simple homogenisation 
:py:func:`maxMaterials` to calculate the composite properties:

.. code-block:: python

    udcomp = mixMaterials(carbon, epoxy, 0.6)

Print the properties of the UD material:

.. code-block:: python

    print(udcomp)

Output:

::

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     1.334e+11 , E2     :     8.784e+09 
    nu12   :          0.26 , G12    :     3.254e+09 

---

**Stiffness Matrix**

The stiffness matrix :math:`\mathbf{Q}` can be calculated as follows:

.. code-block:: python

    Q = udcomp.getQ()
    print("The Q matrix of the composite is:\n\n", Q, "\n")

Output:

::

    [[1.34036479e+11 2.29414891e+09 0.00000000e+00]
     [2.29414891e+09 8.82364964e+09 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 3.25420247e+09]]

**Rotated Matrices**

Calculate the rotated stiffness matrix at 20°:

.. code-block:: python

    Qbar = udcomp.getQbar(20.0)
    print("The Q matrix of the composite under an angle of 20 degrees is:\n\n", Qbar, "\n")

The following output is written to the screen:

::

    [[1.06262605e+11 1.47349731e+10 3.56698011e+10]
     [1.47349731e+10 9.22203891e+09 5.04355063e+09]
     [3.56698011e+10 5.04355063e+09 1.61034402e+10]]

Calculate the rotated stiffness matrix at -20°:

.. code-block:: python

    Qbar = udcomp.getQbar(-20.0)
    print("The Q matrix of the composite under an angle of minus 20 degrees is:\n\n", Qbar, "\n")

This will give the following output:

::

    [[ 1.06262605e+11  1.47349731e+10 -3.56698011e+10]
     [ 1.47349731e+10  9.22203891e+09 -5.04355063e+09]
     [-3.56698011e+10 -5.04355063e+09  1.61034402e+10]]

**Observations**

Note that the :math:`Q_{11}`, :math:`Q_{12}`, :math:`Q_{21}`, :math:`Q_{22}`, and :math:`Q_{66}` 
terms of the stiffness matrices :math:`\mathbf{Q}^{20}` and :math:`\mathbf{Q}^{-20}` matrices 
are identical. The :math:`Q_{16}` and :math:`Q_{26}` terms switch signs between the two matrices.
   

