# Classical Laminate Theory

## Example: Homogenisation

**Problem Statement**

Consider a fibre-reinforced plastic that consists of uni-directional
carbon fibres embedded in an epoxy matrix. The fibre volume fraction
$V_f = 0.6$. The properties of the transversely isotropic fibre are:
$E_{fL} = 220 \, \text{GPa}$, $E_{fT} = 20 \, \text{GPa}$,
$\nu_f = 0.2$, $G_{\text{f}} = 91.7 \, \text{GPa}$.

The properties of the isotropic epoxy matrix are:
$E_{\text{m}} = 3.6 \, \text{GPa}$, $\nu_{\text{m}} = 0.35$,
$G_{\text{m}} = 1.33 \, \text{GPa}$.

Determine the homogenised properties of the composite.

**Solution**

Import the functions `TransverseIsotropic` and `mixMaterials` from the
`pycomposites` package:

``` python
from pycomposites import TransverseIsotropic, mixMaterials
```

Create two materials, `carbon` and
`epoxy`, with the correct properties. Note
that `carbon` is transversely isotropic and
`epoxy` is isotropic.

``` python
carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)

print("The properties of carbon are:\n", carbon)
print("The properties of epoxy are:\n", epoxy)
```

The properties of carbon are:

``` text
Elastic Properties:
-----------------------------------------------------------
E1     :     2.200e+11 , E2     :     2.200e+10 
nu12   :          0.20 , G12    :     9.170e+10 

Thermal expansion coefficients:
-----------------------------------------------------------
alpha1 :     0.000e+00 , alpha2 :     0.000e+00 
```

The properties of epoxy are:

``` text
Elastic Properties:
-----------------------------------------------------------
E1     :     3.600e+09 , E2     :     3.600e+09 
nu12   :          0.35 , G12    :     1.330e+09 

Thermal expansion coefficients:
-----------------------------------------------------------
alpha1 :     0.000e+00 , alpha2 :     0.000e+00 
```

The properties of the composite can be calculated using the
`mixMaterials` function, which takes the two materials and the volume
fraction as input.

``` python
udcomp = mixMaterials(carbon, epoxy, 0.6)

print("Material properties of the composite material:\n\n", udcomp, "\n")
```

Material properties of the composite material:

``` text
Elastic Properties:
-----------------------------------------------------------
E1     :     1.334e+11 , E2     :     7.226e+09 
nu12   :          0.26 , G12    :     3.254e+09 

Thermal expansion coefficients:
-----------------------------------------------------------
alpha1 :     0.000e+00 , alpha2 :     0.000e+00 
```

## Example 2: Classical Laminate Theory

**Problem Statement**

Consider a fibre-reinforced plastic consisting of uni-directional carbon
fibres embedded in an epoxy matrix. The fibre volume fraction is
$V_f = 0.6$. The properties of the transversely isotropic fibre are:

- $E_f^L = 220 \, \text{GPa}$
- $E_f^T = 20 \, \text{GPa}$
- $\nu_f = 0.2$
- $G_f = 91.7 \, \text{GPa}$

The properties of the isotropic epoxy matrix are:

- $E_m = 3.6 \, \text{GPa}$
- $\nu_m = 0.35$
- $G_m = 1.33 \, \text{GPa}$

Determine the $\mathbf{Q}$ matrix of this material, as well as the
rotated matrices $\mathbf{Q}^{20}$ and $\mathbf{Q}^{-20}$. Evaluate the
results.

**Solution**

Import the necessary classes and functions:

``` python
from pycomposites import TransverseIsotropic, mixMaterials
```

For this carbon fibre composite material, the carbon fibres and the
epoxy matrix are modeled as separate transversely isotropic materials.
The T-300 carbon fibres have the following properties:
$E_1 = 220 \, \text{GPa}$, $E_2 = 22 \, \text{GPa}$, $\nu_{12} = 0.2$
and $G_{12} = 91.7 \, \text{GPa}$. The epoxy matrix is isotropic with
the following properties: $E = 3.6 \, \text{GPa}$, $\nu = 0.35$ and
$G = 1.33 \, \text{GPa}$.

Define the materials in Python using the function `TransverseIsotropic`:

``` python
carbon = TransverseIsotropic([220e9, 22e9], 0.2, 91.7e9)
epoxy = TransverseIsotropic(3.6e9, 0.35, 1.33e9)
```

The properties of the carbon fibre and epoxy matrix are:

``` python
print(carbon)
print(epoxy)
```

Output:

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

The composite consists of 60% fibres and 40% epoxy. Use simple
homogenisation `maxMaterials` to calculate the composite properties:

``` python
udcomp = mixMaterials(carbon, epoxy, 0.6)
```

Print the properties of the UD material:

``` python
print(udcomp)
```

Output:

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     1.334e+11 , E2     :     8.784e+09 
    nu12   :          0.26 , G12    :     3.254e+09 

---

**Stiffness Matrix**

The stiffness matrix $\mathbf{Q}$ can be calculated as follows:

``` python
Q = udcomp.getQ()
print("The Q matrix of the composite is:\n\n", Q, "\n")
```

Output:

    [[1.34036479e+11 2.29414891e+09 0.00000000e+00]
     [2.29414891e+09 8.82364964e+09 0.00000000e+00]
     [0.00000000e+00 0.00000000e+00 3.25420247e+09]]

**Rotated Matrices**

Calculate the rotated stiffness matrix at 20°:

``` python
Qbar = udcomp.getQbar(20.0)
print("The Q matrix of the composite under an angle of 20 degrees is:\n\n", Qbar, "\n")
```

The following output is written to the screen:

    [[1.06262605e+11 1.47349731e+10 3.56698011e+10]
     [1.47349731e+10 9.22203891e+09 5.04355063e+09]
     [3.56698011e+10 5.04355063e+09 1.61034402e+10]]

Calculate the rotated stiffness matrix at -20°:

``` python
Qbar = udcomp.getQbar(-20.0)
print("The Q matrix of the composite under an angle of minus 20 degrees is:\n\n", Qbar, "\n")
```

This will give the following output:

    [[ 1.06262605e+11  1.47349731e+10 -3.56698011e+10]
     [ 1.47349731e+10  9.22203891e+09 -5.04355063e+09]
     [-3.56698011e+10 -5.04355063e+09  1.61034402e+10]]

**Observations**

Note that the $Q_{11}$, $Q_{12}$, $Q_{21}$, $Q_{22}$, and $Q_{66}$ terms
of the stiffness matrices $\mathbf{Q}^{20}$ and $\mathbf{Q}^{-20}$
matrices are identical. The $Q_{16}$ and $Q_{26}$ terms switch signs
between the two matrices.
