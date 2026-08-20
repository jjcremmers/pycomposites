# Failure

**Problem Statement**

Calculate the strength of a single layer of E-glass/epoxy composite
loaded by a uniaxial force under an angle $\theta$. The elastic
properties of the material are:

- $E_1 = 39.0 \, \text{GPa}$
- $E_2 = 8.6 \, \text{GPa}$
- $\nu_{12} = 0.28$
- $G_{12} = 3.254 \, \text{GPa}$

The strengths of the material are:

- $X_T = 1080 \, \text{MPa}$
- $X_C = 620 \, \text{MPa}$
- $Y_T = 39 \, \text{MPa}$
- $Y_C = 128 \, \text{MPa}$
- $S = 89 \, \text{MPa}$

Determine the strength for all angles $0 < \theta < 90$ using the
following failure criteria:

1.  Maximum Stress
2.  Maximum Strain
3.  Tsai-Wu failure criterion

For the Tsai-Wu criterion, you may assume:

$$f_{12} = -12\sqrt{f_{11}f_{22}}$$

Plot the results for each failure criterion in a graph.

**Solution**

Create the material object and assign the fracture properties. Print the
properties to verify.

``` python
from pycomposites import TransverseIsotropic

glassepoxy = TransverseIsotropic([39.0e9, 8.6e9], 0.28, 3.254e9)
glassepoxy.setFailureProperties([1080e6, 620e6, 39e6, 128e6, 89e6])

print(glassepoxy)
```

Output:

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     3.900e+10 , E2     :     8.600e+09 
    nu12   :          0.28 , G12    :     3.254e+09 

    Strengths and failure model parameters:
    -----------------------------------------------------------
    Xt     :     1.080e+09 , Xc     :     6.200e+08 
    Yt     :     3.900e+07 , Yc     :     1.280e+08 
    S      :     8.900e+07

To calculate the maximum allowable stress, the process involves:

1.  Creating a unit stress vector for a stress magnitude of 1 Pa under
    an angle $\theta$.
2.  Multiplying the unit stress by a scaling factor $sf$, which is
    iteratively adjusted until the failure index exceeds 1.0.
3.  Storing the corresponding $sf$ as the maximum allowable stress for
    that angle.

The calculations are performed for 90 intervals between angles $0^\circ$
and $90^\circ$.

Initialize the required lists and parameters:

``` python
from numpy import array, dot, linspace, zeros, pi, cos, sin

theta = list(range(0, 91))
maxStress = []
maxStrain = []
TsaiWu = []
sigma = zeros(3)
tiny = 1.0e-6
```

Loop over all angles, calculate the stress components, and evaluate the
failure criteria:

``` python
for t in theta:
    trad = t * pi / 180

    sigma[0] = cos(trad) ** 2
    sigma[1] = sin(trad) ** 2
    sigma[2] = -sin(trad) * cos(trad)

    # Maximum Stress
    fi = 0.0
    sf_low, sf_high = 0.0, 1.0e15

    while fi > 1.0 + tiny or fi < 1.0 - tiny:
        sf = 0.5 * (sf_low + sf_high)
        fi = glassepoxy.getFIMaximumStress(sf * sigma)

        if fi > 1.0:
            sf_high = sf
        else:
            sf_low = sf

    maxStress.append(sf)

    # Maximum Strain
    fi = 0.0
    sf_low, sf_high = 0.0, 1.0e15

    while fi > 1.0 + tiny or fi < 1.0 - tiny:
        sf = 0.5 * (sf_low + sf_high)
        fi = glassepoxy.getFIMaximumStrain(sf * sigma)

        if fi > 1.0:
            sf_high = sf
        else:
            sf_low = sf

    maxStrain.append(sf)

    # Tsai-Wu
    fi = 0.0
    sf_low, sf_high = 0.0, 1.0e15

    while fi > 1.0 + tiny or fi < 1.0 - tiny:
        sf = 0.5 * (sf_low + sf_high)
        fi = glassepoxy.getFITsaiWu(sf * sigma)

        if fi > 1.0:
            sf_high = sf
        else:
            sf_low = sf

    TsaiWu.append(sf)
```

Plot the results for each failure criterion:

``` python
import matplotlib.pyplot as plt

plt.figure()
plt.plot(theta, maxStress, 'r-.', label="Max. Stress")
plt.plot(theta, maxStrain, 'b--', label="Max. Strain")
plt.plot(theta, TsaiWu, 'g', label="Tsai-Wu")
plt.xlabel('Angle (degrees)')
plt.ylabel('Strength [Pa]')
plt.legend()
plt.show()
```

This generates a plot showing the strength values for each failure
criterion as a function of angle.
