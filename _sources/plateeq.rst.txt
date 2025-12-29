Plate Equations
===============

.. currentmodule:: pycomposites.composite

**Problem Statement**

A flat, simply supported panel with dimensions :math:`a = 400 \, \text{mm}` and :math:`b = 600 \, \text{mm}` is loaded by a point load with magnitude :math:`F = 1200 \, \text{N}` at the point :math:`x_0 = 150 \, \text{mm}` and :math:`y_0 = 350 \, \text{mm}`. The panel is made of 6 layers of UCHSC200_SE84 UD material (see Canvas for the properties), stacked in the sequence :math:`[0, 90, 0]_S`.

Calculate the out-of-plane displacement of the panel. Plot the results as a contour plot.

**Solution**

Create the transverse isotropic material model and the laminate. Calculate the :math:`\mathbf{D}` matrix and the parameters :math:`D_1`, :math:`D_2`, and :math:`D_3`.

.. code-block:: python

    from composite import TransverseIsotropic, Laminate

    compUD = TransverseIsotropic([130e9, 7.2e9], 0.337, 4.2e9, [0.57e-6, 35.1e-6], 1514.)
    print(compUD)

    lam = Laminate()
    lam.addMaterial('UCHSC200_SE84', compUD)
    lam.addLayer('UCHSC200_SE84', 0, 0.2e-3)
    lam.addLayer('UCHSC200_SE84', 90, 0.2e-3)
    lam.addLayer('UCHSC200_SE84', 0, 0.2e-3)
    lam.addLayer('UCHSC200_SE84', 0, 0.2e-3)
    lam.addLayer('UCHSC200_SE84', 90, 0.2e-3)
    lam.addLayer('UCHSC200_SE84', 0, 0.2e-3)

    print(lam)

    D = lam.getD()
    D1 = D[0, 0]
    D2 = D[1, 1]
    D3 = D[0, 1] + 2 * D[2, 2]

    print("D1: ", D1, " D2: ", D2, " D3: ", D3)

The output of this snippet of code is:

::

    Elastic Properties:
    -----------------------------------------------------------
    E1     :     1.300e+11 , E2     :     7.200e+09 
    nu12   :          0.34 , G12    :     4.200e+09 
    rho    :       1514.00

    Thermal expansion coefficients:
    -----------------------------------------------------------
    alpha1 :     5.700e-07 , alpha2 :     3.510e-05 

    Laminate properties:
    -----------------------------------------------------------
    layer   thick orient.  material
    -----------------------------------------------------------
        0   0.0002      0   UCHSC200_SE84
        1   0.0002     90   UCHSC200_SE84
        2   0.0002      0   UCHSC200_SE84
        3   0.0002      0   UCHSC200_SE84
        4   0.0002     90   UCHSC200_SE84
        5   0.0002      0   UCHSC200_SE84

    D1:  14.224941196641147  D2:  5.656915190635289  D3:  1.5612132386157986

Next, the dimensions of the panel and the loacation and magnitude of the point load
are entered.

.. code-block:: python

    F = -1200
    a = 0.4
    b = 0.6
    x0 = 0.15
    y0 = 0.35

With these variables, the terms :math:`\mathbb{A}` and :math:`\mathbb{B}` can be 
calculated. In the analysis we will use the first 5 terms in both :math:`x` and :math:`y` 
direction, i.e. :math:`m=1..5` and :math:`n=1..5`.

.. math::

    \mathbb{B}_{mn} = \frac{4F}{ab} \sin \left( \frac{m \pi x_0}{a} \right) \sin \left( \frac{n \pi y_0}{b} \right)

.. math::

    \mathbb{A}_{mn} = \frac{\mathbb{B}_{mn}}{\pi^4 \left( D_1 \frac{m^4}{a^4} + 2 D_3 \frac{m^2 n^2}{a^2 b^2} + D_2 \frac{n^4}{b^4} \right)}

In the python code, this looks as follows:

.. code-block:: python

    import numpy as np
    from math import sin, pi

    nmax = 6
    B = np.zeros((nmax+1, nmax+1))
    A = np.zeros((nmax+1, nmax+1))

    for m in range(1, nmax):
        for n in range(1, nmax):
            B[m, n] = 4 * F * sin(m * pi * x0 / a) * sin(n * pi * y0 / b) / (a * b)
            A[m, n] = B[m, n] / (pi**4 * (D1 * (m / a)**4 + 2 * D3 * (m / a)**2 * (n / b)**2 + D2 * (n / b)**4))

Note that we create numpy matrices of dimensions 6 by 6 for this. The first row and column remain unused.

The displacement field :math:`w(x,y)` and the distributed load can 
.. math::

    w(x, y) = \sum_{m=1}^{5} \sum_{n=1}^{5} \mathbb{A}_{mn} \sin \left( \frac{m \pi x}{a} \right) \sin \left( \frac{n \pi y}{b} \right)

.. math::

    q(x, y) = \sum_{m=1}^{5} \sum_{n=1}^{5} \mathbb{B}_{mn} \sin \left( \frac{m \pi x}{a} \right) \sin \left( \frac{n \pi y}{b} \right)

.. code-block:: python

    nplot = 40
    xplot = np.linspace(0, a, nplot)
    yplot = np.linspace(0, b, nplot)
    q = np.zeros((nplot, nplot))
    w = np.zeros((nplot, nplot))

    for i, y in enumerate(yplot):
        for j, x in enumerate(xplot):
            for m in range(1, nmax+1):
                for n in range(1, nmax+1):
                    q[i, j] += B[m, n] * sin(m * pi * x / a) * sin(n * pi * y / b)
                    w[i, j] += A[m, n] * sin(m * pi * x / a) * sin(n * pi * y / b)


The load :math:`q(x, y)` can be plotted as a 3D surface plot or a contour plot.

3D Surface Plot:

.. code-block:: python

    import matplotlib.pyplot as plt
    from matplotlib import cm

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    x, y = np.meshgrid(xplot, yplot)
    surf = ax.plot_surface(x, y, -q, alpha=0.8, cmap=cm.coolwarm)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('q')
    plt.title('Point Load q(x, y)')
    plt.show()

Contour Plot:

.. code-block:: python

    fig, ax = plt.subplots()
    contourplot = ax.contourf(x, y, q, cmap=cm.coolwarm)
    fig.colorbar(contourplot)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.title('Point Load q(x, y)')
    plt.show()

Out-of-plane Displacement:

.. code-block:: python

    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    surf = ax.plot_surface(x, y, w, alpha=0.8, cmap=cm.coolwarm)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('w')
    plt.title('Vertical Displacement w(x, y)')
    plt.show()

