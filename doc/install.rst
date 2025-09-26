Installation Guide
==================

PyComposites is a Python library for the analysis and design of composite materials.  
It can be installed directly from the
`GitHub source <https://github.com/jjcremmers/pycomposites>`_.

Only the **Python API** is provided (no command-line interface).

Requirements
------------

- Python **3.6 or newer**
- ``pip`` (Python package manager)
- (Optional) a virtual environment for isolating dependencies:

.. code-block:: bash

   python3 -m venv .venv
   source .venv/bin/activate    # Linux / macOS
   .venv\Scripts\activate       # Windows PowerShell

Installation Steps
------------------

1. **Clone the repository**

   .. code-block:: bash

      git clone https://github.com/jjcremmers/pycomposites.git
      cd pycomposites

2. **Install with pip**

   .. code-block:: bash

      pip install .

   For developers who want to make local changes and test immediately, install in editable mode:

   .. code-block:: bash

      pip install -e .

Verifying the Installation
--------------------------

After installation, verify that the API is available by starting Python and importing the package:

.. code-block:: python

   from composite import Layer
   print(Layer)

Quick Start Example
-------------------

The following example shows how to use PyComposites after installation.
Replace the function names with the relevant routines from your project.

.. code-block:: python

   from composites import Laminate, Material

   # Define a unidirectional carbon/epoxy ply
   carbon_epoxy = Material(E1=135e9, E2=10e9, G12=5e9, nu12=0.3)

   # Create a simple cross-ply laminate [0/90/0]
   laminate = Laminate([carbon_epoxy, carbon_epoxy, carbon_epoxy],
                       [0, 90, 0])

   # Compute stiffness matrix
   A, B, D = laminate.abd_matrix()

   print("A-matrix:")
   print(A)

If this example runs without errors, the installation is successful.

Updating PyComposites
---------------------

To update to the latest version from GitHub:

.. code-block:: bash

   cd pycomposites
   git pull origin main
   pip install --upgrade .

