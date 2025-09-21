Installation Guide
==================

PyComposites can be installed directly from the
`GitHub source <https://github.com/jjcremmers/pycomposites>`_.
It provides a **Python API** for composite material analysis.

Requirements
------------

- Python **3.9 or newer**
- ``pip`` (Python package manager)
- (Optional but recommended) a virtual environment:

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

   For developers who want to make local changes and test immediately, 
   install in editable mode:

   .. code-block:: bash

      pip install -e .

Verifying the Installation
--------------------------

After installation, verify that the API is available.

**Example usage**

.. code-block:: python

   import pycomposites

   # Example: access a function or class
   # (replace with actual usage in your project)
   result = pycomposites.__version__
   print("PyComposites version:", result)

If this runs without errors, the installation is successful.

Updating PyComposites
---------------------

To update to the latest version from GitHub:

.. code-block:: bash

   cd pycomposites
   git pull origin main
   pip install --upgrade .
