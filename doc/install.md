# Installation Guide

PyComposites is a Python library for the analysis and design of
composite materials. It can be installed directly from the [GitHub
source](https://github.com/jjcremmers/pycomposites).

Only the **Python API** is provided (no command-line interface).

## Requirements

- Python **3.10 or newer**
- `pip` (Python package manager)
- `git` if you want to download the repository with `git clone`
- NumPy, installed automatically by `pip install .`
- Matplotlib and Jupyter for the plotting and notebook examples,
  installed by `pip install ".[examples]"`

## Installation Steps

These steps are written for users with little Python experience. Use one
of the two paths below:

- **Windows PowerShell** if you installed Python directly on Windows.
- **WSL** if you use Ubuntu or another Linux distribution inside Windows
  Subsystem for Linux.

### Windows PowerShell

Open **PowerShell** and run these commands.

1.  **Check that Python is installed**

    ``` powershell
    python --version
    ```

    You should see Python 3.10 or newer, for example `Python 3.12.3`.

    If `python` is not found, install Python from
    <https://www.python.org/downloads/>. During installation, enable
    **Add python.exe to PATH**.

2.  **Download PyComposites**

    With git:

    ``` powershell
    git clone https://github.com/jjcremmers/pycomposites.git
    cd pycomposites
    ```

    Without git, download the repository as a ZIP file from GitHub,
    unzip it, and open a terminal in the unzipped `pycomposites` folder.

3.  **Create and activate a virtual environment**

    A virtual environment keeps PyComposites and its dependencies
    separate from the rest of your computer.

    ``` powershell
    python -m venv .venv
    .\.venv\Scripts\Activate.ps1
    ```

    If PowerShell blocks activation, run this once and then activate
    again:

    ``` powershell
    Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
    ```

    After activation, your terminal prompt usually starts with
    `(.venv)`.

4.  **Upgrade pip**

    ``` powershell
    python -m pip install --upgrade pip
    ```

5.  **Install PyComposites**

    ``` powershell
    python -m pip install .
    ```

    To run all examples, including plotting scripts and notebooks,
    install the optional example dependencies:

    ``` powershell
    python -m pip install ".[examples]"
    ```

    For developers who want to make local changes and test immediately,
    install in editable mode:

    ``` powershell
    python -m pip install -e ".[examples]"
    ```

### WSL

Open your WSL distribution, for example **Ubuntu**, and run these
commands in the Linux terminal.

1.  **Check that Python is installed**

    ``` bash
    python3 --version
    ```

    You should see Python 3.10 or newer. If Python is missing, install
    it with:

    ``` bash
    sudo apt update
    sudo apt install python3 python3-venv python3-pip git
    ```

2.  **Download PyComposites**

    With git:

    ``` bash
    git clone https://github.com/jjcremmers/pycomposites.git
    cd pycomposites
    ```

    Without git, download the repository as a ZIP file on Windows, unzip
    it, and copy or move it into your WSL home directory. Working inside
    the Linux home directory, for example `/home/yourname/pycomposites`,
    is usually faster and more reliable than working under `/mnt/c/...`.

3.  **Create and activate a virtual environment**

    ``` bash
    python3 -m venv .venv
    source .venv/bin/activate
    ```

    After activation, your terminal prompt usually starts with
    `(.venv)`.

4.  **Upgrade pip**

    ``` bash
    python -m pip install --upgrade pip
    ```

5.  **Install PyComposites**

    ``` bash
    python -m pip install .
    ```

    To run all examples, including plotting scripts and notebooks,
    install the optional example dependencies:

    ``` bash
    python -m pip install ".[examples]"
    ```

    For developers who want to make local changes and test immediately,
    install in editable mode:

    ``` bash
    python -m pip install -e ".[examples]"
    ```

## Verifying the Installation

Complete the installation by running these tests from the repository
root while the virtual environment is activated. The commands are the
same in Windows PowerShell and WSL after the virtual environment is
active.

1.  **Import test**

    This checks that Python can find the installed package.

    ``` bash
    python -c "from pycomposites import TransverseIsotropic, Laminate; material = TransverseIsotropic([135e9, 10e9], 0.3, 5e9); laminate = Laminate(); laminate.addMaterial('UD', material); laminate.addLayer('UD', 0.0, 0.125e-3); laminate.addLayer('UD', 90.0, 0.125e-3); print(laminate.getA())"
    ```

    If this command prints a matrix and no error message, the package
    import is working.

2.  **Example test**

    Run a standalone example:

    ``` bash
    python examples/python/CLT/example1.py
    ```

    This should print material properties for carbon, epoxy, and the
    homogenised composite material.

3.  **Unit test suite**

    Run all automated tests:

    ``` bash
    python -m unittest discover -s test -v
    ```

    The installation is complete when the final lines end with `OK`.

## Troubleshooting

- If `python` is not found, try `python3`.
- If `pip` installs into the wrong Python, always use
  `python -m pip ...` instead of only `pip ...`.
- If imports fail, check that the virtual environment is activated and
  that you installed the package from the repository root.
- If a plotting example fails with a Matplotlib error, install the
  optional example dependencies with
  `python -m pip install ".[examples]"`.
- If PowerShell blocks activation on Windows, run
  `Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser`
  and then run `.\.venv\Scripts\Activate.ps1` again.

## Quick Start Example

The following example shows how to use PyComposites after installation.

``` python
from pycomposites import TransverseIsotropic, Laminate

# Define a unidirectional carbon/epoxy ply material.
carbon_epoxy = TransverseIsotropic([135e9, 10e9], 0.3, 5e9)

# Create a simple cross-ply laminate [0/90/0].
laminate = Laminate()
laminate.addMaterial("UD", carbon_epoxy)
laminate.addLayer("UD", 0.0, 0.125e-3)
laminate.addLayer("UD", 90.0, 0.125e-3)
laminate.addLayer("UD", 0.0, 0.125e-3)

# Compute stiffness matrices.
A = laminate.getA()
B = laminate.getB()
D = laminate.getD()

print("A-matrix:")
print(A)
```

If this example runs without errors, the installation is successful.

## Updating PyComposites

To update to the latest version from GitHub:

``` bash
cd pycomposites
git pull origin main
python -m pip install --upgrade .
```

## Leaving the Virtual Environment

When you are finished, deactivate the virtual environment:

``` bash
deactivate
```
