# PyComposites

PyComposites is a teaching-oriented Python package for thermo-mechanical
analysis of thin-walled composite materials.

The package contains tools for Classical Laminate Theory (CLT), laminate
failure criteria, and simple plate-equation examples. It is intended for
students who need to install the code locally and run the examples during a
course or project.

## Requirements

- Python 3.10 or newer
- pip
- git, if you want to download the repository from GitHub

NumPy is installed automatically with the package. The plotting and notebook
examples additionally use Matplotlib and Jupyter.

## Installation

These steps assume little or no Python experience. Use one of the two paths
below:

- **Windows PowerShell** if you installed Python directly on Windows.
- **WSL** if you use Ubuntu or another Linux distribution inside Windows
  Subsystem for Linux.

### Windows PowerShell

Open **PowerShell** and run these commands.

#### 1. Check that Python is installed

```powershell
python --version
```

You should see Python 3.10 or newer, for example `Python 3.12.3`.

If `python` is not found, install Python from <https://www.python.org/downloads/>.
During installation, enable **Add python.exe to PATH**.

#### 2. Download PyComposites

With git:

```powershell
git clone https://github.com/jjcremmers/pycomposites.git
cd pycomposites
```

Without git, download the repository as a ZIP file from GitHub, unzip it, and
open a terminal in the unzipped `pycomposites` folder.

#### 3. Create and activate a virtual environment

A virtual environment keeps this package and its dependencies separate from
the rest of your computer.

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
```

If PowerShell blocks activation, run this once and then activate again:

```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

After activation, your prompt usually starts with `(.venv)`.

#### 4. Upgrade pip

```powershell
python -m pip install --upgrade pip
```

#### 5. Install PyComposites

```powershell
python -m pip install .
```

To run plotting examples and Jupyter notebooks as well, install the optional
example dependencies:

```powershell
python -m pip install ".[examples]"
```

For local development, use editable mode instead:

```powershell
python -m pip install -e ".[examples]"
```

#### 6. Test the installation

First check that Python can import the package:

```powershell
python -c "from pycomposites import TransverseIsotropic, Laminate; material = TransverseIsotropic([135e9, 10e9], 0.3, 5e9); laminate = Laminate(); laminate.addMaterial('UD', material); laminate.addLayer('UD', 0.0, 0.125e-3); laminate.addLayer('UD', 90.0, 0.125e-3); print(laminate.getA())"
```

Then run one example:

```powershell
python examples/python/CLT/example1.py
```

Finally run the test suite:

```powershell
python -m unittest discover -s test -v
```

The installation is complete when the import check runs, the example prints
material properties, and the test suite ends with `OK`.

#### 7. Deactivate the environment when finished

```powershell
deactivate
```

### WSL

Open your WSL distribution, for example **Ubuntu**, and run these commands in
the Linux terminal.

#### 1. Check that Python is installed

```bash
python3 --version
```

You should see Python 3.10 or newer. If Python is missing, install it with:

```bash
sudo apt update
sudo apt install python3 python3-venv python3-pip git
```

#### 2. Download PyComposites

With git:

```bash
git clone https://github.com/jjcremmers/pycomposites.git
cd pycomposites
```

Without git, download the repository as a ZIP file on Windows, unzip it, and
copy or move it into your WSL home directory. Working inside the Linux home
directory, for example `/home/yourname/pycomposites`, is usually faster and
more reliable than working under `/mnt/c/...`.

#### 3. Create and activate a virtual environment

```bash
python3 -m venv .venv
source .venv/bin/activate
```

After activation, your prompt usually starts with `(.venv)`.

#### 4. Upgrade pip

```bash
python -m pip install --upgrade pip
```

#### 5. Install PyComposites

```bash
python -m pip install .
```

To run plotting examples and Jupyter notebooks as well, install the optional
example dependencies:

```bash
python -m pip install ".[examples]"
```

For local development, use editable mode instead:

```bash
python -m pip install -e ".[examples]"
```

#### 6. Test the installation

```bash
python -c "from pycomposites import TransverseIsotropic, Laminate; material = TransverseIsotropic([135e9, 10e9], 0.3, 5e9); laminate = Laminate(); laminate.addMaterial('UD', material); laminate.addLayer('UD', 0.0, 0.125e-3); laminate.addLayer('UD', 90.0, 0.125e-3); print(laminate.getA())"
python examples/python/CLT/example1.py
python -m unittest discover -s test -v
```

The installation is complete when the import check runs, the example prints
material properties, and the test suite ends with `OK`.

#### 7. Deactivate the environment when finished

```bash
deactivate
```

## Examples

Standalone Python examples are in:

- `examples/python/CLT`
- `examples/python/PlateEquations`

Run one from the repository root after installation:

```bash
python examples/python/CLT/example1.py
```

Notebook examples are in `examples/jupyter`.

## Tests

Run the test suite from the repository root with the virtual environment
activated:

```bash
python -m unittest discover -s test -v
```

## License

[LICENSE.txt](LICENSE.txt)
