# Installing Python, VS Code and Jupyter Notebook on Windows

This guide explains how to set up a Python programming environment on Windows using:

- Python
- Visual Studio Code (VS Code)
- Jupyter Notebook
- Python virtual environments
- LaTeX-style mathematics in Jupyter notebooks
- PDF export

The recommended setup is relatively lightweight and does not require Anaconda.

---

## 1. Install Python

Download a current stable 64-bit version of Python 3 from:

https://www.python.org/downloads/windows/

Run the installer.

> **Important:** During installation, select the option **Add python.exe to PATH**.

After installation, open **PowerShell** or the **Windows Terminal** and check that Python is available:

```powershell
python --version
python -m pip --version
```

Both commands should display version information.

---

## 2. Install Visual Studio Code

Download and install Visual Studio Code from:

https://code.visualstudio.com/

Start VS Code and open the **Extensions** panel.

Install the following Microsoft extensions:

- **Python**
- **Jupyter**

These extensions allow you to run Python programs and Jupyter notebooks directly from VS Code.

---

## 3. Install Git

Download and install Git for Windows from:

https://git-scm.com/download/win

During installation, the default options are usually suitable. After the
installation finishes, open a new **PowerShell** or **Windows Terminal**
window and verify that Git is available:

```powershell
git --version
```

The command should display the installed Git version.

Git is useful for downloading projects from GitHub and for keeping track of
changes to your own work.

---

## 4. Create a Project Folder and Python Environment

It is good practice to create a separate Python environment for each course or project.

Open PowerShell and create a project directory:

```powershell
mkdir my-python-course
cd my-python-course
```

Create a virtual environment:

```powershell
python -m venv .venv
```

Activate the environment:

```powershell
.venv\Scripts\activate
```

After activation, the command prompt should start with something similar to:

```text
(.venv) C:\...\my-python-course>
```

The `.venv` directory now contains a separate Python environment for this project.

---

## 5. Install Jupyter and Scientific Python Packages

With the virtual environment activated, first update `pip`:

```powershell
python -m pip install --upgrade pip
```

Install Jupyter and some commonly used scientific Python packages:

```powershell
python -m pip install jupyter ipykernel numpy scipy matplotlib pandas
```

Additional packages can always be installed later using:

```powershell
python -m pip install package-name
```

---

## 6. Open the Project in VS Code

From the project directory, start VS Code:

```powershell
code .
```

Alternatively, start VS Code normally and use **File → Open Folder** to open the project directory.

In VS Code, press:

```text
Ctrl + Shift + P
```

and select:

```text
Python: Select Interpreter
```

Select the Python interpreter inside the `.venv` directory:

```text
.venv\Scripts\python.exe
```

This tells VS Code to use the Python environment created specifically for this project.

---

## 7. Create a Jupyter Notebook

Create a new file with the extension:

```text
.ipynb
```

For example:

```text
introduction.ipynb
```

VS Code will automatically open the file as a Jupyter notebook.

In the upper-right corner of the notebook, check the selected **kernel**. It should correspond to the `.venv` environment created above.

You can now create and execute Python cells, for example:

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0.0, 2.0 * np.pi, 100)
y = np.sin(x)

plt.plot(x, y)
plt.xlabel("x")
plt.ylabel("sin(x)")
plt.show()
```

---

## 8. Using LaTeX Mathematics in a Notebook

Jupyter notebooks support LaTeX-style mathematical notation through **MathJax**.

You normally do **not** need to install a separate LaTeX distribution such as MiKTeX or TeX Live.

Create a **Markdown cell** in the notebook.

Inline mathematics can be written using `$...$`:

```markdown
The Cauchy stress tensor is denoted by $\boldsymbol{\sigma}$.
```

For equations on a separate line, use `$$...$$`:

```markdown
The linear momentum balance is

$$
\nabla \cdot \boldsymbol{\sigma}
+ \rho \mathbf{b}
= \rho \ddot{\mathbf{u}}.
$$
```

More complicated expressions are also possible:

```markdown
The strain and stress are given by

$$
\begin{aligned}
\boldsymbol{\varepsilon}
    &= \frac{1}{2}
       \left(
       \nabla\mathbf{u}
       + \nabla\mathbf{u}^{T}
       \right),\\
\boldsymbol{\sigma}
    &= \mathbf{C}:\boldsymbol{\varepsilon}.
\end{aligned}
$$
```

This makes Jupyter notebooks particularly useful for engineering and scientific computing, because a notebook can combine:

- explanatory text;
- mathematical equations;
- Python code;
- numerical results;
- tables;
- plots and visualisations.

---

## 9. Exporting a Jupyter Notebook to PDF

A finished Jupyter notebook can be exported to PDF. There are several ways to do this.

### Option 1: Export to HTML and Print to PDF

This is the simplest method and does **not** require a LaTeX installation.

In VS Code, open the notebook and use the **Export** option in the notebook toolbar.

Choose:

```text
HTML
```

Open the generated HTML file in a web browser.

In the browser, select:

```text
Print → Save as PDF
```

or, on Windows:

```text
Ctrl + P → Microsoft Print to PDF
```

This method generally preserves:

- formatted Markdown;
- mathematical equations;
- Python output;
- figures;
- tables.

For most course assignments and reports, this is the recommended approach.

### Option 2: Export from the Command Line

Jupyter includes the `nbconvert` tool. A notebook can be converted to HTML using:

```powershell
jupyter nbconvert --to html introduction.ipynb
```

This creates:

```text
introduction.html
```

Open this file in a browser and print it to PDF.

### Option 3: Direct LaTeX/PDF Export

Jupyter can also generate a PDF through LaTeX:

```powershell
jupyter nbconvert --to pdf introduction.ipynb
```

This approach requires a separate LaTeX installation on Windows, such as **MiKTeX** or **TeX Live**.

It can produce higher-quality typesetting, particularly for documents containing many mathematical equations, but the installation is considerably larger and more complicated.

For beginning Python users, installing a complete LaTeX system solely for PDF export is usually unnecessary.

> **Recommended:** Use **HTML → Print → PDF** unless you specifically need LaTeX-quality document production.

---

## 10. Recommended Project Structure

A simple project or course directory could look like:

```text
my-python-course/
│
├── .venv/
├── notebooks/
│   ├── introduction.ipynb
│   └── examples.ipynb
├── src/
│   └── functions.py
├── data/
└── requirements.txt
```

The `.venv` directory contains the Python environment.

The `notebooks` directory contains Jupyter notebooks, while reusable Python functions and classes can be placed in `src`.

---

## 11. Python Interpreter versus Jupyter Kernel

There are three related components in this setup:

```text
Python installation
        │
        ▼
Virtual environment (.venv)
        │
        ├── Python interpreter
        │
        └── installed packages
                │
                ▼
          Jupyter kernel
                │
                ▼
             VS Code
```

When working in VS Code, it is important that both the **Python interpreter** and the **Jupyter kernel** refer to the correct virtual environment.

If Python code works from PowerShell but not from a notebook, checking the selected Jupyter kernel is therefore a good first troubleshooting step.

---

## 12. Summary

For most scientific and engineering applications, this setup provides a clean and flexible Python environment without requiring Anaconda.

A separate LaTeX installation is **not required** for displaying mathematical equations inside Jupyter notebooks.

For PDF submission, the easiest workflow is:

```text
Jupyter Notebook
      ↓
Export to HTML
      ↓
Open in browser
      ↓
Print / Save as PDF
```

A full LaTeX installation is only required when using Jupyter's direct LaTeX-based PDF export.
