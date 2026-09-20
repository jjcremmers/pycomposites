# Examples

The source code for the examples is provided in the directories
`examples/jupyter` and `examples/python`.

These examples include practical exercises covering:

- **Classical laminate theory**: stiffness matrices, stress-strain relations,
  and laminate response.
- **Failure analysis**: application of common failure criteria to composite
  laminates.
- **Plate equations**: numerical examples of plate bending and deflection
  problems.

The Jupyter notebooks (`examples/jupyter`) are interactive and well suited for
teaching, self-study, and demonstration purposes.

The Python scripts (`examples/python`) provide standalone implementations that
can be executed directly, making them useful as templates for custom projects.

````{note}
To run the Jupyter notebooks, install the optional example dependencies:

```bash
python -m pip install ".[examples]"
jupyter lab
```
````

```{toctree}
:maxdepth: 2

clt
failure
plateeq
```
