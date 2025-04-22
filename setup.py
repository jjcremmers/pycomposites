from setuptools import find_packages, setup

with open("README.md","r") as f:
    long_description = f.read()
    
setup(
    name="pycomposites",
    version="1.1.0",
    description="A collection of classes and functions to anaylse the thermo mechanicalbehaviour of composites.",
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jjcremmers/pycomposites",
    author="Joris Remmers",
    author_email="j.j.c.remmers@tue.nl",    
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ], 
    python_requires='>=3.10',
)

