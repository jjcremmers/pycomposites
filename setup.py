from setuptools import find_packages, setup

with open("README.md","r") as f:
    long_description = f.read()

setup(
    name="pyComposites",
    version="1.0.0",
    description="A collection of class and function for the thermo-mechanical analysis of composite materials",
    package_dir={"": "app"},
    packages=find_packages(where="app"),    
    long_description="pycomposites module",
    long_description_content_type="text/markdown",
    url="https://gitlab.tue.nl/jremmers/pycomposites.git",
    author="Joris Remmers",
    author_email="j.j.c.remmers@tue.nl",    
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],    
    python_requires='>=3.6',
)

