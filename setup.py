from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()
    
setup(
    name="pycomposites",
    version="2026.9",
    description="Tools to analyse the thermo-mechanical behaviour of composite laminates.",
    packages=find_packages(),
    install_requires=[
        "numpy",
    ],
    extras_require={
        "examples": [
            "matplotlib",
            "jupyter",
        ],
        "docs": [
            "myst-parser",
            "sphinx",
            "sphinx-rtd-theme",
        ],
    },
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
    python_requires=">=3.10",
)
