import setuptools

setuptools.setup(
    include_package_data=True,
    name="composites",
    version="1.0.0",
    author="Joris Remmers",
    author_email="j.j.c.remmers@tue.nl",
    description="A collection of class and function for the thermo-mechanical analysis of composite materials",
    long_description="pycomposites module",
    long_description_content_type="text/markdown",
    #url="https://github.com/yourusername/dawn",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],    
    python_requires='>=3.6',
)

