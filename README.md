
# eVTOL BSplines

A general repository containing many useful tools for creating BSplines, which are in turn useful for eVTOL robotics applications


## Linux Setup




### Dependencies Setup


Install the following dependencies libraries using pip

```
pip install numpy
pip install matplotlib
pip install scipy
pip install pyrsistent
pip install setuptools
```



This library requires the use of the bspline_generator class, created by David Christensen. This library is available on Github, and currently must be manually pip installed on your local machine. To do this:


Clone the following repository
```
git clone https://github.com/davidcGIThub/bspline_generator.git
```


Navigate to the working directory and install it using:

```
cd bspline_generator
python3 setup.py sdist bdist_wheel
pip install .
```


### Main Repository Setup

In your home directory, clone the eVTOL_BSplines main repository

```
git clone https://github.com/byu-magicc/eVTOL_BSplines.git
```

Navigate into the eVTOL_BSplines directory

```
cd eVTOL_BSplines
```

Remove old build artifacts and install the package using pip.

```
rm -rf build dist *.egg-info
pip install -e .
```

## Library Use Guide


At present, the limited functionality of the wrapper repository is contained within the `path_generator_simplified` file. To access it, import as:

```
from eVTOL_BSplines.path_generator_simplified import path_generator_simplified
```

The class `path_generator_simplified` currently contains just a waypoint generator function. The vision for `path_generator_simplified` is to contain simplified and wrapped versions of the clunky functions from within the main `path_generator` repository.








