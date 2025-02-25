
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









