# B-Spline Path Generator

This README provides instructions on how to run the examples from scratch.

## Linux Prerequisites

Ensure you have Python, `pip`, git, and cmake, Google Test, Eigen3, and tkinter is installed on your system. Then, install the required packages:

```bash
sudo apt install python3
sudo apt install python3-pip
sudo apt install git
sudo apt install libgtest-dev
sudo apt install libeigen3-dev
sudo apt install python3-tk

```

## Linux Setup

Now install the pyton dependencies through pip
```bash
pip install numpy
pip install matplotlib
pip install scipy
pip install pyrsistent
pip install setuptools
```
Clone the following repository

```bash
git clone https://github.com/davidcGIThub/bspline_generator.git
```
Navigate to the working directory and install

```bash
cd bspline_generator
python3 setup.py sdist bdist_wheel
pip install .
```

Now clone this repository

```bash
git clone https://github.com/davidcGIThub/path_generator.git
```

 Build the c++ code. Navigate to

```bash
/path_generator/PathObjectivesAndConstraints
```
and then run the following commands

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build .
```
