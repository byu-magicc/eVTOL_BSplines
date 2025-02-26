from setuptools import setup, find_packages
from setuptools.command.build_ext import build_ext
import subprocess
import os
import shutil



#creates the function to build the external CMakes
class BuildCMakeExt(build_ext):
    #creates the run function
    def run(self):

        #creates the build_directory
        path_constraints_directory = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "eVTOL_BSplines",
            "submodules",
            "path_generator",
            "PathObjectivesAndConstraints"
        )
        #creates the build directory
        build_dir = os.path.join(path_constraints_directory, "build")


        # Run CMake and build
        self.spawn(["cmake", ".."])
        self.spawn(["make", "-j$(nproc)"])

        # Ensure .so file is copied into package directory
        so_file = os.path.join(build_dir, "libPathObjectivesAndConstraints.so")
        target_dir = os.path.join(build_dir, "..", "python_wrappers")

        if not os.path.exists(target_dir):
            os.makedirs(target_dir)
        shutil.copy(so_file, target_dir)


#creates the setup function
setup(
    name="eVTOL_BSplines",
    version="0.1.0",
    packages=find_packages(include=["eVTOL_BSplines", "eVTOL_BSplines.*"]),
    include_package_data=True,
    install_requires=["numpy"],
    cmdclass={"build_ext": BuildCMakeExt}
)