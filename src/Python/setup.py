from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

ext_modules = [
    Pybind11Extension(
        "wows_shell",
        ["wows_shell_python.cpp"],
        cxx_std="c++17"
    ),
]

setup(
    name="wows_shell",
    version="1.1.2",
    description="Python extension of wows_shell",
    url="https://github.com/jcw780/wows_shell",
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules
)
