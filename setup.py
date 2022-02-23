import sys
import cpufeature
from pathlib import Path
from glob import glob
from setuptools import setup
from pybind11.setup_helpers import Pybind11Extension, build_ext

WIN = sys.platform.startswith("win32")

ext_modules = [
    Pybind11Extension(
        "wows_shell",
        ["src/Python/wows_shell_python.cpp"],
        cxx_std="17",
        include_dirs=["src"]
    ),
]

cxx_compile_args = ext_modules[0].extra_compile_args
if WIN:
    if cpufeature.CPUFeature['AVX2']:
        cxx_compile_args.append('/arch:AVX2')
    elif cpufeature.CPUFeature['AVX']:
        cxx_compile_args.append('/arch:AVX')
else:
    cxx_compile_args.append('-march=native')

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="wows_shell",
    version="1.2.0",
    description="Python extension of wows_shell",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/jcw780/wows_shell",
    cmdclass={"build_ext": build_ext},
    ext_modules=ext_modules,
)
