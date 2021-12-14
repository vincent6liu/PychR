import sys
from setuptools import setup
from warnings import warn

if sys.version_info.major != 3:
    raise RuntimeError("PychR requires Python 3")
if sys.version_info.minor < 9:
    warn("Analysis methods were developed using Python 3.9")

# get version
with open("src/pychr/version.py") as f:
    exec(f.read())

setup(
    name="PychR",
    version=__version__,  # read in from the exec of version.py; ignore error
    description=(
        "Utility functions to extract data from Arrow files and ArchR projects."
    ),
    url="https://github.com/vincent6liu/PychR",
    author="Vincent Liu",
    author_email="liuv@stanford.edu",
    package_dir={"": "src"},
    packages=["pychr"],
    install_requires=[
        "numpy>=1.20.2",
        "pandas>=1.2.4",
        "scipy>=1.6.3",
        "h5py"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires='>=3.9',
)