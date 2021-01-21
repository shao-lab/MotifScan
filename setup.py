import os
import re
import sys

from setuptools import find_packages, setup, Extension

py_version = sys.version_info[:2]
if py_version < (3, 6):
    raise RuntimeError("MotifScan requires Python 3.6+ to install!")

description = "A package for motif discovery and motif enrichment analysis"

pkg_dir = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(pkg_dir, 'motifscan', '__init__.py')) as fin:
    version = re.search(r"__version__ = '(.*?)'", fin.read()).group(1)

with open(os.path.join(pkg_dir, 'README.rst')) as fin:
    long_description = fin.read()

install_requires = [
    "numpy>=1.15",
    "scipy>=1.0",
    "requests",
    "tqdm>=4.42.1",
    "pysam>=0.15.0",
    "matplotlib>=3.0.0"
]

extras_require = {
    "test": ["pytest>=4.0.0",
             "pytest-cov>=2.8.0"],
    "docs": ["sphinx>=2.0.0",
             "sphinx_rtd_theme"]
}

classifiers = [
    "Development Status :: 5 - Production/Stable",
    "Environment :: Console",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: BSD License",
    "Operating System :: Unix",
    "Operating System :: POSIX",
    "Operating System :: MacOS",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python :: 3.6",
    "Programming Language :: Python :: 3.7",
    "Programming Language :: Python :: 3.8",
    "Topic :: Scientific/Engineering :: Bio-Informatics"
]

ext_modules = [
    Extension("motifscan.motif.cscore", ["motifscan/motif/cscore.c"],
              extra_compile_args=["-std=c99"])
]

setup(
    name="motifscan",
    version=version,
    description=description,
    long_description=long_description,
    long_description_content_type='text/x-rst',
    url="https://github.com/shao-lab/MotifScan",
    project_urls={
        "Bug Tracker": "https://github.com/shao-lab/MotifScan/issues",
        "Documentation": "https://motifscan.readthedocs.io",
        "Source Code": "https://github.com/shao-lab/MotifScan",
    },
    author="Hayden Sun",
    author_email="sunhongduo@picb.ac.cn",
    license="BSD",
    packages=find_packages(),
    entry_points={"console_scripts": ["motifscan=motifscan.cli.main:main"]},
    python_requires=">=3.6",
    install_requires=install_requires,
    extras_require=extras_require,
    ext_modules=ext_modules,
    classifiers=classifiers,
    zip_safe=False,
)
