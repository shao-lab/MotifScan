import os
from setuptools import setup, Extension, find_packages

os.environ['CC'] = 'g++'
os.environ['LDSHARED'] = 'g++'

VERSION = "1.1"
INSTALL_REQUIRES = ['numpy',
                    'scipy',
                    'pandas>=0.17.0',
                    'matplotlib'
                    ]

CLASSIFIERS = ['Development Status :: 3 - Alpha',
               'Environment :: Console',
               'Intended Audience :: Education',
               'Intended Audience :: Science/Research',
               'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
               'Operating System :: POSIX :: Linux',
               'Operating System :: Unix',
               'Programming Language :: Python',
               'Programming Language :: Python :: 2',
               'Programming Language :: Python :: 2.7',
               'Topic :: Scientific/Engineering :: Bio-Informatics',
               ]

with open('README.md', 'r') as fin:
    long_description = fin.read()

setup(name="motifscan",
      version=VERSION,
      author="Jiawei Wang",
      author_email="jerryeah@gmail.com",
      maintainer="Hayden Sun",
      maintainer_email="sunhongduo@picb.ac.cn",
      url="http://bioinfo.sibs.ac.cn/shaolab/motifscan",
      description="A motif discovery package developed by Shao lab in SIBS, CAS",
      long_description=long_description,
      license='GPLv3',
      packages=find_packages(exclude=['*.tests']),
      scripts=['bin/motifscan', 'bin/motifcompile', 'bin/genomecompile'],
      ext_modules=[Extension(name='motifscan.score_c', sources=['motifscan/lib/score.cpp'], extra_link_args=["--shared"])],
      install_requires=INSTALL_REQUIRES,
      classifiers=CLASSIFIERS,
      zip_safe=False
      )
