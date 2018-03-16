MotifScan
=========

|pypi| |license|

.. |pypi| image:: https://img.shields.io/pypi/v/motifscan.svg
   :target: https://pypi.python.org/pypi/motifscan

.. |license| image:: https://img.shields.io/pypi/l/MAnorm.svg
   :target: https://github.com/shao-lab/MAnorm/blob/master/LICENSE

Introduction
------------

MotifScan is a precise and easy-use motif discovery tool based on given motifs. To search for candidate motif targets
in given DNA sequences, the program scans them with a window of the same length as the motif, and defined raw motif
score of the sequence S in window as the ratio of the probability to observe target sequence S given the motif's
Position Weight Matrix (PWM) M and the probability to observe S given the genome background B. For each annotated motif,
we modeled the genome background distribution of motif score by randomly sampling the genome for 10^6 times, and defined
the targets of this motif as those candidates whose motif score was higher than the cutoff. And the enrichment of each
motif was represented by the ratio of motif target densities at input regions compared to random control regions,
together with a p-value calculated from hypergeometric distribution. It is noticeable that MotifScan is not a
de novo motif discovery tool.

Documentation
-------------

To see the full documentation of MAnorm, please refer to: http://motifscan.readthedocs.io/en/latest/

Installation
------------

The latest version release of MotifScan is available at
`PyPI <https://pypi.python.org/pypi/motifscan>`__:

::

    $ pip install motifscan

MAnorm uses `setuptools <https://setuptools.readthedocs.io/en/latest/>`__ for installation from source code.
The source code of MAnorm is hosted on GitHub: https://github.com/shao-lab/MotifScan

You can clone the repo and execute the following command under source directory:

::

    $ python setup.py install

Usage
-----

::

    $

**Note:** Using -h/--help for the details of all arguments.

License
-------

`BSD 3-Clause
License <https://github.com/shao-lab/MotifScan/blob/master/LICENSE>`__


