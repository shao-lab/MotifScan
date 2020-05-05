MotifScan
=========

.. image:: https://github.com/shao-lab/MotifScan/workflows/Python%20package/badge.svg
   :alt: Github Actions
   :target: https://github.com/shao-lab/MotifScan/actions
.. image:: https://readthedocs.org/projects/motifscan/badge/?version=latest
   :alt: Documentation Status
   :target: https://motifscan.readthedocs.io/en/latest/?badge=latest
.. image:: https://img.shields.io/pypi/v/motifscan.svg
   :alt: PyPI
   :target: https://pypi.org/project/motifscan/
.. image:: https://img.shields.io/pypi/pyversions/motifscan.svg
   :alt: Python Version
   :target: https://pypi.org/project/motifscan/
.. image:: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square
   :alt: Bioconda
   :target: http://bioconda.github.io/recipes/motifscan/README.html
.. image:: https://codecov.io/gh/shao-lab/MotifScan/branch/master/graph/badge.svg
   :alt: Codecov
   :target: https://codecov.io/gh/shao-lab/MotifScan
.. image:: https://img.shields.io/pypi/l/motifscan.svg
   :alt: License
   :target: https://github.com/shao-lab/MotifScan/blob/master/LICENSE

Introduction
============

**Scan input genomic regions with known DNA motifs**

Given a set of input genomic regions, MotifScan scans the sequences to
detect the occurrences of known motifs. It can also applies a statistical test
on each motif to check whether the motif is significantly over- or under-represented
(enriched or depleted) in the input genomic regions compared to another set of control
regions.

Citation
========

`Sun, H., Wang, J., Gong, Z. et al. Quantitative integration of epigenomic variation and
transcription factor binding using MAmotif toolkit identifies an important role of IRF2
as transcription activator at gene promoters. Cell Discov 4, 38 (2018).`__

.. __: https://doi.org/10.1038/s41421-018-0045-y

Documentation
=============

To see the full documentation of MotifScan, please refer to: https://motifscan.readthedocs.io

Installation
============

The latest version release of MotifScan is available at
`PyPI <https://pypi.python.org/pypi/motifscan>`__:

::

    $ pip install motifscan

Or you can install MotifScan via conda:

::

    $ conda install -c bioconda motifscan

Usage
=====

Install genome assemblies
-------------------------

Install from a remote database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can download genome assemblies from the `UCSC`_ database.

First, display all available genome assemblies:

.. code-block:: shell

    $ motifscan genome --list-remote

Then, install a genome assembly (e.g. hg19):

.. code-block:: shell

    $ motifscan genome --install -n hg19 -r hg19

Install with local files
^^^^^^^^^^^^^^^^^^^^^^^^

To install a genome assembly locally, you have to prepare a FASTA file
containing the genome sequences and a genome annotation file (refGene.txt).

.. code-block:: shell

    $ motifscan genome --install -n hg19 -i <hg19.fa> -a <refGene.txt>


Install and build motif sets
----------------------------

Install from a remote database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Users can install motif PFMs sets in the `JASPAR`_ 2020 database.

First, display all available motif PFMs sets in JASPAR 2020:

.. code-block:: shell

    $ motifscan motif --list-remote

Then, install a JASPAR motif PFMs set (e.g. vertebrates_non-redundant):

.. code-block:: shell

    $ motifscan motif --install -n <motif_set> -r vertebrates_non-redundant -g hg19


Install with local files
^^^^^^^^^^^^^^^^^^^^^^^^

Install a motif set with local PFMs file:

.. code-block:: shell

   $ motifscan motif --install -n <motif_set> -i <pfms.jaspar> -g hg19

Build PFMs for additional genome
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Build the motif PFMs set for another installed genome assembly hg38:

.. code-block:: shell

   $ motifscan motif --build <motif_set> -g hg38

Scanning Motifs
---------------

After the data preparation steps, you can now scan a set of genomic regions to
detect the occurrences of known motifs.

.. code-block:: shell

   $ motifscan scan -i regions.bed -g hg19 -m <motif_set> -o <output_dir>

.. _UCSC: https://genome.ucsc.edu/
.. _JASPAR: http://jaspar.genereg.net/

**Note:** Using -h/--help for the details of all arguments.


License
=======

`BSD 3-Clause
License <https://github.com/shao-lab/MotifScan/blob/master/LICENSE>`__
