.. _quickstart:

==========
Quickstart
==========

.. contents::
   :local:

First, to check whether MotifScan is properly installed, you can inspect the
version of MotifScan by ``-v/--version`` option:

.. code-block:: shell

    $ motifscan --version

Configuration
=============

MotifScan requires basic data files including genome sequences and motif
PFMs (Position Frequency Matrices) to detect the binding sites of known motifs.
Before scanning, users should install genome assemblies and motif sets from
a remote database or with local prepared files.

Default Installation Location
-----------------------------

Newly installed genome assemblies are placed under ``$HOME/.motifscan/genomes/``,
if you want to change it:

.. code-block:: shell

    $ motifscan config --set-default-genome <path>

As for motif sets (PFMs/PWMs), the default path is under ``$HOME/.motifscan/motifs/``,
you can also change it with command:

.. code-block:: shell

    $ motifscan config --set-default-motif <path>

Please check :ref:`Config<userguide config>` for all the details about configurations.


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
===============

After the data preparation steps, you can now scan a set of genomic regions to
detect the occurrences of known motifs.

.. code-block:: shell

   $ motifscan scan -i regions.bed -g hg19 -m <motif_set> -o <output_dir>

.. _UCSC: https://genome.ucsc.edu/
.. _JASPAR: http://jaspar.genereg.net/