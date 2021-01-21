ChangeLog
=========

v1.3.0 (2021-01-21)
-------------------

* Use pthread for concurrent scanning (significant speed improvement)
* Optimize for memory usage
* Support for multiple rounds sampling to calculate score cutoffs
* Add matrix id in output files
* Bug fixes and other performance improvements

v1.2.2 (2020-05-24)
-------------------

* Replace new special characters '/' and '*' to '_' for motif names
* Better type definition in score.C

v1.2.1 (2020-05-18)
-------------------

* Fix a variable definition in C loop (C99)

v1.2.0 (2020-05-05)
-------------------

* Brand-new MotifScan package for Python version 3.6+.
* Use a config module to manage the genome and motif data paths.
* Download genome assembly from UCSC and motif sets from JASPAR2020.
* Support multiple input formats (BED, MACS, MACS2, narrowPeak etc.).
* Enable to run multiple processes when scanning motifs.
* Speed up calculating motif scores with a C extension.
* Bug fixes and performance improvements.
