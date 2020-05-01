.. _faq:

FAQ
===

How can I get help?
-------------------

See :ref:`get help`.

I can not install/run MotifScan on Windows
------------------------------------------
MotifScan does not support the Windows platform. Try `WSL`_ for a workaround if
you only have windows systems.

The Python version I use is not the same as MotifScan requires
--------------------------------------------------------------

MotifScan requires Python 3.6+. If you are not using these versions, it is recommended to
create a separated Python environment with `virtualenv`_ or `conda`_.

Can MotifScan discover *de novo* motifs
---------------------------------------

No, MotifScan can only scan with known motifs given the position matrices.

What data files do MotifScan require?
-------------------------------------

MotifScan requires genome data files and motif data files. The genome data files
include a sequence FASTA file and an optional refGene gene annotation file if you
want to perform the motif enrichment analysis. And motif PFMs (Position frequency
matrices) are required to be scanned with. All these data files have to be configured
(installed) properly before scanning. See :ref:`User Guide<userguide>` for more details.


.. _WSL: https://docs.microsoft.com/en-us/windows/wsl/about
.. _virtualenv: https://pypi.org/project/virtualenv/
.. _conda: https://conda.io
