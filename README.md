endo
====

This software aims to simplify the analysis of genome-wide data
provided by the TCGA (http://cancergenome.nih.gov/). This is
achieved by loading and validating the available data and
providing a range of methods for querying, analyzing, and
exporting it.

Data Access
-----------
Please note that some of the data accessed and analyzed by this
software might not be freely available. One must apply to the
TCGA for access to the data.

Prerequisites
-------------
This software requires a working installation of the SWI-Prolog
environment licensed under the Lesser GNU Public License and
available at www.swi-prolog.org.

License
-------
This software is licensed under the MIT license included in the
LICENSE file in the top-level directory of the distribution and
at http://opensource.org/licenses/MIT.

Using This Software
-------------------
Installation of the software is not necessary. Start the
SWI-Prolog top-level interpreter in the top-level directory of
the project and consult the load.pl file:

[source]
----
$ cd endo
$ swipl
?- [load].
----

The tools provided by this software are documented extensively in
the comments in each source code file. Attempting to use the
software without first reading the commented code is done at your
own risk.
