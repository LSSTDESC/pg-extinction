extinction
===========================================================

An unofficial python port of "dust_getval" by Sogo Mineo for HSC.
The original "dust_getval" was written by D. Schlegel 1998.

This version has been adapted as needed for use with LSST
DC2 data 

Requirements
-----------------------------------------------------------

  * python 3.6
  * numpy
  * astropy

Usage
-----------------------------------------------------------

```
from extinction.dustval import Extinction

extinction = Extinction()

ra  = numpy.array([....]) # in radians
dec = numpy.array([....]) # in radians

ebv = extinction.get_Ebv(ra, dec) # E(B-V)
```
