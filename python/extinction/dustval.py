# Copyright (C) 2018  Sogo Mineo
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# This is an unofficial python port of "dust_getval"
# originally written by D. Schlegel, 19 Jan 1998, Durham.

import numpy
import astropy.coordinates
import astropy.io.fits
import astropy.units
import astropy.wcs

import os


class Extinction(object):
    """
    Whole-sky map of E(B-V)
    """
    def __init__(self):
        self.__dust_ngp = None
        self.__dust_sgp = None

    def get_Ebv(self, ra, dec):
        """
        @param ra  (numpy.array) in radians
        @param dec (numpy.array) in radians
        @return (numpy.array) E(B-V)
        """
        ra  = numpy.asarray(ra , dtype=float)
        dec = numpy.asarray(dec, dtype=float)

        galcoord = astropy.coordinates.SkyCoord(
            ra    = ra  * astropy.units.radian,
            dec   = dec * astropy.units.radian,
            frame = "icrs"
        ).galactic

        gall = galcoord.l.degree
        galb = galcoord.b.degree

        ebv = numpy.full(shape=gall.shape, fill_value=numpy.nan, dtype=numpy.float32)

        if numpy.any(galb >= 0):
            index = (galb >= 0)
            if self.__dust_ngp is None:
                self.__dust_ngp = _DustMap("SFD_dust_4096_ngp.fits")
            ebv[index] = self.__dust_ngp.getvalue(
                gall[index], galb[index]
            )

        if numpy.any(galb < 0):
            index = (galb < 0)
            if self.__dust_sgp is None:
                self.__dust_sgp = _DustMap("SFD_dust_4096_sgp.fits")
            ebv[index] = self.__dust_sgp.getvalue(
                gall[index], galb[index]
            )

        return ebv


class _DustMap(object):
    """
    Hemisphere map of E(B-V)
    """
    def __init__(self, name):
        """
        @param name (str): FITS file name relative to "maps/"
        """
        path = os.path.join(os.path.dirname(__file__), "maps", name)
        self.__hdu = astropy.io.fits.open(path)[0]
        self.__wcs = astropy.wcs.WCS(self.__hdu.header)

    def getvalue(self, gall, galb):
        """
        Get E(B-V) at (gall, galb).
        @param gall (numpy.array) galactic coordinate "l" in degrees
        @param galb (numpy.array) galactic coordinate "b" in degrees
        """
        gall = numpy.asarray(gall, dtype=float)
        galb = numpy.asarray(galb, dtype=float)

        height, width = self.__hdu.data.shape
        x, y = self.__wcs.all_world2pix(gall, galb, 0.0)

        x0 = numpy.floor(x)
        y0 = numpy.floor(y)

        wx1 = (x - x0).astype(numpy.float32) # Weight of x1
        wy1 = (y - y0).astype(numpy.float32) # Weight of y1

        index = (x0 < 0.0)
        x0 [index] = 0.0
        wx1[index] = 0.0
        index = (y0 < 0.0)
        y0 [index] = 0.0
        wy1[index] = 0.0

        index = (x0 > width - 2.0)
        x0 [index] = width - 2.0
        wx1[index] = 1.0
        index = (y0 > height - 2.0)
        y0 [index] = height - 2.0
        wy1[index] = 1.0

        x0 = x0.astype(int)
        y0 = y0.astype(int)
        x1 = x0 + 1
        y1 = y0 + 1

        v00 = self.__hdu.data[y0, x0]
        v01 = self.__hdu.data[y0, x1]
        v10 = self.__hdu.data[y1, x0]
        v11 = self.__hdu.data[y1, x1]

        wx0 = numpy.float32(1.0) - wx1
        wy0 = numpy.float32(1.0) - wy1

        if(numpy.any((v00 == 0) & (wy0*wx0 != 0))
        or numpy.any((v01 == 0) & (wy0*wx1 != 0))
        or numpy.any((v10 == 0) & (wy1*wx0 != 0))
        or numpy.any((v11 == 0) & (wy1*wx1 != 0))
        ):
            raise RuntimeError("No data at the given coordinates.")

        return wy0*(v00*wx0 + v01*wx1) + wy1*(v10*wx0 + v11*wx1)
