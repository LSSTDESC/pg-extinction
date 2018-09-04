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

from extinction import dustval

import numpy

import unittest


class Extinction(unittest.TestCase):
    def test_get_Ebv(self):
        radec = numpy.array([
            [ra, dec] for ra in range(0, 360, 60) for dec in [-60, 0, 60]
        ], dtype=float)

        answer = numpy.array([
            0.013533305, 0.031818148, 0.9793637  ,
            0.019306011, 0.4595568  , 0.6899398  ,
            0.1923654  , 0.035891615, 0.068283126,
            0.9508033  , 0.031325743, 0.015246875,
            0.3727805  , 0.12535407 , 0.011143208,
            0.05354808 , 0.19928998 , 0.1123919  ,
        ], dtype=numpy.float32)

        myanswer = dustval.Extinction().get_Ebv(
            radec[..., 0]*(numpy.pi/180),
            radec[..., 1]*(numpy.pi/180),
        )

        diff = numpy.max(numpy.abs((myanswer - answer) / answer))
        self.assertLess(diff, 1e-4)


class TestDustMap(unittest.TestCase):
    def test_galbzero(self):
        ngp = dustval._DustMap("SFD_dust_4096_ngp.fits")
        sgp = dustval._DustMap("SFD_dust_4096_sgp.fits")

        radius = 2048
        n = int(100*(2*numpy.pi*radius))

        gall = numpy.arange(n).astype(float) * (360.0 / n)
        galb = numpy.full(shape=gall.shape, fill_value=0.0)
        # getvalue() may raise error
        # if (gall, galb) should be outside the circle of zea projection,
        # which should not happen.
        ngp.getvalue(gall, galb)
        sgp.getvalue(gall, galb)


if __name__ == "__main__":
    unittest.main()
