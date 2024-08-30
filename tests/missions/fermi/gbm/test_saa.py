# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#

import os
import copy
import numpy as np
import unittest

from astropy.time import TimeDelta
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.saa import *


class TestGbmSaa(unittest.TestCase):
    
    def setUp(self):
        # list of times to test GbmSaaCollection.at method
        self.t_polygons = [
            148900505, 297842765, 299086681, 300808412, 488816193,
            676400255, 678869577, 681588310, 701846688, 721849115,
            732696665, 744059705, 744716255, 1145534405
        ]


    def test_polygon1(self):
        saa = GbmSaaPolygon1()
        vals = [-30.000, -22.600, 2.500, 5.200, 5.200, 4.600,
                0.700, -8.600, -9.900, -12.500, -21.700, -30.000, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [33.900,  24.500, -18.600, -25.700, -36.000, -42.000,
                -58.800, -93.100, -97.500, -98.500, -92.100, -86.100,  33.900]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 13)

    def test_polygon2(self):
        saa = GbmSaaPolygon2()
        vals = [-30.000, -21.670, -13.330, -5.000, -6.000, -7.000,
                -9.000, -11.000, -14.000, -17.000, -23.500, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [21.000, 0.000, -21.000, -42.000, -50.000, -58.000,
                -65.000, -72.000, -76.000, -80.000, -80.000, -80.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 12)

    def test_polygon3(self):
        saa = GbmSaaPolygon3()
        vals = [-30.000, -14.500, -9.300, 1.000, 2.000, 3.000,
                -2.000, -7.000, -11.000, -15.000, -22.500, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [26.000, -2.000, -11.300, -30.000, -37.500, -45.000,
                -62.500, -80.000, -85.000, -90.000, -88.500, -87.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 12)

    def test_polygon4(self):
        saa = GbmSaaPolygon4()
        vals = [-30.000, -2.000, -0.500, 0.500, 1.000, 1.000,
                0.000, -2.000, -7.500, -11.500, -16.500, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [40.000, -27.000, -35.000, -47.000, -60.000, -65.000,
                -75.000, -82.000, -90.000, -94.000, -97.000, -92.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 12)

    def test_polygon5(self):
        saa = GbmSaaPolygon5()
        vals = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                 -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
                  -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 13)

    def test_polygon6(self):
        saa = GbmSaaPolygon6()
        vals = [-30.000, -19.867, -9.733, -3.400, -2.000, -2.000,
                -5.000, -10.155, -12.880, -16.220, -22.404, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [-5.000, -16.398, -28.103, -35.605, -38.400, -45.000,
                -65.000, -84.000, -89.200, -90.300, -90.300, -86.100]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 12)

    def test_polygon8(self):
        saa = GbmSaaPolygon8()
        vals = [-30.000, -19.867, -9.733, -3.400, 0.000, 0.000,
                -3.000, -8.155, -10.880, -16.220, -22.404, -30.000, -30.000]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [0.000, -11.398, -23.103, -30.605, -36.400, -43.000,
                -65.000, -84.000, -89.200, -90.300, -90.300, -86.100, 0.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 13)

    def test_polygon10(self):
        saa = GbmSaaPolygon10()
        vals = [-24.395, -30.000, -30.000, -30.000, -30.000, -24.060,
                -16.220, -8.638, -6.155, -1.000, 2.000, 2.000,
                -3.400, -19.570802973653127, -24.395]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [22.000, 22.000, 0.000, -2.000, -86.100, -90.300,
                -90.300, -88.738, -84.000, -65.000, -45.000, -38.400,
                -30.605, -11.999457706441582, 22.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 15)

    def test_polygon13(self):
        saa = GbmSaaPolygon13()
        vals = [-24.395, -30.000, -30.000, -30.000, -30.000, -30.000,
                -30.000, -22.646, -16.220, -9.404, -6.155, -1.000,
                2.000, 2.000, -3.650, -20.136067618721196, -24.395]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [22.000, 22.000, 5.500, 0.000, -2.000, -40.000,
                -86.100, -91.300, -91.800, -89.700, -84.000, -65.000,
                -45.000, -38.400, -30.605, -8.015646247668737, 22.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 17)

    def test_polygon14(self):
        saa = GbmSaaPolygon14()
        vals = [-24.395, -30.000, -30.000, -30.000, -30.000, -30.000,
                -30.000, -22.646, -16.220, -9.404, -6.155, -1.000,
                2.000, 2.000, -3.650, -17.867934687794072, -24.395]
        self.assertListEqual(saa.latitude.tolist(), vals)
        vals = [22.000, 22.000, 5.500, 0.000, -2.340, -40.000,
                -86.100, -91.300, -91.800, -89.700, -84.000, -65.000,
                -45.000, -38.400, -30.605, -11.123461787369832, 22.000]
        self.assertListEqual(saa.longitude.tolist(), vals)
        self.assertEqual(saa.num_points, 17)

    def test_collection(self):

        saa_collection = GbmSaaCollection()

        for t, poly in zip(self.t_polygons, saa_collection.polygons):
            t_fermi = Time(t, format='fermi')
            self.assertListEqual(saa_collection.at(t_fermi).latitude.tolist(), poly.latitude.tolist())
            self.assertListEqual(saa_collection.at(t_fermi).longitude.tolist(), poly.longitude.tolist())

        # ensure we fail at out-of-bounds time
        bad_time = Time(-1, format='fermi')
        with self.assertRaises(ValueError):
            saa_collection.at(bad_time)

    def test_collection_sort(self):

        saa_collection = GbmSaaCollection()

        # create copies of existing time ranges before modification
        time_range0 = copy.deepcopy(saa_collection.polygons[0]._time_range)
        time_range1 = copy.deepcopy(saa_collection.polygons[1]._time_range)

        # check that polygons are sorted
        for i in range(saa_collection.num_polygons - 1):
            this_poly = saa_collection.polygons[i]
            next_poly = saa_collection.polygons[i + 1]
            self.assertTrue(next_poly._time_range._low > this_poly._time_range._low)
            self.assertTrue(next_poly._time_range._high > this_poly._time_range._high)

        # ensure we fail on duplicate time range
        saa_collection.polygons[0]._time_range = time_range1
        with self.assertRaises(ValueError):
            saa_collection.sort()

        # ensure we fail on overlapping low range
        saa_collection.polygons[0]._time_range = time_range0
        saa_collection.polygons[1]._time_range._low = time_range0._low
        with self.assertRaises(ValueError):
            saa_collection.sort()

        # ensure we fail on overlapping high range
        saa_collection.polygons[0]._time_range._high = time_range1._high
        saa_collection.polygons[1]._time_range = time_range1
        with self.assertRaises(ValueError):
            saa_collection.sort()

        # ensure we fail on overlapping mid range
        delta = 0.5 * (time_range1._high.value - time_range1._low.value)
        saa_collection.polygons[0]._time_range._high = time_range1._low + TimeDelta(delta, format='sec')
        saa_collection.polygons[1]._time_range = time_range1
        with self.assertRaises(ValueError):
            saa_collection.sort()

    def test_gbmsaa(self):
        self.assertListEqual(GbmSaa().latitude.tolist(), GbmSaaPolygon5().latitude.tolist())
        self.assertListEqual(GbmSaa().longitude.tolist(), GbmSaaPolygon5().longitude.tolist())
