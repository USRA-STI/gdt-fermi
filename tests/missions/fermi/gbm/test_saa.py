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
        self.t_polygon1 = Time(743454904, format='fermi')
        self.t_polygon2 = Time(743454906, format='fermi')

    def test_polygon1_latitude(self):
        vals = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                 -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
        saa = GbmSaaPolygon1()
        self.assertListEqual(saa.latitude.tolist(), vals)

    def test_polygon1_longitude(self):
        vals = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
                  -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]
        saa = GbmSaaPolygon1()
        self.assertListEqual(saa.longitude.tolist(), vals)
    
    def test_polygon1_num_points(self):
        saa = GbmSaaPolygon1()
        self.assertEqual(saa.num_points, 13)

    def test_polygon2_latitude(self):
        vals = [-24.395, -30.000, -30.000, -30.000, -30.000, -24.060,
                -16.220, -8.638, -6.155, -1.000, 2.000, 2.000,
                -3.400, -19.570802973653127, -24.395]
        saa = GbmSaaPolygon2()
        self.assertListEqual(saa.latitude.tolist(), vals)

    def test_polygon2_longitude(self):
        vals = [22.000, 22.000, 0.000, -2.000, -86.100, -90.300,
                -90.300, -88.738, -84.000, -65.000, -45.000, -38.400,
                -30.605, -11.999457706441582, 22.000]
        saa = GbmSaaPolygon2()
        self.assertListEqual(saa.longitude.tolist(), vals)
    
    def test_polygon2_num_points(self):
        saa = GbmSaaPolygon2()
        self.assertEqual(saa.num_points, 15)


    def test_collection(self):

        saa_collection = GbmSaaCollection()

        for t, poly in [(self.t_polygon1, GbmSaaPolygon1()), 
                        (self.t_polygon2, GbmSaaPolygon2())]:
            self.assertListEqual(saa_collection.at(t).latitude.tolist(), poly.latitude.tolist())
            self.assertListEqual(saa_collection.at(t).longitude.tolist(), poly.longitude.tolist())

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
        self.assertListEqual(GbmSaa().latitude.tolist(), GbmSaaPolygon1().latitude.tolist())
        self.assertListEqual(GbmSaa().longitude.tolist(), GbmSaaPolygon1().longitude.tolist())
