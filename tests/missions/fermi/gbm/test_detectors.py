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
import numpy as np
import unittest

from gdt.missions.fermi.gbm.detectors import *


class TestGbmDetectors(unittest.TestCase):
    
    def test_n0(self):
        det = GbmDetectors.n0
        self.assertEqual(det.full_name, 'NAI_00')
        self.assertEqual(det.number, 0)
        self.assertEqual(det.azimuth.value, 45.89)
        self.assertEqual(det.zenith.value, 20.58)

    def test_n1(self):
        det = GbmDetectors.n1
        self.assertEqual(det.full_name, 'NAI_01')
        self.assertEqual(det.number, 1)
        self.assertEqual(det.azimuth.value, 45.11)
        self.assertEqual(det.zenith.value, 45.31)

    def test_n2(self):
        det = GbmDetectors.n2
        self.assertEqual(det.full_name, 'NAI_02')
        self.assertEqual(det.number, 2)
        self.assertEqual(det.azimuth.value, 58.44)
        self.assertEqual(det.zenith.value, 90.21)

    def test_n3(self):
        det = GbmDetectors.n3
        self.assertEqual(det.full_name, 'NAI_03')
        self.assertEqual(det.number, 3)
        self.assertEqual(det.azimuth.value, 314.87)
        self.assertEqual(det.zenith.value, 45.24)

    def test_n4(self):
        det = GbmDetectors.n4
        self.assertEqual(det.full_name, 'NAI_04')
        self.assertEqual(det.number, 4)
        self.assertEqual(det.azimuth.value, 303.15)
        self.assertEqual(det.zenith.value, 90.27)

    def test_n5(self):
        det = GbmDetectors.n5
        self.assertEqual(det.full_name, 'NAI_05')
        self.assertEqual(det.number, 5)
        self.assertEqual(det.azimuth.value, 3.35)
        self.assertEqual(det.zenith.value, 89.79)

    def test_n6(self):
        det = GbmDetectors.n6
        self.assertEqual(det.full_name, 'NAI_06')
        self.assertEqual(det.number, 6)
        self.assertEqual(det.azimuth.value, 224.93)
        self.assertEqual(det.zenith.value, 20.43)

    def test_n7(self):
        det = GbmDetectors.n7
        self.assertEqual(det.full_name, 'NAI_07')
        self.assertEqual(det.number, 7)
        self.assertEqual(det.azimuth.value, 224.62)
        self.assertEqual(det.zenith.value, 46.18)

    def test_n8(self):
        det = GbmDetectors.n8
        self.assertEqual(det.full_name, 'NAI_08')
        self.assertEqual(det.number, 8)
        self.assertEqual(det.azimuth.value, 236.61)
        self.assertEqual(det.zenith.value, 89.97)

    def test_n9(self):
        det = GbmDetectors.n9
        self.assertEqual(det.full_name, 'NAI_09')
        self.assertEqual(det.number, 9)
        self.assertEqual(det.azimuth.value, 135.19)
        self.assertEqual(det.zenith.value, 45.55)

    def test_na(self):
        det = GbmDetectors.na
        self.assertEqual(det.full_name, 'NAI_10')
        self.assertEqual(det.number, 10)
        self.assertEqual(det.azimuth.value, 123.73)
        self.assertEqual(det.zenith.value, 90.42)

    def test_nb(self):
        det = GbmDetectors.nb
        self.assertEqual(det.full_name, 'NAI_11')
        self.assertEqual(det.number, 11)
        self.assertEqual(det.azimuth.value, 183.74)
        self.assertEqual(det.zenith.value, 90.32)

    def test_b0(self):
        det = GbmDetectors.b0
        self.assertEqual(det.full_name, 'BGO_00')
        self.assertEqual(det.number, 12)
        self.assertEqual(det.azimuth.value, 0.00)
        self.assertEqual(det.zenith.value, 90.00)

    def test_b1(self):
        det = GbmDetectors.b1
        self.assertEqual(det.full_name, 'BGO_01')
        self.assertEqual(det.number, 13)
        self.assertEqual(det.azimuth.value, 180.00)
        self.assertEqual(det.zenith.value, 90.00)

    def test_bgo(self):
        bgo = [det.name for det in GbmDetectors.bgo()]
        self.assertListEqual(bgo, ['b0', 'b1'])

    def test_nai(self):
        nai = [det.name for det in GbmDetectors.nai()]
        self.assertListEqual(nai, ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6',
                                   'n7', 'n8', 'n9', 'na', 'nb'])

    def test_is_bgo(self):
        self.assertTrue(GbmDetectors.b0.is_bgo())
        self.assertFalse(GbmDetectors.n0.is_bgo())
        
    def test_is_nai(self):
        self.assertTrue(GbmDetectors.na.is_nai())
        self.assertFalse(GbmDetectors.b1.is_nai())
        
