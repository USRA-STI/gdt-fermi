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

from gdt.missions.fermi.gbm.saa import *


class TestGbmSaa(unittest.TestCase):
    
    def setUp(self):
        self.saa = GbmSaa()

    def test_latitude(self):
        vals = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                 -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
        self.assertListEqual(self.saa.latitude.tolist(), vals)

    def test_longitude(self):
        vals = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
                  -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]
        self.assertListEqual(self.saa.longitude.tolist(), vals)
    
    def test_num_points(self):
        self.assertEqual(self.saa.num_points, 13)
