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

from gdt.missions.fermi.mcilwainl import *


class TestCalcMcilwainL(unittest.TestCase):
    
    def test_one_lat_lon(self):
        mcl = calc_mcilwain_l(-25.62, 78.57)
        self.assertAlmostEqual(mcl, 1.61, places=2)
        
    def test_multi_lat_lon(self):
        lats = [-20.07, 25.52, -18.42, 20.78]
        lons = [116.63, 113.87, 128.83, 147.23]
        test_vals = [1.38, 1.12, 1.33, 1.09]
        mcls = calc_mcilwain_l(lats, lons)
        for i in range(4):
            self.assertAlmostEqual(mcls[i], test_vals[i], places=2)


