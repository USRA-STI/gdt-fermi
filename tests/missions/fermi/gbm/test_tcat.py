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
import unittest
from tempfile import TemporaryDirectory
from gdt.core import data_path
from gdt.missions.fermi.gbm.tcat import *

tcat_file = data_path / 'fermi-gbm/glg_tcat_all_bn190222537_v01.fit'


@unittest.skipIf(not tcat_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestTcat(unittest.TestCase):
    
    def setUp(self):
        self.tcat = Tcat.open(tcat_file)
    
    def tearDown(self):
        self.tcat.close()
    
    def test_fermi_location(self):
        lon, lat = self.tcat.fermi_location
        self.assertEqual(lon.value, 74.2500)
        self.assertEqual(lat.value, -24.8500)
    
    def test_localizing_instrument(self):
        self.assertEqual(self.tcat.localizing_instrument, 'Fermi, GBM')
    
    def test_location(self):
        coord = self.tcat.location
        self.assertEqual(coord.ra.value, 147.32)
        self.assertEqual(coord.dec.value, 60.94)
    
    def test_name(self):
        self.assertEqual(self.tcat.name, 'GRB190222537')
    
    def test_time_range(self):
        self.assertEqual(self.tcat.time_range.tstart, 572532678.644472)
        self.assertEqual(self.tcat.time_range.tstop, 572533293.055776)
    
    def test_trigtime(self):
        self.assertEqual(self.tcat.trigtime, 572532812.150778)    
    
    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.tcat.write(this_path, overwrite=True)
            tcat = Tcat.open(os.path.join(this_path, self.tcat.filename))

            self.assertTupleEqual(tcat.fermi_location, self.tcat.fermi_location)
            self.assertEqual(tcat.localizing_instrument, self.tcat.localizing_instrument)
            self.assertEqual(tcat.name, self.tcat.name)
            self.assertEqual(tcat.time_range, self.tcat.time_range)
            self.assertEqual(tcat.trigtime, self.tcat.trigtime)

            tcat.close()
