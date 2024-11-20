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
from gdt.missions.fermi.gbm.tte import *
from gdt.core.binning.binned import combine_by_factor
from gdt.core.binning.unbinned import bin_by_time

tte_file = data_path / 'fermi-gbm/glg_tte_n9_bn090131090_v01.fit'


@unittest.skipIf(not tte_file.exists(), "test files aren't downloaded. run `gdt-data download fermi-gbm`.")
class TestGbmTte(unittest.TestCase):
    
    def setUp(self):
        self.tte = GbmTte.open(tte_file)
    
    def tearDown(self):
        self.tte.close()
    
    def test_detector(self):
        self.assertEqual(self.tte.detector, 'n9')
    
    def test_energy_range(self):
        self.assertAlmostEqual(self.tte.energy_range[0], 4.380218, places=6)
        self.assertAlmostEqual(self.tte.energy_range[1], 2000.0, places=6)

    def test_event_deadtime(self):
        self.assertEqual(self.tte.event_deadtime, 2.6e-6)

    def test_filename(self):
        self.assertEqual(self.tte.filename, os.path.basename(tte_file))
    
    def test_headers(self):
        self.assertEqual(self.tte.headers.num_headers, 4)
    
    def test_num_chans(self):
        self.assertEqual(self.tte.num_chans, 128)

    def test_overflow_deadtime(self):
        self.assertEqual(self.tte.overflow_deadtime, 1.0e-5)

    def test_time_range(self):
        t0, t1 = self.tte.time_range
        t0 += self.tte.trigtime
        t1 += self.tte.trigtime
        self.assertAlmostEqual(t0, 255060537.657974, places=6)
        self.assertAlmostEqual(t1, 255060863.884596, places=6)
    
    def test_trigtime(self):
        self.assertAlmostEqual(self.tte.trigtime, 255060563.149072, places=6)     
        
    def test_rebin_energy(self):
        tte2 = self.tte.rebin_energy(combine_by_factor, 2)
        self.assertEqual(tte2.num_chans, self.tte.num_chans//2)
    
    def test_to_phaii(self):
        phaii = self.tte.to_phaii(bin_by_time, 1.024)
        self.assertEqual(phaii.data.num_times, 319)
        self.assertEqual(phaii.detector, self.tte.detector)
        self.assertEqual(phaii.headers[0]['FILETYPE'], 'PHAII')
        self.assertEqual(phaii.headers[0]['DATE-OBS'], self.tte.headers[0]['DATE-OBS'])
        self.assertEqual(phaii.headers[0]['DATE-END'], self.tte.headers[0]['DATE-END'])
        self.assertEqual(phaii.headers[0]['TSTART'], self.tte.headers[0]['TSTART'])
        self.assertEqual(phaii.headers[0]['TSTOP'], self.tte.headers[0]['TSTOP'])
        self.assertEqual(phaii.headers[0]['TRIGTIME'], self.tte.headers[0]['TRIGTIME'])
        self.assertEqual(phaii.headers[0]['OBJECT'], self.tte.headers[0]['OBJECT'])
        self.assertEqual(phaii.headers[0]['RA_OBJ'], self.tte.headers[0]['RA_OBJ'])
        self.assertEqual(phaii.headers[0]['DEC_OBJ'], self.tte.headers[0]['DEC_OBJ'])
        self.assertEqual(phaii.headers[0]['ERR_RAD'], self.tte.headers[0]['ERR_RAD'])

    def test_slice_energy(self):
        tte2 = self.tte.slice_energy((50.0, 300.0))
        emin, emax = tte2.energy_range
        self.assertAlmostEqual(emin, 48.38871, places=5)
        self.assertAlmostEqual(emax, 305.1512, places=4)

    def test_slice_time(self):
        tte2 = self.tte.slice_time((10.0, 20.0))
        t0, t1 = tte2.time_range
        t0 += tte2.trigtime
        t1 += tte2.trigtime
        self.assertAlmostEqual(t0, 255060573.1, places=1)
        self.assertAlmostEqual(t1, 255060583.1, places=1)
    
    def test_to_pha(self):
        pha = self.tte.to_pha()
        self.assertEqual(self.tte.num_chans, pha.num_chans)

    def test_to_spectrum(self):
        spec = self.tte.to_spectrum()
        self.assertEqual(self.tte.num_chans, spec.size)

    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.tte.write(this_path, overwrite=True)
            tte = GbmTte.open(os.path.join(this_path, self.tte.filename))
            self.assertListEqual(tte.data.times.tolist(), self.tte.data.times.tolist())
            self.assertListEqual(tte.data.channels.tolist(), self.tte.data.channels.tolist())
            self.assertListEqual(tte.ebounds.low_edges(), self.tte.ebounds.low_edges())
            self.assertListEqual(tte.ebounds.high_edges(), self.tte.ebounds.high_edges())
            self.assertListEqual(tte.gti.low_edges(), self.tte.gti.low_edges())
            self.assertListEqual(tte.gti.high_edges(), self.tte.gti.high_edges())
            self.assertEqual(tte.trigtime, self.tte.trigtime)
            self.assertEqual(tte.detector, self.tte.detector)
            self.assertEqual(tte.headers[1], self.tte.headers[1])
            self.assertEqual(tte.headers[2], self.tte.headers[2])
            self.assertEqual(tte.headers[3], self.tte.headers[3])
            tte.close()


