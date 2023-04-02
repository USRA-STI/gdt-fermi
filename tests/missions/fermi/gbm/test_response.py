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
from gdt.missions.fermi.gbm.response import *
from gdt.core.spectra.functions import PowerLaw

rsp_file = data_path / 'fermi-gbm/glg_cspec_n9_bn090131090_v01.rsp'
rsp2_file = data_path / 'fermi-gbm/glg_cspec_n9_bn090131090_v00.rsp2'


@unittest.skipIf(not rsp_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestGbmRsp(unittest.TestCase):

    def setUp(self):
        self.rsp = GbmRsp.open(rsp_file)

    def tearDown(self):
        self.rsp.close()

    def test_detector(self):
        self.assertEqual(self.rsp.detector, 'n9')

    def test_num_chans(self):
        self.assertEqual(self.rsp.num_chans, 128)

    def test_num_ebins(self):
        self.assertEqual(self.rsp.num_ebins, 140)

    def test_tcent(self):
        self.assertAlmostEqual(self.rsp.tcent, 174.339, places=3)

    def test_trigtime(self):
        self.assertAlmostEqual(self.rsp.trigtime, 255060563.149072, places=6)

    def test_tstart(self):
        self.assertAlmostEqual(self.rsp.tstart, -132.866, places=3)

    def test_tstop(self):
        self.assertAlmostEqual(self.rsp.tstop, 481.545, places=3)

    def test_fold_spectrum(self):
        ebins = self.rsp.fold_spectrum(PowerLaw().fit_eval, (0.01, -2.2), exposure=2.0)
        self.assertEqual(ebins.size, self.rsp.num_chans)
        self.assertEqual(ebins.exposure[0], 2.0)

    def test_rebin(self):
        rsp = self.rsp.rebin(factor=2)
        self.assertEqual(rsp.num_chans, self.rsp.num_chans // 2)
        self.assertEqual(rsp.num_ebins, self.rsp.num_ebins)

    def test_resample(self):
        rsp = self.rsp.resample(num_photon_bins=70)
        self.assertEqual(rsp.num_chans, self.rsp.num_chans)
        self.assertEqual(rsp.num_ebins, self.rsp.num_ebins // 2)

    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.rsp.write(this_path, overwrite=True)
            rsp = GbmRsp.open(os.path.join(this_path, self.rsp.filename))
            self.assertEqual(rsp.detector, self.rsp.detector)
            for i in range(self.rsp.num_chans):
                self.assertListEqual(rsp.drm.matrix[:, i].tolist(),
                                     self.rsp.drm.matrix[:, i].tolist())
            self.assertListEqual(rsp.ebounds.low_edges(), self.rsp.ebounds.low_edges())
            self.assertListEqual(rsp.ebounds.high_edges(), self.rsp.ebounds.high_edges())
            self.assertEqual(rsp.trigtime, self.rsp.trigtime)
            self.assertEqual(rsp.tstart, self.rsp.tstart)
            self.assertEqual(rsp.tstop, self.rsp.tstop)
            rsp.close()


@unittest.skipIf(not rsp2_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestGbmRsp2(unittest.TestCase):

    def setUp(self):
        self.rsp2 = GbmRsp2.open(rsp2_file)

    def tearDown(self):
        self.rsp2.close()

    def test_detector(self):
        self.assertEqual(self.rsp2.detector, 'n9')

    def test_num_chans(self):
        self.assertEqual(self.rsp2.num_chans, 128)

    def test_num_drms(self):
        self.assertEqual(self.rsp2.num_drms, 8)

    def test_num_ebins(self):
        self.assertEqual(self.rsp2.num_ebins, 140)

    def test_tcent(self):
        vals = [-2.304, 13.056, 28.416, 43.777, 61.697, 80.129, 96.514, 108.802]
        for i in range(8):
            self.assertAlmostEqual(self.rsp2.tcent[i], vals[i], places=3)

    def test_trigtime(self):
        self.assertAlmostEqual(self.rsp2.trigtime, 255060563.149072, places=6)

    def test_tstart(self):
        vals = [-9.984, 5.376, 20.736, 36.097, 51.457, 71.937, 88.322, 104.706]
        for i in range(8):
            self.assertAlmostEqual(self.rsp2.tstart[i], vals[i], places=3)

    def test_tstop(self):
        vals = [5.376, 20.736, 36.097, 51.457, 71.937, 88.322, 104.706, 112.8980]
        for i in range(8):
            self.assertAlmostEqual(self.rsp2.tstop[i], vals[i], places=3)

    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.rsp2.write(this_path, overwrite=True)
            rsp2 = GbmRsp2.open(os.path.join(this_path, self.rsp2.filename))
            self.assertEqual(rsp2.detector, self.rsp2.detector)
            self.assertListEqual(rsp2.ebounds.low_edges(), self.rsp2.ebounds.low_edges())
            self.assertListEqual(rsp2.ebounds.high_edges(), self.rsp2.ebounds.high_edges())
            self.assertEqual(rsp2.trigtime, self.rsp2.trigtime)
            self.assertListEqual(rsp2.tstart.tolist(), self.rsp2.tstart.tolist())
            self.assertListEqual(rsp2.tstop.tolist(), self.rsp2.tstop.tolist())

            for i in range(self.rsp2.num_drms):
                for j in range(self.rsp2.num_chans):
                    self.assertListEqual(rsp2[i].drm.matrix[:, j].tolist(),
                                         self.rsp2[i].drm.matrix[:, j].tolist())

            rsp2.close()
