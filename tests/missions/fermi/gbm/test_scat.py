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
from tempfile import TemporaryDirectory
import astropy.io.fits as fits
from gdt.core import data_path
from gdt.missions.fermi.gbm.scat import *

scat_file = data_path / 'fermi-gbm/glg_scat_all_bn170817529_flnc_comp_v00.fit'


@unittest.skipIf(not scat_file.exists(), "test files aren't downloaded. run 'gdt-data download fermi-gbm'.")
class TestGbmModelFit(unittest.TestCase):

    def setUp(self):
        scat = Scat.open(scat_file)
        self.fit = scat.model_fits[0]

    def test_duration_fluence(self):
        self.assertAlmostEqual(self.fit.duration_fluence.value, 0.46573272, places=6)
        self.assertAlmostEqual(self.fit.duration_fluence.uncertainty[0], 0.05764648, places=6)

    def test_energy_fluence_50_300(self):
        self.assertAlmostEqual(self.fit.energy_fluence_50_300.value, 9.85259234e-08, places=6)
        self.assertAlmostEqual(self.fit.energy_fluence_50_300.uncertainty[0], 1.27686235e-08, places=6)

    def test_photon_flux_50_300(self):
        self.assertAlmostEqual(self.fit.photon_flux_50_300.value, 1.8259016, places=6)
        self.assertAlmostEqual(self.fit.photon_flux_50_300.uncertainty[0], 0.2260026, places=6)

    def test_covariance(self):
        vals = [3.8231263e-04, -9.7585416e-01, 1.0429098e-02, 0.0000000e+00,
                -9.7585416e-01, 2.9397126e+03, -2.3200571e+01, 0.0000000e+00,
                1.0429098e-02, -2.3200573e+01, 3.4899354e-01, 0.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00]
        covar = self.fit.covariance.flatten()
        for i in range(covar.size):
            self.assertAlmostEqual(covar[i], vals[i], places=4)

    def test_dof(self):
        self.assertEqual(self.fit.dof, 479)

    def test_energy_fluence(self):
        self.assertAlmostEqual(self.fit.energy_fluence.value, 1.4034558e-07, places=6)
        self.assertAlmostEqual(self.fit.energy_fluence.uncertainty[0], 3.0259528e-08, places=6)

    def test_energy_flux(self):
        self.assertAlmostEqual(self.fit.energy_flux.value, 5.5022377e-07, places=6)
        self.assertAlmostEqual(self.fit.energy_flux.uncertainty[0], 1.1863225e-07, places=6)

    def test_flux_energy_range(self):
        self.assertTupleEqual(self.fit.flux_energy_range, (10.0, 1000.0))

    def test_name(self):
        self.assertEqual(self.fit.name, 'Comptonized, Epeak')

    def test_parameters(self):
        names = ['Amplitude', 'Epeak', 'Index', 'Pivot E =fix']
        vals = [0.03199916, 215.09435, 0.14382371, 100.]
        errs = [0.01955282, 54.219116, 0.5907568, 0.0]

        params = self.fit.parameters
        for i in range(4):
            self.assertEqual(params[i].name, names[i])
            self.assertAlmostEqual(params[i].value, vals[i], places=5)
            self.assertAlmostEqual(params[i].uncertainty[0], errs[i], places=5)

    def test_photon_fluence(self):
        self.assertAlmostEqual(self.fit.photon_fluence.value, 0.71739393, places=6)
        self.assertAlmostEqual(self.fit.photon_fluence.uncertainty[0], 0.11288685, places=6)

    def test_photon_flux(self):
        self.assertAlmostEqual(self.fit.photon_flux.value, 2.8125374, places=6)
        self.assertAlmostEqual(self.fit.photon_flux.uncertainty[0], 0.442572037, places=6)

    def test_stat_name(self):
        self.assertEqual(self.fit.stat_name, 'Castor C-STAT')

    def test_stat_value(self):
        self.assertAlmostEqual(self.fit.stat_value, 516.5009, places=4)

    def test_time_range(self):
        self.assertAlmostEqual(self.fit.time_range[0], -0.192, places=3)
        self.assertAlmostEqual(self.fit.time_range[1], 0.064, places=3)


@unittest.skipIf(not scat_file.exists(), "test files aren't downloaded. run 'gdt-data download fermi-gbm'.")
class TestGbmDetectorData(unittest.TestCase):

    def setUp(self):
        scat = Scat.open(scat_file)
        self.det = scat.detectors[0]
        with fits.open(scat_file) as f:
            self.data = f[1].data

    def test_active(self):
        self.assertTrue(self.det.active)

    def test_channel_mask(self):
        mask = np.zeros(128, dtype=bool)
        mask[3:125] = True
        self.assertListEqual(self.det.channel_mask.tolist(), mask.tolist())

    def test_channel_range(self):
        self.assertTupleEqual(self.det.channel_range, (3, 124))

    def test_datatype(self):
        self.assertEqual(self.det.datatype, 'TTE')

    def test_detector(self):
        self.assertEqual(self.det.detector, 'BGO_00')

    def test_energy_edges(self):
        self.assertListEqual(self.det.energy_edges.tolist(),
                             self.data['E_EDGES'][0].tolist())

    def test_energy_range(self):
        self.assertAlmostEqual(self.det.energy_range[0], 284.65, places=2)
        self.assertAlmostEqual(self.det.energy_range[1], 40108., places=1)

    def test_filename(self):
        self.assertEqual(self.det.filename, 'glg_tte_b0_bn170817529_v00.fit')

    def test_instrument(self):
        self.assertEqual(self.det.instrument, 'GBM')

    def test_num_chans(self):
        self.assertEqual(self.det.num_chans, 128)

    def test_photon_counts(self):
        self.assertListEqual(self.det.photon_counts.tolist(),
                             self.data['PHTCNTS'][0].tolist())

    def test_photon_errors(self):
        self.assertListEqual(self.det.photon_errors.tolist(),
                             self.data['PHTERRS'][0].tolist())

    def test_photon_model(self):
        self.assertListEqual(self.det.photon_model.tolist(),
                             self.data['PHTMODL'][0].tolist())

    def test_response(self):
        self.assertEqual(self.det.response, 'glg_cspec_b0_bn170817529_v04.rsp')

    def test_time_range(self):
        self.assertAlmostEqual(self.det.time_range[0], -0.192, places=3)
        self.assertAlmostEqual(self.det.time_range[1], 0.064, places=3)


@unittest.skipIf(not scat_file.exists(), "test files aren't downloaded. run 'gdt-data download fermi-gbm'.")
class TestScat(unittest.TestCase):

    def setUp(self):
        self.scat = Scat.open(scat_file)

    def test_detectors(self):
        det_names = [det.detector for det in self.scat.detectors]
        self.assertListEqual(det_names, ['BGO_00', 'NAI_01', 'NAI_02', 'NAI_05'])

    def test_model_fits(self):
        model_names = [fit.name for fit in self.scat.model_fits]
        self.assertListEqual(model_names, ['Comptonized, Epeak'])

    def test_num_detectors(self):
        self.assertEqual(self.scat.num_detectors, 4)

    def test_num_fits(self):
        self.assertEqual(self.scat.num_fits, 1)

    def test_add_detector_data(self):
        self.scat.add_detector_data(self.scat.detectors[0])
        self.assertEqual(self.scat.num_detectors, 5)
        self.scat._detectors = self.scat._detectors[:-1]

    def test_add_model_fit(self):
        self.scat.add_model_fit(self.scat.model_fits[0])
        self.assertEqual(self.scat.num_fits, 2)
        self.scat._model_fits = self.scat._model_fits[:-1]

    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.scat.write(this_path, overwrite=True)
            scat = Scat.open(os.path.join(this_path, self.scat.filename))

            self.assertEqual(scat.num_detectors, self.scat.num_detectors)
            self.assertEqual(scat.num_fits, self.scat.num_fits)
            self.assertEqual(scat.model_fits[0].stat_value,
                             self.scat.model_fits[0].stat_value)
            self.assertEqual(scat.model_fits[0].stat_name,
                             self.scat.model_fits[0].stat_name)

            for i in range(4):
                self.assertEqual(scat.detectors[i].photon_counts.tolist(),
                                 self.scat.detectors[i].photon_counts.tolist())
                self.assertEqual(scat.detectors[i].photon_errors.tolist(),
                                 self.scat.detectors[i].photon_errors.tolist())
                self.assertEqual(scat.detectors[i].photon_model.tolist(),
                                 self.scat.detectors[i].photon_model.tolist())

            for i in range(3):
                self.assertEqual(scat.headers[i], self.scat.headers[i])

            scat.close()

    def test_errors(self):
        with self.assertRaises(TypeError):
            self.scat.add_detector_data(self.scat.model_fits[0])

        with self.assertRaises(TypeError):
            self.scat.add_model_fit(self.scat.detectors[0])

        with self.assertRaises(TypeError):
            Scat.from_data(self.scat.detectors[0], self.scat.model_fits)

        with self.assertRaises(TypeError):
            Scat.from_data(self.scat.detectors, self.scat.model_fits[0])

        with self.assertRaises(TypeError):
            Scat.from_data(self.scat.detectors, self.scat.model_fits, headers=1)
