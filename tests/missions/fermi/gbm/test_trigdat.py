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
from gdt.missions.fermi.gbm.detectors import *
from gdt.missions.fermi.gbm.trigdat import *
from gdt.core.coords import SpacecraftFrame

trigdat_file = data_path / 'fermi-gbm/glg_trigdat_all_bn170101116_v01.fit'


@unittest.skipIf(not trigdat_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestTrigdat(unittest.TestCase):

    def setUp(self):
        self.trigdat = Trigdat.open(trigdat_file)

    def tearDown(self):
        self.trigdat.close()

    def test_backrates(self):
        self.assertIsInstance(self.trigdat.backrates, BackRates)

    def test_fsw_locations(self):
        self.assertIsInstance(self.trigdat.fsw_locations[0], FswLocation)

    def test_gti(self):
        gti = self.trigdat.gti
        self.assertAlmostEqual(gti.low_edges()[0], 504931505.00888, places=5)
        self.assertAlmostEqual(gti.high_edges()[0], 504932119.419874, places=5)

    def test_maxrates(self):
        self.assertIsInstance(self.trigdat.maxrates[0], MaxRates)

    def test_num_maxrates(self):
        self.assertEqual(self.trigdat.num_maxrates, 3)

    def test_poshist(self):
        poshist = self.trigdat.poshist
        self.assertIsInstance(poshist, SpacecraftFrame)
        self.assertAlmostEqual(poshist.quaternion[0].x, 0.26009563, places=6)
        self.assertAlmostEqual(poshist.quaternion[0].y, 0.02415343, places=6)
        self.assertAlmostEqual(poshist.quaternion[0].z, 0.50508887, places=6)
        self.assertAlmostEqual(poshist.quaternion[0].w, -0.8225887, places=6)
        self.assertAlmostEqual(poshist.quaternion[-1].x, 0.05106141, places=6)
        self.assertAlmostEqual(poshist.quaternion[-1].y, 0.12925173, places=6)
        self.assertAlmostEqual(poshist.quaternion[-1].z, 0.64235514, places=6)
        self.assertAlmostEqual(poshist.quaternion[-1].w, -0.7537019, places=6)

        self.assertAlmostEqual(poshist.obsgeoloc[0].x.to('km').value, -673.25, places=2)
        self.assertAlmostEqual(poshist.obsgeoloc[0].y.to('km').value, 6624.75, places=2)
        self.assertAlmostEqual(poshist.obsgeoloc[0].z.to('km').value, 1848., places=2)
        self.assertAlmostEqual(poshist.obsgeoloc[-1].x.to('km').value, -4462., places=2)
        self.assertAlmostEqual(poshist.obsgeoloc[-1].y.to('km').value, 4394.5, places=2)
        self.assertAlmostEqual(poshist.obsgeoloc[-1].z.to('km').value, 2904.75, places=2)

        self.assertAlmostEqual(poshist.obstime[0].fermi, 504931505.0088801, places=6)
        self.assertAlmostEqual(poshist.obstime[-1].fermi, 504932111.22787404, places=6)

    def test_time_range(self):
        self.assertAlmostEqual(self.trigdat.time_range.tstart,
                               504931505.00888, places=5)
        self.assertAlmostEqual(self.trigdat.time_range.tstop,
                               504932119.419874, places=5)

    def test_triggered_detectors(self):
        self.assertListEqual(self.trigdat.triggered_detectors, ['n9', 'na', 'nb'])

    def test_trigrates(self):
        self.assertIsInstance(self.trigdat.trigrates, MaxRates)

    def test_trigtime(self):
        self.assertAlmostEqual(self.trigdat.trigtime, 504931642.867272, places=6)

    def test_to_phaii(self):
        phaii = self.trigdat.to_phaii('n0')
        self.assertEqual(phaii.detector, 'n0')
        self.assertEqual(phaii.trigtime, self.trigdat.trigtime)
        self.assertEqual(np.round(phaii.data.time_widths.min(), 3), 1.024)
        self.assertEqual(np.round(phaii.data.time_widths.max(), 3), 8.192)

        phaii = self.trigdat.to_phaii('n9', timescale=256)
        self.assertEqual(phaii.detector, 'n9')
        self.assertEqual(np.round(phaii.data.time_widths.min(), 3), 0.256)
        self.assertEqual(np.round(phaii.data.time_widths.max(), 3), 8.192)

        phaii = self.trigdat.to_phaii('b0', timescale=64)
        self.assertEqual(phaii.detector, 'b0')
        self.assertEqual(np.round(phaii.data.time_widths.min(), 3), 0.064)
        self.assertEqual(np.round(phaii.data.time_widths.max(), 3), 8.192)

        # incorrect detector
        with self.assertRaises(ValueError):
            self.trigdat.to_phaii('yo mama')

        # incorrect timescale
        with self.assertRaises(ValueError):
            self.trigdat.to_phaii('n0', timescale=10)

    def test_sum_detectors(self):
        phaii = self.trigdat.sum_detectors(['n0', 'n1'])
        phaii_n0 = self.trigdat.to_phaii('n0')
        phaii_n1 = self.trigdat.to_phaii('n1')
        self.assertEqual(phaii.detector, 'n0+n1')
        self.assertEqual(phaii.trigtime, self.trigdat.trigtime)
        self.assertEqual(np.round(phaii.data.time_widths.min(), 3), 1.024)
        self.assertEqual(np.round(phaii.data.time_widths.max(), 3), 8.192)
        self.assertListEqual(phaii.data.counts.flatten().tolist(),
                             (phaii_n0.data.counts
                              + phaii_n1.data.counts).flatten().tolist())
        self.assertListEqual(phaii.data.exposure.tolist(),
                             ((phaii_n0.data.exposure
                               + phaii_n1.data.exposure) / 2.0).tolist())

        phaii = self.trigdat.sum_detectors(['n7', 'n8', 'n9'], timescale=64)
        phaii_n7 = self.trigdat.to_phaii('n7', timescale=64)
        phaii_n8 = self.trigdat.to_phaii('n8', timescale=64)
        phaii_n9 = self.trigdat.to_phaii('n9', timescale=64)
        self.assertEqual(phaii.detector, 'n7+n8+n9')
        self.assertEqual(phaii.trigtime, self.trigdat.trigtime)
        self.assertEqual(np.round(phaii.data.time_widths.min(), 3), 0.064)
        self.assertEqual(np.round(phaii.data.time_widths.max(), 3), 8.192)
        self.assertListEqual(phaii.data.counts.flatten().tolist(),
                             (phaii_n7.data.counts
                              + phaii_n8.data.counts
                              + phaii_n9.data.counts).flatten().tolist())
        self.assertListEqual(phaii.data.exposure.tolist(),
                             ((phaii_n7.data.exposure
                               + phaii_n8.data.exposure
                               + phaii_n9.data.exposure) / 3.0).tolist())

        # incorrect detector
        with self.assertRaises(ValueError):
            self.trigdat.sum_detectors(['n0', 'yo mama'])

        # incorrect timescale
        with self.assertRaises(ValueError):
            self.trigdat.sum_detectors(['n0', 'n1'], timescale=10)

    def test_write(self):
        with TemporaryDirectory() as this_path:
            self.trigdat.write(this_path, overwrite=True)
            trigdat = Trigdat.open(os.path.join(this_path, self.trigdat.filename))

            self.assertListEqual(trigdat.gti.low_edges(), self.trigdat.gti.low_edges())
            self.assertListEqual(trigdat.gti.high_edges(), self.trigdat.gti.high_edges())
            self.assertEqual(trigdat.num_maxrates, self.trigdat.num_maxrates)
            self.assertEqual(trigdat.time_range.tstart, self.trigdat.time_range.tstart)
            self.assertEqual(trigdat.time_range.tstop, self.trigdat.time_range.tstop)
            self.assertEqual(trigdat.trigtime, self.trigdat.trigtime)
            self.assertListEqual(trigdat._data['RATE'].flatten().tolist(),
                                 self.trigdat._data['RATE'].flatten().tolist())
            self.assertListEqual(trigdat._data['TIME'].tolist(),
                                 self.trigdat._data['TIME'].tolist())
            self.assertListEqual(trigdat._data['ENDTIME'].tolist(),
                                 self.trigdat._data['ENDTIME'].tolist())
            self.assertListEqual(trigdat._data['SCATTITD'].flatten().tolist(),
                                 self.trigdat._data['SCATTITD'].flatten().tolist())
            self.assertListEqual(trigdat._data['EIC'].flatten().tolist(),
                                 self.trigdat._data['EIC'].flatten().tolist())

            trigdat.close()

    def test_from_data(self):
        phaiis = [self.trigdat.to_phaii(det.name) for det in GbmDetectors]
        poshist = self.trigdat.poshist
        trigtime = self.trigdat.trigtime
        time_range = self.trigdat.time_range.as_tuple()
        trigrate = self.trigdat.trigrates
        maxrates = self.trigdat.maxrates
        backrates = self.trigdat.backrates
        fswlocs = self.trigdat.fsw_locations
        filename = self.trigdat.filename
        headers = self.trigdat.headers
        trigdat = Trigdat.from_data(phaiis, poshist, trigtime,
                                    time_range=time_range, trigrate=trigrate,
                                    maxrates=maxrates, backrates=backrates,
                                    fswlocations=fswlocs, filename=filename,
                                    headers=headers)

        self.assertListEqual(trigdat.gti.low_edges(), self.trigdat.gti.low_edges())
        self.assertListEqual(trigdat.gti.high_edges(), self.trigdat.gti.high_edges())
        self.assertEqual(trigdat.num_maxrates, self.trigdat.num_maxrates)
        self.assertEqual(trigdat.time_range.tstart, self.trigdat.time_range.tstart)
        self.assertEqual(trigdat.time_range.tstop, self.trigdat.time_range.tstop)
        self.assertEqual(trigdat.trigtime, self.trigdat.trigtime)

        numtimes = 134
        self.assertTupleEqual(trigdat._data['TIME'].shape, (numtimes,))
        self.assertTupleEqual(trigdat._data['ENDTIME'].shape, (numtimes,))
        self.assertTupleEqual(trigdat._data['SCATTITD'].shape, (numtimes, 4))
        self.assertTupleEqual(trigdat._data['EIC'].shape, (numtimes, 3))
        self.assertTupleEqual(trigdat._data['RATE'].shape, (numtimes, 8, 14))

        # wrong phaii
        with self.assertRaises(TypeError):
            Trigdat.from_data(maxrates, poshist, trigtime)

        # wrong poshist
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, trigtime, trigtime)

        # wrong value for trigtime
        with self.assertRaises(ValueError):
            Trigdat.from_data(phaiis, poshist, 100.0)

        # wrong trigrate
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, trigrate=1)

        # wrong maxratesrates
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, maxrates=1)
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, maxrates=maxrates[0])

        # wrong backrates
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, backrates=1)

        # wrong fsw locations
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, fswlocations=1)
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, fswlocations=fswlocs[0])

        # wrong headers
        with self.assertRaises(TypeError):
            Trigdat.from_data(phaiis, poshist, trigtime, headers=1)


@unittest.skipIf(not trigdat_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestMaxRates(unittest.TestCase):

    def setUp(self):
        with Trigdat.open(trigdat_file) as trigdat:
            self.maxrates = trigdat.trigrates

    def test_all_rates(self):
        with fits.open(trigdat_file) as f:
            self.assertListEqual(self.maxrates.all_rates.flatten().tolist(),
                                 f[1].data['TRIGRATE'].flatten().tolist())

    def test_eic(self):
        with fits.open(trigdat_file) as f:
            self.assertListEqual(self.maxrates.eic.flatten().tolist(),
                                 f[1].data['EIC'].flatten().tolist())

    def test_num_chans(self):
        self.assertEqual(self.maxrates.num_chans, 8)

    def test_num_dets(self):
        self.assertEqual(self.maxrates.num_dets, 14)

    def test_quaternion(self):
        with fits.open(trigdat_file) as f:
            self.assertListEqual(self.maxrates.quaternion.flatten().tolist(),
                                 f[1].data['SCATTITD'].flatten().tolist())

    def test_time_range(self):
        tstart, tstop = self.maxrates.time_range
        self.assertAlmostEqual(tstart, 504931642.739272)
        self.assertAlmostEqual(tstop, 504931642.867272)

    def test_timescale(self):
        self.assertEqual(self.maxrates.timescale, 128)

    def test_get_detector(self):
        rates = self.maxrates.get_detector('n1')
        vals = [344., 24., 48., 248., 384., 24., 56., 192.]
        self.assertListEqual(rates.tolist(), vals)

    def test_create(self):
        # empty
        maxrates = MaxRates()
        self.assertIsNone(maxrates.all_rates)
        self.assertIsNone(maxrates.eic)
        self.assertEqual(maxrates.num_chans, 0)
        self.assertEqual(maxrates.num_dets, 0)
        self.assertIsNone(maxrates.quaternion)
        self.assertTupleEqual(maxrates.time_range, (None, None))
        self.assertIsNone(maxrates.timescale)
        self.assertIsNone(maxrates.get_detector('n0'))

        tstart, tstop = self.maxrates.time_range
        maxrates = MaxRates.create(tstart=tstart, tstop=tstop,
                                   eic=self.maxrates.eic,
                                   quaternion=self.maxrates.quaternion,
                                   rates=self.maxrates.all_rates)
        self.assertListEqual(maxrates.all_rates.flatten().tolist(),
                             self.maxrates.all_rates.flatten().tolist())
        self.assertListEqual(maxrates.eic.flatten().tolist(),
                             self.maxrates.eic.flatten().tolist())
        self.assertListEqual(maxrates.quaternion.flatten().tolist(),
                             self.maxrates.quaternion.flatten().tolist())
        self.assertTupleEqual(maxrates.time_range, self.maxrates.time_range)
        self.assertEqual(maxrates.num_chans, self.maxrates.num_chans)
        self.assertEqual(maxrates.num_dets, self.maxrates.num_dets)
        self.assertEqual(maxrates.timescale, self.maxrates.timescale)

        # wrong quaternion
        with self.assertRaises(ValueError):
            MaxRates.create(quaternion=self.maxrates.quaternion[1:])

        # wrong eic
        with self.assertRaises(ValueError):
            MaxRates.create(eic=self.maxrates.eic[1:])

    def test_write(self):
        with TemporaryDirectory() as this_path:
            with Trigdat.open(trigdat_file) as trigdat:
                trigdat.write(this_path, overwrite=True)
                filename = trigdat.filename

            trigdat = Trigdat.open(os.path.join(this_path, filename))
            maxrates = trigdat.trigrates

            self.assertListEqual(maxrates.all_rates.flatten().tolist(),
                                 self.maxrates.all_rates.flatten().tolist())
            self.assertListEqual(maxrates.eic.flatten().tolist(),
                                 self.maxrates.eic.flatten().tolist())
            self.assertListEqual(maxrates.quaternion.flatten().tolist(),
                                 self.maxrates.quaternion.flatten().tolist())
            self.assertTupleEqual(maxrates.time_range, self.maxrates.time_range)
            self.assertEqual(maxrates.num_chans, self.maxrates.num_chans)
            self.assertEqual(maxrates.num_dets, self.maxrates.num_dets)
            self.assertEqual(maxrates.timescale, self.maxrates.timescale)

            trigdat.close()


@unittest.skipIf(not trigdat_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestBackRates(unittest.TestCase):

    def setUp(self):
        with Trigdat.open(trigdat_file) as trigdat:
            self.backrates = trigdat.backrates

    def test_all_rates(self):
        with fits.open(trigdat_file) as f:
            self.assertListEqual(self.backrates.all_rates.flatten().tolist(),
                                 f[2].data['BCKRATES'].flatten().tolist())

    def test_num_chans(self):
        self.assertEqual(self.backrates.num_chans, 8)

    def test_num_dets(self):
        self.assertEqual(self.backrates.num_dets, 14)

    def test_quality(self):
        self.assertListEqual(self.backrates.quality.tolist(), [1, 0])

    def test_time_range(self):
        tstart, tstop = self.backrates.time_range
        self.assertAlmostEqual(tstart, 504931609.913434)
        self.assertAlmostEqual(tstop, 504931642.68143404)

    def test_get_detector(self):
        rates = self.backrates.get_detector('n1')
        vals = [276., 54., 42., 138., 387., 33., 31., 240.]
        self.assertListEqual(rates.tolist(), vals)

    def test_create(self):
        # empty
        backrates = BackRates()
        self.assertIsNone(backrates.all_rates)
        self.assertEqual(backrates.num_chans, 0)
        self.assertEqual(backrates.num_dets, 0)
        self.assertListEqual(backrates.quality.tolist(), [0, 0])
        self.assertTupleEqual(backrates.time_range, (None, None))
        self.assertIsNone(backrates.get_detector('n0'))

        tstart, tstop = self.backrates.time_range
        backrates = BackRates.create(tstart=tstart, tstop=tstop,
                                     quality=self.backrates.quality.tolist(),
                                     rates=self.backrates.all_rates)
        self.assertListEqual(backrates.all_rates.flatten().tolist(),
                             self.backrates.all_rates.flatten().tolist())
        self.assertListEqual(backrates.quality.flatten().tolist(),
                             self.backrates.quality.flatten().tolist())
        self.assertTupleEqual(backrates.time_range, self.backrates.time_range)
        self.assertEqual(backrates.num_chans, self.backrates.num_chans)
        self.assertEqual(backrates.num_dets, self.backrates.num_dets)

        # wrong quality
        with self.assertRaises(ValueError):
            BackRates.create(quality=self.backrates.quality[1:].tolist())
        with self.assertRaises(TypeError):
            BackRates.create(quality=1)

    def test_write(self):
        with TemporaryDirectory() as this_path:
            with Trigdat.open(trigdat_file) as trigdat:
                trigdat.write(this_path, overwrite=True)
                filename = trigdat.filename

            trigdat = Trigdat.open(os.path.join(this_path, filename))
            backrates = trigdat.backrates

            self.assertListEqual(backrates.all_rates.flatten().tolist(),
                                 self.backrates.all_rates.flatten().tolist())
            self.assertListEqual(backrates.quality.flatten().tolist(),
                                 self.backrates.quality.flatten().tolist())
            self.assertTupleEqual(backrates.time_range, self.backrates.time_range)
            self.assertEqual(backrates.num_chans, self.backrates.num_chans)
            self.assertEqual(backrates.num_dets, self.backrates.num_dets)

            trigdat.close()


@unittest.skipIf(not trigdat_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestFswLocation(unittest.TestCase):

    def setUp(self):
        with Trigdat.open(trigdat_file) as trigdat:
            self.fswloc = trigdat.fsw_locations[0]

    def test_fluence(self):
        with fits.open(trigdat_file) as f:
            self.assertAlmostEqual(self.fswloc.fluence,
                                   f[3].data['FLUENCE'][0], places=6)

    def test_hardness_ratio(self):
        with fits.open(trigdat_file) as f:
            self.assertAlmostEqual(self.fswloc.hardness_ratio,
                                   f[3].data['HDRATIO'][0], places=6)

    def test_intensity(self):
        with fits.open(trigdat_file) as f:
            self.assertAlmostEqual(self.fswloc.intensity,
                                   f[3].data['INTNSITY'][0], places=6)

    def test_location(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.location[0], f[3].data['RA'][0])
            self.assertEqual(self.fswloc.location[1], f[3].data['DEC'][0])
            self.assertEqual(self.fswloc.location[2], f[3].data['STATERR'][0])

    def test_location_sc(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.location_sc[0], f[3].data['TR_SCAZ'][0])
            self.assertEqual(self.fswloc.location_sc[1], f[3].data['TR_SCZEN'][0])

    def test_next_classification(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.next_classification[0], 'GROJ422')
            self.assertEqual(self.fswloc.next_classification[1],
                             f[3].data['RELIABLT'][0, 1])

    def test_rates(self):
        with fits.open(trigdat_file) as f:
            self.assertListEqual(self.fswloc.rates.flatten().tolist(),
                                 f[3].data['LOCRATES'][0, :].flatten().tolist())

    def test_significance(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.significance, f[3].data['SIGMA'][0])

    def test_spectrum(self):
        self.assertEqual(self.fswloc.spectrum, 'soft')

    def test_time(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.time, f[3].data['TIME'][0])

    def test_timescale(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.timescale, f[3].data['TRIG_TS'][0])

    def test_top_classification(self):
        with fits.open(trigdat_file) as f:
            self.assertEqual(self.fswloc.top_classification[0], 'GRB')
            self.assertEqual(self.fswloc.top_classification[1],
                             f[3].data['RELIABLT'][0, 0])

    def test_create(self):
        # empty
        fswloc = FswLocation()
        self.assertIsNone(fswloc.fluence)
        self.assertIsNone(fswloc.hardness_ratio)
        self.assertIsNone(fswloc.intensity)
        self.assertTupleEqual(fswloc.location, (None, None, None))
        self.assertTupleEqual(fswloc.location_sc, (None, None))
        self.assertTupleEqual(fswloc.top_classification, (20.0, 1.0))
        self.assertTupleEqual(fswloc.next_classification, (20.0, 0.0))
        self.assertIsNone(fswloc.rates)
        self.assertIsNone(fswloc.significance)
        self.assertIsNone(fswloc.spectrum)
        self.assertIsNone(fswloc.time)
        self.assertIsNone(fswloc.timescale)

        fswloc = FswLocation.create(time=self.fswloc.time,
                                    radec=self.fswloc.location[:2],
                                    err=self.fswloc.location[2],
                                    algorithm=self.fswloc.spectrum,
                                    class1=self.fswloc.top_classification[0],
                                    class1_prob=self.fswloc.top_classification[1],
                                    class2=self.fswloc.next_classification[0],
                                    class2_prob=self.fswloc.next_classification[1],
                                    intensity=self.fswloc.intensity,
                                    hardness=self.fswloc.hardness_ratio,
                                    fluence=self.fswloc.fluence,
                                    sigma=self.fswloc.significance,
                                    rates=self.fswloc.rates,
                                    timescale=self.fswloc.timescale,
                                    azzen=self.fswloc.location_sc)
        self.assertEqual(fswloc.fluence, self.fswloc.fluence)
        self.assertEqual(fswloc.hardness_ratio, self.fswloc.hardness_ratio)
        self.assertEqual(fswloc.intensity, self.fswloc.intensity)
        self.assertTupleEqual(fswloc.location, self.fswloc.location)
        self.assertTupleEqual(fswloc.location_sc, self.fswloc.location_sc)
        self.assertTupleEqual(fswloc.next_classification, self.fswloc.next_classification)
        self.assertListEqual(fswloc.rates.tolist(), self.fswloc.rates.tolist())
        self.assertEqual(fswloc.significance, self.fswloc.significance)
        self.assertEqual(fswloc.spectrum, self.fswloc.spectrum)
        self.assertEqual(fswloc.time, self.fswloc.time)
        self.assertEqual(fswloc.timescale, self.fswloc.timescale)
        self.assertTupleEqual(fswloc.top_classification, self.fswloc.top_classification)

        # wrong RA/Dec
        with self.assertRaises(ValueError):
            FswLocation.create(radec=self.fswloc.location[:1])
        with self.assertRaises(TypeError):
            FswLocation.create(radec=self.fswloc.location[0])

        # wrong spectrum
        with self.assertRaises(ValueError):
            FswLocation.create(algorithm='meow')

        # wrong top classification
        with self.assertRaises(ValueError):
            FswLocation.create(class1='cat')
        # wrong top classification prob
        with self.assertRaises(ValueError):
            FswLocation.create(class1_prob=-0.1)

        # wrong next classification
        with self.assertRaises(ValueError):
            FswLocation.create(class2='dog')
        # wrong next classification prob
        with self.assertRaises(ValueError):
            FswLocation.create(class2_prob=9.1)

        # wrong sum of class probs
        with self.assertRaises(ValueError):
            FswLocation.create(class1_prob=0.7, class2_prob=0.5)

        # wrong rates
        with self.assertRaises(ValueError):
            FswLocation.create(rates=self.fswloc.rates[:4])

        # wrong Az/Zen
        with self.assertRaises(ValueError):
            FswLocation.create(azzen=self.fswloc.location_sc[:1])
        with self.assertRaises(TypeError):
            FswLocation.create(azzen=self.fswloc.location_sc[0])

    def test_write(self):
        with TemporaryDirectory() as this_path:
            with Trigdat.open(trigdat_file) as trigdat:
                trigdat.write(this_path, overwrite=True)
                filename = trigdat.filename

            trigdat = Trigdat.open(os.path.join(this_path, filename))
            fswloc = trigdat.fsw_locations[0]

            self.assertEqual(fswloc.fluence, self.fswloc.fluence)
            self.assertEqual(fswloc.hardness_ratio, self.fswloc.hardness_ratio)
            self.assertEqual(fswloc.intensity, self.fswloc.intensity)
            self.assertTupleEqual(fswloc.location, self.fswloc.location)
            self.assertTupleEqual(fswloc.location_sc, self.fswloc.location_sc)
            self.assertTupleEqual(fswloc.next_classification, self.fswloc.next_classification)
            self.assertListEqual(fswloc.rates.tolist(), self.fswloc.rates.tolist())
            self.assertEqual(fswloc.significance, self.fswloc.significance)
            self.assertEqual(fswloc.spectrum, self.fswloc.spectrum)
            self.assertEqual(fswloc.time, self.fswloc.time)
            self.assertEqual(fswloc.timescale, self.fswloc.timescale)
            self.assertTupleEqual(fswloc.top_classification, self.fswloc.top_classification)

            trigdat.close()
