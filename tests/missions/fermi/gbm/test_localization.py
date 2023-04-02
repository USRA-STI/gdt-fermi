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
from pathlib import Path

from gdt.core import data_path
from gdt.missions.fermi.gbm.localization import *
from gdt.core.healpix import HealPixLocalization

chi2grid_file = Path(__file__).parent / 'data/chi2grid_bn190531568_v00.dat'
hpx_file = data_path / 'fermi-gbm/glg_healpix_all_bn190915240_v00.fit'


# RA=48.0, Dec=4.0, 1-deg radius
consistent_loc = HealPixLocalization.from_gaussian(48.0, 4.0, 3.0)
    
# RA=210.0, Dec=0.0, 1-deg radius
inconsistent_loc = HealPixLocalization.from_gaussian(210.0, 0.0, 3.0)

# RA=330.0, Dec=0.0, 1-deg radius
behind_earth_loc = HealPixLocalization.from_gaussian(330.0, 0.0, 3.0)

medium_consistent_loc = HealPixLocalization.from_gaussian(48.0, 4.0, 5.0)

large_consistent_loc = HealPixLocalization.from_gaussian(48.0, 4.0, 10.0)

@unittest.skipIf(not hpx_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestGbmHealPixFromFile(unittest.TestCase):
    
    def setUp(self):
        self.hpx = GbmHealPix.open(hpx_file)
    
    def tearDown(self):
        self.hpx.close()
    
    def test_centroid(self):
        ra, dec = self.hpx.centroid
        self.assertAlmostEqual(ra, 48.87, places=2)
        self.assertAlmostEqual(dec, 4.18, places=2)
    
    def test_detectors(self):
        n0 = self.hpx.n0_pointing
        self.assertAlmostEqual(n0.ra.value, 146.6, places=1)
        self.assertAlmostEqual(n0.dec.value, 37.0, places=1)
        n1 = self.hpx.n1_pointing
        self.assertAlmostEqual(n1.ra.value, 177.9, places=1)
        self.assertAlmostEqual(n1.dec.value, 37.9, places=1)
        n2 = self.hpx.n2_pointing
        self.assertAlmostEqual(n2.ra.value, 235.2, places=1)
        self.assertAlmostEqual(n2.dec.value, 32.5, places=1)
        n3 = self.hpx.n3_pointing
        self.assertAlmostEqual(n3.ra.value, 141.1, places=1)
        self.assertAlmostEqual(n3.dec.value, -11.9, places=1)
        n4 = self.hpx.n4_pointing
        self.assertAlmostEqual(n4.ra.value, 148.7, places=1)
        self.assertAlmostEqual(n4.dec.value, -57.7, places=1)
        n5 = self.hpx.n5_pointing
        self.assertAlmostEqual(n5.ra.value, 204.7, places=1)
        self.assertAlmostEqual(n5.dec.value, -14.2, places=1)
        n6 = self.hpx.n6_pointing
        self.assertAlmostEqual(n6.ra.value, 103.6, places=1)
        self.assertAlmostEqual(n6.dec.value, 19.9, places=1)
        n7 = self.hpx.n7_pointing
        self.assertAlmostEqual(n7.ra.value, 82.2, places=1)
        self.assertAlmostEqual(n7.dec.value, 4.9, places=1)
        n8 = self.hpx.n8_pointing
        self.assertAlmostEqual(n8.ra.value, 53.8, places=1)
        self.assertAlmostEqual(n8.dec.value, -31.2, places=1)
        n9 = self.hpx.n9_pointing
        self.assertAlmostEqual(n9.ra.value, 76.3, places=1)
        self.assertAlmostEqual(n9.dec.value, 65.4, places=1)
        na = self.hpx.na_pointing
        self.assertAlmostEqual(na.ra.value, 329.2, places=1)
        self.assertAlmostEqual(na.dec.value, 56.9, places=1)
        nb = self.hpx.nb_pointing
        self.assertAlmostEqual(nb.ra.value, 24.8, places=1)
        self.assertAlmostEqual(nb.dec.value, 13.8, places=1)
        b0 = self.hpx.b0_pointing
        self.assertAlmostEqual(b0.ra.value, 203.0, places=1)
        self.assertAlmostEqual(b0.dec.value, -17.1, places=1)
        b1 = self.hpx.b1_pointing
        self.assertAlmostEqual(b1.ra.value, 23.1, places=1)
        self.assertAlmostEqual(b1.dec.value, 17.1, places=1)

    def test_filename(self):
        self.assertEqual(self.hpx.filename, hpx_file.name)
    
    def test_frame(self):
        frame = self.hpx.frame
        self.assertListEqual(frame.obsgeoloc.xyz.value.tolist(), 
                             [-5039500., 4254000., -2067500.])
        self.assertListEqual(frame.quaternion.xyz.tolist(), 
                             [-0.223915, 0.447149, 0.860062])
        self.assertEqual(frame.quaternion.w, -0.101055)
        self.assertEqual(frame.obstime.fermi, 590219102.911008)

    def test_geo_location(self):
        loc = self.hpx.geo_location
        self.assertAlmostEqual(loc.ra.value, 319.83, places=2)
        self.assertAlmostEqual(loc.dec.value, 17.41, places=2)

    def test_geo_probability(self):
        self.assertLess(self.hpx.geo_probability, 0.01)
    
    def test_geo_radius(self):
        self.assertLessEqual(self.hpx.geo_radius.value, 68.0)
        self.assertGreaterEqual(self.hpx.geo_radius.value, 66.0)

    def test_headers(self):
        headers = self.hpx.headers
        self.assertEqual(headers.num_headers, 2)
    
    def test_npix(self):
        self.assertEqual(self.hpx.npix, 196608)

    def test_nside(self):
        self.assertEqual(self.hpx.nside, 128)

    def test_quaternion(self):
        self.assertListEqual(self.hpx.quaternion.xyz.tolist(), 
                             [-0.223915, 0.447149, 0.860062])
        self.assertEqual(self.hpx.quaternion.w, -0.101055)

    def test_scpos(self):
        self.assertListEqual(self.hpx.scpos.tolist(), 
                             [-5039500., 4254000., -2067500.])

    def test_sun_location(self):
        sun_loc = self.hpx.sun_location
        self.assertAlmostEqual(sun_loc.ra.value, 172.50, places=2)
        self.assertAlmostEqual(sun_loc.dec.value, 3.24, places=2)

    def test_trigtime(self):
        self.assertEqual(self.hpx.trigtime, 590219102.911008)

    def test_observable_fraction(self):
        frac = self.hpx.observable_fraction(consistent_loc)
        self.assertGreaterEqual(frac, 0.99)

        frac = self.hpx.observable_fraction(behind_earth_loc)
        self.assertLessEqual(frac, 0.01)

    def test_region_probability(self):
        # consistent small localization
        prob = self.hpx.region_probability(consistent_loc)
        self.assertGreaterEqual(prob, 0.99)
        
        # inconsistent small localization
        prob = self.hpx.region_probability(inconsistent_loc)
        self.assertLessEqual(prob, 0.01)

        # small localization behind earth
        prob = self.hpx.region_probability(behind_earth_loc)
        self.assertLessEqual(prob, 1e-10)
        
        # consistent medium-sized localization
        prob = self.hpx.region_probability(medium_consistent_loc)
        self.assertGreaterEqual(prob, 0.99)

        # consistent large localization
        prob = self.hpx.region_probability(large_consistent_loc)
        self.assertGreaterEqual(prob, 0.98)
        
        # incorrect prior values
        with self.assertRaises(ValueError):
            self.hpx.region_probability(consistent_loc, prior=-1.0)
        with self.assertRaises(ValueError):
            self.hpx.region_probability(consistent_loc, prior=1.1)

    def test_remove_earth(self):
        new_hpx = self.hpx.remove_earth()
        self.assertEqual(new_hpx.trigtime, self.hpx.trigtime)
        self.assertListEqual(new_hpx.scpos.tolist(), self.hpx.scpos.tolist())
        self.assertListEqual(new_hpx.quaternion.xyz.tolist(), 
                             self.hpx.quaternion.xyz.tolist())
        self.assertEqual(new_hpx.quaternion.w, self.hpx.quaternion.w)
        self.assertEqual(new_hpx.filename, self.hpx.filename)
        
        ra, dec = new_hpx.centroid
        self.assertAlmostEqual(ra, 48.87, places=2)
        self.assertAlmostEqual(dec, 4.18, places=2)

    def test_source_probability(self):
        # consistent
        prob = self.hpx.source_probability(48.0, 4.0)  
        self.assertGreaterEqual(prob, 0.99)      
                
        # inconsistent
        prob = self.hpx.source_probability(148.0, 0.0)  
        self.assertLessEqual(prob, 0.01)
        
        # behind earth
        prob = self.hpx.source_probability(330.0, 0.0)  
        self.assertEqual(prob, 0.0)

        # incorrect prior values
        with self.assertRaises(ValueError):
            self.hpx.source_probability(0.0, 0.0, prior=-1.0)
        with self.assertRaises(ValueError):
            self.hpx.source_probability(0.0, 0.0, prior=1.1)

    def test_from_data(self):
        new_hpx = GbmHealPix.from_data(self.hpx.prob, trigtime=self.hpx.trigtime,
                                       headers=self.hpx.headers, 
                                       filename=self.hpx.filename,
                                       quaternion=self.hpx.quaternion,
                                       scpos=self.hpx.scpos)
        self.assertListEqual(new_hpx.prob.tolist(), self.hpx.prob.tolist())      
        self.assertEqual(new_hpx.trigtime, self.hpx.trigtime)
        self.assertEqual(new_hpx.headers, self.hpx.headers)    
        self.assertEqual(new_hpx.filename, self.hpx.filename)
    
    def test_multiply(self):
        hpx = GbmHealPix.multiply(self.hpx, self.hpx.remove_earth())
        self.assertEqual(hpx.nside, 128)

        hpx = GbmHealPix.multiply(self.hpx, self.hpx.remove_earth(), output_nside=64)
        self.assertEqual(hpx.nside, 64)

        hpx = GbmHealPix.multiply(self.hpx, self.hpx.remove_earth(), primary=1)
        self.assertEqual(hpx.nside, 128)
        
        # the primary map is not a GbmHealPix
        with self.assertRaises(TypeError):
            GbmHealPix.multiply(consistent_loc, self.hpx)
        with self.assertRaises(TypeError):
            GbmHealPix.multiply(self.hpx, consistent_loc, primary=1)

    def test_from_data_errors(self):
        # wrong headers
        with self.assertRaises(TypeError):
            GbmHealPix.from_data(self.hpx.prob, headers=self.hpx.headers[0])
        
        # wrong quaternion
        with self.assertRaises(TypeError):
            GbmHealPix.from_data(self.hpx.prob, quaternion=1)
        with self.assertRaises(ValueError):
            GbmHealPix.from_data(self.hpx.prob, quaternion=self.hpx.quaternion.xyz)

        # wrong scpos
        with self.assertRaises(TypeError):
            GbmHealPix.from_data(self.hpx.prob, scpos=1)
        with self.assertRaises(ValueError):
            GbmHealPix.from_data(self.hpx.prob, scpos=self.hpx.scpos[1:])

    def test_write(self):

        with TemporaryDirectory() as this_path:
            self.hpx.write(this_path)
            hpx = GbmHealPix.open(os.path.join(this_path, self.hpx.filename))

            self.assertListEqual(hpx.prob.tolist(), self.hpx.prob.tolist())
            self.assertListEqual(hpx.sig.tolist(), self.hpx.sig.tolist())
            self.assertListEqual(hpx.scpos.tolist(), self.hpx.scpos.tolist())
            self.assertListEqual(hpx.quaternion.xyz.tolist(),
                                 self.hpx.quaternion.xyz.tolist())
            self.assertEqual(hpx.quaternion.w, self.hpx.quaternion.w)
            self.assertEqual(hpx.trigtime, self.hpx.trigtime)
            self.assertEqual(hpx.headers[1], self.hpx.headers[1])
            hpx.close()


@unittest.skipIf(not hpx_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestGbmHealPixNoFrame(unittest.TestCase):
    
    def setUp(self):
        hpx = GbmHealPix.open(hpx_file)
        self.hpx = GbmHealPix.from_data(hpx.prob, headers=hpx.headers,
                                        filename=hpx.filename)
    
    def test_detectors(self):
        n0 = self.hpx.n0_pointing
        self.assertAlmostEqual(n0.ra.value, 146.6, places=1)
        self.assertAlmostEqual(n0.dec.value, 37.0, places=1)
        n1 = self.hpx.n1_pointing
        self.assertAlmostEqual(n1.ra.value, 177.9, places=1)
        self.assertAlmostEqual(n1.dec.value, 37.9, places=1)
        n2 = self.hpx.n2_pointing
        self.assertAlmostEqual(n2.ra.value, 235.2, places=1)
        self.assertAlmostEqual(n2.dec.value, 32.5, places=1)
        n3 = self.hpx.n3_pointing
        self.assertAlmostEqual(n3.ra.value, 141.1, places=1)
        self.assertAlmostEqual(n3.dec.value, -11.9, places=1)
        n4 = self.hpx.n4_pointing
        self.assertAlmostEqual(n4.ra.value, 148.7, places=1)
        self.assertAlmostEqual(n4.dec.value, -57.7, places=1)
        n5 = self.hpx.n5_pointing
        self.assertAlmostEqual(n5.ra.value, 204.7, places=1)
        self.assertAlmostEqual(n5.dec.value, -14.2, places=1)
        n6 = self.hpx.n6_pointing
        self.assertAlmostEqual(n6.ra.value, 103.6, places=1)
        self.assertAlmostEqual(n6.dec.value, 19.9, places=1)
        n7 = self.hpx.n7_pointing
        self.assertAlmostEqual(n7.ra.value, 82.2, places=1)
        self.assertAlmostEqual(n7.dec.value, 4.9, places=1)
        n8 = self.hpx.n8_pointing
        self.assertAlmostEqual(n8.ra.value, 53.8, places=1)
        self.assertAlmostEqual(n8.dec.value, -31.2, places=1)
        n9 = self.hpx.n9_pointing
        self.assertAlmostEqual(n9.ra.value, 76.3, places=1)
        self.assertAlmostEqual(n9.dec.value, 65.4, places=1)
        na = self.hpx.na_pointing
        self.assertAlmostEqual(na.ra.value, 329.2, places=1)
        self.assertAlmostEqual(na.dec.value, 56.9, places=1)
        nb = self.hpx.nb_pointing
        self.assertAlmostEqual(nb.ra.value, 24.8, places=1)
        self.assertAlmostEqual(nb.dec.value, 13.8, places=1)
        b0 = self.hpx.b0_pointing
        self.assertAlmostEqual(b0.ra.value, 203.0, places=1)
        self.assertAlmostEqual(b0.dec.value, -17.1, places=1)
        b1 = self.hpx.b1_pointing
        self.assertAlmostEqual(b1.ra.value, 23.0, places=1)
        self.assertAlmostEqual(b1.dec.value, 17.1, places=1)

    def test_frame(self):
        self.assertIsNone(self.hpx.frame)

    def test_geo_location(self):
        loc = self.hpx.geo_location
        self.assertAlmostEqual(loc.ra.value, 319.83, places=2)
        self.assertAlmostEqual(loc.dec.value, 17.41, places=2)

    def test_geo_probability(self):
        self.assertLess(self.hpx.geo_probability, 0.01)

    def test_geo_radius(self):
        self.assertEqual(self.hpx.geo_radius.value, 67.5)

    def test_quaternion(self):
        self.assertIsNone(self.hpx.quaternion)

    def test_scpos(self):
        self.assertIsNone(self.hpx.scpos)

    def test_sun_location(self):
        loc = self.hpx.sun_location
        self.assertAlmostEqual(loc.ra.value, 172.50, places=2)
        self.assertAlmostEqual(loc.dec.value, 3.24, places=2)
    
    def test_observable_fraction(self):
        frac = self.hpx.observable_fraction(consistent_loc)
        self.assertGreaterEqual(frac, 0.99)

    def test_region_probability(self):
        # consistent small localization
        prob = self.hpx.region_probability(consistent_loc)
        self.assertGreaterEqual(prob, 0.99)
    
    def test_source_probability(self):
        # consistent
        prob = self.hpx.source_probability(48.0, 4.0)  
        self.assertGreaterEqual(prob, 0.99)      


@unittest.skipIf(not hpx_file.exists(), "test files aren't downloaded. run gdt-download-data.")
class TestGbmHealPixNoHeaders(unittest.TestCase):
    
    def setUp(self):
        hpx = GbmHealPix.open(hpx_file)
        self.hpx = GbmHealPix.from_data(hpx.prob, filename=hpx.filename)

    def test_detectors(self):
        with self.assertRaises(AttributeError):
            self.hpx.n0_pointing

    def test_geo_location(self):
        self.assertIsNone(self.hpx.geo_location)

    def test_geo_probability(self):
        self.assertIsNone(self.hpx.geo_probability)

    def test_sun_location(self):
        self.assertIsNone(self.hpx.sun_location)

    def test_observable_fraction(self):
        with self.assertRaises(RuntimeError):
            self.hpx.observable_fraction(consistent_loc)
 
    def test_region_probability(self):
        # consistent small localization
        prob = self.hpx.region_probability(consistent_loc)
        self.assertGreaterEqual(prob, 0.99)

    def test_remove_earth(self):
        with self.assertRaises(RuntimeError):
            self.hpx.remove_earth()

    def test_source_probability(self):
        # consistent
        prob = self.hpx.source_probability(48.0, 4.0)  
        self.assertGreaterEqual(prob, 0.99)      
 

class TestGbmHealPixFromChi2Grid(unittest.TestCase):
    
    def setUp(self):
        chi2grid = Chi2Grid.open(str(chi2grid_file))
        chi2grid.quaternion = [-0.223915, 0.447149, 0.860062, -0.101055]
        chi2grid.scpos = [-5039500., 4254000., -2067500.]
        chi2grid.trigtime = 581002688
        self.hpx = GbmHealPix.from_chi2grid(chi2grid)

    def test_centroid(self):
        ra, dec = self.hpx.centroid
        self.assertAlmostEqual(ra, 275.98, places=2) 
        self.assertAlmostEqual(dec, 38.68, places=2) 
    
    def test_quaternion(self):
        self.assertListEqual(self.hpx.quaternion.xyz.tolist(), 
                             [-0.223915, 0.447149, 0.860062])
        self.assertEqual(self.hpx.quaternion.w, -0.101055)
    
    def test_scpos(self):
        self.assertListEqual(self.hpx.scpos.tolist(), 
                             [-5039500., 4254000., -2067500.])
    
    def test_trigtime(self):
        self.assertEqual(self.hpx.trigtime, 581002688)
    
        
class TestChi2Grid(unittest.TestCase):
    
    def setUp(self):
        self.chi2grid = Chi2Grid.open(str(chi2grid_file))
        self.numpts = 41168
    
    def test_azimuth(self):
        az = self.chi2grid.azimuth
        self.assertEqual(az.size, self.numpts)
        self.assertListEqual((az <= 360.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        self.assertListEqual((az >= 0.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())

    def test_chisq(self):
        chisq = self.chi2grid.chisq
        self.assertEqual(chisq.size, self.numpts)
        self.assertListEqual((chisq >= 0.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())

    def test_dec(self):
        dec = self.chi2grid.dec
        self.assertEqual(dec.size, self.numpts)
        self.assertListEqual((dec >= -90.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        self.assertListEqual((dec <= 90.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        
    def test_numpts(self):
        self.assertEqual(self.chi2grid.numpts, self.numpts)

    def test_quaternion(self):
        self.assertIsNone(self.chi2grid.quaternion)
        
        quat = [0.0, 1.0, 0.0, 1.0]
        self.chi2grid.quaternion = quat
        self.assertListEqual(self.chi2grid.quaternion.tolist(), quat)
        
        # wrong number of elements
        with self.assertRaises(ValueError):
            self.chi2grid.quaternion = quat[1:]
    
    def test_ra(self):
        ra = self.chi2grid.ra
        self.assertEqual(ra.size, self.numpts)
        self.assertListEqual((ra >= 0.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        self.assertListEqual((ra <= 360.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())

    def test_scpos(self):
        self.assertIsNone(self.chi2grid.scpos)
        
        scpos = np.array([4., 3., -2.]) * 100.0
        self.chi2grid.scpos = scpos
        self.assertListEqual(self.chi2grid.scpos.tolist(), scpos.tolist())
    
        # wrong number of elements
        with self.assertRaises(ValueError):
            self.chi2grid.scpos = scpos[1:]

    def test_significance(self):
        sig = self.chi2grid.significance
        self.assertEqual(sig.size, self.numpts)
        self.assertListEqual((sig >= 0.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        self.assertListEqual((sig <= 1.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())

    def test_trigtime(self):
        self.assertIsNone(self.chi2grid.trigtime)
        
        trigtime = 696638874.0
        self.chi2grid.trigtime = trigtime
        self.assertEqual(self.chi2grid.trigtime, trigtime)
    
        # wrong value
        with self.assertRaises(ValueError):
            self.chi2grid.trigtime = 'bad input'
         
    def test_zenith(self):
        zen = self.chi2grid.zenith
        self.assertEqual(zen.size, self.numpts)
        self.assertListEqual((zen >= 0.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        self.assertListEqual((zen <= 180.0).tolist(), 
                             np.ones(self.numpts, dtype=bool).tolist())
        
    def test_from_data(self):
        chi2grid2 = Chi2Grid.from_data(self.chi2grid.azimuth, self.chi2grid.zenith,  
                                       self.chi2grid.ra, self.chi2grid.dec,
                                       self.chi2grid.chisq)
        self.assertListEqual(chi2grid2.azimuth.tolist(), 
                             self.chi2grid.azimuth.tolist())
        self.assertListEqual(chi2grid2.zenith.tolist(), 
                             self.chi2grid.zenith.tolist())
        self.assertListEqual(chi2grid2.ra.tolist(), self.chi2grid.ra.tolist())
        self.assertListEqual(chi2grid2.dec.tolist(), self.chi2grid.dec.tolist())
        self.assertListEqual(chi2grid2.chisq.tolist(), self.chi2grid.chisq.tolist())
        
        # not all array lengths match
        with self.assertRaises(ValueError):
            Chi2Grid.from_data(self.chi2grid.azimuth[1:], self.chi2grid.zenith,
                               self.chi2grid.ra, self.chi2grid.dec, 
                               self.chi2grid.chisq)


# systematic model tests

class TestGbutsO3Model(unittest.TestCase):
    
    def test_sigma(self):
        sigma, _ = gbuts_o3_model()
        self.assertAlmostEqual(np.rad2deg(sigma[0]), 2.7, places=1)

    def test_frac(self):
        _, frac = gbuts_o3_model()
        self.assertListEqual(frac, [1.0])
        
    
class TestHitlModel(unittest.TestCase):
    
    def test_az_gt_292_5(self):
        sigma, frac = hitl_model(293.0)
        self.assertListEqual(sigma, [np.deg2rad(4.17), np.deg2rad(15.3)])
        self.assertListEqual(frac, [0.918])        

    def test_az_le_67_5(self):
        sigma, frac = hitl_model(60.0)
        self.assertListEqual(sigma, [np.deg2rad(4.17), np.deg2rad(15.3)])
        self.assertListEqual(frac, [0.918])        

    def test_az_between_112_5_and_247_5(self):
        sigma, frac = hitl_model(60.0)
        self.assertListEqual(sigma, [np.deg2rad(4.17), np.deg2rad(15.3)])
        self.assertListEqual(frac, [0.918])        

    def test_az_between_67_5_and_112_5(self):
        sigma, frac = hitl_model(75.0)
        self.assertListEqual(sigma, [np.deg2rad(2.31), np.deg2rad(13.2)])
        self.assertListEqual(frac, [0.884])        

    def test_az_between_247_5_and_292_5(self):
        sigma, frac = hitl_model(250.0)
        self.assertListEqual(sigma, [np.deg2rad(2.31), np.deg2rad(13.2)])
        self.assertListEqual(frac, [0.884])        


class TestGaModel(unittest.TestCase):
    
    def test_sigma(self):
        sigma, _ = ga_model()
        self.assertListEqual(sigma, [np.deg2rad(3.72), np.deg2rad(13.7)])

    def test_frac(self):
        _, frac = ga_model()
        self.assertListEqual(frac, [0.804])
    

class TestRoboBaModel(unittest.TestCase):
    
    def test_long(self):
        sigma, frac = robo_ba_model('long')
        self.assertListEqual(sigma, [np.deg2rad(1.86), np.deg2rad(4.14)])
        self.assertListEqual(frac, [0.579])

    def test_short(self):
        sigma, frac = robo_ba_model('short')
        self.assertListEqual(sigma, [np.deg2rad(2.55), np.deg2rad(4.43)])
        self.assertListEqual(frac, [0.390])

    def test_error(self):
        with self.assertRaises(ValueError):
            robo_ba_model('bad input')


class TestUntargetedSearchModel(unittest.TestCase):
    
    def test_sigma(self):
        sigma, _ = untargeted_search_model()
        self.assertListEqual(sigma, [np.deg2rad(5.53)])

    def test_frac(self):
        _, frac = untargeted_search_model()
        self.assertListEqual(frac, [1.0])
    
