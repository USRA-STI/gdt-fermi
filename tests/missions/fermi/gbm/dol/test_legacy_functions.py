# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Very closely based on the program DoL (Daughter of Locburst).
# Written by:
#               Valerie Connaughton
#               University of Alabama in Huntsville (UAH)
#
# Included in the Gamma-ray Data Tools with permission from UAH
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
import numpy as np
from tests import tests_path
from pathlib import Path
from gdt.missions.fermi.gbm.localization.dol import legacy_spectral_models, __data_dir__
from gdt.missions.fermi.gbm.localization.dol.legacy_functions import *

from . import __test_ref_dir__, __test_args_dir__


class TestLegacyFunctions(unittest.TestCase):
    r""" Tests legacy functions for consistency """

    @classmethod
    def setUpClass(cls):
        r""" Function to setup any class globals that you want
        to use across individual member functions """
        # load dictionary with variables needed as args to legacy functions
        cls.args = cls.loadz(__test_args_dir__ + "legacy_functions.npz")

        # load dictionary with reference for function outputs
        cls.ref = cls.loadz(__test_ref_dir__ + "legacy_functions.npz")

        # number of points to test for ra/dec, az/zen
        n = 100

        # ra and dec values to test
        ra = np.float64(np.radians(np.linspace(-1, 361., n)))
        dec = np.float64(np.radians(np.linspace(-91., 91., n)))
        cls.ra, cls.dec = np.array(np.meshgrid(ra, dec)).T.reshape(-1, 2).T

        # az and zen values to test
        az = np.linspace(0, 2 * np.pi, n, dtype=np.float32)
        zen = np.linspace(-0.5 * np.pi, 0.5 * np.pi, n, dtype=np.float32)
        cls.az, cls.zen = np.array(np.meshgrid(az, zen)).T.reshape(-1, 2).T

        # set this to True to generate new reference file
        cls.update_ref = False

    @classmethod
    def tearDownClass(cls):
        if cls.update_ref:
            path = __test_ref_dir__ + "legacy_functions.npz"
            np.savez_compressed(path, cls.ref)

    @staticmethod
    def fix_keys(d):
        r""" Fix dictionary keys to be regular strings instead of byte strings

        Parameters
        ----------
        d : dict
            Input dictionary with byte string keys
        """
        for key in list(d.keys()):
            try:
                new_key = key.decode("utf-8")
                d[new_key] = d.pop(key)
            except BaseException:
                pass

        # END for (key)

    @classmethod
    def loadz(cls, path):
        r""" Load compressed numpy file

        Parameters
        ----------
        path : str
            Path to compressed numpy file (.npz)

        Returns
        -------
        d : dict
            Uncompressed dictionary object
        """
        # load dict from comptressed numpy file
        npz = np.load(path, allow_pickle=True, encoding='bytes')
        d = npz.f.arr_0.item()

        # fix formatting issues between python2/3 strings
        cls.fix_keys(d)
        for key, val in d.items():
            if isinstance(val, dict):
                cls.fix_keys(d[key])
            # END if (dict))
            try:
                val = val.decode("utf-8")
                d[key] = val
            except BaseException:
                pass
                # END for (key, val)

        return d

    def get_ref(self, key, value):
        r""" Function to retrieve reference value(s) for a given test

        Paramaters
        ----------
        key : str 
            Key name for reference values of the test
        value : multiple types
            Value to use when updating the reference value. The value type
            depends on the key though most are np.ndarray or dict objects

        Returns
        -------
        ref : multiple types
            Reference value(s) for this test. The value type depends on the
            key though most are np.ndarray or dict objects
        """
        # update reference values when requested
        if self.update_ref:
            self.ref[key] = value

        return self.ref[key]

    def test_arctan2(self):
        r""" Test behavior of arctan2 function """
        # input values
        x = np.linspace(-1, 1, 10)
        y = np.linspace(-1, 1, 10)
        x, y = np.array(np.meshgrid(x, y)).T.reshape(-1, 2).T

        test_results = []
        for i in range(x.size):
            test_results.append(arctan2(x[i], y[i]))
        test_results.append(arctan2(1e-8, 4e-10))
        test_results.append(arctan2(2.0, -3.8))
        test_results = np.array(test_results, np.float64)

        # compare test results to reference
        ref = self.get_ref("test_arctan2", test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_choose_energy_data(self):
        r""" Test behavior of data selected within a given energy range """
        c_brates, c_mrates, cenergies = choose_energy_data(
            self.args['ndet'], self.args['nen'],
            self.args['crange'], self.args['mrates'],
            self.args['brates'], self.args['energies'])

        test_results = {'c_brates': c_brates,
                        'c_mrates': c_mrates,
                        'cenergies': cenergies}

        # ensure all outputs match
        ref = self.get_ref("test_choose_energy_data", test_results)
        self.assertTrue(np.allclose(c_brates, ref['c_brates'], rtol=1e-02))
        self.assertTrue(np.allclose(c_mrates, ref['c_mrates'], rtol=1e-02))
        self.assertTrue(np.allclose(cenergies, ref['cenergies'], rtol=1e-02))

    def test_get_geocenter(self):
        r""" Test calculation of the direction of Earth's center """
        sc_quat = self.args['sc_quat']
        sc_pos1 = self.args['sc_pos']
        sc_pos2 = np.array([0, 0, -0.0000000000050002])

        geo = []
        geo.append(get_geocenter(sc_quat, sc_pos1))
        geo.append(get_geocenter(sc_quat, sc_pos2))

        # compare to reference
        ref = self.get_ref('test_get_geocenter', np.array(geo, dtype=object))
        for i, g in enumerate(geo):
            for j, val in enumerate(g):
                self.assertTrue(np.allclose(val, ref[i][j], rtol=1e-02))

    def test_read_table(self):
        r""" Test the method for reading table of expected rates
             for each location on the sky """

        path = __data_dir__ + "band_1deg_50_300_norm.npy"
        ndet = self.args['ndet']
        npoints = self.args['npoints']

        table, idbno = read_table(path, ndet, npoints)
        test_results = {'table': table, 'idbno': idbno}

        # compare to reference values 
        ref = self.get_ref('test_read_table', test_results)
        self.assertTrue(np.allclose(table, ref['table'], rtol=1e-02))
        self.assertEqual(idbno, ref['idbno'])

        # check ability to read text file
        path = __test_args_dir__ + "locrates_1deg_50_300_normal_n.dat"
        table, idbno = read_table(path, ndet, npoints)
        test_results = {'table': table, 'idbno': idbno}

        # compare to reference values 
        ref = self.get_ref('test_read_table', test_results)
        self.assertTrue(np.allclose(table, ref['table'], rtol=1e-02))
        self.assertEqual(idbno, ref['idbno'])

    def test_ang_to_cart_zen(self):
        r""" Test conversion of az/zen to position vector """
        # array to store test results
        dtype = [('az', np.float32), ('zen', np.float32),
                 ('x', np.float32), ('y', np.float32), ('z', np.float32)]
        test_results = np.zeros(self.az.size + 1, dtype)

        # test scalar input behavior
        scalar_ret = ang_to_cart_zen(self.az[0], self.zen[0])
        test_results[0] = (self.az[0], self.zen[0],) + tuple(scalar_ret)

        # test vector input behavior
        vector_ret = ang_to_cart_zen(self.az, self.zen)
        for i in range(self.az.size):
            test_results[i + 1] = (self.az[i], self.zen[i],) \
                                  + tuple(vector_ret[i])

        # compare to reference values
        ref = self.get_ref('test_ang_to_cart_zen', test_results)
        self.assertTrue(
            np.allclose(test_results.view(np.float32), ref.view(np.float32), rtol=1e-02))

    def test_ang_to_cart_dec(self):
        r""" Test conversion of ra/dec to position vector """
        # array to store test results
        dtype = [('ra', np.float64), ('dec', np.float64),
                 ('x', np.float32), ('y', np.float32), ('z', np.float32)]
        test_results = np.zeros(self.ra.size, dtype)

        # loop to generate test results
        for i in range(self.ra.size):
            pos = ang_to_cart_dec(self.ra[i], self.dec[i])
            test_results[i] = (self.ra[i], self.dec[i],) + tuple(pos)

        # update reference values when requested
        ref = self.ref['test_ang_to_cart_dec']
        if self.update_ref:
            ref = test_results

        # compare to reference values
        for key in test_results.dtype.names:
            self.assertTrue(np.allclose(test_results[key], ref[key], rtol=1e-02))

    def test_j2000_to_sc(self):
        r""" Test conversion of j2000 direction to local space craft az/zen """
        # loop to gather test output
        test_results = []
        for i in range(self.ra.size):
            test_results.append(
                j2000_to_sc(self.args['scx'],
                            self.args['scy'],
                            self.args['scz'],
                            ang_to_cart_dec(self.ra[i], self.dec[i])))
        # END for (i)
        test_results = np.array(test_results)

        # compare test output to reference
        ref = self.get_ref('test_j2000_to_sc', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_sc_to_j2000(self):
        r""" Test conversion of local space craft az/zen to j2000 direction """
        # loop to gather test output
        test_results = []
        for i in range(self.ra.size):
            test_results.append(
                sc_to_j2000(self.args['scx'],
                            self.args['scy'],
                            self.args['scz'],
                            ang_to_cart_zen(self.az[i], self.zen[i])))
        # END for (i)
        test_results = np.array(test_results)

        # compare test output to reference
        ref = self.get_ref('test_sc_to_j2000', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-2))

    def test_get_occult(self):
        r""" Test function which returns locations occulted by the Earth """
        test_result = get_occult(self.args['sc_pos'],
                                 self.args['npoints'],
                                 self.args['j2000grid'])

        # compare test output to reference
        ref = self.get_ref('test_get_occult', test_result)
        self.assertTrue(np.allclose(test_result, ref, rtol=1e-02))

    def test_get_good_angle(self):
        r""" Test function to get angle between directions """
        input_xy1 = np.array([1, 3], np.float32)
        input_xy2 = np.array([5, 2], np.float32)
        input_xy3 = np.array([-5, 3], np.float32)
        input_xy4 = np.array([5, -1], np.float32)

        test_results = np.array(
            [get_good_angle(input_xy1, input_xy2),
             get_good_angle(input_xy3, input_xy4)])

        # compare test output to reference
        ref = self.get_ref('test_get_good_angle', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_compute_scat(self):
        r""" Test the computation of scatter rates """
        front, back = compute_scat(self.args['npoints'],
                                   self.args['rgrid'],
                                   self.args['cenergies'],
                                   self.args['geodir'],
                                   self.args['nai_az'],
                                   self.args['nai_zen'],
                                   self.args['nai_unit_vec'],
                                   self.args['back_unit_vec'])

        test_result = {'front': front, 'back': back}

        # compare test output to reference
        ref = self.get_ref('test_compute_scat', test_result)
        rtol = 1e-2
        mask = front > rtol
        self.assertTrue(np.allclose(front[mask], ref['front'][mask], rtol))
        mask = back > rtol
        self.assertTrue(np.allclose(back[mask], ref['back'][mask], rtol))

    def test_read_earthpoints(self):
        r""" Test ability to read earthpoints from file """
        earthpoints, gperp, gpmag, elev = read_earthpoints(
            self.args['geodir'], self.args['nai_unit_vec'],
            self.args['back_unit_vec'], 0)
        earthpoints, back_gperp, back_gpmag, back_elev = read_earthpoints(
            self.args['geodir'], self.args['nai_unit_vec'],
            self.args['back_unit_vec'], 1)

        test_results = {'earthpoints': earthpoints,
                        'gperp': gperp, 'gpmag': gpmag, 'elev': elev,
                        'back_gperp': back_gperp, 'back_gpmag': back_gpmag,
                        'back_elev': back_elev}
        # compare to reference values
        ref = self.get_ref('test_read_earthpoints', test_results)
        rtol=1e-02
        self.assertTrue(np.allclose(elev, ref['elev'], rtol=rtol))
        self.assertTrue(np.allclose(gperp, ref['gperp'], rtol=rtol))
        self.assertTrue(np.allclose(gpmag, ref['gpmag'], rtol=rtol))
        self.assertTrue(np.allclose(back_elev, ref['back_elev'], rtol=rtol))
        self.assertTrue(np.allclose(back_gperp, ref['back_gperp'], rtol=rtol))
        self.assertTrue(np.allclose(back_gpmag, ref['back_gpmag'], rtol=rtol))
        self.assertTrue(np.allclose(earthpoints, ref['earthpoints'], rtol=rtol))

        # test reading legacy text file
        earthpoints, gperp, gpmag, elev = read_earthpoints(
            self.args['geodir'], self.args['nai_unit_vec'],
            self.args['back_unit_vec'], 0,
            path=__test_args_dir__ + "earth_points.dat")
        # compare to reference values
        self.assertTrue(np.allclose(earthpoints, ref['earthpoints'], rtol=rtol))

    def test_read_scatter_data(self):
        r""" Test ability to read scattering data """
        # try default reading from .npy file
        test_result = read_scatter_data()

        # compare to reference
        ref = self.get_ref('test_read_scatter_data', test_result)
        self.assertTrue(np.allclose(test_result, ref, rtol=1e-02))

        # try reading legacy text file
        test_result = read_scatter_data(path=__test_args_dir__ + "alocdat_comp.dat")

        # compare to reference
        ref = self.get_ref('test_read_scatter_data', test_result)
        self.assertTrue(np.allclose(test_result, ref, rtol=1e-02))

    def test_get_scattered_rates(self):
        r""" Test function used to obtain rate in each detector
             with scattering affects included. """
        fscattered_rates = get_scattered_rates(
            rcart=self.args['rcart'],
            gperp=self.ref['test_read_earthpoints']['gperp'],
            gpmag=self.ref['test_read_earthpoints']['gpmag'],
            elev=self.ref['test_read_earthpoints']['elev'],
            earthpoints=self.ref['test_read_earthpoints']['earthpoints'],
            response_res=self.args['response_res'],
            grid_points=self.args['grid_points'],
            scatterdata=self.ref['test_read_scatter_data'],
            scat_geom=0,
            nai_az=self.args['nai_az'],
            nai_zen=self.args['nai_zen'],
            nai_unit_vec=self.args['nai_unit_vec'],
            back_unit_vec=self.args['back_unit_vec'])

        bscattered_rates = get_scattered_rates(
            rcart=self.args['rcart'],
            gperp=self.ref['test_read_earthpoints']['back_gperp'],
            gpmag=self.ref['test_read_earthpoints']['back_gpmag'],
            elev=self.ref['test_read_earthpoints']['back_elev'],
            earthpoints=self.ref['test_read_earthpoints']['earthpoints'],
            response_res=self.args['response_res'],
            grid_points=self.args['grid_points'],
            scatterdata=self.ref['test_read_scatter_data'],
            scat_geom=1,
            nai_az=self.args['nai_az'],
            nai_zen=self.args['nai_zen'],
            nai_unit_vec=self.args['nai_unit_vec'],
            back_unit_vec=self.args['back_unit_vec'])

        test_results = {'fscattered_rates': fscattered_rates,
                        'bscattered_rates': bscattered_rates}

        # compare to reference
        ref = self.get_ref('test_get_scattered_rates', test_results)
        rtol = 1e-2
        mask = fscattered_rates > rtol
        self.assertTrue(np.allclose(fscattered_rates[mask], ref['fscattered_rates'][mask], rtol=rtol))
        mask = bscattered_rates > rtol
        self.assertTrue(np.allclose(bscattered_rates[mask], ref['bscattered_rates'][mask], rtol=rtol))

    def test_crossprod(self):
        r""" Test computation of cross product between vectors """
        x = np.linspace(1e-2, 1, 10, np.float32)
        y = np.linspace(1e-2, 1, 10, np.float32)
        z = np.linspace(1e-2, 1, 10, np.float32)
        v = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)

        test_results = []
        for i in range(len(v)):
            test_results.append(crossprod(v[i], v[-i - 1]))
        test_results = np.array(test_results)

        # compare to reference
        ref = self.get_ref('test_crossprod', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_all_idx(self):
        r""" Test function which returns array of indices """
        idx = np.indices((100, 100), np.int32)
        test_result = all_idx(idx, axis=0)

        # compare to reference
        ref = self.get_ref('test_all_idx', test_result)
        for result, r in zip(test_result, ref):
            self.assertTrue(np.allclose(result, r, rtol=1e-02))

    def test_initialize_det_geometry(self):
        r""" Test detector geometry initialization function """
        nai_az, nai_zen, nai_unit_vec, back_unit_vec = \
            initialize_det_geometry()

        self.assertTrue(np.allclose(nai_az, self.args['nai_az'], rtol=1e-02))
        self.assertTrue(np.allclose(nai_zen, self.args['nai_zen'], rtol=1e-02))
        self.assertTrue(np.allclose(nai_unit_vec, self.args['nai_unit_vec'], rtol=1e-02))
        self.assertTrue(np.allclose(back_unit_vec, self.args['back_unit_vec'], rtol=1e-02))

    def test_get_spec(self):
        r""" Test function used to retrieve spectrum """
        test_results = []
        for spec in [legacy_spectral_models.band_hard,
                     legacy_spectral_models.band_norm,
                     legacy_spectral_models.band_soft,
                     legacy_spectral_models.comp_hard,
                     legacy_spectral_models.comp_norm,
                     legacy_spectral_models.comp_soft,
                     "pl,index=-2.00,amp=10.0"]:
            test_results.append(
                get_spec(spec=spec,
                         energies=legacy_tenergies,
                         mid_energies=self.args['mid_energies'],
                         erange=self.args['erange']))
        # END for (spec)
        test_results = np.array(test_results)

        # compare to reference
        ref = self.get_ref('test_get_spec', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_add_scat(self):
        r""" Test function which adds scattering effects to rates """
        test_result = add_scat(
            npoints=self.args['npoints'],
            rgrid=self.args['rgrid'],
            front_scattered_rates=self.args['front_scattered_rates'],
            back_scattered_rates=self.args['back_scattered_rates'],
            cenergies=self.args['cenergies'],
            spec=legacy_spectral_models.band_norm)

        # compare to reference
        ref = self.get_ref('test_add_scat', test_result)
        self.assertTrue(np.allclose(test_result, ref, rtol=1e-02))

    def test_find_best_location(self):
        r""" Test function which scans to find best source location on the sky """
        test_result = find_best_location(
            ndet=self.args['ndet'],
            vpoints=self.args['npoints'],
            udet=self.args['udet'],
            usedet=self.args['usedet'],
            loctable_entries=self.args['loctable_entries'],
            sduration=self.args['sduration'],
            c_mrates=self.args['c_mrates'],
            c_brates=self.args['c_brates'],
            visible=self.args['visible'])

        # compare to reference
        ref = self.get_ref('test_find_best_location', test_result)
        for i in range(len(test_result)):
            self.assertTrue(np.allclose(test_result[i], ref[i], rtol=1e-02))

    def test_get_det_geometry(self):
        r""" Test function which computes detector geometry """
        test_results = []
        for i in range(self.args['nai_zen'].size):
            test_results.append(
                get_det_geometry(self.az, self.zen,
                                 nai_az=self.args['nai_az'][i],
                                 nai_zen=self.args['nai_zen'][i])
            )
        # END for (i)

        # compare to reference
        ref = self.get_ref('test_get_det_geometry', test_results)
        for i in range(len(test_results)):
            mask = ~np.isclose(test_results[i], ref[i])
            self.assertTrue(np.allclose(test_results[i], ref[i], rtol=1e-2))

    def test_eq2000_to_gal_r(self):
        r""" Test conversion from J2000 equatorial coord to galactic coord """
        dtype = [('ra', np.float64), ('dec', np.float64),
                 ('l', np.float64), ('b', np.float64)]
        test_results = np.empty(self.ra.size, dtype)

        for i in range(self.ra.size):
            test_results = (self.ra[i], self.dec[i],) + \
                           tuple(eq2000_to_gal_r(self.ra[i], self.dec[i]))
        # END for (i)

        # compare to reference
        ref = self.get_ref('test_eq2000_to_gal_r', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))

    def test_get_etog_matrix(self):
        r""" Test calculation of equatorial coord
             to galactic coord rotation matrix """
        etog = ETOG().value

        # compare to reference
        ref = self.get_ref('test_get_etog_matrix', etog)
        self.assertTrue(np.allclose(etog, ref, rtol=1e-02))

    def test_sun_loc(self):
        r""" Test calculation of sun location """
        test_results = np.array(
            [sun_loc(self.args['sc_time']),
             sun_loc(132), sun_loc(-1322), sun_loc(1)])

        # compare to reference 
        ref = self.get_ref('test_sun_loc', test_results)
        self.assertTrue(np.allclose(test_results, ref, rtol=1e-02))


# END class TestLegacyFunctions

if __name__ == '__main__':
    unittest.main()
