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
import sys
import copy
import runpy
import shutil
import unittest
import tempfile
import numpy as np

from gdt.core.coords import Quaternion
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.frame import FermiFrame
from astropy.coordinates.representation import CartesianRepresentation

from . import __test_ref_dir__, __test_args_dir__

class TestLegacyDol(unittest.TestCase):
    r""" Tests legacy functions for consistency

    Note: These unit tests use allclose() with a 1% tolerance because
    legacy float32 support is not consistent across systems. In particular,
    this impacts sine and cosine operations throughout the program.
    """

    @classmethod
    def setUpClass(cls):
        r""" Function to setup any class globals that you want
        to use across individual member functions """
        cls.tmp_dir = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        r""" Function to clean up at end of tests """
        # remove the temp directory
        shutil.rmtree(cls.tmp_dir)

    @staticmethod
    def diff_dolfile(path1, path2):
        r""" Check for differences between two dol files

        Parameters
        ----------
        path1 : str
            Path to file #1
        path1 : str
            Path to file #1

        Returns
        -------
        diff : bool
            True when the files are different
        """
        content1 = (" ".join(open(path1).readlines())).split()
        content2 = (" ".join(open(path2).readlines())).split()

        # remove version string
        content1.pop(7)
        content2.pop(7)

        return not np.allclose(np.array(content1, np.float64), np.array(content2, np.float64), rtol=0.01)

    @staticmethod
    def diff_chi2grid(path1, path2):
        r""" Check for differences between two chi2grid files

        Parameters
        ----------
        path1 : str
            Path to file #1
        path2 : str
            Path to file #2

        Returns
        -------
        diff : bool
            True when the files are different
        """
        header1 = open(path1).readline()
        header2 = open(path2).readline()

        if len(header1) != len(header2):
            return True

        arr1 = np.loadtxt(path1, skiprows=1)
        arr2 = np.loadtxt(path2, skiprows=1)

        return not np.allclose(arr1[:, 2], arr2[:, 2], rtol=0.01)

    def diff(self, path1, path2):
        r""" Check for differences between two files

        Parameters
        ----------
        path1 : str
            Path to file #1
        path2 : str
            Path to file #2

        Returns
        -------
        diff : bool
            True when the files are different
        """
        if not os.path.exists(path1) or not os.path.exists(path2):
            return True

        if "chi2grid" in path1:
            return self.diff_chi2grid(path1, path2)
        return self.diff_dolfile(path1, path2)

    def compare(self, args, output, reference=None):
        r""" Compare the output of legacy_dol() for given arguments

        Parameters
        ----------
        args : str
            Path to numpy file with argument dictionary for legacy_dol()
        dol_file : str
            Name of output dol file
        chi2_file : str
            Name of output chi2 grid file

        Returns
        -------
        comparison : bool
            True when outputs match output stored in test data directory
        """
        args = np.load(__test_args_dir__ + args, allow_pickle=True, encoding='bytes').item()

        # construct argv to pass to legacy_dol __main__
        argv = ["dol.legacy_dol"]
        for k, v in args.items():
            argv.append(str("--" + k))
            if k in ['mrates', 'brates']:
                for i in v:
                    argv.extend([str(j) for j in i])
            else:
                argv.extend([str(i) for i in np.atleast_1d(v)])
        # ensure we output test files to temporary dir
        argv.extend(["--dir", self.tmp_dir])

        # save original argv so we can restore it after
        orig_argv = copy.copy(sys.argv)

        # run main program in legacy_dol.
        # Note: we can also run legacy_dol() function but doing it
        # this way captures a test of the command line args too.
        sys.argv = argv
        runpy.run_module("gdt.missions.fermi.gbm.localization.dol.legacy_dol", args, '__main__')

        # restore original argv
        sys.argv = orig_argv

        # default uses same file names as current output for comparison
        if reference is None:
            reference = output

        # check for differences and missing files
        bad = []
        for i in range(len(output)):
            path1 = os.path.join(self.tmp_dir, output[i])
            path2 = os.path.join(__test_ref_dir__, reference[i])
            bad.append(self.diff(path1, path2))

        # return False if bad things are found
        if np.any(bad):
            return False
        return True

    def test_243216766(self):
        args = "243216766.000.npy"
        output = ["dol_243216766.000.txt", "chi2grid_bn243216766.000_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_387829736(self):
        args = "387829736.000.npy"
        output = ["dol_387829736.000.txt", "chi2grid_bn387829736.000_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_448700543(self):
        args = "448700543.392.npy"
        output = ["dol_448700543.392.txt", "chi2grid_bn448700543.392_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_448700544(self):
        r""" This tests a case where we know there is a difference """
        args = "448700544.392.npy"
        output = ["dol_448700544.392.txt", "chi2grid_bn448700544.392_v00.dat"]
        self.assertFalse(self.compare(args, output))

    def test_469700577(self):
        args = "469700577.000.npy"
        output = ["dol_469700577.000.txt", "chi2grid_bn469700577.000_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_566723806(self):
        args = "566723806.604.npy"
        output = ["dol_566723806.604.txt", "chi2grid_bn566723806.604_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_573748874(self):
        args = "573748874.438.npy"
        output = ["dol_573748874.438.txt", "chi2grid_bn573748874.438_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_575442156(self):
        args = "575442156.911.npy"
        output = ["dol_575442156.911.txt", "chi2grid_bn575442156.911_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_576076018(self):
        args = "576076018.283.npy"
        output = ["dol_576076018.283.txt", "chi2grid_bn576076018.283_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_576160056(self):
        args = "576160056.977.npy"
        output = ["dol_576160056.977.txt", "chi2grid_bn576160056.977_v00.dat"]
        self.assertTrue(self.compare(args, output))

    def test_to_GbmHealPix(self):

        args_path = os.path.join(__test_args_dir__, "576160056.977.npy")
        args = args = np.load(args_path, allow_pickle=True, encoding='bytes').item()

        from gdt.missions.fermi.gbm.localization.dol.legacy_dol import legacy_DoL

        dol = legacy_DoL()
        loc = dol.eval(
            args['crange'], args['mrates'], args['brates'],
            args['sduration'], args['bgduration'], args['sc_pos'], args['sc_quat'],
            args['energies'], args['fra'], args['fdec'], args['sc_time'], scat_opt=1)

        frame = FermiFrame(
            obstime=Time(args['sc_time'], format='fermi'),
            obsgeoloc=CartesianRepresentation(*args['sc_pos'], unit='km'),
            quaternion=Quaternion(args['sc_quat'], scalar_first=False))

        healpix = dol.to_GbmHealPix(loc, frame)

        for i, p in [(178067, 3.4297927e-06),
                     ( 32104, 5.1333986e-06),
                     ( 77421, 4.1037010e-06),
                     (193711, 3.3427739e-06),
                     (113747, 3.4839798e-06)]:
            self.assertTrue(np.isclose(healpix.prob[i], p, rtol=0.01))

# END class TestLegacyDol

if __name__ == "__main__":
    unittest.main()
