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

from gdt.missions.fermi.time import *
from gdt.missions.fermi.gbm.finders import *

download_dir = data_dir = os.path.dirname(os.path.abspath(__file__))


class TestTriggerFinder(unittest.TestCase):
    
    def setUp(self):
        self.finder = TriggerFinder()
    
    def test_set_trigger(self):
        self.finder.cd('080916009')
        self.assertEqual(self.finder.num_files, 109)
        self.finder.cd('170817529')
        self.assertEqual(self.finder.num_files, 128)
    
    def test_ls(self):
        self.finder.cd('170817529')
        [self.assertTrue('ctime' in file) for file in self.finder.ls_ctime()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_cspec()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_rsp(ctime=False)]
        [self.assertTrue('ctime' in file) for file in self.finder.ls_rsp(cspec=False)]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_rsp2(ctime=False)]
        [self.assertTrue('ctime' in file) for file in self.finder.ls_rsp2(cspec=False)]
        [self.assertTrue('lc' in file) for file in self.finder.ls_lightcurve()]
        [self.assertTrue('.fit' in file) for file in self.finder.ls_cat_files()]
        self.assertTrue('trigdat' in self.finder.ls_trigdat()[0])
        [self.assertTrue(('healpix' in file) or ('skymap' in file) or 
                         ('loclist' in file) or ('locprob' in file) or
                         ('locplot' in file)) for file in self.finder.ls_localization()]
            
    def test_get(self):
        self.finder.cd('170817529')
        self.finder.get_cat_files(download_dir)
        cat_files = self.finder.ls_cat_files()
        [os.remove(os.path.join(download_dir, file)) for file in cat_files]


class TestContinuousFinder(unittest.TestCase):
    def setUp(self):
        self.finder = ContinuousFinder()
    
    def test_set_time(self):
        self.finder.cd(Time(604741251.0, format='fermi'))
        self.assertEqual(self.finder.num_files, 379)
    
    def test_ls(self):
        self.finder.cd(Time(604741251.0, format='fermi'))
        [self.assertTrue('ctime' in file) for file in self.finder.ls_ctime()]
        [self.assertTrue('cspec' in file) for file in self.finder.ls_cspec()]
        [self.assertTrue('poshist' in file) for file in self.finder.ls_poshist()]
        [self.assertTrue('spechist' in file) for file in self.finder.ls_spechist()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte()]
        [self.assertTrue('tte' in file) for file in self.finder.ls_tte(full_day=True)]

    def test_get(self):
        self.finder.cd(Time(604741251.0, format='fermi'))
        self.finder.get_poshist(download_dir)
        self.finder.get_ctime(download_dir, dets=('n0', 'n1', 'n2'))
        
        files = self.finder.ls_poshist()
        files.extend(self.finder.ls_ctime())
        for file in files:
            try:
                os.remove(os.path.join(download_dir, file))
            except:
                pass
    
    def test_reconnect(self):
        finder = ContinuousFinder()
        finder = ContinuousFinder()
        self.finder.cd(Time(604741251.0, format='fermi'))


@unittest.skipIf(
    os.environ.get('SKIP_HEASARC_FTP_TESTS', False), 'Skipping HEASARC FTP tests'
)
class TestFtpFinders(unittest.TestCase):
    def test_trigger(self):
        finder = TriggerFtp()
        finder.cd('080916009')
        self.assertEqual(finder.num_files, 109)
        finder.cd('170817529')
        self.assertEqual(finder.num_files, 128)

    def test_continuous(self):
        finder = ContinuousFtp()
        finder.cd(Time(604741251.0, format='fermi'))
        self.assertEqual(finder.num_files, 379)
