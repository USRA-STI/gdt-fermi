# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Developed by: Jacob Smith
#               University of Alabama in Huntsville
#               Center for Space Plasma and Aeronomic Research
#
# This software is not subject to EAR.
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

from astropy import units as u

from gdt.missions.fermi.gbm.tte import GbmTte
from gdt.missions.fermi.gbm.finders import TriggerFinder
from gdt.missions.fermi.gbm.trigger import GbmOnboardTrigger

download_dir = os.path.dirname(os.path.abspath(__file__))


class TestGbmTrigger(unittest.TestCase):

    @classmethod 
    def setUpClass(cls):

        # download TTE data
        finder = TriggerFinder('241212390')
        dets = ['n%d' % i for i in range(10)] + ['na', 'nb']
        files = finder.get_tte(dets=dets, download_dir=download_dir)
        ttes = [GbmTte.open(file) for file in files]

        # initialize onboard trigger class
        cls.onboard_trig = GbmOnboardTrigger()
        cls.onboard_trig.prepare_data(ttes)

        # cleanup TTE data
        for file in files:
            file.unlink(missing_ok=True)

        # create a set of triggers for use with plotting methods
        cls.triggers = cls.onboard_trig.apply_trigger()

    def test_apply_trigger(self):

        triggers = self.onboard_trig.apply_trigger()
        self.assertEqual(len(triggers), 30)
        self.assertEqual(triggers[0].alg_num, 14)
        self.assertEqual(triggers[0].alg.timescale, 2048)
        self.assertEqual(triggers[0].time, 0.0)
        self.assertEqual(triggers[0].triggered_det_names, ['n9', 'na'])

        # check again with holdoff and debug options
        triggers_with_holdoff = self.onboard_trig.apply_trigger(holdoff=300, debug=True)
        self.assertEqual(len(triggers_with_holdoff), 1)
        self.assertEqual(triggers[0].alg_num, triggers_with_holdoff[0].alg_num)
        self.assertEqual(triggers[0].alg.timescale, triggers_with_holdoff[0].alg.timescale)
        self.assertEqual(triggers[0].time, triggers_with_holdoff[0].time)
        self.assertEqual(triggers[0].triggered_det_names, triggers_with_holdoff[0].triggered_det_names)

    def test_exec_errors(self):

        finder = TriggerFinder('241212390')
        dets = ['n0', 'b0']
        files = finder.get_tte(dets=dets, download_dir=download_dir)
        ttes = [GbmTte.open(file) for file in files]

        onboard_trig = GbmOnboardTrigger()
        # tte data includes BGO
        with self.assertRaises(ValueError):
            onboard_trig.prepare_data(ttes)
        # tte data includes duplicate NAI
        with self.assertRaises(ValueError):
            onboard_trig.prepare_data([ttes[1], ttes[1]])
        # tte data without all detectors present
        with self.assertRaises(ValueError):
            onboard_trig.prepare_data([ttes[1]])

        # cleanup TTE files
        for file in files:
            file.unlink(missing_ok=True)
