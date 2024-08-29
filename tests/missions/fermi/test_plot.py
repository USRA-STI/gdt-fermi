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
import pytest
import numpy as np
import unittest
import matplotlib.colors as colors

from gdt.core import data_path
from gdt.missions.fermi.plot import *
from gdt.missions.fermi.gbm.saa import GbmSaaPolygon1
from gdt.missions.fermi.gbm.poshist import GbmPosHist


class TestFermiEarthPlot(unittest.TestCase):

    def test_mcilwainl(self):

        earth_plot = FermiEarthPlot()
        mcilwainl = earth_plot.mcilwainl

        self.assertTrue(mcilwainl.colorbar.vmax, 1.7)

        mcilwainl.alpha = 0.8
        self.assertEqual(mcilwainl.alpha, 0.8)
        self.assertTrue('viridis' in str(mcilwainl))

        norm = colors.Normalize(vmin=-100, vmax=100)
        mcilwainl.norm = norm
        self.assertEqual(mcilwainl.norm.vmax, norm.vmax)

        with self.assertRaises(TypeError):
            mcilwainl.color = 'bad'

    def test_saa_plot(self):

        saa = GbmSaaPolygon1()
        earth_plot = FermiEarthPlot(saa)

    def test_spacecraft_plot(self):
        test_file = data_path / 'fermi-gbm/glg_poshist_all_170101_v01.fit'
        if not test_file.exists():
            pytest.skip("test files aren't downloaded. run gdt-data download fermi-gbm.")

        frames = GbmPosHist.open(test_file).get_spacecraft_frame()
        start, center, stop = [frames[i].obstime for i in [0, 2500, 5000]]
        earth_plot = FermiEarthPlot()
        earth_plot.add_spacecraft_frame(frames, start, stop, trigtime=center)
        earth_plot.standard_title()
        self.assertEqual(
            earth_plot._m.title.get_text(), 
            'Latitude, East Longitude: (24.45S, 120.63W)\nMcIlwain L: 1.21')
