# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
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

from gdt.core import data_path
from gdt.missions.fermi.gbm.finders import TriggerFinder
from gdt.missions.fermi.gbm.trigdat import Trigdat
from gdt.missions.fermi.gbm.classification import Classification

trigdat_file = data_path / 'fermi-gbm/glg_trigdat_all_bn170101116_v01.fit'

@unittest.skipIf(not trigdat_file.exists(), "test files aren't downloaded. run 'gdt-data download fermi-gbm'.")
class TestClassifier(unittest.TestCase):
    
    def setUp(self):
        self.trigdat = Trigdat.open(trigdat_file)
   
    def tearDown(self):
        self.trigdat.close()

    def test_with_loc(self):

        c = Classification()
        det_nums = c.classify_trigdat(self.trigdat, verbose=True)

        self.assertTrue(np.all(det_nums == [9, 10]))
        for i, (name, prob) in enumerate([("GRB", 0.863), ("GROJ_0422_32", 0.096)]):
            self.assertTrue(c.rankings()[i][0].name == name)
            self.assertTrue(np.isclose(c.rankings()[i][1], prob, atol=0.01))
        
    def test_without_loc(self):

        c = Classification()
        det_nums = c.classify_trigdat(self.trigdat, use_loc=False)

        self.assertTrue(np.all(det_nums == [10, 9]))
        for i, (name, prob) in enumerate([("GRB", 0.863), ("GROJ_0422_32", 0.096)]):
            self.assertTrue(c.rankings()[i][0].name == name)
            self.assertTrue(np.isclose(c.rankings()[i][1], prob, atol=0.01))

    def test_tgf(self):

        path = TriggerFinder('241220239').get_trigdat(".")[0]
        trigdat = Trigdat.open(path)

        c = Classification()
        c.classify_trigdat(trigdat)

        self.assertTrue(c.rankings()[0][0].name == 'TGF')
        self.assertTrue(c.rankings()[0][1], 1.0)

        trigdat.close()
        try:
            os.remove(path)
        except:
            pass

    def test_below_horizon(self):

        path = TriggerFinder('240831759').get_trigdat(".")[0]
        trigdat = Trigdat.open(path)

        c = Classification()
        c.classify_trigdat(trigdat)

        self.assertTrue(c.rankings()[0][0].name == 'BELOW_HORIZON')
        self.assertTrue(c.rankings()[0][1], 1.0)

        trigdat.close()
        try:
            os.remove(path)
        except:
            pass

    def test_distant_particles(self):

        path = TriggerFinder('240828632').get_trigdat(".")[0]
        trigdat = Trigdat.open(path)

        c = Classification()
        c.classify_trigdat(trigdat)

        self.assertTrue(c.rankings()[0][0].name == 'DISTANT_PARTICLES')
        self.assertTrue(c.rankings()[0][1], 0.892)

        trigdat.close()
        try:
            os.remove(path)
        except:
            pass

    def test_local_particles(self):

        path = TriggerFinder('211104116').get_trigdat(".")[0]
        trigdat = Trigdat.open(path)

        c = Classification()
        c.classify_trigdat(trigdat)

        self.assertTrue(c.rankings()[0][0].name == 'LOCAL_PARTICLES')
        self.assertTrue(c.rankings()[0][1], 1.0)

        trigdat.close()
        try:
            os.remove(path)
        except:
            pass

    def test_repr(self):

        c = Classification()
        c.classify_trigdat(self.trigdat)
        self.assertEqual(
            str(c),
            "Classification(\n"
            "               ERROR               0.0000\n"
            "               UNRELIABLE_LOCATION 0.0000\n"
            "               LOCAL_PARTICLES     0.0000\n"
            "               BELOW_HORIZON       0.0000\n"
            "               GRB                 0.8633\n"
            "               GENERIC_SGR         0.0044\n"
            "               GENERIC_TRANSIENT   0.0359\n"
            "               DISTANT_PARTICLES   0.0000\n"
            "               SOLAR_FLARE         0.0000\n"
            "               CYG_X1              0.0000\n"
            "               SGR_1806_20         0.0000\n"
            "               GROJ_0422_32        0.0964\n"
            "               RESERVED_1          0.0000\n"
            "               RESERVED_2          0.0000\n"
            "               RESERVED_3          0.0000\n"
            "               RESERVED_4          0.0000\n"
            "               RESERVED_5          0.0000\n"
            "               RESERVED_6          0.0000\n"
            "               RESERVED_7          0.0000\n"
            "               TGF                 0.0000)")

    def test_repr_html(self):

        c = Classification()
        c.classify_trigdat(self.trigdat)
        self.assertEqual(
            c._repr_html_(),
            "<p>Classification:</p><table>"
            "<tr><th>classification</th><th>probability</th></tr>"
            "<tr><td>ERROR</td><td>0.0000</td></tr>"
            "<tr><td>UNRELIABLE_LOCATION</td><td>0.0000</td></tr>"
            "<tr><td>LOCAL_PARTICLES</td><td>0.0000</td></tr>"
            "<tr><td>BELOW_HORIZON</td><td>0.0000</td></tr>"
            "<tr><td>GRB</td><td>0.8633</td></tr>"
            "<tr><td>GENERIC_SGR</td><td>0.0044</td></tr>"
            "<tr><td>GENERIC_TRANSIENT</td><td>0.0359</td></tr>"
            "<tr><td>DISTANT_PARTICLES</td><td>0.0000</td></tr>"
            "<tr><td>SOLAR_FLARE</td><td>0.0000</td></tr>"
            "<tr><td>CYG_X1</td><td>0.0000</td></tr>"
            "<tr><td>SGR_1806_20</td><td>0.0000</td></tr>"
            "<tr><td>GROJ_0422_32</td><td>0.0964</td></tr>"
            "<tr><td>RESERVED_1</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_2</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_3</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_4</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_5</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_6</td><td>0.0000</td></tr>"
            "<tr><td>RESERVED_7</td><td>0.0000</td></tr>"
            "<tr><td>TGF</td><td>0.0000</td></tr></table>")
