#  CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
#  Contract No.: CA 80MSFC17M0022
#  Contractor Name: Universities Space Research Association
#  Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
#  Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
#  Developed by: William Cleveland and Adam Goldstein
#                Universities Space Research Association
#                Science and Technology Institute
#                https://sti.usra.edu
#
#  Developed by: Daniel Kocevski
#                National Aeronautics and Space Administration (NASA)
#                Marshall Space Flight Center
#                Astrophysics Branch (ST-12)
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#   in compliance with the License. You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software distributed under the License
#  is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#  implied. See the License for the specific language governing permissions and limitations under the
#  License.
#
import pytest
import numpy as np
from gdt.core import data_path
from gdt.core.coords import Quaternion
from gdt.missions.fermi.gbm.poshist import GbmPosHist


@pytest.fixture
def test_file():
    return data_path / 'fermi-gbm/glg_poshist_all_170101_v01.fit'


def test_get_spacecraft_frame(test_file):
    if not test_file.exists():
        pytest.skip("test files aren't downloaded. run gdt-download-data.")

    with GbmPosHist.open(test_file) as poshist:
        frame = poshist.get_spacecraft_frame()

        # This file has 86520 rows
        assert frame.obstime.size == 86520

        # Let's verify the first row
        pos = frame[0]
        assert pos.obstime.fermi == 504921545.740106
        assert pos.quaternion == Quaternion((-0.21889445774812064, 0.00818763828216734, 0.6516071490237402,
                                             -0.7262412149469831))

        assert np.all(pos.obsgeoloc.xyz.value == np.array((-6315714.0, -1549390.6, 2302684.5), dtype=np.float32))
        assert str(pos.obsgeoloc.xyz.unit) == 'm'

        assert np.all(pos.obsgeovel.xyz.value == np.array((1011.3916, -7244.9014, -2100.9446), dtype=np.float32))
        assert str(pos.obsgeovel.xyz.unit) == 'm / s'

        # Let's verify the last row
        pos = frame[-1]
        assert pos.obstime.fermi == 505008064.340076
        assert pos.quaternion == Quaternion((0.8262526368132503, -0.4260856438475033, -0.18342245964476556,
                                             0.31955250830508913))

        assert np.all(pos.obsgeoloc.xyz.value == np.array((-3306917.0, -6046107.0, -393885.53), dtype=np.float32))
        assert str(pos.obsgeoloc.xyz.unit) == 'm'

        assert np.all(pos.obsgeovel.xyz.value == np.array((6111.002, -3138.5576, -3268.467), dtype=np.float32))
        assert str(pos.obsgeovel.xyz.unit) == 'm / s'


def test_get_spacecraft_states(test_file):
    if not test_file.exists():
        pytest.skip("test files aren't downloaded. run gdt-download-data.")

    with GbmPosHist.open(test_file) as poshist:
        states = poshist.get_spacecraft_states()

        # This file has 86520 rows
        assert len(states) == 86520

        # First row
        # Value of flag = 1 which means "in sun" and "not in SAA"
        state = states[0]
        assert state['saa'] == False
        assert state['sun'] == True

        # Last row
        # Value of flag = 1 which means "in sun" and "not in SAA"
        state = states[-1]
        assert state['saa'] == False
        assert state['sun'] == True

        # Row 3683
        # Value of flag = 0 which means "not in sun" and "not in SAA"
        state = states[3682]
        assert state['saa'] == False
        assert state['sun'] == False

        # Row 3309
        # Value of flag = 2 which means "not in sun" and "in SAA"
        state = states[3308]
        assert state['saa'] == True
        assert state['sun'] == False

        # Row 2930
        # Value of flag = 3 which means "in sun" and "in SAA"
        state = states[2929]
        assert state['saa'] == True
        assert state['sun'] == True
