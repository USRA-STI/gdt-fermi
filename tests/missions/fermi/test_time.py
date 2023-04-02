"""test_time.py: Tests the Fermi conversion routines"""

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

import pytest
from pathlib import Path
from gdt.missions.fermi.time import Time


@pytest.fixture
def trigger_times():
    retval = []
    test_file = Path(__file__).parent.joinpath('data', 'trigger_data.csv')
    with test_file.open() as fp:
        while True:
            line = fp.readline().strip()
            if not line:
                break
            met, bn, iso = line.split(',')
            retval.append({
                'met': float(met),
                'bn': bn.strip(),
                'iso': iso.strip()
            })
    return retval


def test_convert_met_to_utc_datetime(trigger_times):
    for trigger in trigger_times:
        t = Time(trigger['met'], format='fermi')
        # Have to convert to datetime, since astropy only supports fractional seconds to 3 decimal places
        assert t.utc.datetime.strftime('%Y-%m-%d %H:%M:%S.%f') == trigger['iso']


def test_convert_met_to_burst_number(trigger_times):
    for trigger in trigger_times:
        t = Time(trigger['met'], format='fermi')
        assert t.gbm_bn == trigger['bn']


def test_convert_utc_to_met(trigger_times):
    for trigger in trigger_times:
        t = Time(trigger['iso'], format='iso', scale='utc')
        assert t.fermi == trigger['met']
