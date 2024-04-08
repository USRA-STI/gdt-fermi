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

import datetime
import re

import erfa
import numpy as np
from astropy.time import TimeFromEpoch, TimeUnique, ScaleValueError, Time
from astropy.time.utils import day_frac

__all__ = ['FermiSecTime', 'GbmBurstNumber', 'Time']


class FermiSecTime(TimeFromEpoch):
    """Represents the number of seconds elapsed since Jan 1, 2001, 00:00:00 UTC,
    including leap seconds
    """
    name = 'fermi'
    """(str): Name of the mission"""

    unit = 1.0 / 86400
    """(float): unit in days"""

    epoch_val = '2001-01-01 00:01:04.184'
    """(str): The epoch in Terrestrial Time"""

    epoch_val2 = None

    epoch_scale = 'tt'
    """(str): The scale of :attr:`epoch_val`"""

    epoch_format = 'iso'
    """(str): Format of :attr:`epoch_val`"""


class GbmBurstNumber(TimeUnique):
    """Represent date as Fermi GBM burst number"""

    name = 'gbm_bn'
    """(str): The name of the time format"""

    _bn_pattern = re.compile(r'^(?P<year>\d\d)(?P<month>\d\d)(?P<day>\d\d)(?P<frac>\d\d\d)$', re.I | re.S)

    _second_corrections = [
        # 2008-01-01 00:00:00.000 UTC < x < 2009-01-01 00:00:00.000 UTC adjust seconds by 1
        (datetime.datetime(2008, 1, 1), datetime.datetime(2009, 1, 1), 1),
        # 2009-01-01 00:00:00.000 UTC < met <= 2009-01-13 00:00:00.000 UTC adjust seconds by 2
        (datetime.datetime(2009, 1, 1), datetime.datetime(2009, 1, 13), 2)
    ]

    def _check_val_type(self, val1, val2):
        if not all(self._bn_pattern.match(val) is not None for val in val1.flat):
            raise TypeError('Input values for {} class must be '
                            'burst numbers'.format(self.name))
        if val2 is not None:
            raise ValueError(
                f'{self.name} objects do not accept a val2 but you provided {val2}')
        return val1, None

    def set_jds(self, val1, val2):
        """Convert burst number contained in val1 to jd1, jd2"""
        # Iterate through the datetime objects, getting year, month, etc.
        iterator = np.nditer([val1, None, None, None, None, None, None],
                             flags=['refs_ok', 'zerosize_ok'],
                             op_dtypes=[None] + 5 * [np.intc] + [np.double])

        for val, iy, im, iday, ihr, imin, dsec in iterator:
            m = self._bn_pattern.match(str(val))
            if not m:
                raise ValueError('burst number was not the correct format.')

            bn = m.groupdict()
            iy[...] = int(bn['year']) + 2000
            im[...] = int(bn['month'])
            iday[...] = int(bn['day'])

            secs_of_day = int(bn['frac']) * 86.400

            ihr[...] = secs_of_day // 3600
            secs_of_day %= 3600
            imin[...] = secs_of_day // 60
            dsec[...] = secs_of_day % 60

        jd1, jd2 = erfa.dtf2d(self.scale.upper().encode('ascii'),
                              *iterator.operands[1:])
        self.jd1, self.jd2 = day_frac(jd1, jd2)

    def to_value(self, parent=None, out_subfmt=None):
        """Convert to string representation of GBM burst number"""
        if out_subfmt is not None:
            # Out_subfmt not allowed for this format, so raise the standard
            # exception by trying to validate the value.
            self._select_subfmts(out_subfmt)

        # GBM Burst Numbers are in the UTC scale, so make sure that we use UTC
        if self.scale != 'utc':
            if parent is None:
                raise ValueError('cannot compute value without parent Time object')
            try:
                tm = getattr(parent, 'utc')
            except Exception as err:
                raise ScaleValueError("Cannot convert from '{}' to UTC, got error:\n{}"
                                      .format(self.name, self.scale, err)) from err

            jd1, jd2 = tm._time.jd1, tm._time.jd2
        else:
            jd1, jd2 = self.jd1, self.jd2

        # Rather than define a value property directly, we have a function,
        # since we want to be able to pass in timezone information.
        scale = self.scale.upper().encode('ascii')
        iys, ims, ids, ihmsfs = erfa.d2dtf(scale, 6,  # 6 for microsec
                                           jd1, jd2)
        ihrs = ihmsfs['h']
        imins = ihmsfs['m']
        isecs = ihmsfs['s']
        ifracs = ihmsfs['f']
        iterator = np.nditer([iys, ims, ids, ihrs, imins, isecs, ifracs, None],
                             flags=['refs_ok', 'zerosize_ok'],
                             op_dtypes=7 * [None] + [object])

        for iy, im, iday, ihr, imin, isec, ifracsec, out in iterator:
            day_secs = ihr * 3600 + imin * 60 + isec

            # Adjust results to match the burst numbers issued early in the mission
            dt = datetime.datetime(iy, im, iday, ihr, imin, isec)
            for beg_dt, end_dt, adj_secs in self._second_corrections:
                if beg_dt < dt <= end_dt:
                    # print(f'{beg_dt} < {dt} <= {end_dt}')
                    day_secs += adj_secs
            frac = min(int(round(day_secs / 86400 * 1000)), 999)
            out[...] = f'{iy - 2000:02d}{im:02d}{iday:02d}{frac:03d}'

        return iterator.operands[-1]

    value = property(to_value)
