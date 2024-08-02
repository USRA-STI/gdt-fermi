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
import numpy as np

from ..time import Time
from gdt.core.geomagnetic import SouthAtlanticAnomaly
from gdt.core.data_primitives import Range

__all__ = ['GbmSaaPolygon1', 'GbmSaaPolygon2', 'GbmSaa', 'GbmSaaCollection']

class GbmSaaPolygon1(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
     used from launch (0 seconds) until 19:15 UTC on July 23, 2024.
    """
    _latitude = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                 -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
    _longitude = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
                  -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]

    _time_range = Range(Time(0, format='fermi'),
                        Time(743454905, format='fermi'))


class GbmSaaPolygon2(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 19:15 UTC on July 23, 2024 until mission end which uses a
    placeholder value of Jan 1, 2050 for now.
    """
    _latitude = [-24.395, -30.000, -30.000, -30.000, -30.000, -24.060,
                 -16.220, -8.638, -6.155, -1.000, 2.000, 2.000,
                 -3.400, -19.570802973653127, -24.395]
    _longitude = [22.000, 22.000, 0.000, -2.000, -86.100, -90.300,
                  -90.300, -88.738, -84.000, -65.000, -45.000, -38.400,
                  -30.605, -11.999457706441582, 22.000]

    _time_range = Range(Time(743454905, format='fermi'),
                        Time(1546300805, format='fermi'))


class GbmSaa(GbmSaaPolygon1):
    """Class providing backwards compatibility for code written prior to v2.1.1
    where the GbmSaa class defined GbmSaaPolygon1"""
    pass


class GbmSaaCollection():
    """Collection of SAA boundary definitions"""

    def __init__(self):
        """Constructor"""
        self._polygons = [GbmSaaPolygon1(), GbmSaaPolygon2()]
        self.sort()

    def at(self, time: Time):
        """Return the active SAA polygon for a given time.

        Args:
            time (astropy.time.Time): time for determining the set
                                      of SAA points in latitude/longitude

        Returns:
            (:class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`)
        """
        for polygon in self.polygons:
            if polygon._time_range.contains(time):
                # return the first polygon that contains this time
                return polygon

        raise ValueError(f"This collection of SAA polygons does not contain time {time}")

    @property
    def num_polygons(self):
        """Number of polygons included in the collection

        Returns:
            (int)
        """
        return len(self.polygons)

    @property
    def polygons(self):
        """Returns a list of polygons included in the collection

        Returns:
            (list)
        """
        return self._polygons

    def sort(self):
        """Sorts SAA polygons according to time periods"""
        high = [p._time_range._high for p in self.polygons]

        sorted_polygons = []
        for i in np.argsort(high):
            low, high = self.polygons[i]._time_range.as_tuple()

            # sanity check to ensure the time range of this polygon
            # does not overlap with another polygon. Note that the
            # low end of the time range can equal the high end of
            # the previous time range
            for p in sorted_polygons:
                if (p._time_range.contains(low) and low != p._time_range._high) \
                    or p._time_range.contains(high):
                    raise ValueError("Polygon time ranges overlap.")

            sorted_polygons.append(self.polygons[i])

        self._polygons = sorted_polygons
