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

__all__ = ['GbmSaaPolygon1', 'GbmSaaPolygon2', 'GbmSaaPolygon3',
           'GbmSaaPolygon4', 'GbmSaaPolygon5', 'GbmSaaPolygon6',
           'GbmSaaPolygon7', 'GbmSaaPolygon8', 'GbmSaaPolygon9',
           'GbmSaaPolygon10', 'GbmSaaPolygon11', 'GbmSaaPolygon12',
           'GbmSaaPolygon13', 'GbmSaaPolygon14', 'GbmSaaPolygon15',
           'GbmSaa', 'GbmSaaCollection']


class GbmSaaPolygon1(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from launch (0 seconds) until 18:30:09 UTC on June 9, 2010. This
    is the original SAA polygon.
    """
    _latitude = [-30.000, -22.600, 2.500, 5.200, 5.200, 4.600,
                 0.700, -8.600, -9.900, -12.500, -21.700, -30.000, -30.000]
    _longitude = [33.900,  24.500, -18.600, -25.700, -36.000, -42.000,
                  -58.800, -93.100, -97.500, -98.500, -92.100, -86.100,  33.900]

    _time_range = Range(Time(0, format='fermi'),
                        Time(297801011, format='fermi'))


class GbmSaaPolygon2(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 18:30:09 UTC on June 9, 2010 until 17:41:58 UTC on June 10, 2010.

    Note: This smaller region was proposed by Jerry Fishman to increase
    observing time. It resulted in too many local particle triggers so it
    was adjusted via GbmSaaPolygon3 and GbmSaaPolygon4.
    """
    _latitude = [-30.000, -21.670, -13.330, -5.000, -6.000, -7.000,
                 -9.000, -11.000, -14.000, -17.000, -23.500, -30.000]
    _longitude = [21.000, 0.000, -21.000, -42.000, -50.000, -58.000,
                  -65.000, -72.000, -76.000, -80.000, -80.000, -80.000]

    _time_range = Range(Time(297801011, format='fermi'),
                        Time(297884520, format='fermi'))


class GbmSaaPolygon3(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 17:41:58 UTC on June 10, 2010 until 13:34:00 UTC on July 8, 2010.

    Note: This is a slightly larger area version of GbmSaaPolygon2 to
    remove local particle triggers occurring near SAA entry/exit.
    """
    _latitude = [-30.000, -14.500, -9.300, 1.000, 2.000, 3.000,
                 -2.000, -7.000, -11.000, -15.000, -22.500, -30.000]
    _longitude = [26.000, -2.000, -11.300, -30.000, -37.500, -45.000,
                  -62.500, -80.000, -85.000, -90.000, -88.500, -87.000]

    _time_range = Range(Time(297884520, format='fermi'),
                        Time(300288842, format='fermi'))


class GbmSaaPolygon4(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 13:34:00 UTC on July 8, 2010 until 14:13:00 UTC on July 20, 2010.

    Note: This is a slightly larger area version of GbmSaaPolygon3 to
    remove local particle triggers occurring near SAA entry/exit.
    """
    _latitude = [-30.000, -2.000, -0.500, 0.500, 1.000, 1.000,
                 0.000, -2.000, -7.500, -11.500, -16.500, -30.000]
    _longitude = [40.000, -27.000, -35.000, -47.000, -60.000, -65.000,
                  -75.000, -82.000, -90.000, -94.000, -97.000, -92.000]

    _time_range = Range(Time(300288842, format='fermi'),
                        Time(301327982, format='fermi'))


class GbmSaaPolygon5(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 14:13:00 UTC on July 20, 2010 until 14:20:00 UTC on June 7, 2022.

    Note: Revised polygon proposed by Jerry Fishman, Michael Briggs, and
    Vandiver Chaplin. This is smaller area than the original GbmSaaPolygon1
    but avoids issues with local particle triggers seen in GbmSaaPolygon2-4.
    """
    _latitude = [-30.000, -19.867, -9.733, 0.400, 2.000, 2.000,
                 -1.000, -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
    _longitude = [33.900, 12.398, -9.103, -30.605, -38.400, -45.000,
                  -65.000, -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]

    _time_range = Range(Time(301327982, format='fermi'),
                        Time(676304405, format='fermi'))


class GbmSaaPolygon6(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 14:20:00 UTC on June 7, 2022 until 19:35:00 UTC on June 9, 2022.

    Note: Several day test of a new, smaller area polygon proposed by
    Michael Briggs.
    """
    _latitude = [-30.000, -19.867, -9.733, -3.400, -2.000, -2.000,
                 -5.000, -10.155, -12.880, -16.220, -22.404, -30.000]
    _longitude = [-5.000, -16.398, -28.103, -35.605, -38.400, -45.000,
                  -65.000, -84.000, -89.200, -90.300, -90.300, -86.100]

    _time_range = Range(Time(676304405, format='fermi'),
                        Time(676496105, format='fermi'))


class GbmSaaPolygon7(GbmSaaPolygon5):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 19:35:00 UTC on June 9, 2022 until 18:10:45 UTC on August 3, 2022.

    Note: Reversion to the same shape as GbmSaaPolygon5 following
    the test of GbmSaaPolygon6.
    """
    _time_range = Range(Time(676496105, format='fermi'),
                        Time(681243050, format='fermi'))


class GbmSaaPolygon8(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 18:10:45 UTC on August 3, 2022 until 17:59:26 UTC on August 11, 2022.

    Note: Several day test of a new, smaller area polygon proposed by
    Michael Briggs.
    """
    _latitude = [-30.000, -19.867, -9.733, -3.400, 0.000, 0.000,
                 -3.000, -8.155, -10.880, -16.220, -22.404, -30.000, -30.000]
    _longitude = [0.000, -11.398, -23.103, -30.605, -36.400, -43.000,
                  -65.000, -84.000, -89.200, -90.300, -90.300, -86.100, 0.000]

    _time_range = Range(Time(681243050, format='fermi'),
                        Time(681933571, format='fermi'))


class GbmSaaPolygon9(GbmSaaPolygon5):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 17:59:26 UTC on August 11, 2022 until 16:50:00 UTC on November 15, 2023.

    Note: Reversion to the same shape as GbmSaaPolygon5 following
    the test of GbmSaaPolygon8.
    """
    _time_range = Range(Time(681933571, format='fermi'),
                        Time(721759805, format='fermi'))


class GbmSaaPolygon10(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 16:50:00 UTC on November 15, 2023 until 18:27:00 UTC on November 17, 2023.

    Note: Several day test of a new, smaller area polygon proposed by
    Michael Briggs.
    """
    _latitude = [-24.395, -30.000, -30.000, -30.000, -30.000, -24.060,
                 -16.220, -8.638, -6.155, -1.000, 2.000, 2.000,
                 -3.400, -19.570802973653127, -24.395]
    _longitude = [22.000, 22.000, 0.000, -2.000, -86.100, -90.300,
                  -90.300, -88.738, -84.000, -65.000, -45.000, -38.400,
                  -30.605, -11.999457706441582, 22.000]

    _time_range = Range(Time(721759805, format='fermi'),
                        Time(721938425, format='fermi'))


class GbmSaaPolygon11(GbmSaaPolygon5):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 18:27:00 UTC on November 17, 2023 until 19:15:00 UTC on July 23, 2024.

    Note: Reversion to the same shape as GbmSaaPolygon5 following
    the test of GbmSaaPolygon10.
    """
    _time_range = Range(Time(721938425, format='fermi'),
                        Time(743454905, format='fermi'))


class GbmSaaPolygon12(GbmSaaPolygon10):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 19:15 UTC on July 23, 2024 until 19:15:00 UTC on August 6, 2024.

    Note: Initial implementation of the new, smaller area polygon proposed by
    Michael Briggs.
    """
    _time_range = Range(Time(743454905, format='fermi'),
                        Time(744664505, format='fermi'))


class GbmSaaPolygon13(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 19:15:00 UTC on August 6, 2024 until 00:00:00 UTC on August 8, 2024.

    Note: Minor adjustment to GbmSaaPolygon12 which reduces the number of local
    particle triggers occurring during SAA entry (western side of polygon).
    """
    _latitude = [-24.395, -30.000, -30.000, -30.000, -30.000, -30.000,
                 -30.000, -22.646, -16.220, -9.404, -6.155, -1.000,
                 2.000, 2.000, -3.650, -20.136067618721196, -24.395]
    _longitude = [22.000, 22.000, 5.500, 0.000, -2.000, -40.000,
                  -86.100, -91.300, -91.800, -89.700, -84.000, -65.000,
                  -45.000, -38.400, -30.605, -8.015646247668737, 22.000]

    _time_range = Range(Time(744664505, format='fermi'),
                        Time(744768005, format='fermi'))


class GbmSaaPolygon14(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 00:00:00 UTC on August 8, 2024 until 19:20:45 UTC on September 30, 2024.

    Note: Minor adjustment to GbmSaaPolygon13 which reduces the number of local
    particle triggers occurring during SAA exit (eastern side of polygon).
    """
    _latitude = [-24.395, -30.000, -30.000, -30.000, -30.000, -30.000,
                 -30.000, -22.646, -16.220, -9.404, -6.155, -1.000,
                 2.000, 2.000, -3.650, -17.867934687794072, -24.395]
    _longitude = [22.000, 22.000, 5.500, 0.000, -2.340, -40.000,
                  -86.100, -91.300, -91.800, -89.700, -84.000, -65.000,
                  -45.000, -38.400, -30.605, -11.123461787369832, 22.000]

    _time_range = Range(Time(744768005, format='fermi'),
                        Time(749416625, format='fermi'))


class GbmSaaPolygon15(GbmSaaPolygon5):
    """The coordinates of the GBM SAA boundary in latitude and East longitude
    used from 19:20:45 UTC on September 30, 2024 until mission end which uses a
    placeholder value of Jan 1, 2050 for now.

    Note: Reversion to the same shape as GbmSaaPolygon5 due to too many
    local particle triggers.
    """
    _time_range = Range(Time(749416850, format='fermi'),
                        Time(1546300805, format='fermi'))


class GbmSaa(GbmSaaPolygon5):
    """Class providing backwards compatibility for code written prior to v2.1.1
    where the GbmSaa class defined GbmSaaPolygon5"""
    pass


class GbmSaaCollection():
    """Collection of SAA boundary definitions"""

    def __init__(self):
        """Constructor"""
        self._polygons = [
            GbmSaaPolygon1(), GbmSaaPolygon2(), GbmSaaPolygon3(),
            GbmSaaPolygon4(), GbmSaaPolygon5(), GbmSaaPolygon6(),
            GbmSaaPolygon7(), GbmSaaPolygon8(), GbmSaaPolygon9(),
            GbmSaaPolygon10(), GbmSaaPolygon11(), GbmSaaPolygon12(),
            GbmSaaPolygon13(), GbmSaaPolygon14(), GbmSaaPolygon15()]
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
