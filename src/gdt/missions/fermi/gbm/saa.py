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
from ..time import Time
from gdt.core.geomagnetic import SouthAtlanticAnomaly

__all__ = ['GbmSaa']

class GbmSaa(SouthAtlanticAnomaly):
    """The coordinates of the GBM SAA boundary in latitude and East longitude.
    Note that the boundary was updated on July 23, 2024 to support slightly
    more GRB detections during O4b.
    """

    def __init__(self, time: Time):
        """Constructor

        Args:
            time (astropy.time.Time): time for determining the set
                                      of SAA points in latitude/longitude
        """
        self.update(time)
        super().__init__()

    def update(self, time: Time):
        """Update the SAA polygon points to the version used
        at the given time.

        Args:
            time (astropy.time.Time): time for determining the set
                                      of SAA points in latitude/longitude
        """
        if time < Time(743454905, format='fermi'):
            # original SAA region used for the first 16 years of the mission
            self._latitude = [
                -30.000, -19.867, -9.733, 0.400, 2.000, 2.000, -1.000,
                -6.155, -8.880, -14.220, -18.404, -30.000, -30.000]
            self._longitude = [
                33.900, 12.398, -9.103, -30.605, -38.400, -45.000, -65.000,
                -84.000, -89.200, -94.300, -94.300, -86.100, 33.900]
        else:
            # SAA region made smaller on 2024-07-23T19:15:00 UTC to slightly
            # increase the number of GRB detections per year
            self._latitude = [
                -24.395, -30.000, -30.000, -30.000, -30.000, -24.060,
                -16.220, -8.638, -6.155, -1.000, 2.000, 2.000,
                -3.400, -19.570802973653127, -24.395]
            self._longitude = [
                22.000, 22.000, 0.000, -2.000, -86.100, -90.300,
                -90.300, -88.738, -84.000, -65.000, -45.000, -38.400,
                -30.605, -11.999457706441582, 22.000]

