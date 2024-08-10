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
import astropy.units as u
from astropy.timeseries import TimeSeries
import astropy.coordinates.representation as r
from gdt.core.coords import Quaternion
from gdt.core.file import FitsFileContextManager
from gdt.core.coords.spacecraft import SpacecraftFrameModelMixin, SpacecraftStatesModelMixin
from gdt.core.coords.spacecraft import SpacecraftFrame
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.headers import PosHistHeaders
from .detectors import GbmDetectors

__all__ = ['GbmPosHist']

FERMI_TO_UNIX_OFFSET = 978307200.0

class GbmPosHist(SpacecraftFrameModelMixin, SpacecraftStatesModelMixin, FitsFileContextManager):
    """Class for reading a GBM Position history file.
    """
    def _reorder_bytes(self, arr):
        """Method to reorder bytes according to old and new numpy API"""
        if np.__version__ >= '2.0.0':
            return arr.view(arr.dtype.newbyteorder()).byteswap()
        return arr.byteswap().newbyteorder()

    def get_spacecraft_frame(self) -> SpacecraftFrame:
        sc_frame = SpacecraftFrame(
            obsgeoloc=r.CartesianRepresentation(
                x=self._reorder_bytes(self.column(1, 'POS_X')),
                y=self._reorder_bytes(self.column(1, 'POS_Y')),
                z=self._reorder_bytes(self.column(1, 'POS_Z')),
                unit=u.m
            ),
            obsgeovel=r.CartesianRepresentation(
                x=self._reorder_bytes(self.column(1, 'VEL_X')) * u.m / u.s,
                y=self._reorder_bytes(self.column(1, 'VEL_Y')) * u.m / u.s,
                z=self._reorder_bytes(self.column(1, 'VEL_Z')) * u.m / u.s,
                unit=u.m / u.s
            ),
            quaternion=Quaternion(
                self._reorder_bytes(
                    self.columns_as_array(1, ['QSJ_1', 'QSJ_2', 'QSJ_3', 'QSJ_4']))),
            obstime=Time(self.column(1, 'SCLK_UTC'), format='fermi'),
            detectors=GbmDetectors
        )
        return sc_frame

    def get_spacecraft_states(self) -> TimeSeries:
        saa = self._in_saa(self.column(1, 'FLAGS'))
        series = TimeSeries(
            time=Time(self.column(1, 'SCLK_UTC'), format='fermi'),
            data={
                'sun': self._in_sun(self.column(1, 'FLAGS')),
                'saa': saa,
                'good': ~saa
            }
        )
        return series

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a GBM POSHIST FITS file.
        
        Args:
            file_path (str): The file path of the FITS file
        
        Returns:        
            (:class:`GbmPosHist`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = PosHistHeaders.from_headers(hdrs)
        return obj

    @staticmethod
    def _in_sun(flags: np.array) -> np.array:
        return (flags & 0x01) != 0

    @staticmethod
    def _in_saa(flags: np.array) -> np.array:
        return (flags & 0x02) != 0
