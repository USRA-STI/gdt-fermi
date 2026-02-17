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
import astropy.coordinates.representation as r
import astropy.io.fits as fits
import astropy.table as table
from astropy.timeseries import TimeSeries
import astropy.units as u
from gdt.core.coords import Quaternion
from gdt.core.file import FitsFileContextManager
from gdt.core.coords.spacecraft import SpacecraftFrameModelMixin, SpacecraftStatesModelMixin
from ..frame import FermiFrame
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.headers import PosHistHeaders
from .detectors import GbmDetectors

__all__ = ['GbmPosHist']

FERMI_TO_UNIX_OFFSET = 978307200.0

class GbmPosHist(SpacecraftFrameModelMixin, SpacecraftStatesModelMixin, FitsFileContextManager):
    """Class for reading a GBM Position history file, which contains the 
    spacecraft position, attitude (orientation), velocity, and various state
    information.
    """
    def _reorder_bytes(self, arr):
        """Method to reorder bytes according to old and new numpy API"""
        if np.__version__ >= '2.0.0':
            return arr.view(arr.dtype.newbyteorder()).byteswap()
        return arr.byteswap().newbyteorder()

    def get_spacecraft_frame(self) -> FermiFrame:
        """Retrieve the spacecraft frame information (position, attitude, and
           velocity).
        
        Returns:        
            (:class:`FermiFrame`)
        """
        sc_frame = FermiFrame(
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
        """Retrieve the spacecraft state information:
          *  `sun`: True if the sun is visible (not behind Earth)
          *  `saa`: True if Fermi is inside the GBM definition of the SAA
          *  `good`: True if spacecraft is in a good time interval.  For GBM, this
                     is essentially the opposite of the `saa` flag. 
        
        Returns:        
            (:class:`astropy.timeseries.TimeSeries`)
        """
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
    def merge(cls, poshist1, poshist2):
        """Merge two GbmPosHist objects into a single object. The typical use
        case is to merge two consecutive objects into one.  The order of input
        for the files does not matter, and any duplicate time entries will be
        removed.
        
        Args:
            poshist1 (:class:`GbmPosHist`): The first position history object
            poshist2 (:class:`GbmPosHist`): The second position history object
        
        Returns:        
            (:class:`GbmPosHist`)
        """
        # convert FITS_rec object to tables so to merge
        table1 = table.Table(poshist1._hdulist[1].data)
        table2 = table.Table(poshist2._hdulist[1].data)
        
        # determine which is the reference based on the start time
        # and merge data
        if poshist1.headers[0]['TSTART'] < poshist2.headers[0]['TSTART']:
            new_table = table.vstack([table1, table2])
            ref_obj = poshist1
            other_obj = poshist2
        else:
            new_table = table.vstack([table2, table1])
            ref_obj = poshist2
            other_obj = poshist1
        
        # remove duplicate time entries
        new_table = table.unique(new_table, keys='SCLK_UTC')
        
        # update the primary header with the date/time range info
        prihdr = ref_obj._hdulist[0].header.copy()
        prihdr['DATE-END'] = other_obj._hdulist[0].header['DATE-END']
        prihdr['TSTOP'] = other_obj._hdulist[0].header['TSTOP']
        prihdr['FILENAME'] = ''
        
        last_infile = list(prihdr['INFILE*'].keys())[-1]
        last_infile = int(last_infile.split('INFILE')[1])
        other_infiles = list(other_obj._hdulist[0].header['INFILE*'].values())
        num_infiles = len(other_infiles)
        for i, infile in enumerate(other_infiles):
            prihdr[f'INFILE{last_infile + 1 + i:02}'] = infile
        
        # update the table header with the date/time range info
        tblhdr = ref_obj._hdulist[1].header.copy()
        tblhdr['DATE-END'] = other_obj._hdulist[0].header['DATE-END']
        tblhdr['TSTOP'] = other_obj._hdulist[0].header['TSTOP']
        
        # create the merged data HDU list
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=prihdr)
        hdulist.append(primary_hdu)
        hdulist.append(fits.BinTableHDU(data=new_table, header=tblhdr))
        
        # create the new object
        obj = cls()
        obj._hdulist = hdulist
        obj._headers = PosHistHeaders.from_headers([hdu.header \
                                                    for hdu in obj.hdulist])
        return obj            

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
