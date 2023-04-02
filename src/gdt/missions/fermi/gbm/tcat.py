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
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord, Latitude, Longitude
from gdt.core.data_primitives import TimeRange
from gdt.core.file import FitsFileContextManager
from .headers import TcatHeaders

__all__ = ['Tcat']

class Tcat(FitsFileContextManager):
    """TCAT (Trigger CATalog) file class.
    """
    def __init__(self):
        super().__init__()
    
    @property
    def fermi_location(self):
        """(astropy.coordinates.Longitude, astropy.coordinates.Latitude): Fermi's 
        orbital longitude and latitude"""
        return (Longitude(self.headers[0]['GEO_LONG'], unit='deg'), 
                Latitude(self.headers[0]['GEO_LAT'], unit='deg'))
    
    @property
    def localizing_instrument(self):
        """(str): The localizing instrument"""
        return self.headers[0]['LOC_SRC']

    @property
    def location(self):
        """(astropy.coordinates.SkyCoord): Location in equatorial coordinates"""
        return SkyCoord(self.headers[0]['RA_OBJ'], self.headers[0]['DEC_OBJ'], 
                        frame='icrs', unit='deg')

    @property
    def name(self):
        """(str): Name of the trigger"""
        return self.headers[0]['OBJECT']

    @property
    def time_range(self):
        """(:class:`~gdt.core.data_primitives.TimeRange`): The time range"""
        return TimeRange(self.headers[0]['TSTART'], self.headers[0]['TSTOP'])

    @property
    def trigtime(self):
        """(float): The trigger time"""
        return self.headers[0]['TRIGTIME']
    
    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a TCAT FITS file and return the Tcat object
        
        Args:
            file_path (str): The file path of the FITS file
        
        Returns:        
            (:class:`Tcat`)
        """
        obj = super().open(file_path, **kwargs)
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = TcatHeaders.from_headers(hdrs)        
        return obj
    
    def _build_hdulist(self):
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self['PRIMARY'])
        hflag = 0
        for key, val in self['PRIMARY'].items():
            if key == 'HISTORY':
                primary_hdu.header['HISTORY'][hflag] = val
                hflag += 1
            else:
                primary_hdu.header[key] = val
        hdulist.append(primary_hdu)

    def __repr__(self):
        return '<Tcat: {}>'.format(self.name)