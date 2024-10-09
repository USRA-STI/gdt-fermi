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
import astropy.io.fits as fits

from gdt.core.phaii import Phaii
from gdt.core.data_primitives import Ebounds, Gti, TimeEnergyBins
from .detectors import GbmDetectors
from .headers import PhaiiHeaders, PhaiiTriggerHeaders
from ..time import Time

__all__ = ['GbmPhaii', 'Ctime', 'Cspec']


class GbmPhaii(Phaii):
    """PHAII class for GBM time history spectra.
    """
    @property
    def detector(self):
        """(str): The detector name"""
        try:
            return GbmDetectors.from_full_name(self.headers[0]['DETNAM']).name
        except:
            return self.headers[0]['DETNAM']
        
    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a GBM PHAII FITS file and return the GbmPhaii object
        
        Args:
            file_path (str): The file path of the FITS file
        
        Returns:        
            (:class:`GbmPhaii`)
        """
        obj = super().open(file_path, **kwargs)
        trigtime = None
        
        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]
        if 'TRIGTIME' in hdrs[0].keys():
            headers = PhaiiTriggerHeaders.from_headers(hdrs)
            trigtime = float(headers['PRIMARY']['TRIGTIME'])
        else:
            headers = PhaiiHeaders.from_headers(hdrs)
        
        # the channel energy bounds
        ebounds = Ebounds.from_bounds(obj.column(1, 'E_MIN'), 
                                      obj.column(1, 'E_MAX'))
        
        # the 2D time-channel counts data 
        time = obj.column(2, 'TIME')
        endtime = obj.column(2, 'ENDTIME')
        exposure = obj._assert_exposure(obj.column(2, 'EXPOSURE'))
        if trigtime is not None:
            time -= trigtime
            endtime -= trigtime
        data = TimeEnergyBins(obj.column(2, 'COUNTS'), time, endtime, exposure,
                              obj.column(1, 'E_MIN'), obj.column(1, 'E_MAX'),
                              quality=obj.column(2, 'QUALITY'))
        
        # the good time intervals
        gti_start = obj.column(3, 'START')
        gti_stop = obj.column(3, 'STOP')
        if trigtime is not None:
            gti_start -= trigtime
            gti_stop -= trigtime
        gti = Gti.from_bounds(gti_start, gti_stop)
            
        
        if headers[0]['DATATYPE'] == 'CSPEC':
            class_ = Cspec
        elif headers[0]['DATATYPE'] == 'CTIME':
            class_ = Ctime
        else:
            class_ = cls
          
        obj.close()
        
        return class_.from_data(data, gti=gti, trigger_time=trigtime, 
                                filename=obj.filename, headers=headers)
    
    def _build_hdulist(self):

        # create FITS and primary header
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)
        
        # the ebounds extension
        ebounds_hdu = self._ebounds_table()
        hdulist.append(ebounds_hdu)
        
        # the spectrum extension
        spectrum_hdu = self._spectrum_table()
        hdulist.append(spectrum_hdu)        
        
        # the GTI extension
        gti_hdu = self._gti_table()
        hdulist.append(gti_hdu)
        
        return hdulist
        
    def _build_headers(self, trigtime, tstart, tstop, num_chans):
        
        headers = self.headers.copy()
        for hdu in headers:
            hdu['TSTART'] = tstart
            hdu['TSTOP'] = tstop
            try:
                hdu['DETCHANS'] = num_chans
            except:
                pass
            if trigtime is not None:
                hdu['TRIGTIME'] = trigtime
        
        return headers
    
    def _ebounds_table(self):
        chan_col = fits.Column(name='CHANNEL', format='1I', 
                               array=np.arange(self.num_chans, dtype=int))
        emin_col = fits.Column(name='E_MIN', format='1E', unit='keV', 
                               array=self.ebounds.low_edges())
        emax_col = fits.Column(name='E_MAX', format='1E', unit='keV', 
                               array=self.ebounds.high_edges())
        
        hdu = fits.BinTableHDU.from_columns([chan_col, emin_col, emax_col], 
                                            header=self.headers['EBOUNDS'])
        for key, val in self.headers['EBOUNDS'].items():
            hdu.header[key] = val

        return hdu

    def _spectrum_table(self):
        tstart = np.copy(self.data.tstart)
        tstop = np.copy(self.data.tstop)
        if self.trigtime is not None:
            tstart += self.trigtime
            tstop += self.trigtime
        
        counts_col = fits.Column(name='COUNTS', 
                                 format='{}I'.format(self.num_chans), 
                                 bzero=32768, bscale=1, unit='count',
                                 array=self.data.counts)
        expos_col = fits.Column(name='EXPOSURE', format='1E', unit='s', 
                                array=self.data.exposure)
        qual_col = fits.Column(name='QUALITY', format='1I', 
                               array=self.data.quality)
        time_col = fits.Column(name='TIME', format='1D', unit='s', 
                               bzero=self.trigtime, array=tstart)
        endtime_col = fits.Column(name='ENDTIME', format='1D', unit='s', 
                                  bzero=self.trigtime, array=tstop)
        hdu = fits.BinTableHDU.from_columns([counts_col, expos_col, qual_col, 
                                             time_col, endtime_col], 
                                            header=self.headers['SPECTRUM'])

        for key, val in self.headers['SPECTRUM'].items():
            hdu.header[key] = val
        hdu.header.comments['TZERO1'] = 'offset for unsigned integers'
        hdu.header.comments['TSCAL1'] = 'data are not scaled'
        if self.trigtime is not None:
            hdu.header.comments['TZERO4'] = 'Offset, equal to TRIGTIME'
            hdu.header.comments['TZERO5'] = 'Offset, equal to TRIGTIME'
        return hdu

    def _gti_table(self):
        tstart = np.array(self.gti.low_edges())
        tstop = np.array(self.gti.high_edges())
        if self.trigtime is not None:
            tstart += self.trigtime
            tstop += self.trigtime

        start_col = fits.Column(name='START', format='1D', unit='s', 
                                bzero=self.trigtime, array=tstart)
        stop_col = fits.Column(name='STOP', format='1D', unit='s', 
                                bzero=self.trigtime, array=tstop)
        hdu = fits.BinTableHDU.from_columns([start_col, stop_col], 
                                            header=self.headers['GTI'])
        
        for key, val in self.headers['GTI'].items():
            hdu.header[key] = val        
        if self.trigtime is not None:
            hdu.header.comments['TZERO1'] = 'Offset, equal to TRIGTIME'
            hdu.header.comments['TZERO2'] = 'Offset, equal to TRIGTIME'
        return hdu
        

class Cspec(GbmPhaii):
    """Class for GBM CSPEC data.
    """


class Ctime(GbmPhaii):
    """Class for GBM CTIME data.
    """

        
