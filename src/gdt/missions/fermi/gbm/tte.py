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

from gdt.core.tte import PhotonList
from gdt.core.data_primitives import Ebounds, Gti, EventList
from .detectors import GbmDetectors
from .headers import TteHeaders, TteTriggerHeaders, PhaiiHeaders, \
                     PhaiiTriggerHeaders
from ..time import Time
from .phaii import GbmPhaii

__all__ = ['GbmTte']


class GbmTte(PhotonList):
    """Class for Time-Tagged Event data.
    """ 
    @property
    def detector(self):
        """(str): The detector name"""
        try:
            return GbmDetectors.from_full_name(self.headers[0]['DETNAM']).name
        except:
            return self.headers[0]['DETNAM']

    def to_phaii(self, bin_method, *args, time_range=None, energy_range=None,
                 channel_range=None, **kwargs):
        if self.trigtime is not None:
            headers = PhaiiTriggerHeaders()
        else:
            headers = PhaiiHeaders()
        
        # do not copy the value of these keys
        exceptions = ['CREATOR', 'DATATYPE', 'EXTNAME', 'FILENAME', 'FILETYPE',
                      'HDUCLAS1']
        # copy over the key values for each header
        for i in range(self.headers.num_headers):
            for key, val in self.headers[i].items():
                if key in exceptions:
                    continue
                try:
                    headers[i][key] = val        
                except:
                    # header key is present in TTE but not in PHAII
                    pass
        headers['PRIMARY']['FILETYPE'] = 'PHAII'
        headers['PRIMARY']['DATATYPE'] = 'PHAII'
        
        return super().to_phaii(bin_method, *args, time_range=time_range, 
                                energy_range=energy_range, 
                                channel_range=channel_range, 
                                phaii_class=GbmPhaii, headers=headers, **kwargs)
    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a TTE FITS file and return the TTE object

        Args:
            file_path (str): The file path of the FITS file
        
        Returns:        
            (:class:`GbmTte`)
        """
        obj = super().open(file_path, **kwargs)
        trigtime = None
        
        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]
        if 'TRIGTIME' in hdrs[0].keys():
            headers = TteTriggerHeaders.from_headers(hdrs)
            trigtime = float(headers['PRIMARY']['TRIGTIME'])
        else:
            headers = TteHeaders.from_headers(hdrs)
        
        # the channel energy bounds
        ebounds = Ebounds.from_bounds(obj.column(1, 'E_MIN'), 
                                      obj.column(1, 'E_MAX'))
        
        # data
        times = obj.column(2, 'TIME')
        if trigtime is not None:
            times -= trigtime
        data = EventList(times=times, channels=obj.column(2, 'PHA'), 
                         ebounds=ebounds)
        
        # the good time intervals
        gti_start = obj.column(3, 'START')
        gti_stop = obj.column(3, 'STOP')
        if trigtime is not None:
            gti_start -= trigtime
            gti_stop -= trigtime
        gti = Gti.from_bounds(gti_start, gti_stop)

        obj.close()
        
        return cls.from_data(data, gti=gti, trigger_time=trigtime, 
                             filename=obj.filename,
                             headers=headers, 
                             event_deadtime=headers[2]['EVT_DEAD'],
                             overflow_deadtime=1e-5)

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
        
        # the events extension
        events_hdu = self._events_table()
        hdulist.append(events_hdu)        
        
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

    def _events_table(self):
        times = np.copy(self.data.times)
        if self.trigtime is not None:
            times += self.trigtime
        
        time_col = fits.Column(name='TIME', format='1D', bzero=self.trigtime,
                               unit='s', array=times)
        pha_col = fits.Column(name='PHA', format='1I', array=self.data.channels)
        
        hdu = fits.BinTableHDU.from_columns([time_col, pha_col], 
                                            header=self.headers['EVENTS'])
        for key, val in self.headers['EVENTS'].items():
            hdu.header[key] = val
        if self.trigtime is not None:
            hdu.header.comments['TZERO1'] = 'Offset, equal to TRIGTIME'
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
