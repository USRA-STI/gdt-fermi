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
import os
import astropy.io.fits as fits
import astropy.units as u
import astropy.coordinates.representation as r
import numpy as np
import warnings

from ..time import *
from .detectors import GbmDetectors
from .headers import TrigdatHeaders, PhaiiTriggerHeaders
from .phaii import GbmPhaii
from gdt.core.coords.spacecraft import SpacecraftFrame
from gdt.core.coords.quaternion import Quaternion
from gdt.core.data_primitives import TimeEnergyBins, TimeRange, Gti
from gdt.core.file import FitsFileContextManager

__all__ = ['Trigdat', 'MaxRates', 'BackRates', 'FswLocation']

# Map the classification numbers to the string names
classifications = {0: 'ERROR', 1: 'UNRELOC', 2: 'LOCLPAR', 3: 'BELOWHZ',
                   4: 'GRB', 5: 'SGR', 6: 'TRANSNT', 7: 'DISTPAR',
                   8: 'SFL', 9: 'CYGX1', 10: 'SGR1806', 11: 'GROJ422',
                   19: 'TGF', 20: 'UNCERT', 21: 'GALBIN'}
# localization spectra
spectrum = ['hard', 'normal', 'soft']


class Trigdat(FitsFileContextManager):
    """Class for the GBM Trigger Data
    """
    def __init__(self):
        super().__init__()
        self._data = None
        self._trigrates = None
        self._maxrates = None
        self._backrates = None
        self._fsw_locations = None
        self._emin = np.array([3.4, 10.0, 22.0, 44.0, 95.0, 300., 500., 800.])
        self._emax = np.array([10., 22.0, 44.0, 95.0, 300., 500., 800., 2000.])
        self._time_range = None
        self._detectors = [det.name for det in GbmDetectors]
        self._trigtime = None

    @property
    def backrates(self):
        """(:class:`BackRates`): A BackRates object containing the info from 
                                 the on-board background estimate"""
        return self._backrates

    @property
    def fsw_locations(self):
        """(list of :class:`FswLocation`): A list of flight-software-determined 
                                           locations for the trigger"""
        return self._fsw_locations

    @property
    def gti(self):
        """(:class:`~gdt.core.data_primitives.Gti`): The Good Time Intervals"""
        return self._gti

    @property
    def maxrates(self):
        """(list of :class:`MaxRates`): A list of MaxRates objects, each 
                                        containing maxrates info"""
        return self._maxrates

    @property
    def num_maxrates(self):
        """(int): The number of MaxRates issued by the flight software"""
        return len(self._maxrates)

    @property
    def poshist(self):
        """(:class:`~.poshist.GbmPosHist`): The position/attitude history"""
        return self._poshist

    @property
    def time_range(self):
        """(:class:`~gdt.core.data_primitives.TimeRange`): The time range of 
                                                           the data"""
        return self._time_range

    @property
    def triggered_detectors(self):
        """(list of str): The detectors that were triggered"""
        try:
            detmask = self.headers['PRIMARY']['DET_MASK']
            detmask = np.array(list(detmask)).astype(int) == 1
        except:
            return None
        
        if detmask.size == 14:
            return (np.array(self._detectors)[detmask]).tolist()
        elif detmask.size == 12:
            return (np.array(self._detectors[:-2])[detmask]).tolist()
        else:
            return None

    @property
    def trigrates(self):
        """(:class:`MaxRates`): The trigger information and rates"""
        return self._trigrates

    @property
    def trigtime(self):
        """(float): The trigger time"""
        return self._trigtime

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open and read a Trigdat file
        
        Args:
            file_path (str):  The file path of the trigdat file
        
        Returns:
            (:class:`Trigdat`)
        """
        obj = super().open(file_path, **kwargs)

        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]
        obj._headers = TrigdatHeaders.from_headers(hdrs)
        obj._trigtime = obj._headers[0]['TRIGTIME']
        
        # store trigrate, maxrates, backrates, and fsw location
        trigrate_data =  obj.hdulist['TRIGRATE'].data
        if trigrate_data.size > 0:
            obj._trigrates = MaxRates.from_recarray(trigrate_data[0])
        obj._maxrates = [MaxRates.from_recarray(maxrate) for maxrate in
                         obj.hdulist['MAXRATES'].data]
        
        try:
            obj._backrates = BackRates.from_recarray(obj.hdulist['BCKRATES'].data[0])
        except:
            warnings.warn('No BCKRATES data in this file.', RuntimeWarning)      
        obj._fsw_locations = [FswLocation.from_recarray(ob_calc) \
                              for ob_calc in obj.hdulist['OB_CALC'].data]
        
        # store the data            
        obj._data = obj.hdulist['EVNTRATE'].data
        obj._data.sort(order='TIME')
        # incorrectly shaped in the file, so we have to fix that here
        obj._rates = obj._data['RATE'].reshape(-1, 14, 8)
                    
        # store the position history
        idx, dt = obj._time_indices(1024)
        eic = obj._data['EIC'][idx]
        quat = obj._data['SCATTITD'][idx]
        times = obj._data['TIME'][idx]
        obj._poshist = SpacecraftFrame(
            obsgeoloc=r.CartesianRepresentation(x=eic[:,0], y=eic[:,1], 
                                                z=eic[:,2], unit=u.km),
            quaternion=Quaternion(quat),
            obstime=Time(times, format='fermi'),
            detectors=GbmDetectors)
                

        obj._time_range = TimeRange(obj._data['ENDTIME'][0] - dt[0], 
                                    obj._data['ENDTIME'][-1])
        
        obj._gti = Gti.from_list(obj._gti_from_times(obj._data['TIME'],
                                                     obj._data['ENDTIME']))
        return obj

    @classmethod
    def from_data(cls, phaii_list, poshist, trigtime, time_range=None, 
                  trigrate=None, maxrates=None, backrates=None, 
                  fswlocations=None, filename=None, headers=None):
        """Create a Trigdat object from a list of GbmPhaii objects and a 
        SpacecraftFrame object.
        
        Note:
            The list of PHAIIs does not need to include every detector, as 
            missing detectors will be assumed to have zero rates.  However, 
            do not include more than one file from the same detector in the 
            list. The PHAII list can be in any order.
        
        Args:
            phaii_list (list of :class:`~.phaii.GbmPhaii`): A list of GbmPhaii 
                                                           objects
            poshist (:class:`~gdt.core.coords.SpacecraftFrame`): 
                The spacecraft frame
            trigtime (float): The trigger time
            time_range (tuple, optional): The time range of the Trigdat. If not
                                         specified, uses the time range of the
                                         GbmPhaiis
            trigrate (:class:`MaxRates`, optional): The trigger rates object
            maxrates (list of :class:`MaxRates`, optional): The maxrates objects
            backrates (:class:`BackRates`, optional): The background rates
                                                      object
            fswlocations (list of :class:`FswLocation`, optional): The flight
                         software location objects
        
        Returns:
            (:class:`Trigdat`)
        """
        # make sure phaii_list contains the correct objects and desired time
        # is contained within all files
        for phaii in phaii_list:
            if not isinstance(phaii, GbmPhaii):
                raise TypeError('phaii_list must be a list of GbmPhaii objects')
            phaii_time_range = TimeRange(*phaii.time_range)
        
        # make sure the poshist is a PosHist object and desired time is 
        # contained within the file
        if not isinstance(poshist, SpacecraftFrame):
            raise TypeError('poshist must be a SpacecraftFrame object')
        poshist_time_range = TimeRange(poshist.obstime[0].fermi, 
                                       poshist.obstime[-1].fermi)
        if not poshist_time_range.contains(trigtime):
            raise ValueError('{} does not contain the specified ' \
                             'trigtime'.format(poshist))
        
        # slice the files by time, if requested
        trange = (time_range[0]-trigtime, time_range[1]-trigtime)
        if time_range is not None:
            phaii_list = [phaii.slice_time(trange) for phaii in phaii_list]
        
        
        obj = cls()
        obj._poshist = poshist
        obj._trigtime = trigtime
        
        # create the data record array
        num_chans = phaii_list[0].num_chans
        num_dets = len(obj._detectors)
        num_bins = phaii_list[0].data.num_times
        dtypes = [('TIME', '>f8'), ('ENDTIME', '>f8'), 
                  ('SCATTITD', '>f4', (4,)), ('EIC', '>f4', (3,)), 
                  ('RATE', '>f4', (num_chans, num_dets))]
        obj._data = np.recarray((num_bins,), dtype=dtypes)
        tstart = phaii_list[0].data.tstart
        tstop = phaii_list[0].data.tstop
        tcent = phaii_list[0].data.time_centroids
        
        # interpolate poshist to match the bin edges as is standard
        fermi_times = Time(tstart+trigtime, format='fermi')
        obj._poshist = poshist.at(fermi_times)
        
        # to preserve the detector ordering assumed by the trigdat 
        # datatype, we must pull the phaii file from the list that
        # corresponds with the right detector in our sequence.  if the 
        # detector isn't represented in the phaii_list, then the rates will
        # be zeros.
        det_mask = np.zeros(14, dtype=int)
        for i in range(num_bins):
            rates = np.zeros((14,8))
            for j in range(len(obj._detectors)):
                try:
                    det_idx = [phaii.detector for \
                               phaii in phaii_list].index(obj._detectors[j])
                    det_mask[j] = 1
                except:
                    continue
                rates[j,:] = phaii_list[det_idx].data.counts[i,:] * \
                      (1.024/phaii_list[det_idx].data.time_widths[i,np.newaxis])
            
            
            quat = np.append(obj._poshist[i].quaternion.xyz, 
                             obj._poshist[i].quaternion.w)
            eic = obj._poshist.obsgeoloc[0].xyz.to('km').value
            obj._data[i] = (tstart[i]+trigtime, tstop[i]+trigtime, quat,
                           eic, rates.T)

        obj._time_range = TimeRange(obj._data['TIME'][0], 
                                    obj._data['ENDTIME'][-1])
        obj._gti = Gti.from_list(obj._gti_from_times(obj._data['TIME'],
                                                     obj._data['ENDTIME']))
    
        # assign/create trigrates, maxrates, backrates, and fsw_locations              
        if trigrate is not None:
            if not isinstance(trigrate, MaxRates):
                raise TypeError('trigrate must be a MaxRates object')
            obj._trigrates = trigrate
        else:
            obj._trigrates = MaxRates.create()
        
        if maxrates is not None:
            try:
                iter(maxrates)
            except:
                raise TypeError('maxrates must be a list of MaxRate objects')
            if any([not isinstance(m, MaxRates) for m in maxrates]):
                raise TypeError('maxrates must be a list of MaxRate objects')
            obj._maxrates = maxrates
        else:
            obj._maxrates = [MaxRates.create()]
        
        if backrates is not None:
            if not isinstance(backrates, BackRates):
                raise TypeError('backrates must be a BackRats object')
            obj._backrates = backrates
        else:
            obj._backrates = BackRates.create()
        
        if fswlocations is not None:
            try:
                iter(fswlocations)
            except:
                raise TypeError('fsw_locations must be a list of FswLocation ' \
                                'objects')
            if any([not isinstance(f, FswLocation) for f in fswlocations]):
                raise TypeError('fsw_locations must be a list of FswLocation ' \
                                'objects')
            obj._fsw_locations = fswlocations
        else:
            obj._fsw_locations = [FswLocation.create()]
                 
        if filename is not None:
            obj._filename = str(filename)
        
        if headers is not None:
            if not isinstance(headers, TrigdatHeaders):
                raise ValueError('headers must be a TrigdatHeaders object')
            obj._headers = headers
        else:
            headers = TrigdatHeaders()
            pos = obj._poshist.scpos.interpolate(obj.trigtime)
            headers[0]['FILENAME'] = obj._filename
            for hdu in headers:
                try:
                    hdu['TSTART'] = obj._data['TIME'][0]
                    hdu['TSTOP'] = obj._data['ENDTIME'][-1]
                except:
                    pass
                try:
                    hdu['TRIGTIME'] = obj._trigtime
                except:
                    pass
                try:
                    hdu['GEO_LONG'] = pos.longitude.value
                    hdu['GEO_LAT'] = pos.latitude.value
                except:
                    pass
                try:
                    hdu['DETTYPE'] = 'BOTH'
                except:
                    pass
            obj._headers = headers
            
        obj._headers[0]['DET_MASK'] = ''.join(det_mask.astype(str))
        
        return obj
    
    def to_phaii(self, detector, timescale=1024):
        """Convert the data for a detector to a :class:`~.phaii.GbmPhaii`
        object.
        
        Note:
            A standard Trigdat will contain data on 8192, 1024, 256, and 62 ms
            timescales.  A non-standard Trigdat can be created from other 
            data types (see :meth:`Trigdat.from_data`) and therefore the 
            timescales will likely be different and arbitrary. For that reason,
            the ``timescale`` keyword will be ignored for non-standard Trigdat.
        
        Args:
            detector (str): The detector to convert
            timescale (int, optional): 
                The minimum timescale in ms of the data to return. Available 
                options are 1024, 256, and 64. This is ignored if the trigdat
                has non-standard timescales.
        
        Returns:
            (:class:`~.phaii.GbmPhaii`)
        """
        # check for valid detectors and timescale
        detector = detector.lower()
        if detector not in self._detectors:
            raise ValueError('Illegal detector name')

        # grab the correct detector rates (stored as 14 dets x 8 channels)
        det_idx = self._detectors.index(detector)

        if (timescale != 1024) and (timescale != 256) and (timescale != 64):
            raise ValueError('Illegal Trigdat resolution. Available resolutions: \
                              1024, 256, 64')

        # return the requested timescales
        try:
            time_idx, dt = self._time_indices(timescale)
        except:
            warnings.warn('\nTrigdat.to_ctime(): ' \
                          'Exporting from non-standard file. ' \
                          'Timescale is ignored.')
            time_idx = np.arange(self._rates.shape[0])
            dt = self._data['ENDTIME']-self._data['TIME']
            
        # calculate counts and exposure
        counts = np.round(self._rates[time_idx, det_idx, :] * \
                          (dt[:, np.newaxis] / 1.024))
        exposure = self._calc_exposure(counts, dt)

        # the 'TIME' and 'ENDTIME' edges are not aligned before being written to
        # file, so we must fix this to prevent TimeEnergyBins from thinking 
        # every single bin is a single discontiguous segment
        tstop = self._data['ENDTIME'][time_idx] - self.trigtime
        tstart = self._fix_tstart(tstop, dt)

        # create the Time-Energy histogram
        bins = TimeEnergyBins(counts, tstart, tstop, exposure,
                              self._emin, self._emax)
        
        headers = self._set_phaii_headers(detector, bins.num_chans)
        
        obj = GbmPhaii.from_data(bins, gti=self.gti, trigger_time=self.trigtime,
                                 headers=headers)
        return obj

    def sum_detectors(self, detectors, timescale=1024):
        """Sum the data from a list of detectors and convert to a GbmPhaii 
        object. The exposures from different detectors are averaged.
        
        Args:
            detectors (list of str): The detectors to sum
            timescale (int, optional): 
                The minimum timescale in ms of the data to return.  Available 
                options are 1024, 256, and 64.

        Returns:
            (:class:`~.phaii.GbmPhaii`)
        """
        # check for valid detectors and timescale
        for det in detectors:
            det = det.lower()
            if det not in self._detectors:
                raise ValueError('Illegal detector name')

        if (timescale != 1024) and (timescale != 256) and (timescale != 64):
            raise ValueError('Illegal Trigdat resolution. Available resolutions: \
                              1024, 256, 64')
        
        try:
            time_idx, dt = self._time_indices(timescale)
        except:
            warnings.warn('\nTrigdat.to_ctime(): ' \
                          'Exporting from non-standard file. ' \
                          'Timescale is ignored.')
            time_idx = np.arange(self._rates.shape[0])
            dt = self._data['ENDTIME']-self._data['TIME']

        counts = None
        exposure = None
        for det in detectors:
            # grab the correct detector rates (stored as 14 dets x 8 channels)
            det_idx = self._detectors.index(det)

            # calculate counts
            c = np.round(self._rates[time_idx, det_idx, :] * \
                         (dt[:, np.newaxis] / 1.024))
            if counts is None:
                counts = c        
            else:
                counts += c
            
            # calculate exposure
            e = self._calc_exposure(c, dt)
            if exposure is None:
                exposure = e
            else:
                exposure += e
        
        exposure /= float(len(detectors))
            

        # the 'TIME' and 'ENDTIME' edges are not aligned before being written to
        # file, so we must fix this to prevent TimeEnergyBins from thinking 
        # every single bin is a single discontiguous segment
        tstop = self._data['ENDTIME'][time_idx] - self.trigtime
        tstart = self._fix_tstart(tstop, dt)

        # create the Time-Energy histogram
        bins = TimeEnergyBins(counts, tstart, tstop, exposure,
                              self._emin, self._emax)

        det_str = '+'.join(detectors)
        headers = self._set_phaii_headers(det_str, bins.num_chans)
        
        obj = GbmPhaii.from_data(bins, gti=self.gti, trigger_time=self.trigtime,
                                 headers=headers)
        return obj

    def _build_hdulist(self):
        
        # create FITS and primary header
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)
        
        # TRIGRATE extension
        trigrate_hdu = fits.BinTableHDU(data=self.trigrates.to_recarray(),
                                        name='TRIGRATE', 
                                        header=self.headers['TRIGRATE'])
        for key, val in self.headers['TRIGRATE'].items():
            trigrate_hdu.header[key] = val
        trigrate_hdu.header.comments['TTYPE1'] = 'Beginning of accumulation, ' \
                                                 'calculated value'
        trigrate_hdu.header.comments['TTYPE2'] = 'End of accumulation, same as'\
                                                 ' PCKTTIME'
        trigrate_hdu.header.comments['TTYPE3'] = 'Spacecraft attitude ' \
                                                 'quaternions'
        trigrate_hdu.header.comments['TTYPE4'] = 'Spacecraft position: Earth ' \
                                                 'X, Y, & Z'
        trigrate_hdu.header.comments['TTYPE5'] = 'Rates-14 detectors, 8 channels'
        hdulist.append(trigrate_hdu)
        

        # BCKRATES extension
        backrates_hdu = fits.BinTableHDU(data=self.backrates.to_recarray(),
                                         name='BCKRATES', 
                                         header=self.headers['BCKRATES'])
        for key, val in self.headers['BCKRATES'].items():
            backrates_hdu.header[key] = val
        backrates_hdu.header.comments['TTYPE1'] = 'Start time of the ' \
                                                  'background accumulation'
        backrates_hdu.header.comments['TTYPE2'] = 'Creation time of the ' \
                                                  'background model'
        backrates_hdu.header.comments['TTYPE3'] = 'Quality Flag for the ' \
                                                  'background model'
        backrates_hdu.header.comments['TTYPE4'] = 'Rates-14 detectors, 8 channels'                            
        hdulist.append(backrates_hdu)

        # OB_CALC extension
        obcalc_hdu = fits.BinTableHDU(data=np.array([loc.to_recarray() for \
                                                    loc in self.fsw_locations]),
                                      name='OB_CALC', 
                                      header=self.headers['OB_CALC'])
        for key, val in self.headers['OB_CALC'].items():
            obcalc_hdu.header[key] = val
        obcalc_hdu.header.comments['TTYPE1'] = 'End time of location rates ' \
                                               'accumulation'
        obcalc_hdu.header.comments['TTYPE2'] = 'Calculated Right Ascension of '\
                                               'source'
        obcalc_hdu.header.comments['TTYPE3'] = 'Calculated Declination of source'
        obcalc_hdu.header.comments['TTYPE4'] = 'Error radius of calculated' \
                                               'location'
        obcalc_hdu.header.comments['TTYPE5'] = 'Location algorithm used in ' \
                                               'calculation'
        obcalc_hdu.header.comments['TTYPE6'] = 'Two highest probable event ' \
                                               'classes'
        obcalc_hdu.header.comments['TTYPE7'] = 'Reliabilities of probable ' \
                                               'event classes'
        obcalc_hdu.header.comments['TTYPE8'] = 'Standardized count intensity'
        obcalc_hdu.header.comments['TTYPE9'] = 'Standard hardness ratio ' \
                                               '(channels 3 / 4)'
        obcalc_hdu.header.comments['TTYPE10'] = 'Fluence counts in 3 + 4 ' \
                                                'energy band'
        obcalc_hdu.header.comments['TTYPE11'] = 'Significance of localization '\
                                                'in sigma'
        obcalc_hdu.header.comments['TTYPE12'] = 'Localization rates-12 NaI ' \
                                                'detectors'
        obcalc_hdu.header.comments['TTYPE13'] = 'Location data timescale'
        obcalc_hdu.header.comments['TTYPE14'] = 'Source azimuth with respect ' \
                                                'to Spacecraft'
        obcalc_hdu.header.comments['TTYPE15'] = 'Source zenith with respect to'\
                                                ' Spacecraft'
        hdulist.append(obcalc_hdu)

        # MAXRATES extension
        maxrates_hdu = fits.BinTableHDU(data=np.array([mr.to_recarray() for \
                                                       mr in self.maxrates]),
                                      name='MAXRATES', 
                                      header=self.headers['MAXRATES'])
        for key, val in self.headers['MAXRATES'].items():
            maxrates_hdu.header[key] = val
        maxrates_hdu.header.comments['TTYPE1'] = 'Beginning of accumulation, ' \
                                                 'calculated value'
        maxrates_hdu.header.comments['TTYPE2'] = 'End of accumulation, same as'\
                                                 ' PCKTTIME'
        maxrates_hdu.header.comments['TTYPE3'] = 'Spacecraft attitude ' \
                                                 'quaternions'
        maxrates_hdu.header.comments['TTYPE4'] = 'Spacecraft position: Earth ' \
                                                 'X, Y, & Z'
        maxrates_hdu.header.comments['TTYPE5'] = 'Rates-14 detectors, 8 channels'
        hdulist.append(maxrates_hdu)

        # EVNTRATE extension
        evntrate_hdu = fits.BinTableHDU(data=self._data, name='EVNTRATE', 
                                        header=self.headers['EVNTRATE'])
        for key, val in self.headers['EVNTRATE'].items():
            evntrate_hdu.header[key] = val
        evntrate_hdu.header.comments['TTYPE1'] = 'Beginning of accumulation, ' \
                                                 'calculated value'
        evntrate_hdu.header.comments['TTYPE2'] = 'End of accumulation, same as'\
                                                 ' PCKTTIME'
        evntrate_hdu.header.comments['TTYPE3'] = 'Spacecraft attitude ' \
                                                 'quaternions'
        evntrate_hdu.header.comments['TTYPE4'] = 'Spacecraft position: Earth ' \
                                                 'X, Y, & Z'
        evntrate_hdu.header.comments['TTYPE5'] = 'Rates-14 detectors, 8 channels'
        hdulist.append(evntrate_hdu)
        
        return hdulist

    def _fix_tstart(self, tstop, dt):
        # this ensures that edge differences < 1 ms get fixed
        tstart = tstop - dt
        mask = (np.abs(tstart[1:] - tstop[:-1]) < 0.001)
        tstart[1:][mask] = tstop[:-1][mask]
        return tstart

    def _calc_exposure(self, counts, dt):
        """Calculate the exposure
        
        Args:
            counts (np.array): The observed counts in each bin
            dt (np.array): The time bin widths 
        
        Returns:
            (np.array)
        """
        deadtime = np.copy(counts)
        deadtime[:, :-1] *= 2.6e-6  # 2.6 microsec for each count
        deadtime[:, -1] *= 1e-5  # 10 microsec for each count in overflow
        total_deadtime = np.sum(deadtime, axis=1)
        exposure = (1.0 - total_deadtime) * dt
        return exposure

    def _time_indices(self, time_res):
        """Indices into the Trigdat arrays corresponding to the desired time 
        resolution(s) 
        
        Args:
            time_res (int): The time resolution in ms of the data
        
        Returns:
            (np.array, np.array): Indices into the trigdat arrays and the \
                                  bin widths in seconds
        """
        # bin widths
        dt = np.round((self._data['ENDTIME'] - self._data['TIME']) * 1000)
        # background bins
        back_idx = np.where(dt == 8192)[0]
        # 1 s scale bins - this is the minimum amount returned
        idx = np.where(dt == 1024)[0]
        cnt = len(idx)
        # reconcile 8 s and 1 s data
        idx = self._reconcile_timescales(back_idx, idx)

        # reconcile 8 s + 1 s and 256 ms data
        if time_res <= 256:
            tidx = np.where(dt == 256)[0]
            idx = self._reconcile_timescales(idx, tidx)

        # reconcile 8 s + 1 s + 256 ms and 64 ms data
        if time_res == 64:
            tidx = np.where(dt == 64)[0]
            idx = self._reconcile_timescales(idx, tidx)

        # return reconciled indices
        return idx, np.reshape(dt[idx] / 1000.0, len(idx))

    def _reconcile_timescales(self, idx1, idx2):
        """Reconcile indices representing different timescales and glue them 
           together to form a complete (mostly) continuous set of indices
        
        Args: 
            idx1 (np.array): Indices of the "bracketing" timescale
            idx2 (np.array): Indices of the "inserted" timescale 
        
        Returns:
            np.array: Indices of idx2 spliced into idx1
		"""
        # if there is no idx1, then we only have to return idx2 and v.v.
        if idx1.size == 0:
            return idx2
        if idx2.size == 0:
            return idx1

        # bin edges for both selections
        start_times1 = self._data['TIME'][idx1]
        end_times1 = self._data['ENDTIME'][idx1]
        start_times2 = self._data['TIME'][idx2]
        end_times2 = self._data['ENDTIME'][idx2]

        # find where bracketing timescale ends and inserted timescale begins
        mask = end_times1 >= start_times2[0]
        if mask.sum():
            start_idx = (np.where(mask))[0][0]
            idx = np.concatenate((idx1[0:start_idx], idx2))
        else:
            idx = np.concatenate((idx1, idx2))

        # find where inserted timescale ends and bracketing timescale begins again
        mask = start_times1 >= end_times2[-1]
        if mask.sum():
            end_idx = (np.where(mask))[0][0]
            idx = np.concatenate((idx, idx1[end_idx:]))
                
        return idx

    def _gti_from_times(self, tstarts, tstops):
        """Estimate the GTI from the bin start and stop times.
        This may return multiple GTIs if several background packets are missing
		
		Args:
            tstarts (np.array): The start times of the bins
            tstops (np.array): The end times of the bins
        
        Returns:
            [(float, float), ...]: A list of time ranges
        """
        tstart = tstarts[0]
        tstop = tstops[-1]
        dt = tstarts[1:] - tstops[:-1]
        idx = np.where(np.abs(dt) > 10.0)[0]
        if idx.sum() > 0:
            gti = [(tstart, tstops[idx[0] - 1]), (tstarts[idx[0]], tstop)]
        else:
            gti = [(tstart, tstop)]
        return np.array(gti)

    def _set_phaii_headers(self, det, num_chans):
        headers = PhaiiTriggerHeaders()
        try:
            detname = GbmDetectors.from_str(det).full_name
        except:
            detname = det
        
        headers[0]['DETNAM'] = detname
        headers[0]['DATE-OBS'] = self.headers[0]['DATE-OBS']
        headers[0]['DATE-END'] = self.headers[0]['DATE-END']
        headers[0]['TSTART'] = self.headers[0]['TSTART']
        headers[0]['TSTOP'] = self.headers[0]['TSTOP']
        headers[0]['TRIGTIME'] = self.trigtime
        headers[0]['OBJECT'] = self.headers[0]['OBJECT']
        headers[0]['RA_OBJ'] = self.headers[0]['RA_OBJ']
        headers[0]['DEC_OBJ'] = self.headers[0]['DEC_OBJ']
        headers[0]['ERR_RAD'] = self.headers[0]['ERR_RAD']

        headers[1]['DETNAM'] = detname
        headers[1]['DATE-OBS'] = self.headers[0]['DATE-OBS']
        headers[1]['DATE-END'] = self.headers[0]['DATE-END']
        headers[1]['TSTART'] = self.headers[0]['TSTART']
        headers[1]['TSTOP'] = self.headers[0]['TSTOP']
        headers[1]['TRIGTIME'] = self.trigtime
        headers[1]['OBJECT'] = self.headers[0]['OBJECT']
        headers[1]['RA_OBJ'] = self.headers[0]['RA_OBJ']
        headers[1]['DEC_OBJ'] = self.headers[0]['DEC_OBJ']
        headers[1]['ERR_RAD'] = self.headers[0]['ERR_RAD']
        headers[1]['DETCHANS'] = num_chans

        headers[2]['DETNAM'] = detname
        headers[2]['DATE-OBS'] = self.headers[0]['DATE-OBS']
        headers[2]['DATE-END'] = self.headers[0]['DATE-END']
        headers[2]['TSTART'] = self.headers[0]['TSTART']
        headers[2]['TSTOP'] = self.headers[0]['TSTOP']
        headers[2]['TRIGTIME'] = self.trigtime
        headers[2]['OBJECT'] = self.headers[0]['OBJECT']
        headers[2]['RA_OBJ'] = self.headers[0]['RA_OBJ']
        headers[2]['DEC_OBJ'] = self.headers[0]['DEC_OBJ']
        headers[2]['ERR_RAD'] = self.headers[0]['ERR_RAD']
        headers[2]['DETCHANS'] = num_chans

        headers[3]['DETNAM'] = detname
        headers[3]['DATE-OBS'] = self.headers[0]['DATE-OBS']
        headers[3]['DATE-END'] = self.headers[0]['DATE-END']
        headers[3]['TSTART'] = self.headers[0]['TSTART']
        headers[3]['TSTOP'] = self.headers[0]['TSTOP']
        headers[3]['TRIGTIME'] = self.trigtime
        headers[3]['OBJECT'] = self.headers[0]['OBJECT']
        headers[3]['RA_OBJ'] = self.headers[0]['RA_OBJ']
        headers[3]['DEC_OBJ'] = self.headers[0]['DEC_OBJ']
        headers[3]['ERR_RAD'] = self.headers[0]['ERR_RAD']
        return headers        
    
    def __repr__(self):
        s = '<Trigdat: {};\n'.format(self.filename)
        s += ' trigtime={:.2f};'.format(self.trigtime)
        s += ' triggered detectors={}>'.format(self.triggered_detectors)
        return s


class MaxRates:
    """Class for the MAXRATES and TRIGRATE data in Trigdat.
    """
    _dtypes = [('TIME', '>f8'), ('ENDTIME', '>f8'), ('SCATTITD', '>f4', (4,)), 
               ('EIC', '>f4', (3,)), ('TRIGRATE', '>f4', (8,14))]
    
    def __init__(self):
        self._detectors = [det.name for det in GbmDetectors]
        self._time_range = (None, None)
        self._quats = None
        self._eic = None
        self._rates = None

    @property
    def all_rates(self):
        """(np.array): An array (:attr:`num_chans`, :attr:`num_dets`) of 
                       the maxrates"""
        return self._rates

    @property
    def eic(self):
        """(np.array): The position of Fermi in Earth inertial coordinates"""
        return self._eic

    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        if self._rates is not None:
            return self._rates.shape[0]
        else:
            return 0

    @property
    def num_dets(self):
        """(int): The number of detectors"""
        if self._rates is not None:
            return self._rates.shape[1]
        else:
            return 0

    @property
    def quaternion(self):
        """ (np.array): The quaternions at the maxrates time"""
        return self._quats

    @property
    def time_range(self):
        """(float, float): The time range of the maxrates"""
        return self._time_range

    @property
    def timescale(self):
        """(float): The timescale of the maxrates"""
        if self.time_range[0] is not None and self.time_range[1] is not None:
            return round((self.time_range[1] - self.time_range[0]) * 1000.0)
        else:
            return None

    @classmethod
    def create(cls, tstart=None, tstop=None, quaternion=None, eic=None, 
               rates=None):
        """Create a MaxRates object from keywords.
        
        Args:
            tstart (float, optional): The start time of the MaxRates interval
            tstop (float, optional): The stop time of the MaxRates interval
            quaternion (np.array, optional): The attitude quaternion
            eic (np.array, optional): The position in earth inertial coordinates
            rates (np.array): The (num channels, num detectors) rates array
        
        Returns:
            (:class:`MaxRates`)
        """
        obj = cls()
        obj._time_range = (tstart, tstop)
        
        if isinstance(quaternion, (list, tuple)): 
            quaternion = np.array(quaternion)
        if quaternion is not None:
            if quaternion.size != 4:
                raise ValueError('quaternion must have 4 elements')
        obj._quats = quaternion

        if isinstance(eic, (list, tuple)): 
            eic = np.array(eic)
        if eic is not None:
            if eic.size != 3:
                raise ValueError('eic must have 3 elements')
        obj._eic = eic

        if isinstance(rates, (list, tuple)):
            rates = np.array(rates)
        obj._rates = rates
        
        return obj

    @classmethod
    def from_recarray(cls, rec_array):
        """Create from a FITS or numpy record array.
        
        Args:
            rec_array (np.recarray): The FITS TRIGRATE or MAXRATES record array 
                                     from the trigdat file
        
        Returns:
            (:class:`MaxRates`)
        """
        obj = cls()
        obj._time_range = (rec_array['TIME'], rec_array['ENDTIME'])
        obj._quats = rec_array['SCATTITD']
        obj._eic = rec_array['EIC']
        try:
            obj._rates = rec_array['TRIGRATE']
        except:
            obj._rates = rec_array['MAXRATES']
        return obj
    
    def get_detector(self, det):
        """Retrieve the rates for a detector
        
        Args:
            det (str): The detector
        
        Returns:
            (np.array)
        """
        if self._rates is None:
            return None
        mask = (np.array(self._detectors) == det)
        return self._rates[:, mask].squeeze()

    def to_recarray(self):
        """Convert contents of the object to a numpy record array.
        
        Returns:
            (np.recarray)
        """
        dtypes = self._dtypes
        dtypes[-1] = (*dtypes[-1][:-1], (self.num_chans, self.num_dets))
        record = np.array([(*self._time_range, self._quats, self._eic, 
                            self._rates)], dtype=dtypes)
        return record
                            
    def __repr__(self):
        return '<MaxRates: {0} ms>'.format(self.timescale)


class BackRates:
    """Class for the background rates data in Trigdat.
    """
    _dtypes = [('TIME', '>f8'), ('ENDTIME', '>f8'), ('QUALITY', 'u1', (2,)), 
               ('BCKRATES', '>f4', (8,14))]

    def __init__(self):
        self._detectors = [det.name for det in GbmDetectors]
        self._time_range = (None, None)
        self._quality = np.array((0, 0))
        self._rates = None

    @property
    def all_rates(self):
        """(np.array): An array (:attr:`num_chans`, :attr:`num_dets`) of 
                       the background rates"""
        return self._rates
        
    @property
    def num_chans(self):
        """(int): The number of energy channels"""
        if self._rates is not None:
            return self._rates.shape[0]
        else:
            return 0

    @property
    def num_dets(self):
        """(int): The number of detectors"""
        if self._rates is not None:
            return self._rates.shape[1]
        else:
            return 0

    @property
    def quality(self):
        """(int, int): The quality flags for the background"""
        return self._quality

    @property
    def time_range(self):
        """(float, float): The time range of the background rates"""
        return self._time_range

    @classmethod
    def create(cls, tstart=None, tstop=None, quality=None, rates=None):
        """Create a BackRates object from keywords.
        
        Args:
            tstart (float, optional): The start time of the MaxRates interval
            tstop (float, optional): The stop time of the MaxRates interval
            quality (int, int): Quality flags
            rates (np.array): The (num channels, num detectors) rates array
        
        Returns:
            (:class:`BackRates`)
        """
        obj = cls()
        obj._time_range = (tstart, tstop)
        
        if quality is not None:
            if not isinstance(quality, (list, tuple)):
                raise TypeError('quality must be a 2-tuple integer')
            if len(quality) != 2:
                raise ValueError('quality must be a 2-tuple integer')
            obj._quality = np.asarray(quality)

        if isinstance(rates, (list, tuple)):
            rates = np.array(rates)
        obj._rates = rates
        
        return obj

    @classmethod
    def from_recarray(cls, rec_array):
        """Create from a FITS or numpy record array.
        
        Args:
            rec_array (np.recarray): The FITS BCKRATES record array from the 
                                     trigdat file
        
        Returns:
            (:class:`BackRates`)
        """
        obj = cls()
        obj._time_range = (rec_array['TIME'], rec_array['ENDTIME'])
        obj._quality = rec_array['QUALITY']
        obj._rates = rec_array['BCKRATES']
        return obj
    
    def get_detector(self, det):
        """Retrieve the background rates for a detector

        Args:
            det (str): The detector
        
        Returns:
            np.array: An array of size (:attr:`num_chans`,) of background rates 
                      for the detector 
        """
        if self._rates is None:
            return None
        mask = (np.array(self._detectors) == det)
        return self._rates[:, mask].squeeze()

    def to_recarray(self):
        """Convert contents of the object to a numpy record array
        
        Returns:
            (np.recarray)
        """
        dtypes = self._dtypes
        dtypes[-1] = (*dtypes[-1][:-1], (self.num_chans, self.num_dets))
        record = np.array([(*self._time_range, self._quality, self._rates)], 
                          dtype=dtypes)
        return record

    def __repr__(self):
        return '<BackRates: time range={0}>'.format(self.time_range)


class FswLocation:
    """Class for the Flight Software localization information.
    """
    _dtypes = [('TIME', '>f8'), ('RA', '>f4'), ('DEC', '>f4'), 
               ('STATERR', '>f4'), ('LOCALG', '>i2'), ('EVTCLASS', '>i2', (2,)), 
               ('RELIABLT', '>f4', (2,)), ('INTNSITY', '>f4'), 
               ('HDRATIO', '>f4'), ('FLUENCE', '>f4'), ('SIGMA', '>f4'),
               ('LOCRATES', '>i2', (12,)), ('TRIG_TS', '>f4'), 
               ('TR_SCAZ', '>f4'), ('TR_SCZEN', '>f4')]
    def __init__(self):
        self._time = None
        self._location = (None, None, None)
        self._algorithm = None
        self._class1 = (20, 1.0)
        self._class2 = (20, 0.0)
        self._intensity = None
        self._hardness = None
        self._fluence = None
        self._sigma = None
        self._rates = None
        self._timescale = None
        self._azzen = (None, None)
    
    @property
    def fluence(self):
        """(float): The fluence of the localization interval"""
        return self._fluence

    @property
    def hardness_ratio(self):
        """(float): The hardness ratio for the localization interval"""
        return self._hardness

    @property
    def intensity(self):
        """(float): The brightness of the signal"""
        return self._intensity

    @property
    def location(self):
        """(float, float, float): The RA, Dec, and statistical error of the
                                  onboard localization"""
        return self._location

    @property
    def location_sc(self):
        """(float, float): The localization in spacecraft coordinates: 
                           Azimuth, Zenith"""
        return self._azzen

    @property
    def next_classification(self):
        """(str, float): The next most likely classification of the trigger and 
                         the probability"""
        return self._class2

    @property
    def rates(self):
        """(np.array): The rates in each NaI detector used for localization"""
        return self._rates

    @property
    def significance(self):
        """(float): The S/N ratio of the localization interval"""
        return self._sigma

    @property
    def spectrum(self):
        """(str): The spectrum used in the localization"""
        return self._algorithm

    @property
    def time(self):
        """(float): Time at which the localization was calculated"""
        return self._time

    @property
    def timescale(self):
        """(float): The localization interval timescale"""
        return self._timescale

    @property
    def top_classification(self):
        """(str, float): The most likely classification of the trigger and the 
                         probability"""
        return self._class1

    @classmethod
    def create(cls, time=None, radec=(None, None), err=None, algorithm='normal', 
               class1='UNCERT', class1_prob=1.0, class2='UNCERT', 
               class2_prob=0.0, intensity=None, hardness=None, fluence=None, 
               sigma=None, rates=None, timescale=None, azzen=(None, None)):
        """Create a FswLocation object from keywords
        
        Args:
            time (float, optional): The time of the localization
            radec ((float, float), optional): RA and Dec
            err (float, optional): The localization error
            algorithm (str, optional): The algorithm spectrum; one of 'hard',
                                       'normal', or 'soft'
            class1 (str, optional): The primary classification
            class1_prob (float, optional): The primary classification 
                                           probability
            class2 (str, optional): The secondary classification
            class2_prob (float, optional): The secondary classification 
                                           probability
            intensity: (float, optional): Intensity of the signal in the 
                                          localization interval
            hardness: (float, optional): The hardness ratio in the 
                                         localization interval
            fluence: (float, optional): The fluence in the localization interval
            sigma: (float, optional): The S/N ratio of the localization interval
            rates: (np.array, optional): A 12-element array containing the rate
                                         in each NaI used for localization
            timescale (float, optional): The localization timescale
            azzen ((float, float), optional): Spacecraft azimuth and zenith of
                                              the localization
        
        Returns:
            (:class:`FswLocation`)
        """
        
        obj = cls()
        obj._time = time
        
        if not isinstance(radec, (list, tuple)):
            raise TypeError('radec must be a 2-tuple float')
        if len(radec) != 2:
            raise ValueError('radec must be a 2-tuple float')
        obj._location = (*radec, err)
        
        if algorithm not in spectrum:
            raise ValueError('algorithm is not valid')
        obj._algorithm = algorithm
        
        if class1 not in classifications.values():
            raise ValueError('class1 is not a valid class.')
        if class1_prob < 0.0 or class1_prob > 1.0:
            raise ValueError('class1_prob must be in range [0,1]')
        obj._class1 = (class1, class1_prob)
        
        if class2 not in classifications.values():
            raise ValueError('class2 is not a valid class.')
        if class2_prob < 0.0 or class2_prob > 1.0 or \
           class1_prob+class2_prob > 1.0:
            raise ValueError('class2_prob must be in range [0,1] and '\
                             'class1_prob + class2_prob must be <= 1')
        obj._class2 = (class2, class2_prob)
        
        obj._intensity = intensity
        obj._hardness = hardness
        obj._fluence = fluence
        obj._sigma = sigma
        obj._timescale = timescale
        
        if rates is None:
            rates = np.zeros(12)
        else:
            rates = np.asarray(rates)
            if rates.size != 12:
                raise ValueError('rates must be a 12-element array')
        obj._rates = rates
        
        if not isinstance(azzen, (list, tuple)):
            raise TypeError('azzen must be a 2-tuple float')
        if len(azzen) != 2:
            raise ValueError('azzen must be a 2-tuple float')
        obj._azzen = azzen
        
        return obj
        
    @classmethod
    def from_recarray(cls, rec_array):
        """Create from a FITS or numpy record array.
        
        Args:
            rec_array (np.recarray): The FITS OB_CALC record array from the 
                                     trigdat file
        
        Returns:
            (:class:`FswLocation`)
        """
        obj = cls()
        obj._time = rec_array['TIME']
        obj._location = (rec_array['RA'], rec_array['DEC'], rec_array['STATERR'])
        obj._algorithm = spectrum[rec_array['LOCALG'] - 1]
        obj._class1 = (classifications[rec_array['EVTCLASS'][0]],
                       rec_array['RELIABLT'][0])
        obj._class2 = (classifications[rec_array['EVTCLASS'][1]],
                        rec_array['RELIABLT'][1])
        obj._intensity = rec_array['INTNSITY']
        obj._hardness = rec_array['HDRATIO']
        obj._fluence = rec_array['FLUENCE']
        obj._sigma = rec_array['SIGMA']
        obj._rates = rec_array['LOCRATES']
        try:
            obj._timescale = rec_array['TRIG_TS']
        except KeyError:
            pass
        try:
            obj._azzen = (rec_array['TR_SCAZ'], rec_array['TR_SCZEN'])
        except KeyError:
            pass
        return obj
    
    def to_recarray(self):
        """Convert contents of the object to a numpy record array
        
        Returns:
            (np.recarray)
        """
        dtypes = self._dtypes
        algorithm = spectrum.index(self.spectrum)
        class_nums = {v: k for k, v in classifications.items()}
        record = np.array([(self.time, *self.location, algorithm, 
                           (class_nums[self._class1[0]], 
                            class_nums[self._class2[0]]),
                           (self._class1[1], self._class2[1]), self.intensity,
                           self.hardness_ratio, self.fluence, self.significance,
                           self.rates, self.timescale, *self.location_sc)], 
                           dtype=dtypes)
        return record
 
    def __repr__(self):
        s = '<FswLocation: {0}: {1:.1f}%;\n'.format(self.top_classification[0],
                                                self.top_classification[1]*100.)
        try:
            s += ' RA={0:.1f}, Dec={1:.1f}, err={2:.1f}>'.format(*self.location) 
        except:
            s += ' RA=None, Dec=None, err=None>'
        return s
