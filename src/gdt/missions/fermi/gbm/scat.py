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
from gdt.core.data_primitives import Parameter
from gdt.core.file import FitsFileContextManager
from gdt.core.spectra.parameters import *
from .headers import ScatHeaders

__all__ = ['Scat', 'GbmDetectorData', 'GbmModelFit']


class GbmModelFit(ModelFit):
    """A container for the info from a model fit, with values used in the 
    GBM SCAT files.
    """

    def __init__(self):
        super().__init__()
        self._photon_flux_50_300 = None
        self._energy_fluence_50_300 = None
        self._duration_fluence = None

    @property
    def duration_fluence(self):
        """(:class:`~gdt.spectra.parameters.EnergyFluence`): The energy fluence 
        over the duration energy range, nominally 50-300 keV."""
        return self._duration_fluence

    @duration_fluence.setter
    def duration_fluence(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('duration_fluence must be of Parameter type')
        self._duration_fluence = val

    @property
    def energy_fluence_50_300(self):
        """(:class:`~gdt.spectra.parameters.EnergyFluence`): The energy fluence 
        over 50-300 keV"""
        return self._energy_fluence_50_300

    @energy_fluence_50_300.setter
    def energy_fluence_50_300(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('energy_fluence_50_300 must be of Parameter type')
        self._energy_fluence_50_300 = val

    @property
    def photon_flux_50_300(self):
        """(:class:`~gdt.spectra.parameters.PhotonFlux`): The photon flux over 
        50-300 keV"""
        return self._photon_flux_50_300

    @photon_flux_50_300.setter
    def photon_flux_50_300(self, val):
        if not isinstance(val, Parameter):
            raise TypeError('photon_flux_50_300 must be of Parameter type')
        self._photon_flux_50_300 = val

    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
            (astropy.io.fits.BinTableHDU): The FITS table
        """
        numparams = len(self.parameters)
        cols = []
        cols.append(
            fits.Column(name='TIMEBIN', format='2D', array=[self.time_range]))
        i = 0
        for param in self.parameters:
            col = fits.Column(name='PARAM{0}'.format(i), format='3E',
                              array=[param.to_fits_value()])
            cols.append(col)
            i += 1

        cols.append(fits.Column(name='PHTFLUX', format='3E',
                                array=[self.photon_flux.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLNC', format='3E',
                                array=[self.photon_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLUX', format='3E',
                                array=[self.energy_flux.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNC', format='3E',
                                array=[self.energy_fluence.to_fits_value()]))
        cols.append(fits.Column(name='PHTFLUXB', format='3E',
                                array=[
                                    self.photon_flux_50_300.to_fits_value()]))
        cols.append(fits.Column(name='DURFLNC', format='3E',
                                array=[self.duration_fluence.to_fits_value()]))
        cols.append(fits.Column(name='NRGFLNCB', format='3E',
                                array=[
                                    self.energy_fluence_50_300.to_fits_value()]))
        cols.append(fits.Column(name='REDCHSQ', format='2E',
                                array=[[self.stat_value / self.dof] * 2]))
        cols.append(fits.Column(name='CHSQDOF', format='1I', array=[self.dof]))
        cols.append(fits.Column(name='COVARMAT',
                                format='{0}E'.format(int(numparams * numparams)),
                                dim='({0},{0})'.format(numparams),
                                array=[self.covariance]))

        hdu = fits.BinTableHDU.from_columns(cols, name='FIT PARAMS')
        return hdu

    @classmethod
    def from_fits_row(cls, fits_row, model_name, param_names=None,
                      flux_range=(10.0, 1000.0), dur_range=(50.0, 300.0)):
        """Read a FITS row and return a :class:`GbmModelFit` object
        
        Args:
            fits_row (np.recarray): The FITS row
            model_name (str): The model name
            param_names (list, optional): The list of parameter names
            flux_range (tuple, optional): The flux energy range (low, high)
            dur_range (tuple, optional): The duration energy range (low, high)
        
        Returns:
            (:class:`GbmModelFit`)
        """
        time_range = tuple(fits_row['TIMEBIN'])
        nparams = sum([1 for name in fits_row.array.dtype.names if 'PARAM' in name])
        if param_names is None:
            param_names = [''] * nparams

        params = []
        for i in range(nparams):
            param = fits_row['PARAM' + str(i)]
            params.append(Parameter(param[0], tuple(param[1:]),
                                    name=param_names[i]))

        pflux = PhotonFlux(fits_row['PHTFLUX'][0],
                           tuple(fits_row['PHTFLUX'][1:]), flux_range)
        pflnc = PhotonFluence(fits_row['PHTFLNC'][0],
                              tuple(fits_row['PHTFLNC'][1:]), flux_range)
        eflux = EnergyFlux(fits_row['NRGFLUX'][0],
                           tuple(fits_row['NRGFLUX'][1:]), flux_range)
        eflnc = EnergyFluence(fits_row['NRGFLNC'][0],
                              tuple(fits_row['NRGFLNC'][1:]), flux_range)
        pfluxb = PhotonFlux(fits_row['PHTFLUXB'][0],
                            tuple(fits_row['PHTFLUXB'][1:]), (50.0, 300.0))
        eflncb = EnergyFluence(fits_row['NRGFLNCB'][0],
                               tuple(fits_row['NRGFLNCB'][1:]), (50.0, 300.0))
        durflnc = PhotonFluence(fits_row['DURFLNC'][0],
                                tuple(fits_row['DURFLNC'][1:]), dur_range)
        dof = fits_row['CHSQDOF']
        # scat provides the [fit stat, chisq], while bcat is only the fit stat
        try:
            stat_val = fits_row['REDCHSQ'][1] * dof
        except:
            stat_val = fits_row['REDCHSQ'] * dof

        num_params = len(params)
        covar = fits_row['COVARMAT'].reshape(num_params, num_params)

        obj = cls.from_data(model_name, time_range, parameters=params,
                            photon_flux=pflux, photon_fluence=pflnc,
                            energy_flux=eflux, energy_fluence=eflnc,
                            flux_energy_range=flux_range, stat_value=stat_val,
                            dof=dof, covariance=covar,
                            photon_flux_50_300=pfluxb,
                            energy_fluence_50_300=eflncb,
                            duration_fluence=durflnc)
        return obj


class GbmDetectorData(DetectorData):
    """A container for the detector info, with values used in the GBM SCAT files.
    """

    def to_fits_row(self):
        """Return the contained data as a FITS table row
        
        Returns:
            (astropy.io.fits.BinTableHDU)
        """
        numchans = len(self.energy_edges)
        e_dim = str(numchans) + 'E'
        p_dim = str(numchans - 1) + 'E'
        p_unit = 'Photon cm^-2 s^-1 keV^-1'
        fit_int = '{0}: {1} s, '.format(self.time_range[0], self.time_range[1])
        fit_int += '{0}: {1} keV, '.format(self.energy_range[0],
                                           self.energy_range[1])
        fit_int += 'channels {0}: {1}'.format(self.channel_range[0],
                                              self.channel_range[1])

        col1 = fits.Column(name='INSTRUME', format='20A',
                           array=[self.instrument])
        col2 = fits.Column(name='DETNAM', format='20A', array=[self.detector])
        col3 = fits.Column(name='DATATYPE', format='20A',
                           array=[self.datatype])
        col4 = fits.Column(name='DETSTAT', format='20A',
                           array=['INCLUDED' if self.active else 'OMITTED'])
        col5 = fits.Column(name='DATAFILE', format='60A',
                           array=[self.filename])
        col6 = fits.Column(name='RSPFILE', format='60A', array=[self.response])
        col7 = fits.Column(name='FIT_INT', format='60A', array=[fit_int])
        col8 = fits.Column(name='CHANNUM', format='1J', array=[numchans - 1])
        col9 = fits.Column(name='FITCHAN', format='{}J'.format(numchans),
                           array=[self.channel_mask])
        col10 = fits.Column(name='E_EDGES', format=e_dim, unit='keV',
                            array=[self.energy_edges])
        col11 = fits.Column(name='PHTCNTS', format=p_dim, unit=p_unit,
                            array=[self.photon_counts])
        col12 = fits.Column(name='PHTMODL', format=p_dim, unit=p_unit,
                            array=[self.photon_model])
        col13 = fits.Column(name='PHTERRS', format=p_dim, unit=p_unit,
                            array=[self.photon_errors])
        hdu = fits.BinTableHDU.from_columns(
            [col1, col2, col3, col4, col5, col6,
             col7, col8, col9, col10, col11,
             col12, col13], name='DETECTOR DATA')
        return hdu

    @classmethod
    def from_fits_row(cls, fits_row):
        """Read a FITS row and return a GbmDetectorData object.
        
        Returns:
            (:class:`GbmDetectorData`)
        """
        instrument = fits_row['INSTRUME']
        det = fits_row['DETNAM']
        datatype = fits_row['DATATYPE']
        if fits_row['DETSTAT'] == 'INCLUDED':
            active = True
        else:
            active = False
        filename = fits_row['DATAFILE']
        rspfile = fits_row['RSPFILE']
        fit_ints = fits_row['FIT_INT'].split()
        time_range = (float(fit_ints[0][:-1]), float(fit_ints[1]))
        energy_range = (float(fit_ints[3][:-1]), float(fit_ints[4]))
        channel_range = (int(fit_ints[7][:-1]), int(fit_ints[8]))
        numchans = fits_row['CHANNUM']
        channel_mask = np.zeros(numchans, dtype=bool)
        channel_mask[fits_row['FITCHAN'][0]:fits_row['FITCHAN'][1] + 1] = True
        energy_edges = fits_row['E_EDGES']
        photon_counts = fits_row['PHTCNTS']
        photon_model = fits_row['PHTMODL']
        photon_errors = fits_row['PHTERRS']
        obj = cls.from_data(instrument, det, datatype, numchans,
                            filename=filename, active=active, response=rspfile,
                            time_range=time_range, energy_range=energy_range,
                            channel_range=channel_range,
                            channel_mask=channel_mask,
                            energy_edges=energy_edges,
                            photon_counts=photon_counts,
                            photon_model=photon_model,
                            photon_errors=photon_errors)
        return obj


class Scat(FitsFileContextManager):
    """A container class for the spectral fit data in an SCAT (Spectroscopy 
    CATalog) file.
    """

    def __init__(self):
        super().__init__()
        self._detectors = []
        self._model_fits = []

    @property
    def detectors(self):
        """(list): The :class:`GbmDetectorData` objects used in the analysis"""
        return self._detectors

    @property
    def model_fits(self):
        """(list): The :class:`GbmModelFit` objects, one for each model fit"""
        return self._model_fits

    @property
    def num_detectors(self):
        """(int): The number of detectors in the SCAT file"""
        return len(self.detectors)

    @property
    def num_fits(self):
        """(int): The number of model fits"""
        return len(self.model_fits)

    def add_detector_data(self, detector_data):
        """Add a new detector to the Scat
        
        Args:
            detector_data (:class:`GbmDetectorData`): The detector data
        """
        if not isinstance(detector_data, GbmDetectorData):
            raise TypeError("Can only add GbmDetectorData objects")
        self._detectors.append(detector_data)

    def add_model_fit(self, model_fit):
        """Add a new model fit to the Scat
        
        Args:
            model_fit (:class:`GbmModelFit`): The model fit data
        """
        if not isinstance(model_fit, GbmModelFit):
            raise TypeError("Can only add GbmModelFit objects")
        self._model_fits.append(model_fit)

    @classmethod
    def from_data(cls, detector_list, model_fits, headers=None, filename=None):
        """Create a Scat object from data.
        
        Args:
            detector_list (list of :class:`GbmDetectorData`): The detector info 
                                                              used in the fit
            model_fits (list of :class:`GbmModelFit`): The model fits
            headers (:class:`~.headers.ScatHeaders`, optional): The file headers
            filename (str, optional): The filename
        
        Returns:
            (:class:`Scat`)
        """
        obj = cls()
        obj._filename = filename

        try:
            iter(detector_list)
        except:
            raise TypeError('detector_list must be a list of DetectorData objects')
        for det in detector_list:
            obj.add_detector_data(det)

        try:
            iter(model_fits)
        except:
            raise TypeError('model_fits must be a list of ModelFit objects')
        for model_fit in model_fits:
            obj.add_model_fit(model_fit)

        if not isinstance(headers, ScatHeaders):
            raise TypeError('headers must be a ScatHeaders object')
        obj._headers = headers

        return obj

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a SCAT FITS file and create a Scat object.
        
        Args:
            file_path (str): The file path of the FITS file
        
        Returns:
            (:class:`Scat`)
        """
        with super().open(file_path, **kwargs) as obj:
            hdrs = [hdu.header for hdu in obj.hdulist]
            headers = ScatHeaders.from_headers(hdrs)

            # read the detector data HDU
            det_list = []
            for row in obj.hdulist['DETECTOR DATA'].data:
                det_list.append(GbmDetectorData.from_fits_row(row))

            # read the fit params HDU
            model_fits = cls._from_fitparam_hdu(obj.hdulist['FIT PARAMS'])

        obj = cls.from_data(det_list, model_fits, headers=headers,
                            filename=obj.filename)

        return obj

    @staticmethod
    def _from_fitparam_hdu(fit_hdu):
        """format the data from the fitparam HDU so the GbmModelFit list can be
        created."""

        fit_data = fit_hdu.data
        fit_hdr = fit_hdu.header
        nparams = fit_hdr['N_PARAM']

        # find unique models in the event there are multiple components
        models = [fit_hdr.comments['TTYPE' + str(2 + i)].split(':')[0]
                  for i in range(nparams)]
        model = '+'.join(list(set(models)))

        # the parameter names
        param_names = [
            fit_hdr.comments['TTYPE' + str(2 + i)].split(':')[1].strip()
            for i in range(nparams)]

        # populate each model fit
        model_fits = []
        for row in fit_data:
            modelfit = GbmModelFit.from_fits_row(row, model,
                                                 param_names=param_names)
            modelfit.stat_name = fit_hdr['STATISTC']
            model_fits.append(modelfit)

        return model_fits

    def _build_hdulist(self):

        # create FITS and primary header
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)

        # the detector data extension
        det_hdu = self._det_table()
        hdulist.append(det_hdu)

        # the model fit extension
        model_hdu = self._model_table()
        hdulist.append(model_hdu)

        return hdulist

    def _det_table(self):
        det_data = np.concatenate([det.to_fits_row().data \
                                   for det in self.detectors])
        adet = self.detectors[0].to_fits_row()
        cols = []
        for column in adet.columns:
            col = fits.Column(name=column.name, format=column.format,
                              unit=column.unit, array=det_data[column.name])
            cols.append(col)

        det_hdu = fits.BinTableHDU.from_columns(cols, header=self.headers[1])
        for key, val in self.headers[1].items():
            det_hdu.header[key] = val

        det_hdu.header.comments['TTYPE1'] = 'Instrument name for this detector'
        det_hdu.header.comments['TTYPE2'] = \
            'Detector number; if one of several available'
        det_hdu.header.comments['TTYPE3'] = 'Data type used for this analysis'
        det_hdu.header.comments['TTYPE4'] = \
            'Was this detector INCLUDED or OMITTED'
        det_hdu.header.comments['TTYPE5'] = 'Data file name for this dataset'
        det_hdu.header.comments['TTYPE6'] = \
            'Response file name for this dataset'
        det_hdu.header.comments['TTYPE7'] = 'Fit intervals'
        det_hdu.header.comments['TTYPE8'] = \
            'Total number of energy channels for this detector'
        det_hdu.header.comments['TTYPE9'] = \
            'Channels selected in fitting this detector'
        det_hdu.header.comments['TTYPE10'] = \
            'Energy edges for each selected detector'
        det_hdu.header.comments['TTYPE11'] = 'Array of photon counts data'
        det_hdu.header.comments['TTYPE12'] = 'Array of photon model data'
        det_hdu.header.comments['TTYPE13'] = \
            'Array of errors in photon counts data'
        det_hdu.header['NUMFITS'] = (self.num_fits,
                                     'Number of spectral fits in the data')

        return det_hdu

    def _model_table(self):

        model_data = np.concatenate([model.to_fits_row().data \
                                     for model in self.model_fits])
        amodel = self.model_fits[0].to_fits_row()
        cols = []
        for column in amodel.columns:
            col = fits.Column(name=column.name, format=column.format,
                              unit=column.unit, array=model_data[column.name])
            cols.append(col)

        model_hdu = fits.BinTableHDU.from_columns(cols, header=self.headers[2])
        for key, val in self.headers[2].items():
            model_hdu.header[key] = val

        e_range = self.model_fits[0].flux_energy_range
        model_name = self.model_fits[0].name
        param_names = self.model_fits[0].parameter_list()
        numparams = len(param_names)
        statistic = self.model_fits[0].stat_name

        g_range = '({0}-{1} keV)'.format(e_range[0], e_range[1])
        b_range = '(50-300 keV)'

        model_hdu.header.comments['TTYPE1'] = \
            'Start and stop times relative to trigger'
        for i in range(numparams):
            colname = 'TTYPE{0}'.format(i + 2)
            model_hdu.header.comments[colname] = '{0}: {1}'.format(model_name,
                                                                   param_names[i])

        ttypes = ['TTYPE' + str(numparams + 2 + i) for i in range(10)]
        model_hdu.header.comments[ttypes[0]] = \
            'Photon Flux (ph/s-cm^2) std energy ' + g_range
        model_hdu.header.comments[ttypes[1]] = \
            'Photon Fluence (ph/cm^2) std energy ' + g_range
        model_hdu.header.comments[ttypes[2]] = \
            'Energy Flux (erg/s-cm^2) std energy ' + g_range
        model_hdu.header.comments[ttypes[3]] = \
            'Energy Fluence (erg/cm^2) std energy ' + g_range
        model_hdu.header.comments[ttypes[4]] = \
            'Photon Flux (ph/s-cm^2) BATSE energy ' + b_range
        model_hdu.header.comments[ttypes[5]] = \
            'Photon Fluence (ph/cm^2) BATSE energy ' + b_range
        model_hdu.header.comments[ttypes[6]] = \
            'Energy Fluence (erg/cm^2) BATSE energy ' + b_range
        model_hdu.header.comments[ttypes[7]] = \
            'Reduced Chi^2 (1) and fitting statistic (2)'
        model_hdu.header.comments[ttypes[8]] = 'Degrees of Freedom'
        model_hdu.header.comments[ttypes[9]] = \
            'Covariance matrix for the fir (N_PARAM^2)'

        model_hdu.header['N_PARAM'] = (numparams,
                                       'Total number of fit parameters (PARAMn)')
        model_hdu.header['FLU_LOW'] = (e_range[0],
                                       'Lower limit of flux/fluence integration (keV)')
        model_hdu.header['FLU_HIGH'] = (e_range[1],
                                        'Upper limit of flux/fluence integration (keV)')
        model_hdu.header['STATISTC'] = (statistic,
                                        'Indicates merit function used for fitting')
        model_hdu.header['NUMFITS'] = (self.num_fits,
                                       'Number of spectral fits in the data')

        return model_hdu

    def __repr__(self):
        s = '<{0}: {1};\n'.format(self.__class__.__name__, self.filename)
        s += ' {0} detectors; {1} fits>'.format(self.num_detectors,
                                                self.num_fits)
        return s
