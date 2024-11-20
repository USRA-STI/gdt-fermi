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
from gdt.core.headers import Header, FileHeaders
from ..time import Time

__all__ = ['HealpixHeaders', 'PhaHeaders', 'PhaiiHeaders', 'PhaiiTriggerHeaders',
           'PosHistHeaders', 'RspHeaders', 'ScatHeaders', 'TcatHeaders',
           'TrigdatHeaders', 'TteHeaders', 'TteTriggerHeaders']

# mission definitions
_telescope = 'GLAST'
_instrument = 'GBM'
_observer = 'Meegan'
_origin = 'GIOC'
_timesys = 'TT'
_timeunit = 's'
_radecsys = 'FK5'
_equinox = 2000.0
_mjdrefi = 51910
_mjdreff = '7.428703703703703e-4'

# common keyword cards
_ancrfile_card = ('ANCRFILE', 'none', 'Name of corresponding ARF file (if any)')
_areascale_card = ('AREASCAL', 1.0, 'No special scaling of effective area by channel')
_backfile_card = ('BACKFILE', 'none', 'Name of corresponding background file (if any)')
_backscale_card = ('BACKSCAL', 1.0, 'No scaling of background')
_chantype_card = ('CHANTYPE', 'PHA', 'No corrections have been applied')
_corrfile_card = ('CORRFILE', 'none', 'Name of corresponding correction file (if any)')
_corrscale_card = ('CORRSCAL', 1., 'Correction scaling file')
_datatype_card = ('DATATYPE', '', 'GBM datatype used for this file')
_date_card = ('DATE', '', 'file creation date (YYYY-MM-DDThh:mm:ss UT)')
_date_end_card = ('DATE-END', '', 'Date of end of observation')
_date_obs_card = ('DATE-OBS', '', 'Date of start of observation')
_dec_obj_card = ('DEC_OBJ', 0.0, 'Calculated Dec of burst')
_detchans_card = ('DETCHANS', 0, 'Total number of channels in each rate')
_detnam_card = ('DETNAM', '', 'Individual detector name')
_equinox_card = ('EQUINOX', _equinox, 'Equinox for RA and Dec')
_err_rad_card = ('ERR_RAD', 0.0, 'Calculated Location Error Radius')
_extname_card = ('EXTNAME', '', 'name of this binary table extension')
_extver_card = ('EXTVER', 1, 'Version of this extension format')
_filename_card = ('FILENAME', '', 'Name of this file')
_filetype_card = ('FILETYPE', '', 'Name for this type of FITS file')
_filever_card = ('FILE-VER', '1.0.0', 'Version of the format for this filetype')
_filter_card = ('FILTER', 'none', 'The instrument filter in use (if any)')
_grouping_card = ('GROUPING', 0, 'No special grouping has been applied')
_hduclass_card = ('HDUCLASS', 'OGIP', 'Conforms to OGIP standard indicated in HDUCLAS1')
_hduvers_card = ('HDUVERS', '1.2.1', 'Version of HDUCLAS1 format in use')
_infile_card = ('INFILE01', '', 'Level 1 input lookup table file')
_instrument_card = ('INSTRUME', _instrument, 'Specific instrument used for observation')
_mjdrefi_card = ('MJDREFI', _mjdrefi, 'MJD of GLAST reference epoch, integer part')
_mjdreff_card = ('MJDREFF', _mjdreff, 'MJD of GLAST reference epoch, fractional part')
_object_card = ('OBJECT', '', 'Burst name in standard format, yymmddfff')
_observer_card = ('OBSERVER', _observer, 'GLAST Burst Monitor P.I.')
_origin_card = ('ORIGIN', _origin, 'Name of organization making file')
_poiserr_card = ('POISSERR', True, 'Assume Poisson Errors')
_radecsys_card = ('RADECSYS', _radecsys, 'Stellar reference frame')
_ra_obj_card = ('RA_OBJ', 0.0, 'Calculated RA of burst')
_respfile_card = ('RESPFILE', 'none', 'Name of corresponding RMF file (if any)')
_syserr_card = ('SYS_ERR', 0., 'No systematic errors')
_telescope_card = ('TELESCOP', _telescope, 'Name of mission/satellite')
_timesys_card = ('TIMESYS', _timesys, 'Time system used in time keywords')
_timeunit_card = ('TIMEUNIT', _timeunit, 'Time since MJDREF, used in TSTART and TSTOP')
_tstart_card = ('TSTART', 0.0, '[GLAST MET] Observation start time')
_tstop_card = ('TSTOP', 0.0, '[GLAST MET] Observation stop time')
_trigtime_card = ('TRIGTIME', 0.0, 'Trigger time relative to MJDREF, double precision')

#----------------
class GbmHeader(Header):

    def __setitem__(self, key, val):
        if not isinstance(key, tuple) and not isinstance(val, tuple):
            if key.upper() == 'TSTART':
                self['DATE-OBS'] = Time(val, format='fermi').isot
            elif key.upper() == 'TSTOP':
                self['DATE-END'] = Time(val, format='fermi').isot
            else:
                pass

            if 'INFILE' in key.upper():
                super(Header, self).__setitem__(key, val)
                return

        super().__setitem__(key, val)


class DataPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(), _filetype_card, _filever_card, 
                _telescope_card, _instrument_card, _detnam_card, _observer_card,
                _origin_card, _date_card, _date_obs_card, _date_end_card,
                _timesys_card, _timeunit_card,_mjdrefi_card,_mjdreff_card,
                _tstart_card, _tstop_card, _filename_card, _datatype_card]


class DataTriggerPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = DataPrimaryHeader.keywords + [_trigtime_card, _object_card, 
               _radecsys_card, _equinox_card, _ra_obj_card, _dec_obj_card,
               _err_rad_card]


class EboundsHeader(GbmHeader):
    name = 'EBOUNDS'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _hduclass_card,
                ('HDUCLAS1', 'RESPONSE', 'These are typically found in RMF ' \
                                         'files'),
                ('HDUCLAS2', 'EBOUNDS', 'From CAL/GEN/92-002'), _hduvers_card,
                _chantype_card, _filter_card, _detchans_card,
                ('CH2E_VER', '', 'Channel to energy conversion scheme used'),
                ('GAIN_COR', 0.0, 'Gain correction factor applied to energy ' \
                                  'edges')]
    

class EboundsTriggerHeader(GbmHeader):
    name = 'EBOUNDS'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _trigtime_card,
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card,
                _dec_obj_card, _err_rad_card, _hduclass_card,
                ('HDUCLAS1', 'RESPONSE', 'These are typically found in RMF ' \
                                         'files'),
                ('HDUCLAS2', 'EBOUNDS', 'From CAL/GEN/92-002'), _hduvers_card,
                _chantype_card, _filter_card, _detchans_card,
                ('CH2E_VER', '', 'Channel to energy conversion scheme used'),
                ('GAIN_COR', 0.0, 'Gain correction factor applied to energy ' \
                                  'edges')]


class SpectrumHeader(GbmHeader):
    name = 'SPECTRUM'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _filter_card, 
                _areascale_card, _backfile_card, _backscale_card, 
                _corrfile_card, _corrscale_card, _respfile_card, _ancrfile_card,
                _syserr_card, _poiserr_card, _grouping_card, _hduclass_card,
                ('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
                ('HDUCLAS2', 'TOTAL', 'Indicates gross data (source + background)'),
                ('HDUCLAS3', 'COUNT', 'Indicates data stored as counts'),
                ('HDUCLAS4', 'TYPEII', 'Indicates PHA Type II file format'),
                _hduvers_card, _chantype_card, _detchans_card, _extver_card]                   


class SpectrumTriggerHeader(GbmHeader):
    name = 'SPECTRUM'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _trigtime_card,
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card,
                _dec_obj_card, _err_rad_card, _filter_card, _areascale_card,
                _backfile_card, _backscale_card, _corrfile_card, _corrscale_card,
                _respfile_card, _ancrfile_card, _syserr_card, _poiserr_card,
                _grouping_card, _hduclass_card,
                ('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
                ('HDUCLAS2', 'TOTAL', 'Indicates gross data (source + background)'),
                ('HDUCLAS3', 'COUNT', 'Indicates data stored as counts'),
                ('HDUCLAS4', 'TYPEII', 'Indicates PHA Type II file format'),
                _hduvers_card, _chantype_card, _detchans_card, _extver_card]                   


class EventsHeader(GbmHeader):
    name = 'EVENTS'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card,
                ('EVT_DEAD', 2.6e-6, '[s] Deadtime per event'),
                ('EVTDEDHI', 1.0417e-5, '[s] Deadtime per overflow channel event'),
                _detchans_card, _hduclass_card,
                ('HDUCLAS1', 'EVENTS', 'Extension contains Events'),
                _extver_card]


class EventsTriggerHeader(GbmHeader):
    name = 'EVENTS'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _trigtime_card,
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card,
                _dec_obj_card, _err_rad_card, _respfile_card,
                ('EVT_DEAD', 2.6e-6, 'Deadtime per event (s)'),
                _detchans_card, _hduclass_card,
                ('HDUCLAS1', 'EVENTS', 'Extension contains Events'),
                _extver_card]


class GtiHeader(GbmHeader):
    name = 'GTI'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _hduclass_card,
                ('HDUCLAS1', 'GTI', 'Indicates good time intervals'),
                _hduvers_card, _extver_card]


class GtiTriggerHeader(GbmHeader):
    name = 'GTI'
    keywords = GtiHeader.keywords + [_trigtime_card, _object_card, 
               _radecsys_card, _equinox_card, _ra_obj_card, _dec_obj_card,
               _err_rad_card]


class RspPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(), _filetype_card, _filever_card, _date_card,
                _filename_card, _date_obs_card, _date_end_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _trigtime_card,
                _timesys_card, _timeunit_card, _telescope_card, 
                _instrument_card, _detnam_card, _observer_card, _origin_card,
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card, 
                _dec_obj_card, 
                ('DRM_NUM', 1, 'Number of DRMs stored in this file'),
                ('DRM_TYPE', '', 'Data type for which DRM is intended'),
                ('DIRDRMDB', '', 'Path to direct detector response'),
                ('DIRSCTDB', '', 'Path to atmospheric scatter')]


class RspEboundsHeader(EboundsTriggerHeader):
    keywords = EboundsTriggerHeader.keywords.copy()
    keywords.remove(_err_rad_card)


class SpecRespHeader(GbmHeader):
    name = 'SPECRESP MATRIX'
    keywords = [_extname_card, _extver_card, _date_card, _date_obs_card,
                _date_end_card, _mjdrefi_card, _mjdreff_card, _tstart_card,
                _tstop_card, _trigtime_card, _timesys_card, _timeunit_card,
                _telescope_card, _instrument_card, _detnam_card,
                ('MAT_TYPE', '', 'Response Matrix Type'), _observer_card,
                _origin_card,
                ('RSP_NUM', 1, 'Response matrix index number'), 
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card,
                _dec_obj_card,
                ('SRC_AZ', 0.0, 'Azimuth of source in spacecraft coordinates'),                         
                ('SRC_EL', 0.0, 'Elevation of source in spacecraft coordinates'),                           
                ('GEO_AZ', 0.0, 'Azimuth of geocenter in spacecraft coordinates'),
                ('GEO_EL', 0.0, 'Elevation of geocenter in spacecraft coordinates'),
                ('DET_ANG', 0.0, 'Angle between source and detector normal'),
                ('GEO_ANG', 0.0, 'Angle between geocenter and detector normal'),
                _filter_card, _chantype_card,
                ('NUMEBINS', 0, 'Number of true energy bins of the MATRIX'),
                _detchans_card,
                ('INFILE04', '', 'Atmospheric scattering database'),
                _hduclass_card, _hduvers_card,
                ('HDUCLAS1', 'RESPONSE', 'Typically found in RMF files'),
                ('HDUCLAS2', 'RSP_MATRIX', 'From CAL/GEN/92-002')]


class PhaSpectrumHeader(GbmHeader):
    name = 'SPECTRUM'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card, _trigtime_card,
                ('DATATYPE', '', 'GBM datatype used for this file'),
                _object_card, _radecsys_card, _equinox_card, _ra_obj_card,
                _dec_obj_card, _err_rad_card, _filter_card, _areascale_card,
                _backfile_card, _backscale_card, _corrfile_card,
                _corrscale_card, _respfile_card, _ancrfile_card, _syserr_card,
                _poiserr_card, _grouping_card, _hduclass_card,
                ('HDUCLAS1', 'SPECTRUM', 'PHA dataset (OGIP memo OGIP-92-007)'),
                ('HDUCLAS2', 'TOTAL', 'Indicates gross data (source + background)'),
                ('HDUCLAS3', 'RATE', 'Indicates data stored as rates'),
                ('HDUCLAS4', 'TYPEI', 'Indicates PHA Type I file format'),
                _hduvers_card, _chantype_card, _detchans_card,
                ('EXPOSURE', 0.0, 'Accumulation time - deadtime'), _extver_card]


class HealpixPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(),
                ('FILETYPE', 'IMAGE', 'Name for this type of FITS file'),
                _telescope_card, _instrument_card, _observer_card, _origin_card,
                _date_card, _date_obs_card, _date_end_card, _timesys_card,
                _timeunit_card, _mjdrefi_card, _mjdreff_card, _tstart_card,
                _tstop_card, _filename_card, _trigtime_card, _object_card,
                _radecsys_card, _equinox_card, _ra_obj_card, _dec_obj_card,
                _err_rad_card,
                ('THETA', 0.0, '[deg] Angle from spacecraft zenith'),
                ('PHI', 0.0, '[deg] Angle from spacecraft +X axis toward +Y'),
                ('LOC_SRC', 'Fermi, GBM', 'Mission/Instrument providing the localization'),
                ('CLASS', 'GRB', 'Classification of trigger'),
                ('OBJ_CLAS', 'GRB', 'Classification of trigger'),
                ('GEO_LONG', 0.0, '[deg] Spacecraft geographical east longitude'),
                ('GEO_LAT', 0.0, '[deg] Spacecraft geographical north latitude'),
                ('RA_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: RA'),
                ('DEC_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: Dec'),
                ('RA_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: RA'),
                ('DEC_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: Dec'),
                ('LOC_VER', '', 'Version string of localizing software'),
                ('LOC_ENRG', '(50, 300)', 'Energy range used for localization'),
                ('LMETHOD', 'Interactive', 'Method of localization')]


class HealpixImageHeader(GbmHeader):
    name = 'HEALPIX'
    keywords = [('PIXTYPE', 'HEALPIX', 'HEALPIX pixelization'),
                ('ORDERING', 'NESTED', 'Pixel ordering scheme, either RING or NESTED'),
                ('COORDSYS', 'C', 'Ecliptic, Galactic or Celestial (equatorial)'),
                _extname_card,
                ('NSIDE', 0, 'Resolution parameter of HEALPIX'),
                ('FIRSTPIX', 0, 'First pixel # (0 based)'),
                ('LASTPIX', 0, 'Last pixel # (0 based)'),
                ('INDXSCHM', 'IMPLICIT', 'Indexing: IMPLICIT or EXPLICIT'),
                ('OBJECT', 'FULLSKY', 'Sky coverage, either FULLSKY or PARTIAL'),
                ('SUN_RA', 0.0, 'RA of Sun'),
                ('SUN_DEC', 0.0, 'Dec of Sun'),
                ('GEO_RA', 0.0, 'RA of Geocenter relative to Fermi'),
                ('GEO_DEC', 0.0, 'Dec of Geocenter relative to Fermi'),
                ('GEO_RAD', 0.0, 'Radius of the Earth'),
                ('N0_RA', 0.0, 'RA pointing for detector n0'),
                ('N0_DEC', 0.0, 'Dec pointing for detector n0'),
                ('N1_RA', 0.0, 'RA pointing for detector n1'),
                ('N1_DEC', 0.0, 'Dec pointing for detector n1'),
                ('N2_RA', 0.0, 'RA pointing for detector n2'),
                ('N2_DEC', 0.0, 'Dec pointing for detector n2'),
                ('N3_RA', 0.0, 'RA pointing for detector n3'),
                ('N3_DEC', 0.0, 'Dec pointing for detector n3'),
                ('N4_RA', 0.0, 'RA pointing for detector n4'),
                ('N4_DEC', 0.0, 'Dec pointing for detector n4'),
                ('N5_RA', 0.0, 'RA pointing for detector n5'),
                ('N5_DEC', 0.0, 'Dec pointing for detector n5'),
                ('N6_RA', 0.0, 'RA pointing for detector n6'),
                ('N6_DEC', 0.0, 'Dec pointing for detector n6'),
                ('N7_RA', 0.0, 'RA pointing for detector n7'),
                ('N7_DEC', 0.0, 'Dec pointing for detector n7'),
                ('N8_RA', 0.0, 'RA pointing for detector n8'),
                ('N8_DEC', 0.0, 'Dec pointing for detector n8'),
                ('N9_RA', 0.0, 'RA pointing for detector n9'),
                ('N9_DEC', 0.0, 'Dec pointing for detector n9'),
                ('NA_RA', 0.0, 'RA pointing for detector na'),
                ('NA_DEC', 0.0, 'Dec pointing for detector na'),
                ('NB_RA', 0.0, 'RA pointing for detector nb'),
                ('NB_DEC', 0.0, 'Dec pointing for detector nb'),
                ('B0_RA', 0.0, 'RA pointing for detector b0'),
                ('B0_DEC', 0.0, 'Dec pointing for detector b0'),
                ('B1_RA', 0.0, 'RA pointing for detector b1'),
                ('B1_DEC', 0.0, 'Dec pointing for detector b1'),
                ('COMMENT', 'SCPOS: []'), ('COMMENT', 'QUAT: []')]


class ScatPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [_date_card, _filetype_card, Header.creator(), _origin_card,
                _telescope_card, _instrument_card, _observer_card, 
                _mjdrefi_card, _mjdreff_card, _timesys_card, _timeunit_card,
                _date_obs_card, _date_end_card, _tstart_card, _tstop_card,
                _trigtime_card, _filename_card, 
                ('COMMENT', 
                 'This file consists of time-sequenced spectral fit parameters')]


class ScatDetectorHeader(GbmHeader):
    name = 'DETECTOR DATA'
    keywords = [_extname_card, _date_card, _origin_card, _telescope_card, 
                _instrument_card, _observer_card, _mjdrefi_card, _mjdreff_card, 
                _timesys_card, _timeunit_card, _date_obs_card, _date_end_card, 
                _tstart_card, _tstop_card, _trigtime_card, 
                ('NUMFITS', 0, 'Number of spectral fits in the data')]


class ScatFitParamsHeader(GbmHeader):
    name = 'FIT PARAMS'
    keywords = [_extname_card, _date_card, _origin_card, _telescope_card, 
                _instrument_card, _observer_card, _mjdrefi_card, _mjdreff_card, 
                _timesys_card, _timeunit_card, _date_obs_card, _date_end_card, 
                _tstart_card, _tstop_card, _trigtime_card, 
                ('NUMFITS', 0, 'Number of spectral fits in the data'),
                ('N_PARAM', 0, 'Total number of fit parameters (PARAMn)'),
                ('FLU_LOW', 0.0, 'Lower limit of flux/fluence integration (keV)'),
                ('FLU_HIGH', 0.0, 'Upper limit of flux/fluence integration (keV)'),
                ('STATISTC', '', 'Indicates merit function used for fitting')]     


class TcatHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(), _filetype_card, _telescope_card, 
                _instrument_card, _observer_card, _origin_card, _date_card, 
                _date_obs_card, _date_end_card, _timesys_card, _timeunit_card,
                _mjdrefi_card, _mjdreff_card, _tstart_card, _tstop_card, 
                _filename_card, _trigtime_card, _object_card, _radecsys_card, 
                _equinox_card, _ra_obj_card, _dec_obj_card, _err_rad_card,
                ('THETA', 0.0, '[deg] Angle from spacecraft zenith'),
                ('PHI', 0.0, '[deg] Angle from spacecraft +X axis toward +Y'),
                ('LOC_SRC', 'Fermi, GBM', 
                 'Mission/Instrument providing the localization'),
                ('CLASS', '', 'Classification of trigger'),
                ('OBJ_CLAS', '', 'Classification of trigger'),
                ('TRIGSCAL', 0, '[ms] Triggered timescale'),
                ('TRIG_ALG', 0, 'Triggered algorithm number'),
                ('CHAN_LO', 0, 'Trigger channel: low'),
                ('CHAN_HI', 0, 'Trigger channel: high'),
                ('ADC_LO', 0, 'Trigger channel: low (ADC: 0 - 4095)'),
                ('ADC_HI', 0, 'Trigger channel: high (ADC: 0 - 4095)'),
                ('TRIG_SIG', 0.0, 'Trigger significance (sigma)'),
                ('GEO_LONG', 0.0, '[deg] Spacecraft geographical east longitude'),
                ('GEO_LAT', 0.0, '[deg] Spacecraft geographical north latitude'),
                ('DET_MASK', '', 'Triggered detectors: (0-13)'),
                ('RA_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: RA'),
                ('DEC_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: Dec'),
                ('RA_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: RA'),
                ('DEC_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: Dec'),
                ('INFILE01', '', 'Level 0 input data file'),
                ('LOC_VER', '', 'Version string of localizing software'),
                ('LOC_ENRG', '', 'Energy range used for localization'),
                ('RELIABLT', 0.0, 'Reliability of classification'),
                ('GCN_FLAG', '', ''),
                ('HISTORY', '()'),
                ('HISTORY', 'BKG_POLY_ORDER='),
                ('HISTORY', 'BKGINT1 = ()'),
                ('HISTORY', 'BKGINT2 = ()')]


class TrigdatPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(), _filetype_card, _filever_card, 
                _telescope_card, _instrument_card, _observer_card, _origin_card,
                _date_card, _date_obs_card, _date_end_card, _timesys_card,
                _timeunit_card, _mjdrefi_card, _mjdreff_card, _tstart_card,
                _tstop_card, ('DETTYPE', 'NAI', 'Detector type: NAI or BGO'),
                ('DATATYPE', 'TRIGDAT', 'Type of lookup table: CTIME or CSPEC'),
                _filename_card, _trigtime_card, _object_card, _radecsys_card,
                _equinox_card, _ra_obj_card, _dec_obj_card, _err_rad_card,
                ('TRIGSCAL', 0, '[ms] Triggered timescale'),
                ('TRIG_ALG', 0, 'Triggered algorithm number'),
                ('CHAN_LO', 0, 'Trigger channel: low'),
                ('CHAN_HI', 0, 'Trigger channel: high'),
                ('ADC_LO', 0, 'Trigger channel: low (ADC: 0 - 4095)'),
                ('ADC_HI', 0, 'Trigger channel: high (ADC: 0 - 4095)'),
                ('TRIG_SIG', 0.0, 'Trigger significance (sigma)'),
                ('GEO_LONG', 0.0, '[deg] Spacecraft geographical east longitude'),
                ('GEO_LAT', 0.0, '[deg] Spacecraft geographical north latitude'),
                ('DET_MASK', '', 'Triggered detectors: (0-13)'),
                ('RA_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: RA'),
                ('DEC_SCX', 0.0, '[deg] Pointing of spacecraft x-axis: Dec'),
                ('RA_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: RA'),
                ('DEC_SCZ', 0.0, '[deg] Pointing of spacecraft z-axis: Dec'),
                ('INFILE01', '', 'Level 0 input data file')]


class TrigdatSecondaryHeader(GbmHeader):
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card,
                _observer_card, _origin_card, _date_card, _date_obs_card,
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card,
                _mjdreff_card, _tstart_card, _tstop_card,
                ('DETTYPE', '', 'Detector type: NAI or BGO'),
                ('DATATYPE', 'TRIGDAT', 'Type of lookup table: CTIME or CSPEC'),
                _trigtime_card, _object_card, _radecsys_card, _equinox_card,
                _ra_obj_card, _dec_obj_card, _err_rad_card]                
    

class TrigdatTrigrateHeader(TrigdatSecondaryHeader):
    name = 'TRIGRATE'
    keywords = TrigdatSecondaryHeader.keywords.copy()
    keywords.remove(_detnam_card)


class TrigdatBackratesHeader(TrigdatSecondaryHeader):
    name = 'BCKRATES'


class TrigdatObCalcHeader(TrigdatSecondaryHeader):
    name = 'OB_CALC'


class TrigdatMaxratesHeader(TrigdatSecondaryHeader):
    name = 'MAXRATES'


class TrigdatEventrateHeader(TrigdatSecondaryHeader):
    name = 'EVNTRATE'


class PosHistPrimaryHeader(GbmHeader):
    name = 'PRIMARY'
    keywords = [Header.creator(), _filetype_card, _filever_card, 
                _telescope_card, _instrument_card, _detnam_card, _observer_card,
                _origin_card, _date_card, _date_obs_card, _date_end_card, 
                _timesys_card, _timeunit_card, _mjdrefi_card, _mjdreff_card,
                _tstart_card, _tstop_card, _filename_card]


class PosHistDataHeader(GbmHeader):
    name = 'GLAST POS HIST'
    keywords = [_extname_card, _telescope_card, _instrument_card, _detnam_card, 
                _observer_card, _origin_card, _date_card, _date_obs_card, 
                _date_end_card, _timesys_card, _timeunit_card, _mjdrefi_card, 
                _mjdreff_card, _tstart_card, _tstop_card, _extver_card]

#-------------------------------------

class PhaiiHeaders(FileHeaders):
    """FITS headers for continuous CTIME and CSPEC files"""
    _header_templates = [DataPrimaryHeader(), EboundsHeader(), SpectrumHeader(), 
                         GtiHeader()]


class PhaiiTriggerHeaders(FileHeaders):
    """FITS headers for trigger CTIME and CSPEC files"""
    _header_templates = [DataTriggerPrimaryHeader(), EboundsTriggerHeader(), 
                         SpectrumTriggerHeader(), GtiTriggerHeader()]


class TteHeaders(FileHeaders):
    """FITS headers for continuous TTE files"""
    _header_templates = [DataPrimaryHeader(), EboundsHeader(), EventsHeader(), 
                         GtiHeader()]


class TteTriggerHeaders(FileHeaders):
    """FITS headers for trigger TTE files"""
    _header_templates = [DataTriggerPrimaryHeader(), EboundsTriggerHeader(), 
                         EventsTriggerHeader(), GtiTriggerHeader()]


class PhaHeaders(FileHeaders):
    """FITS headers for PHA files"""
    _header_templates = [DataTriggerPrimaryHeader(), EboundsTriggerHeader(), 
                         PhaSpectrumHeader(), GtiTriggerHeader()]


class HealpixHeaders(FileHeaders):
    """FITS headers for localization HEALPix files"""
    _header_templates = [HealpixPrimaryHeader(), HealpixImageHeader()]


class ScatHeaders(FileHeaders):
    """FITS headers for spectral fits files"""
    _header_templates = [ScatPrimaryHeader(), ScatDetectorHeader(), 
                         ScatFitParamsHeader()]


class TcatHeaders(FileHeaders):
    """FITS headers for trigger catalog files"""
    _header_templates = [TcatHeader()]


class TrigdatHeaders(FileHeaders):
    """FITS headers for TRIGDAT files"""
    _header_templates = [TrigdatPrimaryHeader(), TrigdatTrigrateHeader(),
                         TrigdatBackratesHeader(), TrigdatObCalcHeader(),
                         TrigdatMaxratesHeader(), TrigdatEventrateHeader()]


class RspHeaders(FileHeaders):
    """FITS headers for detector response files"""
    _header_templates = [RspPrimaryHeader(), RspEboundsHeader(), 
                         SpecRespHeader()]


class PosHistHeaders(FileHeaders):
    """FITS headers for position history files"""
    _header_templates = [PosHistPrimaryHeader(), PosHistDataHeader()]
        
    
