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
import numpy as np

from gdt.core.response import *
from gdt.core.data_primitives import Ebounds, ResponseMatrix
from .detectors import GbmDetectors
from .headers import RspHeaders

__all__ = ['GbmRsp', 'GbmRsp2']

class GbmRsp(Rsp):
    """Class for GBM single-DRM response files.
    """
    @classmethod
    def open(cls, file_path, **kwargs):
        """Read a single-DRM response file from disk.

        Args:
            file_path (str): The file path 

        Returns:
           (:class:`GbmRsp`)
        """
        obj = super().open(file_path, **kwargs)
        if len(obj.hdulist) > 3:
            raise RuntimeError('{} is not a RSP file; it may be a RSP2 ' \
                               'file'.format(filename))
        trigtime = None

        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]
        headers = RspHeaders.from_headers(hdrs)
        tstart = headers['SPECRESP MATRIX']['TSTART']
        tstop = headers['SPECRESP MATRIX']['TSTOP']
        try:
            trigtime = headers['SPECRESP MATRIX']['TRIGTIME']
        except:
            pass
        
        ebounds = Ebounds.from_bounds(obj.column(1, 'E_MIN'), 
                                      obj.column(1, 'E_MAX'))

        num_ebins = obj.hdulist['SPECRESP MATRIX'].header['NUMEBINS']
        num_chans = obj.hdulist['SPECRESP MATRIX'].header['DETCHANS']
        fchan = np.copy(obj.column(2, 'F_CHAN'))
        nchan = np.copy(obj.column(2, 'N_CHAN'))
        ngrp = np.copy(obj.column(2, 'N_GRP'))
        matrix = cls._decompress_drm(obj.column(2, 'MATRIX'), num_ebins,
                                      num_chans, fchan, nchan)
            
        drm = ResponseMatrix(matrix, obj.column(2, 'ENERG_LO'),
                             obj.column(2, 'ENERG_HI'), obj.column(1, 'E_MIN'),
                             obj.column(1, 'E_MAX'))
        
        det = GbmDetectors.from_full_name(headers[0]['DETNAM'])
        
        obj.close()
        obj = cls.from_data(drm, filename=obj.filename,
                            start_time=tstart, stop_time=tstop, 
                            trigger_time=trigtime, headers=headers,
                            detector=det.name)

        obj._fchan = fchan
        obj._nchan = nchan
        obj._ngrp = ngrp

        return obj

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
        
        # the drm extension
        drm_hdu = self._drm_table()
        hdulist.append(drm_hdu)        
        
        return hdulist

    def _build_headers(self, num_chans, num_ebins):
        
        headers = self.headers.copy()
        headers['EBOUNDS']['DETCHANS'] = num_chans
        headers['SPECRESP MATRIX']['DETCHANS'] = num_chans
        headers['SPECRESP MATRIX']['NUMEBINS'] = num_ebins
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

    def _drm_table(self):
        elo_col = fits.Column(name='ENERG_LO', format='1E', 
                              array=self.drm.photon_bins.low_edges(),
                              unit='keV')
        ehi_col = fits.Column(name='ENERG_HI', format='1E', 
                              array=self.drm.photon_bins.high_edges(),
                              unit='keV')
        ngrp_col = fits.Column(name='N_GRP', format='1I', array=self._ngrp)
        fchan_col = fits.Column(name='F_CHAN', format='PI(1)', 
                                array=self._fchan)
        nchan_col = fits.Column(name='N_CHAN', format='PI(1)', 
                                array=self._nchan)
        matrix_col = fits.Column(name='MATRIX', array=self.drm.matrix,
                                 format='PE({})'.format(self.num_chans),
                                 unit='cm^2')
        
        hdu = fits.BinTableHDU.from_columns([elo_col, ehi_col, ngrp_col, 
                                             fchan_col, nchan_col, matrix_col],
                                         header=self.headers['SPECRESP MATRIX'])
        for key, val in self.headers['SPECRESP MATRIX'].items():
            hdu.header[key] = val

        return hdu

    # mark FIXME: Currently not used
    def _compress_drm(self, spec_index):
        """Compress a DRM by removing all the zeros.  This can result in a 
        file that is ~50% the size of an uncompressed DRM because the DRM is
        largely a triangle matrix.
        
        Also modifies the helper FITS columns that are used to decompress the
        matrix
        
        Args:
            spec_index (int): The index of the DRM in the event of multiple DRMs
        
        Returns:        
            np.array: An array of variable length arrays containing the compressed matrix
        """

        # function to split channels into contiguous groups
        def group_consecutive_indices(indices):
            # change points where consecutive indices > 1
            diff = indices[1:] - indices[0:-1]
            diff_idx = np.where(diff > 1)[0] + 1
            # add first, last indices
            diff_idx = np.concatenate(([0], diff_idx, [idx.size]))
            # group into consecutive indices
            groups = [idx[diff_idx[i]:diff_idx[i + 1]] \
                      for i in range(diff_idx.size - 1)]
            return groups

        drm = self._drm_list[spec_index]

        matrix = np.zeros((self.numebins,), dtype=np.object_)
        first_chans = []
        num_chans = []
        num_groups = []
        for ibin in range(self.numebins):
            idx = np.where(drm[ibin, :] > 0.0)[0]
            # if all zeros
            if idx.size == 0:
                num_groups.append(1)
                first_chans.append(np.array([self.numchans]))
                num_chans.append(np.array([1]))
                matrix[ibin] = np.array([0.0])
            else:
                groups = group_consecutive_indices(idx)
                num_groups.append(len(groups))
                first_chans.append(
                    np.array([group[0] + 1 for group in groups]))
                num_chans.append(np.array([len(group) for group in groups]))
                matrix[ibin] = np.concatenate(
                    [drm[ibin, group] for group in groups])

        self._ngrp_list[spec_index] = np.array(num_groups)
        self._fchan_list[spec_index] = first_chans
        self._nchan_list[spec_index] = num_chans

        return matrix

    # mark FIXME: This assumes NGRP=1, which for GBM is ok, not for general use
    @staticmethod
    def _decompress_drm(matrix, num_photon_bins, num_channels, _fchan, _nchan):
        """Decompresses a DRM using the standard F_CHAN, N_CHAN, and N_GRP
        keywords.
        
        Args:
            drm_data (np.recarray): The DRM data
        
        Returns:        
            (np.array)
        """
        # The format of the compress matrix is a series of groups, for each
        # energy bin, of channels with non-zero values.
        # fchan stands for the first channel of each of these groups
        # and nchan for the number of channels in the group group.
        # Each row in the matrix is a 1D list consisting on the contatenated
        # values of all groups for a given energy bin 
        # Note that in FITS the first index is 1
        drm = np.zeros((num_photon_bins, num_channels))
        for fchans, nchans, effective_area, drm_row \
            in zip(_fchan, _nchan, matrix, drm):

            channel_offset = 0

            for fchan, nchan in zip(fchans, nchans):
                
                start_idx = fchan - 1
                end_idx = start_idx + nchan

                drm_row[start_idx:end_idx] = \
                    effective_area[channel_offset:channel_offset + nchan]
                
                channel_offset += nchan

        return drm


# mark FIXME: interpolate method needs to update header
class GbmRsp2(Rsp2):
    """Class for GBM multiple-DRM response files.
    """    
    @classmethod
    def open(cls, file_path, **kwargs):
        """Read a multiple-DRM response file from disk.

        Args:
            file_path (str): The file path

        Returns:
           (:class:`GbmRsp2`)
        """
        obj = super().open(file_path, **kwargs)

        # get headers
        hdrs = [hdu.header for hdu in obj.hdulist]

        # extract responses into a list
        rsp_list = []
        num_drm = hdrs[0]['DRM_NUM']
        for i in range(num_drm):
            hdr = RspHeaders.from_headers([hdrs[0], hdrs[1], hdrs[i+2]])
            num_ebins = hdrs[i+2]['NUMEBINS']
            num_chans = hdrs[i+2]['DETCHANS']
            drm_data = obj.hdulist[i+2].data
            matrix = drm_data['MATRIX']
            fchan = drm_data['F_CHAN']
            nchan = drm_data['N_CHAN']
            ngrp = drm_data['N_GRP']
            matrix = GbmRsp._decompress_drm(matrix, num_ebins, num_chans,
                                            fchan, nchan)
            
            drm = ResponseMatrix(matrix, drm_data['ENERG_LO'],
                                 drm_data['ENERG_HI'], obj.column(1, 'E_MIN'),
                                 obj.column(1, 'E_MAX'))
            
            det = GbmDetectors.from_full_name(hdrs[i+2]['DETNAM'])
            
            rsp = GbmRsp.from_data(drm, start_time=hdrs[i+2]['TSTART'],
                                   stop_time=hdrs[i+2]['TSTOP'],
                                   trigger_time=hdrs[i+2]['TRIGTIME'],
                                   headers=hdr, detector=det.name)
            rsp._fchan = fchan
            rsp._nchan = nchan
            rsp._ngrp = ngrp
            rsp_list.append(rsp)
        
        obj.close()
        
        return cls.from_rsps(rsp_list, filename=obj.filename)                    
