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
from gdt.core.heasarc import BaseFinder

__all__ = ['TriggerFinder', 'ContinuousFinder', 'TriggerFtp', 'ContinuousFtp']

class GbmFinder(BaseFinder):
    """Subclassing BaseFinder to enable _file_filter() to take a list of
    GBM detectors.
    """    
    def _file_filter(self, file_list, filetype, extension, dets=None):
        """Filters the directory for the requested filetype, extension, and 
        detectors
        
        Args:
            filetype (str): The type of file, e.g. 'cspec'
            extension (str): The file extension, e.g. '.pha'
            dets (list, optional): The detectors. If omitted, then files for 
                                   all detectors are returned

        Returns:
            list: The filtered file list
        """
        files = super()._file_filter(file_list, filetype, extension)
        
        if dets is not None:
            if type(dets) == str:
                dets = [dets]
            files = [f for f in files if
                     any('_' + det + '_' in f for det in dets)]

        return files
    

class TriggerFinder(GbmFinder):
    """A class that interfaces with the HEASARC trigger directories.
    An instance of this class will represent the available files associated
    with a single trigger.
    
    An instance can be created without a trigger number, however a trigger
    number will need to be set by :meth:`cd(tnum) <cd>` to query and download files.
    An instance can also be changed from one trigger number to another without
    having to create a new instance.  If multiple instances are created with
    the keyword protocol='FTP' and exist simultaneously, they will all use a
    single FTP connection.
        
    Parameters:
        tnum (str, optional): A valid trigger number
    """
    _root = '/fermi/data/gbm/triggers'
    
    def cd(self, tnum):
        """Change directory to new trigger.
        
        Args:
            tnum (str): The trigger number
        """
        super().cd(tnum)
        
    def get_cat_files(self, download_dir, **kwargs):
        """Download all catalog files for the trigger.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = []
        files.extend(self._file_filter(self.files, 'bcat', 'fit'))
        files.extend(self._file_filter(self.files, 'scat', 'fit'))
        files.extend(self._file_filter(self.files, 'tcat', 'fit'))
        return self.get(download_dir, files, **kwargs)

    def get_cspec(self, download_dir, dets=None, **kwargs):
        """Download the CSPEC files for the trigger.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'cspec', 'pha', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def get_ctime(self, download_dir, dets=None, **kwargs):
        """Download the CTIME files for the trigger.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'ctime', 'pha', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def get_healpix(self, download_dir, **kwargs):
        """Download the healpix localization file for the trigger.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'healpix', 'fit')
        return self.get(download_dir, files, **kwargs)

    def get_localization(self, download_dir, **kwargs):
        """Download all localization files for the trigger.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = []
        files.extend(self._file_filter(self.files, 'healpix', 'fit'))
        files.extend(self._file_filter(self.files, 'skymap', 'png'))
        files.extend(self._file_filter(self.files, 'loclist', 'txt'))
        files.extend(self._file_filter(self.files, 'locprob', 'fit'))
        files.extend(self._file_filter(self.files, 'locplot', 'png'))
        return self.get(download_dir, files, **kwargs)

    def get_lightcurve(self, download_dir, **kwargs):
        """Download the lightcurve plots for the trigger.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'lc', 'pdf')
        return self.get(download_dir, files, **kwargs)

    def get_rsp(self, download_dir, ctime=True, cspec=True, dets=None,
                **kwargs):
        """Download the response Type-I files for the trigger.
        
        Args:
            download_dir (str): The download directory
            ctime (bool, optional): If True, download the ctime responses. 
                                    Default is True.
            cspec (bool, optional): If True, download the cspec responses. 
                                    Default is True.
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = []
        if cspec:
            files.extend(
                self._file_filter(self.files, 'cspec', 'rsp', dets=dets))
        if ctime:
            files.extend(
                self._file_filter(self.files, 'ctime', 'rsp', dets=dets))
        return self.get(download_dir, files, **kwargs)

    def get_rsp2(self, download_dir, ctime=True, cspec=True, dets=None,
                 **kwargs):
        """Download the response Type-I files for the trigger.
        
        Args:
            download_dir (str): The download directory
            ctime (bool, optional): If True, download the ctime responses. 
                                    Default is True.
            cspec (bool, optional): If True, download the cspec responses. 
                                    Default is True.
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = []
        if cspec:
            files.extend(
                self._file_filter(self.files, 'cspec', 'rsp2', dets=dets))
        if ctime:
            files.extend(
                self._file_filter(self.files, 'ctime', 'rsp2', dets=dets))
        return self.get(download_dir, files, **kwargs)

    def get_trigdat(self, download_dir, **kwargs):
        """Download the trigger data (trigdat) file for the trigger.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'trigdat', 'fit')
        return self.get(download_dir, files, **kwargs)

    def get_tte(self, download_dir, dets=None, **kwargs):
        """Download the TTE files for the trigger.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'tte', 'fit', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def ls_cat_files(self):
        """List all catalog files for the trigger.

        Returns:
            (list of str)
        """
        files = []
        files.extend(self._file_filter(self.files, 'bcat', 'fit'))
        files.extend(self._file_filter(self.files, 'scat', 'fit'))
        files.extend(self._file_filter(self.files, 'tcat', 'fit'))
        return files

    def ls_cspec(self):
        """List all CSPEC files for the trigger.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'cspec', 'pha')

    def ls_ctime(self):
        """List all CTIME files for the trigger.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'ctime', 'pha')

    def ls_lightcurve(self):
        """List all lightcurve plots for the trigger.

        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'lc', 'pdf')

    def ls_localization(self):
        """List all localization files for the trigger.

        Returns:
            (list of str)
        """
        files = []
        files.extend(self._file_filter(self.files, 'healpix', 'fit'))
        files.extend(self._file_filter(self.files, 'skymap', 'png'))
        files.extend(self._file_filter(self.files, 'loclist', 'txt'))
        files.extend(self._file_filter(self.files, 'locprob', 'fit'))
        files.extend(self._file_filter(self.files, 'locplot', 'png'))
        return files

    def ls_rsp(self, ctime=True, cspec=True):
        """List all response Type-I files for the trigger.

        Args:
            ctime (bool, optional): If True, list the ctime responses. 
                                    Default is True.
            cspec (bool, optional): If True, list the cspec responses. 
                                    Default is True.
        
        Returns:
            (list of str)
        """
        files = []
        if cspec:
            files.extend(self._file_filter(self.files, 'cspec', 'rsp'))
        if ctime:
            files.extend(self._file_filter(self.files, 'ctime', 'rsp'))
        return files

    def ls_rsp2(self, ctime=True, cspec=True):
        """List all response Type-II files for the trigger.

        Args:
            ctime (bool, optional): If True, list the ctime responses. 
                                    Default is True.
            cspec (bool, optional): If True, list the cspec responses. 
                                    Default is True.

        Returns:
            (list of str)
        """
        files = []
        if cspec:
            files.extend(self._file_filter(self.files, 'cspec', 'rsp2'))
        if ctime:
            files.extend(self._file_filter(self.files, 'ctime', 'rsp2'))
        return files

    def ls_trigdat(self):
        """List the trigger data (trigdat) file for the trigger.

        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'trigdat', 'fit')

    def ls_tte(self):
        """List all TTE files for the trigger.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'tte', 'fit')

    def _construct_path(self, str_trigger_num):
        """Constructs the path for a trigger
        
        Args:
            str_trigger_num (str): The trigger number

        Returns:
            str: The path of the directory for the trigger
        """
        year = '20' + str_trigger_num[0:2]
        path = os.path.join(self._root, year, 'bn' + str_trigger_num,
                            'current')
        return path


class TriggerFtp(TriggerFinder):
    """Class providing backwards compatibility for code written prior to v2.0.5
    where the TriggerFtp handled interactions with triggered data.

    Parameters:
        args: The set of parameters needed to define the data path
        **kwargs: Options passed to :class:`TriggerFinder` class
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        super().__init__(*args, protocol='FTP', **kwargs)


class ContinuousFinder(GbmFinder):
    """A class that interfaces with the HEASARC continuous daily data
    directories. An instance of this class will represent the available files 
    associated with a single day.
    
    An instance can be created without a time, however a time will need to be 
    set by :meth:`cd(time) <cd>` to query and download files. An instance can also be 
    changed from one time to another without having to create a new instance.  
    If multiple instances are created with the keyword protocol='FTP' and exist 
    simultaneously, they will all use a single FTP connection.
    
 
    Parameters:
        time (astropy.time.Time, optional): The time object
    """
    _root = '/fermi/data/gbm/daily'

    def cd(self, time):
        """Change directory to new date.
        
        Args:
            time (astropy.time.Time): The time
        """
        super().cd(time)

    def get_all(self, download_dir, **kwargs):
        """Download all files within a daily directory.

        Note:
            Use at your own risk. Unless you have a high-bandwidth connection 
            and can handle downloading several GBs, this function is not
            recommended for use.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        return self.get(download_dir, self._file_list, **kwargs)

    def get_cspec(self, download_dir, dets=None, **kwargs):
        """Download the CSPEC files.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'cspec', 'pha', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def get_ctime(self, download_dir, dets=None, **kwargs):
        """Download the CTIME files.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'ctime', 'pha', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def get_poshist(self, download_dir, **kwargs):
        """Download the poshist file.
        
        Args:
            download_dir (str): The download directory
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self.ls_poshist()
        return self.get(download_dir, files, **kwargs)

    def get_spechist(self, download_dir, dets=None, **kwargs):
        """Download the spechist files.
        
        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = self._file_filter(self.files, 'spechist', 'fit', dets=dets)
        return self.get(download_dir, files, **kwargs)

    def get_tte(self, download_dir, dets=None, full_day=False, **kwargs):
        """Download all TTE files associated with a time.
        
        Note:
            Unless you have a high-bandwidth connection and can handle
            downloading several GBs, it is not recommended to download the 
            full day of TTE data.

        Args:
            download_dir (str): The download directory
            dets (list, optional): The detectors' data to download. 
                                   If omitted, will download all.
            full_day (bool, optional): 
                If True, will download the TTE files for the full day.  If False,
                will return the TTE files for the covering the specified time.
                Default is False.
            verbose (bool, optional): If True, will output the download status. 
                                      Default is True.

        Returns:
            (list of Path)
        """
        files = []
        files.extend(self._file_filter(self.files, 'tte', 'fit.gz', dets=dets))
        files.extend(self._file_filter(self.files, 'tte', 'fit', dets=dets))

        if not full_day:
            files = self._filter_tte(files)

        return self.get(download_dir, files, **kwargs)

    def ls_cspec(self):
        """List all CSPEC files.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'cspec', 'pha')

    def ls_ctime(self):
        """List all CTIME files.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'ctime', 'pha')

    def ls_poshist(self):
        """List the poshist file
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'poshist', 'fit')

    def ls_spechist(self):
        """List all spechist files.
        
        Returns:
            (list of str)
        """
        return self._file_filter(self.files, 'spechist', 'fit')

    def ls_tte(self, full_day=False):
        """List all TTE files.

        Args:
            full_day (bool, optional):
                If True, will return the TTE files for the full day.  If False,
                will return the TTE files for the hour covering the specified 
                time. Default is False.
        
        Returns:
            (list of str)
        """
        files = []
        files.extend(self._file_filter(self.files, 'tte', 'fit.gz'))
        files.extend(self._file_filter(self.files, 'tte', 'fit'))

        if not full_day:
            files = self._filter_tte(files)

        return files

    def _construct_path(self, met_obj):
        """Constructs the path for a time
        
        Args:
            met_obj (:class:`.time.Met`): The MET time object

        Returns:
            str: The path of the directory for the time
        """
        path = os.path.join(self._root, met_obj.datetime.strftime('%Y/%m/%d'),
                            'current')
        return path

    def _filter_tte(self, files):
        """Filters a list of TTE files for only the files that contain the
        desired time
        
        Args:
            files (list of str): The list of TTE files

        Returns:
            list of str: The filtered list of files
        """
        id = self._args[0].strftime('%Y%m%d_%Hz')[2:]
        files = [f for f in files if id in f]
        return files


class ContinuousFtp(ContinuousFinder):
    """Class providing backwards compatibility for code written prior to v2.0.5
    where the ContinuousFtp handled interactions with continuous data.

    Parameters:
        args: The set of parameters needed to define the data path
        **kwargs: Options passed to :class:`ContinuousFinder` class
    """
    def __init__(self, *args, **kwargs):
        """Constructor"""
        super().__init__(*args, protocol='FTP', **kwargs)
