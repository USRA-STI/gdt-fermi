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
import astropy.units as u
from gdt.core.detector import Detectors

__all__ = ['GbmDetectors']

# unfortunately Sphinx has a major bug that prevents the autodoc of Enums, 
# so we have to define all of this in the docstring...

class GbmDetectors(Detectors):
    """The GBM Detector name and orientation definitions.    

    .. rubric:: Attributes Summary
    .. autosummary::

      azimuth
      elevation
      full_name
      number
      zenith

    .. rubric:: Methods Summary

    .. autosummary::

      bgo
      from_full_name
      from_num
      from_str
      is_bgo
      is_nai
      nai
      pointing
      skycoord
  
    .. rubric:: Attributes Documentation

    .. autoattribute:: azimuth
    .. autoattribute:: elevation
    .. autoattribute:: full_name
    .. autoattribute:: number
    .. autoattribute:: zenith

    .. rubric:: Methods Documentation

    .. autoattribute:: bgo
    .. automethod:: from_full_name
    .. automethod:: from_num
    .. automethod:: from_str
    .. automethod:: is_bgo
    .. automethod:: is_nai
    .. autoattribute:: nai
    .. automethod:: pointing
    .. automethod:: skycoord
    """
    n0 = ('NAI_00', 0, 45.89 * u.deg, 20.58 * u.deg)
    n1 = ('NAI_01', 1, 45.11 * u.deg, 45.31 * u.deg)
    n2 = ('NAI_02', 2, 58.44 * u.deg, 90.21 * u.deg)
    n3 = ('NAI_03', 3, 314.87 * u.deg, 45.24 * u.deg)
    n4 = ('NAI_04', 4, 303.15 * u.deg, 90.27 * u.deg)
    n5 = ('NAI_05', 5, 3.35 * u.deg, 89.79 * u.deg)
    n6 = ('NAI_06', 6, 224.93 * u.deg, 20.43 * u.deg)
    n7 = ('NAI_07', 7, 224.62 * u.deg, 46.18 * u.deg)
    n8 = ('NAI_08', 8, 236.61 * u.deg, 89.97 * u.deg)
    n9 = ('NAI_09', 9, 135.19 * u.deg, 45.55 * u.deg)
    na = ('NAI_10', 10, 123.73 * u.deg, 90.42 * u.deg)
    nb = ('NAI_11', 11, 183.74 * u.deg, 90.32 * u.deg)
    b0 = ('BGO_00', 12, 0.00 * u.deg, 90.00 * u.deg)
    b1 = ('BGO_01', 13, 180.00 * u.deg, 90.00 * u.deg)

    @classmethod
    def bgo(cls):
        """Get all detectors that are BGOs
        
        Returns:
            (list of :class:`GbmDetectors`)    
        """
        return [x for x in cls if x.is_bgo()]

    @classmethod
    def nai(cls):
        """Get all detectors that are NaIs
    
        Returns:
            (list of :class:`GbmDetectors`)
        """
        return [x for x in cls if x.is_nai()]

    def is_bgo(self):
        """Check if detector is a BGO.
    
        Returns:
            (bool)
        """
        return self.name[0] == 'b'

    def is_nai(self):
        """Check if detector is an NaI.
    
        Returns:
            (bool)
        """
        return self.name[0] == 'n'

