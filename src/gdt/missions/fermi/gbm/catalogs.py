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

from gdt.core import cache_path
from gdt.core.heasarc import BrowseCatalog

__all__ = ['TriggerCatalog', 'BurstCatalog']

gbm_cache_path = os.path.join(cache_path, 'fermigbm')

class TriggerCatalog(BrowseCatalog):
    """Class that interfaces with the GBM Trigger Catalog via HEASARC Browse.

    Parameters:
        cache_path (str): The path where the cached catalog will live.
        cached (bool, optional): Set to True to read from the cached file
                                 instead of querying HEASARC. Default is False.
        verbose (bool, optional): Default is True
    """
    def __init__(self, cache_path=gbm_cache_path, **kwargs):
        super().__init__(cache_path, table='fermigtrig', **kwargs)

class BurstCatalog(BrowseCatalog):
    """Class that interfaces with the GBM Burst Catalog via HEASARC Browse.

    Parameters:
        cache_path (str): The path where the cached catalog will live.
        cached (bool, optional): Set to True to read from the cached file
                                 instead of querying HEASARC. Default is False.
        verbose (bool, optional): Default is True
    """
    def __init__(self, cache_path=gbm_cache_path, **kwargs):
        super().__init__(cache_path, table='fermigbrst', **kwargs)
