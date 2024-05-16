# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Very closely based on the program DoL (Daughter of Locburst).
# Written by:
#               Valerie Connaughton
#               University of Alabama in Huntsville (UAH)
#
# Included in the Gamma-ray Data Tools with permission from UAH
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
from pathlib import Path
# import required packages here to ensure we fail with
# "ImportError: No module named..." when one is missing
import numpy as np

__all__ = ["legacy_spectral_models", "legacy_functions", "legacy_dol", "__data_dir__"]
__author__ = ("Fermi GBM Team")
__version__ = "4.15a0"
__maintainer__ = "Josh Wood"
__email__ = "joshua.r.wood@nasa.gov"
__status__ = "Development"
__data_dir__ = os.path.join(os.path.dirname(__file__), "data/")

