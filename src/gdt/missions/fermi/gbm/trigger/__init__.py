# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Developed by: Jacob Smith
#               University of Alabama in Huntsville
#               Center for Space Plasma and Aeronomic Research
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
import matplotlib.pyplot as plt
from astropy.utils import isiterable

from astropy import units as u

from gdt.core.trigger import Trigger, TriggerAlgorithm, SlidingWindowMethod
from gdt.missions.fermi.gbm.detectors import GbmDetectors

from .algorithms import current_onboard_algorithms

__all__ = ['GbmOnboardTrigger']


class GbmOnboardTrigger(SlidingWindowMethod):
    """Class for applying GBM's on-board trigger algorithm"""

    def __init__(self, algorithms: dict = current_onboard_algorithms, background_window: int = 17024,
                 background_offset: int = 4096, resolution: int = 16,
                 channel_edges: list = None,
                 det_names: list = None,
                 n: int = 2,
                 verbose: bool = True):

        if channel_edges is None:
            channel_edges = [0, 8, 20, 33, 51, 85, 106, 127, 128]

        if det_names is None:
            det_names = [det.name for det in GbmDetectors if det.is_nai()]

        # Note: current algorithms only use nai detectors. 5 additional BGO
        # trigger algorithms are not given in the first GBM burst catalog:
        #       William S. Paciesas et al 2012 ApJS 199 18, DOI 10.1088/0067-0049/199/1/18
        super().__init__(algorithms, background_window, background_offset,
                         resolution, channel_edges, det_names, n, verbose)

    def prepare_data(self, ttes: list, time_range: list = None):
        """Method to prepare TTE data into a PHAII format for use with
        trigger algorithms. During preparation, TTE files are sorted
        according to detector name and we check that there is one 
        file for each NaI detector.

        Args:
            ttes (list): List of GbmTte objects

        Returns:
            (list)
        """
        found_names = []

        def sort_func(tte):
            if tte.detector not in self.det_names:
                raise ValueError("Detector {} not required by trigger".format(tte.detector))
            if tte.detector in found_names:
                raise ValueError("Duplicate TTE file for detector {}".format(tte.detector))
            found_names.append(tte.detector)
            return self.det_names.index(tte.detector)

        ttes.sort(key=sort_func)

        if len(ttes) != len(self.det_names):
            missing = [name for name in self.det_names if name not in found_names]
            raise ValueError("Missing TTE files for detectors {}".format(missing))

        super().prepare_data(ttes, time_range)
