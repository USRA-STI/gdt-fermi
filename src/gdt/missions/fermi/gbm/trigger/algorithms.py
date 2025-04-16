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
# This software is not subject to EAR.
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
from gdt.core.trigger import TriggerAlgorithm

__all__ = ['onboard_algorithms_2024_03_14', 'current_onboard_algorithms']

# The dictionary objects below this point define sets of Fermi-GBM trigger
# algorithms used onboard the Fermi spacecraft at different times during
# its science mission. Each dictionary is labeled according to the date
# it was implemented. Comments in-line with each trigger algorithm note
# the reference where it is defined.

# References:
#   [1] William S. Paciesas et al 2012 ApJS 199 18, DOI 10.1088/0067-0049/199/1/18
#   [2] https://gcn.nasa.gov/circulars/35929 2024
#   [3] Meegan et al 2009 ApJ 702 791, DOI 10.1088/0004-637X/702/1/791

onboard_algorithms_2024_03_14 = {
    1: TriggerAlgorithm(timescale=16, offset=0, channels=(3, 4), threshold=7.5),  # Ref [1]
    2: TriggerAlgorithm(timescale=32, offset=0, channels=(3, 4), threshold=7.5),  # Ref [1]
    3: TriggerAlgorithm(timescale=32, offset=16, channels=(3, 4), threshold=7.5),  # Ref [1]
    4: TriggerAlgorithm(timescale=64, offset=0, channels=(3, 4), threshold=5.0),  # Ref [1]
    5: TriggerAlgorithm(timescale=64, offset=32, channels=(3, 4), threshold=5.0),  # Ref [1]
    6: TriggerAlgorithm(timescale=128, offset=0, channels=(3, 4), threshold=5.0),  # Ref [1]
    7: TriggerAlgorithm(timescale=128, offset=64, channels=(3, 4), threshold=5.0),  # Ref [1]
    8: TriggerAlgorithm(timescale=256, offset=0, channels=(3, 4), threshold=4.0),  # Ref [2]
    9: TriggerAlgorithm(timescale=256, offset=128, channels=(3, 4), threshold=4.0),  # Ref [2]
    10: TriggerAlgorithm(timescale=512, offset=0, channels=(3, 4), threshold=4.0),  # Ref [2]
    11: TriggerAlgorithm(timescale=512, offset=256, channels=(3, 4), threshold=4.0),  # Ref [2]
    12: TriggerAlgorithm(timescale=1024, offset=0, channels=(3, 4), threshold=4.5),  # Ref [1]
    13: TriggerAlgorithm(timescale=1024, offset=512, channels=(3, 4), threshold=4.5),  # Ref [1]
    14: TriggerAlgorithm(timescale=2048, offset=0, channels=(3, 4), threshold=4.5),  # Ref [1]
    15: TriggerAlgorithm(timescale=2048, offset=1024, channels=(3, 4), threshold=4.5),  # Ref [1]
    16: TriggerAlgorithm(timescale=4096, offset=0, channels=(3, 4), threshold=4.5),  # Ref [1]
    17: TriggerAlgorithm(timescale=4096, offset=2048, channels=(3, 4), threshold=4.5),  # Ref [1]
    22: TriggerAlgorithm(timescale=16, offset=0, channels=(2, 2), threshold=8.0),  # Ref [1]
    23: TriggerAlgorithm(timescale=32, offset=0, channels=(2, 2), threshold=8.0),  # Ref [1]
    24: TriggerAlgorithm(timescale=32, offset=16, channels=(2, 2), threshold=8.0),  # Ref [1]
    25: TriggerAlgorithm(timescale=64, offset=0, channels=(2, 2), threshold=5.5),  # Ref [1]
    26: TriggerAlgorithm(timescale=64, offset=32, channels=(2, 2), threshold=5.5),  # Ref [1]
    43: TriggerAlgorithm(timescale=16, offset=0, channels=(5, 7), threshold=8.0),  # Ref [1]
    50: TriggerAlgorithm(timescale=16, offset=0, channels=(4, 7), threshold=8.0),  # Ref [1]
}

# a dictionary that always points to the most recent set of trigger algorithms
current_onboard_algorithms = onboard_algorithms_2024_03_14
