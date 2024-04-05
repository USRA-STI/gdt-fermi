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
import numpy as np

from . import __data_dir__

###############################################################################

###########
# GLOBALS #
###########

# standard band spectra
band_hard = "band,alpha=0.0,beta=-1.5,epeak=1000.0,amp=10.0"
band_norm = "band,alpha=-1.0,beta=-2.3,epeak=230.0,amp=10.0"
band_soft = "band,alpha=-2.0,beta=-3.4,epeak=70.0,amp=10.0"

# rates files for each band spectrum
band_hard_50_300 = os.path.join(__data_dir__, "band_1deg_50_300_hard.npy")
band_norm_50_300 = os.path.join(__data_dir__, "band_1deg_50_300_norm.npy")
band_soft_50_300 = os.path.join(__data_dir__, "band_1deg_50_300_soft.npy")
band_soft_5_50 = os.path.join(__data_dir__, "band_1deg_5_50_soft.npy")

# standard comp spectra
comp_hard = "comp,index=-0.25,epeak=1000.0,amp=10.0"
comp_norm = "comp,index=-1.15,epeak=350.0,amp=10.0"
comp_soft = "comp,index=-1.95,epeak=50.0,amp=10.0"

# rates files for each comp spectrum
comp_hard_50_300 = os.path.join(__data_dir__, "comp_1deg_50_300_hard.npy")
comp_norm_50_300 = os.path.join(__data_dir__, "comp_1deg_50_300_norm.npy")
comp_soft_50_300 = os.path.join(__data_dir__, "comp_1deg_50_300_soft.npy")

###########
# GLOBALS #
###########

###############################################################################

def legacy_band(e, alpha, beta, epeak):
    r""" Band spectrum definition with legacy float32 support

    Parameters
    ----------
    e : np.ndarray
        Energy points to evaluate Band spectrum at
    alpha : float
        Low energy index
    beta : float
        High energy index
    epeak : float
        Energy of spectral peak in keV

    Returns
    -------
    spec : np.ndarray
        Band spectrum for (alpha, beta, epeak) evaluated at e
    """
    if alpha < -1.9:
        alpha = np.float32(-1.9)
    e0 = epeak / (np.float32(2.) + alpha)
    eb = (alpha - beta) * e0
    e100 = np.float32(100.)

    func1 = ((e / e100) ** alpha) * np.exp(-e / e0)
    func2 = ((alpha - beta) * e0 / e100) ** (alpha - beta) \
            * np.exp(beta - alpha) * (e / e100) ** beta

    return np.where(e < eb, func1, func2)


# END legacy_band()

def legacy_comp(e, index, epeak):
    r""" Comptonized spectrum definition with float32 support
 
    Parameters
    ----------
    e : np.ndarray
        Energy points to evaluate comptonized spectrum at
    index : float
        Spectral index for power law component
    epeak : float
        Peak energy in the spectrum dictating exponential cutoff


    Returns
    -------
    spec : np.ndarray
        Comptonized spectrum for (index, epeak) evaluated at e
    """
    return (e / np.float32(100.0)) ** index \
        * np.exp(-e * (np.float32(2.0) + index) / epeak)


# END legacy_comp()

def legacy_pl(e, index):
    r""" Power law spectrum definition with legacy float32 support

    Parameters
    ----------
    e : np.ndarray
        Energy points to evaluate power law spectrum at
    index : float
        Spectral index of the power law

    Returns
    -------
    spec : np.ndarray
        Power law spectrum for (index) evaluated at E)
    """
    return (e / np.float32(100.0)) ** index

# END legacy_pl()
