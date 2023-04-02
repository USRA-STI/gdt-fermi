#  CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
#  Contract No.: CA 80MSFC17M0022
#  Contractor Name: Universities Space Research Association
#  Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
#  Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
#  Developed by: William Cleveland and Adam Goldstein
#                Universities Space Research Association
#                Science and Technology Institute
#                https://sti.usra.edu
#
#  Developed by: Daniel Kocevski
#                National Aeronautics and Space Administration (NASA)
#                Marshall Space Flight Center
#                Astrophysics Branch (ST-12)
#
#   Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
#   in compliance with the License. You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software distributed under the License
#  is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
#  implied. See the License for the specific language governing permissions and limitations under the
#  License.
#
import os
import numpy as np

from gdt.core.types import Numbers

__all__ = ['calc_mcilwain_l']

coeffs_file = os.path.join(os.path.dirname(__file__), 'data', 'McIlwainL_Coeffs.npy')
poly_coeffs = None


def calc_mcilwain_l(latitude: Numbers, longitude: Numbers) -> np.ndarray:
    """Estimate the McIlwain L value given the latitude (-30, +30) and
    East Longitude.  This uses a cubic polynomial approximation to the full
    calculation.

    Args:
        latitude: Latitude in degrees from -30 to +30
        longitude: East longitude in degrees from 0 to 360

    Returns:
        (float or np.array)
    """
    global coeffs_file, poly_coeffs

    latitude = np.asarray([latitude])
    longitude = np.asarray([longitude])
    orig_shape = latitude.shape
    latitude = latitude.flatten()
    longitude = longitude.flatten()

    # Load the poly coeffs data if it isn't already in memory.
    if poly_coeffs is None:
        poly_coeffs = np.load(coeffs_file)

    longitude[longitude < 0.0] += 360.0
    longitude[longitude == 360.0] = 0.0

    bad_idx = (latitude < -30.0) | (latitude > 30.0) | (longitude < 0.0) | (
            longitude >= 360.0)
    if np.sum(bad_idx) != 0:
        raise ValueError(
            'Out of range coordinates for McIlwain L for {0} locations'.format(
                np.sum(bad_idx)))

    idx = np.asarray((longitude / 10.0).astype(int))
    idx2 = np.asarray(idx + 1)
    idx2[idx2 >= 36] = 0
    idx2 = idx2.astype(int)

    longitude_left = 10.0 * idx
    f = (longitude - longitude_left) / 10.0  # interpolation weight, 0 to 1

    mc_l = (1.0 - f) * (
            poly_coeffs[idx, 0] + poly_coeffs[idx, 1] * latitude + poly_coeffs[idx, 2] * latitude ** 2 +
            poly_coeffs[idx, 3] * latitude ** 3) + f * (
                   poly_coeffs[idx2, 0] + poly_coeffs[idx2, 1] * latitude + poly_coeffs[idx2, 2] * latitude ** 2 +
                   poly_coeffs[idx2, 3] * latitude ** 3)
    mc_l = np.squeeze(mc_l.reshape(orig_shape))
    return mc_l[()] if mc_l.shape == () else mc_l
