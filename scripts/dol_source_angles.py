#!/usr/bin/env python
import os
import argparse

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

from gdt.missions.fermi.gbm.localization.dol.legacy_functions import *
from gdt.missions.fermi.gbm.localization.dol.legacy_dol import legacy_DoL


def main():
    nen = legacy_DoL.default_nen

    ########
    # ARGS #
    ########

    p = argparse.ArgumentParser(
        description="Compute spacecraft and earth center angles to a source.",
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--ra", type=np.float32,
                   help="J2000 right ascension of source in degrees")
    p.add_argument("--dec", type=np.float32,
                   help="J2000 declination of source in degrees")
    p.add_argument("--sc_pos", type=np.float32, nargs=3, default=None,
                   help="Spacecraft position (X,Y,Z)")
    p.add_argument("--sc_quat", type=np.float64, nargs=4, default=None,
                   help="Spacecraft quaternion (S,X,Y,Z)")
    p.add_argument("--path", type=str, default="",
                   help="Path to dol_rates_xxx.txt file")
    args = p.parse_args()

    valid_sc = (args.sc_pos is not None) and (args.sc_quat is not None)
    valid_path = os.path.exists(args.path)

    if valid_path and valid_sc:
        raise ValueError("User must supply path to dol_rates_xxx.txt file OR "
                         "sc_pos and sc_quat, not both.")
    if (not valid_path) and (not valid_sc):
        raise ValueError("User must supply path to dol_rates_xxx.txt file OR "
                         "sc_pos and sc_quat")

    if valid_path:
        f = open(args.path)
        info = np.array(f.readline().split())
        args.sc_pos = info[range(-nen - 13, -nen - 10)].astype(np.float32)
        args.sc_quat = info[range(-nen - 10, -nen - 6)].astype(np.float64)
        f.close()
    else:
        args.sc_pos = np.array(args.sc_pos, np.float32)
        args.sc_quat = np.array(args.sc_quat, np.float64)

    # convert to radians
    args.ra /= legacy_dtorad
    args.dec /= legacy_dtorad

    ########
    # ARGS #
    ########

    ###########################################################################

    #################
    # SOURCE ANGLES #
    #################

    # angles of NaI detectors given in spacecraft azimuth, zenith
    nai_az, nai_zen, nai_unit_vec, back_unit_vec = \
        initialize_det_geometry(verbose=False)

    # direction of the earth relative to spacecraft
    geodir, geo_az, geo_zen, scx, scy, scz = get_geocenter(
        args.sc_quat, args.sc_pos, verbose=False)

    # source position vector as well as direction in spacecraft azimuth, zenith
    source_pos = ang_to_cart_dec(args.ra, args.dec)
    source_az, source_zen = j2000_to_sc(scx, scy, scz, source_pos)

    good_ang = get_good_angle(-args.sc_pos, source_pos)
    print(" J2000: dist from source to geocenter (initial) {0:.9}".format(
        good_ang))
    print(" J2000: dist from source to z-axis (initial) {0:.17}".format(
        legacy_dtorad * np.arccos(
            dot(scz, source_pos / np.sqrt(dot(source_pos, source_pos))))))

    # angles of detectors to source given source_az and source_zen in sc frame
    print(" Source angles (initial pos):")
    det_ang_source = get_det_geometry(source_az, source_zen, nai_az, nai_zen)
    for j in range(12):
        print("{0:9} {1:0<18.17} {2:0<10.9} {3:0<10.9}".format(
            j, legacy_dtorad * det_ang_source[j], nai_az[j], nai_zen[j]))

    #################
    # SOURCE ANGLES #
    #################

    ###########################################################################


# END main()

if __name__ == "__main__":
    main()
