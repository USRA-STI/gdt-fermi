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
from astropy.coordinates import FunctionTransform, ICRS, frame_transform_graph
from gdt.core.coords import *
from gdt.core.coords.spacecraft.frame import spacecraft_to_icrs, icrs_to_spacecraft

__all__ = ['FermiFrame', 'fermi_to_icrs', 'icrs_to_fermi']

class FermiFrame(SpacecraftFrame):
    """
    The Fermi spacecraft frame in azimuth and elevation.  The frame is defined 
    as a quaternion that represents a rotation from the Fermi frame to the ICRS 
    frame. This class is a wholesale inheritance of SpacecraftFrame
    
    Example use:
    
        >>> from gdt.core.coords import Quaternion
        >>> quat = Quaternion([-0.218,  0.009,  0.652, -0.726], scalar_first=False)
        >>> fermi_frame = FermiFrame(quaternion=quat)
        >>> coord = SkyCoord(100.0, -30.0, unit='deg')
        >>> az_el = SkyCoord.transform_to(fermi_frame)
    """
    pass
    
@frame_transform_graph.transform(FunctionTransform, FermiFrame, ICRS)
def fermi_to_icrs(fermi_frame, icrs_frame):
    """Convert from the Fermi frame to the ICRS frame.
    
    Args:
        fermi_frame (:class:`FermiFrame`): The Fermi frame
        icrs_frame (:class:`astropy.coordinates.ICRS`)
    
    Returns:
        (:class:`astropy.coordinates.ICRS`)
    """
    return spacecraft_to_icrs(fermi_frame, icrs_frame)

@frame_transform_graph.transform(FunctionTransform, ICRS, FermiFrame)
def icrs_to_fermi(icrs_frame, fermi_frame):
    """Convert from the ICRS frame to the Fermi frame.
    
    Args:
        icrs_frame (:class:`astropy.coordinates.ICRS`)
        fermi_frame (:class:`FermiFrame`): The Fermi frame
    
    Returns:
        (:class:`FermiFrame`)
    """
    return icrs_to_spacecraft(icrs_frame, fermi_frame)

