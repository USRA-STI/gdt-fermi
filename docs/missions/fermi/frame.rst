.. _fermi-frame:
.. |FermiFrame| replace:: :class:`~gdt.missions.fermi.frame.FermiFrame`
.. |GbmPosHist| replace:: :class:`~gdt.missions.fermi.frame.gbm.poshist.GbmPosHist`
.. |Quaternion| replace:: :class:`~gdt.core.coords.Quaternion`

**********************
Fermi Spacecraft Frame
**********************

The Fermi spacecraft frame, |FermiFrame|, is the frame that is aligned
with the Fermi spacecraft coordinate frame, and is represented by a 
quaternion that defines the rotation from spacecraft coordinates to the ICRS
coordinate frame.  This frame takes advantage of the Astropy coordinate frame
design, so we can use the FermiFrame to convert Astropy SkyCoord objects 
between the FermiFrame and any celestial frame.

While the FermiFrame is typically initialized when reading from a mission 
position history file (e.g. |GbmPosHist|) instead of manually by a user, we
can manually define the frame with a |Quaternion|:

    >>> from gdt.core.coords import Quaternion
    >>> from gdt.missions.fermi.frame import *
    >>> quat = Quaternion([-0.218,  0.009,  0.652, -0.726], scalar_first=False)
    >>> fermi_frame = FermiFrame(quaternion=quat)
    >>> fermi_frame
    <FermiFrame: 1 frames;
     obstime=[J2000.000]
     obsgeoloc=[(0., 0., 0.) m]
     obsgeovel=[(0., 0., 0.) m / s]
     quaternion=[(x, y, z, w) [-0.218,  0.009,  0.652, -0.726]]>

Notice that we can also define the frame with an ``obstime``, which is useful
for transforming between the FermiFrame and a non-inertial time-dependent frame; 
an ``obsgeoloc``, which can define the spacecraft location in orbit; and
``obsgeovel``, which defines the spacecraft orbital velocity.

Now let us define a SkyCoord in RA and Dec:

    >>> from astropy.coordinates import SkyCoord
    >>> coord = SkyCoord(100.0, -30.0, unit='deg')
    
And we can simply rotate this into the Fermi frame with the following:
    
    >>> fermi_coord = coord.transform_to(fermi_frame)
    >>> (fermi_coord.az, fermi_coord.el)
    (<Longitude [200.39733555] deg>, <Latitude [-41.88750942] deg>)

We can also transform from the Fermi frame to other frames.  For example, we
define a coordinate in the Fermi frame this way:

    >>> fermi_coord = SkyCoord(50.0, 25.0, frame=fermi_frame, unit='deg')
    
Now we can tranform to ICRS coordinates:

    >>> fermi_coord.icrs
    <SkyCoord (ICRS): (ra, dec) in deg
        [(313.69000519, 26.89158349)]>

or Galactic coordinates:

    >>> fermi_coord.galactic
    <SkyCoord (Galactic): (l, b) in deg
        [(71.5141302, -11.56931006)]>
        
or any other coordinate frames provided by Astropy.



Reference/API
=============

.. automodapi:: gdt.missions.fermi.frame
   :inherited-members:


