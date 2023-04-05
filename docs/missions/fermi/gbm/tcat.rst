.. _gbm-tcat:
.. |Tcat| replace:: :class:`~gdt.missions.fermi.gbm.tcat.Tcat`

******************************
Fermi GBM Trigger Catalog Data
******************************
The Fermi Trigger CATalog (TCAT) file is a FITS file containing high-level
information about an on-board GBM trigger and is ingested at HEASARC to 
include basic information in the trigger catalog.  The TCAT is purely metadata,
having only a single FITS header, but the |Tcat| class contains some convenience
functions for accessing and converting the metadata into GDT objects.

We can open a TCAT file:

    >>> from gdt import test_data
    >>> from gdt.missions.fermi.gbm.tcat import Tcat
    >>> filepath = test_data['fermi-gbm'].joinpath('glg_tcat_all_bn190222537_v01.fit')
    >>> tcat = Tcat.open(filepath)
    >>> tcat
    <Tcat: GRB190222537>
    
We can retrieve some basic properties:

    >>> tcat.name
    'GRB190222537'
    >>> tcat.trigtime
    572532812.150778
    >>> tcat.time_range
    <TimeRange: (572532678.644472, 572533293.055776)>

We can also retrieve the instrument that localized this trigger and what that
localization is:

    >>> tcat.localizing_instrument
    'Fermi, GBM'
    >>> tcat.location
    <SkyCoord (ICRS): (ra, dec) in deg
        (147.32, 60.94)>

Note that while this is obviously localized by GBM, there are cases where a
transient is localized much more precisely by another instrument, and that will
be recorded in the TCAT.

We can also retrieve the position of Fermi in orbit at the time of the detection:

    >>> tcat.fermi_location
    (<Longitude 74.25 deg>, <Latitude -24.85 deg>)
    
Of course the full set of metadata is available in the header:

    >>> tcat.headers
    <TcatHeaders: 1 headers>
    
Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.tcat
   :inherited-members:


