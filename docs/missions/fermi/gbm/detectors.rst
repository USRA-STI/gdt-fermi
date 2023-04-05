.. _gbm-detectors:
.. |GbmDetectors| replace:: :class:`~gdt.missions.fermi.gbm.detectors.GbmDetectors`
.. |Detectors| replace:: :class:`~gdt.core.detector.Detectors`

******************************
Fermi GBM Detector Definitions
******************************
The |GbmDetectors| class contains the naming and orientation definitions of the
GBM detectors.

The GBM detectors have three different naming/indexing conventions, although the
one that is used most is the ``'n0', 'n1',...,'b0', 'b1'`` naming convention.

We can easily retrieve a detector definition by using standard "dot" notation:

    >>> from gdt.missions.fermi.gbm.detectors import GbmDetectors
    >>> GbmDetectors.n0
    <GbmDetectors: n0>

We can retrive the full name of the detector, which is what is mostly used in
the FITS headers of the GBM data files:

    >>> GbmDetectors.nb.full_name
    'NAI_11'

There is also a standard detector indexing scheme that is used for all GBM 
detectors:

    >>> GbmDetectors.b0.number
    12

Since the |GbmDetectors| class inherits from the |Detectors| base class, we 
can also retrieve the pointing information of a GBM detector:

    >>> # detector azimuth, zenith
    >>> GbmDetectors.from_str('n2').pointing()
    (<Quantity 58.44 deg>, <Quantity 90.21 deg>)
    
    >>> # detector elevation
    GbmDetectors.from_full_name('NAI_02').elevation
    <Quantity -0.21 deg>

We can also iterate over all GBM detectors:
    
    >>> # the list of detector names
    >>> print([det.name for det in GbmDetectors])
    ['n0', 'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7', 'n8', 'n9', 'na', 'nb', 
    'b0', 'b1']

We can also get the list of BGO (or NaI) detectors:

    >>> GbmDetectors.bgo()
    [<GbmDetectors: b0>, <GbmDetectors: b1>]

And we can test if a particular detector is an NaI or BGO detector:

    >>> print([det.is_nai() for det in GbmDetectors])
    [True, True, True, True, True, True, True, True, True, True, True, True, 
    False, False]


Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.detectors
   :inherited-members:

