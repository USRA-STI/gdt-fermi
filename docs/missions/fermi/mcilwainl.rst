.. _fermi-mcilwainl:
.. |calc_mcilwain_l| replace:: :func:`~gdt.missions.fermi.mcilwainl.calc_mcilwain_l`

********************************************************************
McIlwain L Parameter for Fermi (:mod:`gdt.missions.fermi.mcilwainl`)
********************************************************************
The function |calc_mcilwain_l| is used to estimate the McIlwain L parameter, 
which describes a subset of the geomagnetic field lines and is correlated with
trapped charged particle density.

This is a specialized function for Fermi, assuming Fermi's approximate altitude,
and which implements a cubic polynomial approximation to the full calculation 
of the McIlwain L parameter as a function of latitude/longitude.  It is only 
valid within the East latitude range of :math:`(-30^\circ, +30^\circ)`, which 
covers the extent of Fermi's orbit.

As an example, we can calculate the McIlwain L parameter if Fermi were at
:math:`20^\circ N, 100^\circ W` (over Central America, just outside of Mexico 
City):

    >>> from gdt.missions.fermi.mcilwainl import *
    >>> # note that we must specify East Longitude
    >>> calc_mcilwain_l(-100.0, 20.0)
    1.44546112

Reference/API
=============

.. automodapi:: gdt.missions.fermi.mcilwainl
   :inherited-members:


