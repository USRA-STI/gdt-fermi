.. _gbm-collection:
.. |GbmDetectorCollection| replace:: :class:`~gdt.missions.fermi.gbm.collection.GbmDetectorCollection`
.. |DataCollection| replace:: :class:`~gdt.core.collection.DataCollection`

**********************************************************************
Fermi GBM Data Collections (:mod:`gdt.missions.fermi.gbm.collections`)
**********************************************************************
The |GbmDetectorCollection| class extends the functionality of |DataCollection|
so that it can be used with GBM's heterogenous detector types: NaI and BGO 
detectors.  This specialized collection class allows us to apply the same 
action over a collection of detectors, but also to specify different parameters
based on if a detector is an NaI or BGO.

As an example, let's assume we have two NaI CSPEC files open and one BGO file
open, and we read them into a collection:

    >>> from gdt.missions.fermi.gbm.collection import GbmDetectorCollection
    >>> cspecs = GbmDetectorCollection.from_list([cspec_n0, cspec_n1, cspec_b0])
    >>> cspecs
    <GbmDetectorCollection: 3 Cspec objects>
    
Just like the base |DataCollection|, we can retrieve properties of each item in
the collection with a single call:

    >>> # get the energy range for each item
    >>> cspec.energy_range()
    [(4.5702357, 2000.0), (4.089358, 2000.0), (113.00731, 50000.0)]

Or call a function that applies the same arguments/keywords to each item:

    >>> # generate the lightcurve over from -10 to +10 s
    >>> cspecs.to_lightcurve(time_range=[-10.0, 10.0])
    [<TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>,
     <TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>,
     <TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>]

But what if we wanted to integrate over a different energy range for the NaI 
and BGO detectors to generate the lightcurves?  For a situation like that,
we specify NaI-specific arguments with ``nai_args`` and BGO-specific arguments
with ``bgo_args``.  Similarly, if we need to pass different keywords, we would
specify ``nai_kwargs`` and ``bgo_kwargs``.  Note that the xxx_args and xxx_kwargs
follow the same convention as normal Python args and kwargs: args are passed
as tuples containing the arguments, and kwargs are passed as dicts.

For example:

    >>> nai_kwargs = {'energy_range': (8.0, 900.0)}
    >>> bgo_kwargs = {'energy_range': (325, 35000.0)}
    >>> cspecs.to_lightcurve(time_range=[-10.0, 10.0], nai_kwargs=nai_kwargs, bgo_kwargs=bgo_kwargs)
    [<TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>,
     <TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>,
     <TimeBins: 14 bins;
      range (-13.824249982833862, 10.752128005027771);
      1 contiguous segments>]
      
Or similarly slice over different energies:
    >>> cspecs.slice_energy(nai_args=((8.0, 900.0),), bgo_args=((325, 35000.0),))
    [<Cspec: 
      trigger time: 484477130.219298;
      time range (-4003.394767999649, 4001.3853880167007);
      energy range (7.3262835, 923.9744)>,
     <Cspec: 
      trigger time: 484477130.219298;
      time range (-4003.394767999649, 4001.3853880167007);
      energy range (7.7042193, 912.946)>,
     <Cspec: 
      trigger time: 484477130.219298;
      time range (-4003.394767999649, 4001.3853880167007);
      energy range (317.86954, 35982.97)>]

For other details about using collections, please see the documentation on 
:external:ref:`Data Collections<core-collection>`.

Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.collection
   :inherited-members:


