.. _gbm-saa:
.. |GbmSaa| replace:: :class:`~gdt.missions.fermi.gbm.saa.GbmSaa`
.. |FermiEarthPlot| replace:: :class:`~gdt.missions.fermi.plot.FermiEarthPlot`

*********************************************************************
Fermi GBM SAA Boundary Definition (:mod:`gdt.missions.fermi.gbm.saa`)
*********************************************************************
The South Atlantic Anomaly (SAA) region is an area of high particle 
particle flux in the Earth's Van Allen belts. Whenever Fermi is near
this region, the detectors are turned off and not collecting data.

The boundary of the SAA region is defined using a polygon in latitude
and longitude. The definition of this boundary was changed once during
GBM operations, resulting in two distinct periods:

* **Period 1 (Launch - July 23, 2024):** Original SAA boundary
  used during the first 16 years of operations. This definition
  is conservative and includes a large portion of the southern
  Atlantic ocean.
* **Period 2 (July 23, 2024 - present):** Updated SAA boundary
  implemented during the fourth gravitational wave observing
  run (O4). This definition is designed to increase detector uptime
  by reducing the area of the SAA boundary to more closely
  follow the limits of the Van Allen belts.

We can retrieve the latitude and longitude values for this
boundary at any time by initializing the |GbmSaa| class
with a time object based on the Astropy ``Time`` class.
For example, we can obtain the SAA boundary for **Period 1**
by using zero seconds in the :ref:`Fermi Mission Elapsed Time<fermi-time>` format:

    >>> from gdt.missions.fermi.time import Time
    >>> from gdt.missions.fermi.gbm.saa import GbmSaa
    >>> saa = GbmSaa(Time(0, format='fermi'))
    >>> saa.latitude
    array([-30.   , -19.867,  -9.733,   0.4  ,   2.   ,   2.   ,  -1.   ,
            -6.155,  -8.88 , -14.22 , -18.404, -30.   , -30.   ])
    >>> saa.longitude
    array([ 33.9  ,  12.398,  -9.103, -30.605, -38.4  , -45.   , -65.   ,
           -84.   , -89.2  , -94.3  , -94.3  , -86.1  ,  33.9  ])

We can also plot the SAA boundary using |FermiEarthPlot|:

    >>> import matplotlib.pyplot as plt
    >>> from gdt.missions.fermi.plot import FermiEarthPlot
    >>> plot = FermiEarthPlot(saa=saa, mcilwain=False)
    >>> plt.show()
    
.. image:: saa_figs/saafig1.png

The SAA region in this plot is marked in red.  For more details on customizing
these plots, including the SAA region, see 
:external:ref:`Plotting Spacecraft in Earth Orbit<plot-earthplot>`.

To get the region at a different time, you can either create a new |GbmSaa| object 
or simply update the existing object using the :meth:`~gdt.missions.fermi.gbm.saa.GbmSaa.update`
method:

    >>> import datetime
    >>> saa.update(Time(datetime.datetime.now(), format='datetime'))

The overlay below shows a comparison of the SAA regions from **Period 1** and **Period 2**.
Note that **Period 2** is essentially just a trimmed down version of **Period 1**.

.. image:: saa_figs/saafig2.png

    
Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.saa
   :inherited-members:


