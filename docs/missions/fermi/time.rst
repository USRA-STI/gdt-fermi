.. _fermi-time:

*******************
Fermi Mission Epoch
*******************

The Fermi Mission epoch, also called the Fermi Mission Elapsed Time (MET) is 
the number of seconds elapsed since January 1, 2001, 00:00:00 UTC, including 
leap seconds.  We have defined a specialized epoch to work with Astropy ``Time``
objects so that Fermi MET can be easily converted to/from other formats and time
scales.

To use this, we simply import and create an astropy Time object with a `'fermi'`
format:

    >>> from gdt.missions.fermi.time import *
    >>> fermi_met = Time(697422649, format='fermi')
    >>> fermi_met
    <Time object: scale='tt' format='fermi' value=697422649.0>
    
Now, say we want to retrieve the GPS timestamp:

    >>> fermi_met.gps
    1359765062.0

The Astropy ``Time`` object readily converts it for us. We can also do the 
reverse conversion:

    >>> gps_time = Time(fermi_met.gps, format='gps')
    >>> gps_time
    <Time object: scale='tai' format='gps' value=1359765062.0>
    
    >>> gps_time.fermi
    697422649.0

And we should, of course, get back the Fermi MET we started with.  This enables
you do do any time conversions already provided by Astropy, as well as time
conversions between other missions within the GDT.

In addition to time conversions, all time formatting available in Astropy is 
also available here.  For example, we can format the Fermi MET in ISO format:

    >>> fermi_met.iso
    '2023-02-07 00:31:53.184'
    
Finally, there is a specialized format associated with the Fermi epoch, which
allows us to output the format for Fermi GBM burst numbers, which is based on
time:

    >>> fermi_met.gbm_bn
    '230207021'

    
Reference/API
=============

.. automodapi:: gdt.missions.fermi.time
   :inherited-members:


