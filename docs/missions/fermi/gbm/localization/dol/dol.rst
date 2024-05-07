.. _dol:

**************************************************************************************
Fermi GBM Localizations using the DoL (:mod:`gdt.missions.fermi.gbm.localization.dol`)
**************************************************************************************

The ``dol`` package provides modules that perform the same legacy DoL localizations used
to create GBM ground localizations by the GBM Science Team. These python modules exactly
replicate the operations of the original Fortran code.

The package is divided into the following modules:

.. toctree::
   :maxdepth: 1

   legacy_dol
   legacy_functions
   legacy_spectral_models

Users can check the links above for additional information on the subroutines of each module.
The remainder of this page will focus on creating an example ground localization
for GRB 170817A.

Preparing Data to Perform a Localization
========================================
The first step to creating a ground localization is obtaining data near the time of a burst.
For simplicity, we will work with a triggered data (trigdat) file in this example because
it conventiently contains binned data for all detectors as well as position history information
about the spacecraft. However, any :ref:`Fermi GBM PHAII formatted data<gbm-phaii>`
and :ref:`position history file<gbm-poshist>` will work as long as the data are binned to
allow selection of counts within the 50-300 keV energy range used for GRB localizations.

We begin by downloading the trigdat file ``glg_trigdat_all_bn170817529_v01.fit`` to our local 
directory using the :class:`~gdt.missions.fermi.gbm.finders.TriggerFtp` class initialized with
GBM burst number 170817529, which corresponds to GRB 170817A.

    >>> from gdt.missions.fermi.gbm.finders import TriggerFtp
    >>> ftp = TriggerFtp("170817529")
    >>> ftp.get_trigdat(".")

We then plot a summed lightcurve for all detectors identified as contributing
to the identification of this trigger. This is done over the 50-300 keV range
near the time of the trigger. We specify a binning timescale of
64 milliseconds when creating this lightcurve because we are examining a
known short GRB. Longer binnings are more appropriate for long GRBs.

    >>> from gdt.core.plot.lightcurve import Lightcurve
    >>> loc_erange = (50.0, 300.0)
    >>> trigdet = trigdat.triggered_detectors
    >>> summed_phaii = trigdat.sum_detectors(trigdet, timescale=64)
    >>> summed_phaii = summed_phaii.slice_energy(loc_erange)
    >>> summed_phaii = summed_phaii.slice_time((-8, 14))
    >>> lcplot1 = Lightcurve(summed_phaii.to_lightcurve())

Upon visual inspection, we see that the GRB appears in the lightcurve as a single pulse
from approximately -0.256 sec to +0.448 sec around the trigger time. This will be our 
source window that we use to compute the observed counts in each detector while the
GRB is active. Let us do that now while also creating a new lightcurve with the source
period highlighted.

    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>>
    >>> src_time = (-0.256, 0.448)
    >>> src_counts = []
    >>> src_exposure = []
    >>> print("\nSource Counts")
    >>> for det in trigdat._detectors:
    ...    # get data for one detector
    ...    phaii = trigdat.to_phaii(det, timescale=64)
    ...    # select data for the source period
    ...    bin = phaii.data.integrate_time(*src_time)
    ...    # get counts and exposure
    ...    src_counts.append(bin.counts.astype(np.int32))
    ...    src_exposure.append(bin.exposure)
    ...    print(f" - {det} {src_counts[-1]}")
    Source Counts
     - n0 [ 48 267 176 127 135  25  18  77]
     - n1 [ 50 301 188 153 155  27  30  66]
     - n2 [ 55 285 190 162 141  31  54  33]
     - n3 [ 64 316 171 131 126  32  22  55]
     - n4 [ 51 293 188 147 113  26  52  27]
     - n5 [ 65 367 217 177 149  30  46  23]
     - n6 [ 54 217 135 114 103  23  53  18]
     - n7 [ 70 275 182 135 107  22  26  58]
     - n8 [ 59 252 155 121 116  21  43  32]
     - n9 [ 48 179 152 115 116  30  82  13]
     - na [ 34  89 110 146 104  25  59  27]
     - nb [ 33 120 101 155 122  44  29  74]
     - b0 [374 204 318 128  26  10  15  95]
     - b1 [358 219 258 103  37  26  15  85]
    >>> avg_src_exposure = np.sum(src_exposure) / np.array(src_exposure).size
    >>> print(" Exposure %.3f sec" % avg_src_exposure)
    Exposure 0.768 sec
    >>> lcplot2 = Lightcurve(summed_phaii.to_lightcurve())
    >>> ax = plt.gca()
    >>> src_span = ax.axvspan(*src_time, color='C1', alpha=0.2)
    >>> l2 = plt.legend([src_span], ["Source Period"], framealpha=1.0,)

Creating Your Own Ground Localization
=====================================
