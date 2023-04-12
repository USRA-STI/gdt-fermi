.. _gbm-scat:
.. |Scat| replace:: :class:`~gdt.missions.fermi.gbm.scat.Scat`
.. |GbmModelFit| replace:: :class:`~gdt.missions.fermi.gbm.scat.GbmModelFit`
.. |GbmDetectorData| replace:: :class:`~gdt.missions.fermi.gbm.scat.GbmDetectorData`
.. |Parameter| replace:: :class:`~gdt.core.data_primitives.Parameter`

*****************************************************************
Fermi GBM Spectral Fits Data (:mod:`gdt.missions.fermi.gbm.scat`)
*****************************************************************
GBM provides standard spectral fits for each GRB triggered on-board. These 
spectral fits are hosted in FITS format at HEASARC as "SCat" (Spectral Catalog) 
files. The SCat files contain the spectral fit data, including fit parameters, 
uncertainties, fluxes/fluences, fit statistic, covariance. The files also 
contain metadata about the detectors and data that were used in the fit. All of 
this information can be easily accessed using the |Scat| class.

We can read a SCat FITS file:

    >>> from gdt.core import data_path
    >>> from gdt.missions.fermi.gbm.scat import Scat
    >>> filepath = data_path.joinpath('fermi-gbm/glg_scat_all_bn170817529_flnc_comp_v00.fit')
    >>> 
    >>> scat
    <Scat: glg_scat_all_bn170817529_flnc_comp_v00.fit;
     4 detectors; 1 fits>

We can acces the list of model fits by:
    
    >>> scat.model_fits
    [<GbmModelFit: Comptonized, Epeak>]
    
The |GbmModelFit| class contains the all of the relevant info for a fit. Let's 
look at a few properties contained in this class:    

    >>> one_fit = scat.model_fits[0]
    >>> one_fit.time_range
    (-0.192, 0.064)
    >>> one_fit.stat_name
    'Castor C-STAT'
    >>> one_fit.stat_value
    516.5009155273438
    >>> one_fit.dof
    479
    >>> one_fit.name
    'Comptonized, Epeak'
    >>> one_fit.parameter_list()
    ['Amplitude', 'Epeak', 'Index', 'Pivot E =fix']

For the actual fit data, they are stored as |Parameter| objects, which group 
together the fit value, uncertainty, and metadata such as the parameter name 
and units. If we want to quickly see the fit information, we can simply print 
the GbmModelFit object:

    >>> print(one_fit)
    Comptonized, Epeak
       Amplitude: 0.03 +/- 0.02
       Epeak: 215.09 +/- 54.22
       Index: 0.14 +/- 0.59
       Pivot E =fix: 1.00e+02 +/- 0.00e+00

Otherwise, we can directly access the parameter fit info by attributes and 
methods of the Parameter class:

    >>> epeak = one_fit.parameters[1]
    >>> epeak.name
    'Epeak'
    >>> epeak.value
    215.09434509277344
    >>> epeak.uncertainty
    (54.219116, 54.219116)
    >>> epeak.one_sigma_range()
    (160.87522888183594, 269.31346130371094)

For each fit, the photon and energy fluxes and fluences are also recorded, and 
those are stored as special Parameter objects as well:

    >>> print(one_fit.photon_flux)
    Photon Flux: 2.81 +/- 0.44 ph/cm^2/s
    
    >>> print(one_fit.energy_fluence)
    Energy Fluence: 1.40e-07 +/- 3.03e-08 erg/cm^2

As for the detector metadata used in the fits, they are stored in 
|GbmDetectorData| objects:

    >>> scat.detectors
    [<GbmDetectorData: BGO_00; TTE;
      time range: (-0.192, 0.064) s;
      energy range: (284.65, 40108.0) keV>,
     <GbmDetectorData: NAI_01; TTE;
      time range: (-0.192, 0.064) s;
      energy range: (8.501, 903.45) keV>,
     <GbmDetectorData: NAI_02; TTE;
      time range: (-0.192, 0.064) s;
      energy range: (9.119, 902.57) keV>,
     <GbmDetectorData: NAI_05; TTE;
      time range: (-0.192, 0.064) s;
      energy range: (8.698, 909.58) keV>]

We can access various properties stored in these objects. For example:

    >>> [det.detector for det in scat.detectors]
    ['BGO_00', 'NAI_01', 'NAI_02', 'NAI_05']
    >>> [det.channel_range for det in scat.detectors]
    [(3, 124), (5, 124), (5, 124), (5, 124)]
    >>> # was the detector used in the fit?
    >>> [det.active for det in scat.detectors]
    [True, True, True, True]

In addition to these properties, the deconvolved photon counts, errors, and 
detector-convolved photon model are stored as well. These can be used for 
making the model count rate plots. For example:

    >>> # the detector-convolved photon model values for the second detector
    >>> scat.detectors[1].photon_model
    array([1.96157284e-02, 1.99605655e-02, 2.02203058e-02, 2.04198491e-02,
           2.05748081e-02, 2.06954200e-02, 2.07886994e-02, 2.08651833e-02,
           ...
           5.15691772e-06, 3.84230407e-06, 2.66911161e-06, 3.89260578e-07],
         dtype=float32)


Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.scat
   :inherited-members:


