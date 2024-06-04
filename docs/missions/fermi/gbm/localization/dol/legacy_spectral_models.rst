.. _gbm-dol-legacy-spectral-models:

.. |band_method| replace:: :meth:`~gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models.legacy_band`

.. |comp_method| replace:: :meth:`~gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models.legacy_comp`

.. |pl_method| replace:: :meth:`~gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models.legacy_pl`

********************************************************************************************************************************
Legacy DoL Spectral Models
********************************************************************************************************************************
(:mod:`gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models`)

This module provides legacy definitions of Band, Power law, and Cutoff power law spectral models
with support for 32-bit floating point operations. The functional definitions
of these models and their corresponding API methods are shown in Table 1.

.. table:: Table 1. Spectral Models

    +------------------+---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Model            | Method        | Functional Form                                                                                                                                                                  |
    +==================+===============+==================================================================================================================================================================================+
    | Band             | |band_method| | :math:`\frac{dN}{dE} = ( \frac{E}{100 \, keV} )^\alpha \, e^{-(\alpha + 2) \, E / E_{peak}} \quad E < E_b, \, \, E_b \equiv (\alpha - \beta) \, E_{peak} \, / \, (\alpha + 2)`   |                                                                                                                                                                                                                                                          
    +                  +               +                                                                                                                                                                                  +
    |                  |               | :math:`\quad \quad \, ( \frac{E}{100 \, keV} )^\beta \,  e^{\alpha - \beta} \, ( \frac{(\alpha - \beta) \, E_{peak}}{100 \, keV (\alpha + 2)} )^{\alpha - \beta} \quad E \ge E_b`|
    +------------------+---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Cutoff power law | |comp_method| | :math:`\frac{dN}{dE} = ( \frac{E}{100 \, keV} )^{index} \, e^{-E / E_{peak} (2+ index)}`                                                                                         |
    +------------------+---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
    | Power law        | |pl_method|   | :math:`\frac{dN}{dE} = ( \frac{E}{100 \, keV} )^{index}`                                                                                                                         |
    +------------------+---------------+----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Additionally, this module provides access to pre-computed detector response matrices
for the six spectra shown in Table 2.


.. table:: Table 2. Pre-computed response matrices

    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | Name           | Model            | Parameters                                                 | Response File(s)     | Energy Range [keV] |
    +================+==================+============================================================+======================+====================+
    | ``band_hard``  | Band             | :math:`\alpha = 0.0, \, \beta = -1.5, \, E_{peak} = 1000.0`| ``band_hard_50_300`` | 50-300             |
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | ``band_norm``  | Band             | :math:`\alpha = -1.0, \, \beta = -2.3, \, E_{peak} = 230.0`| ``band_norm_50_300`` | 50-300             |
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | ``band_soft``  | Band             | :math:`\alpha = -2.0, \, \beta = -3.4, \, E_{peak} = 70.0` | ``band_soft_50_300`` | 50-300             |
    +                +                  +                                                            +                      +                    +
    |                |                  |                                                            | ``band_soft_5_50``   | 5-50               |
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | ``comp_hard``  | Cutoff power law | :math:`index = -0.25, \, E_{peak} = 1000.0`                | ``comp_hard_50_300`` | 50-300             | 
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | ``comp_norm``  | Cutoff power law | :math:`index = -1.15, \, E_{peak} = 350.0`                 | ``comp_norm_50_300`` | 50-300             | 
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+
    | ``comp_soft``  | Cutoff power law | :math:`index = -1.95, \, E_{peak} = 50.0`                  | ``comp_soft_50_300`` | 50-300             | 
    +----------------+------------------+------------------------------------------------------------+----------------------+--------------------+

The spectral parameter string needed to initialize the :class:`~gdt.missions.fermi.gbm.localization.dol.legacy_dol.legacy_DoL`
class with a pre-computed response matrix can be obtained by invoking the spectrum name shown in Table 2 from the
:mod:`~gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models` module:

    >>> from gdt.missions.fermi.gbm.localization.dol import legacy_spectral_models
    >>> legacy_spectral_models.band_hard
    'band,alpha=0.0,beta=-1.5,epeak=1000.0,amp=10.0'

Similarly, the path to the pre-computed response file can be obtained by invoking the listed response
file name:

    >>> from gdt.missions.fermi.gbm.localization.dol import legacy_spectral_models
    >>> legacy_spectral_models.band_hard_50_300
    '/path/to/site-packages/gdt/missions/fermi/gbm/localization/dol/data/band_1deg_50_300_hard.npy'

Each response file consists of the simulated response in all 14 GBM detectors over the specified
energy range computed for source locations spanning the entire sky using a 1-degree grid.
The shape of the response matrix is (14, 41168).

Users looking to work directly with the response grid can load it with

    >>> import numpy as np
    >>> rsp = np.load(path, allow_pickle=True, encoding='bytes').item()[b"table"]

where ``path`` is one of the response files listed in Table 2. See [1]_ for more
details on the creation of the response grid.

References:
"""""""""""

.. [1] `Connaughton, V. et al. 2015, ApJ, 216, 32 <https://iopscience.iop.org/article/10.1088/0067-0049/216/2/32>`_


Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.localization.dol.legacy_spectral_models
   :inherited-members:

