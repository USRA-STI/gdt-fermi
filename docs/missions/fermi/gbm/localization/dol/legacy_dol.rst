.. _gbm-dol-legacy:

*********************************************************************************************
Fermi GBM Localizations using DoL
*********************************************************************************************
(:mod:`gdt.missions.fermi.gbm.localization.dol.legacy_dol`)

DoL related documentation goes here.

Standard behavior is a localization fit over three spectral templates from 50-300 keV.

flux norm guess is based on 

chi2 is computed as follows

S = norm * DRM
chi2 = (M - (S + B))**2 / (S + B)

where var = S + B is the model variance.

M = measured counts, B = background counts model, S = signal model


Another typical use for SGRs is the soft fit from 5-50 keV.

spec = [("Soft_5_50:", legacy_spectral_models.band_soft)]
locrates = [legacy_spectral_models.band_soft_5_50]
crange = [1, 2]


Should probably define shape of locrates

Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.localization.dol.legacy_dol
   :inherited-members:

