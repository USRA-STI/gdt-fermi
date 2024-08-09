.. _notebooks:


Jupyter Notebook Tutorials
==========================

* :download:`PHAII Data Tutorial <notebooks/Fermi_GBM_PHAII_Data_Tutorial.ipynb.tar>`
  
  PHAII data is one of the primary types of science data provided by GBM, and 
  is temporally pre-binned. The two types of PHAII data are CSPEC and CTIME, 
  where CSPEC has 128 energy channels and CTIME has 8 energy channels. Learn 
  how to plot the lightcurves and count spectra of gamma-ray bursts using CSPEC 
  and CTIME data.

----

* :download:`TTE Data Tutorial <notebooks/Fermi_GBM_TTE_Data_Tutorial.ipynb.tar>`
  
  TTE (Time-Tagged Event) data is one of the primary types of science data 
  provided by GBM, and is temporally unbinned. Learn how to plot the lightcurves
  and count spectra of gamma-ray bursts using TTE data.

----

* :download:`Detector Responses Tutorial <notebooks/Fermi_GBM_Detector_Responses_Tutorial.ipynb.tar>`
  
  Detector response files allow you to compare a theoretical photon spectrum to 
  an observed count spectrum. The two types of detector response files are rsp 
  and rsp2 files, where rsp contain a single DRM (Detector Response Matrix), and 
  rsp2 contain more than one DRM. Learn how to read, manipulate, and plot 
  detector response 'rsp2' data.

----

* :download:`Position History Data Tutorial <notebooks/Fermi_GBM_Position_History_Data_Tutorial.ipynb.tar>`
  
  Position history (POSHIST) data contains the spacecraft location in orbit and 
  pointing information for an entire day. Learn how to open, read, and plot 
  position history data.

----

* :download:`Trigger Data Tutorial <notebooks/Fermi_GBM_Trigger_Data_Tutorial.ipynb.tar>`
  
  Trigger data (TRIGDAT) was designed to contain the minimum amount of data 
  required for rapid on-ground characterization and localization of triggers. 
  It contains pre-binned lightcurve data for each detector along with spacecraft
  position and attitude information. Learn how to plot the lightcurves of 
  gamma-ray bursts and the positional data of the Fermi instrument using trigger 
  data.

----

* :download:`Localizations and Skymaps Tutorial <notebooks/Fermi_GBM_Localizations_and_Skymaps_Tutorial.ipynb.tar>`
  
  GBM produces localizations for GRBs that contain the best-modeled systematic 
  uncertainty in the localization, individual detector pointings, and the 
  geocenter location as observed by Fermi. Localization information is contained
  in HEALPix FITS files. Learn how to access localizations for gamma-ray bursts 
  and to plot skymaps, both of the localizations and of other known points, 
  using HEALPix files.

----
   
* :download:`Finding GBM Data Tutorial <notebooks/Finding_GBM_Data_Tutorial.ipynb.tar>`
  
  GBM Data is hosted publicly on the HEASARC FTP server via the Fermi Science 
  Support Center, and is stored in a consistent directory structure. Learn how 
  to earch the HEASARC FTP server for trigger data and continuous data, as well 
  as how to search the GBM catalogs.

----
    
* :download:`Reduction and Export Tutorial <notebooks/Reduction_And_Export_Tutorial.ipynb.tar>`
  
  Often, we would like to reduce GBM data and model the background so that we 
  can examine the source spectrum. The following workflow will show you how you 
  can do this and export the relevant data. Learn how to plot a lightcurve, 
  complete a background fit, and write and export the results.

----

* :download:`Spectral Analysis Tutorial <notebooks/Spectral_Analysis_Tutorial.ipynb.tar>`
  
  Often, one would like to perform a spectral fit on GBM data. The following 
  workflow will guide you through a simple example of this process. Learn how 
  to fit GBM spectral data.

----

* :download:`Spectral Catalog Data Tutorial <notebooks/Fermi_GBM_Spectral_Catalog_Data_Tutorial.ipynb.tar>`
  
  GBM provides standard spectral fits for each GRB triggered on-board, which 
  are hosted as Spectral Catalog (SCat) files that contain data such as fit 
  parameters, uncertainties, fluxes/fluences, fit statistic, and covariance. 
  Learn how to download, access, and manipulate SCat data.
