.. _gdt-fermi:

****************************************************
Welcome to Fermi Gamma-ray Data Tools Documentation!
****************************************************

.. figure:: images/gdt-fermi_logo.png

The Fermi Gamma-ray Data Tools (GDT) is a toolkit for Fermi data built on the 
:external:ref:`GDT Core Package<gdt-core>` and is the next iteration of the 
`Fermi GBM Data Tools <https://fermi.gsfc.nasa.gov/ssc/data/analysis/gbm/gbm_data_tools/gdt-docs>`_. 

The Fermi Gamma-ray Space Telescope was launched on June 11, 2008 and contains
two intruments: the Large Area Telescope (LAT) surveying the sky in the 
~20 MeV--300 GeV energy range, and the Gamma-ray Burst Monitor (GBM) observing 
the full unocculted sky with 14 individual detectors from ~8 keV--40 MeV. 
Currently the toolkit services all Fermi GBM public data, and a future version
will incorporate support for some Fermi LAT data.

.. rubric:: Citing

If you use the Fermi Gamma-ray Data Tools in your research and publications, 
we would definitely appreciate an appropriate acknowledgment and citation! We 
suggest the following BibTex:

::

 @misc{GDT-Fermi,
       author = {Adam Goldstein and William H. Cleveland and Daniel Kocevski},
       title = {Fermi Gamma-ray Data Tools: v2.0.0},
       year = 2023,
       url = {https://github.com/USRA-STI/gdt-fermi}
 }
 

.. rubric:: Additional Resources
 
The Fermi Science Support Center is a fantastic resource for all things Fermi.
Specifically, for GBM, a lot of useful information about the data products can 
be found `here <https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/>`_.  For 
questions, bug reports, and comments, please visit the 
`Fermi Help Desk <https://fermi.gsfc.nasa.gov/ssc/help/>`_.

.. rubric:: Acknowledgments

The Fermi Gamma-ray Data Tools were partially funded by the Fermi Guest Investigator 
program (NNH18ZDA001N) and by Cooperative Agreement 80MSFC17M0022.

***************
Getting Started
***************
.. toctree::
   :maxdepth: 1

   install

Jupyter Notebook Tutorials
==========================

Jupyter notebooks can be found in the /notebooks directory.

* **Fermi GBM PHAII Data Tutorial:** PHAII data is one of the primary types of science data provided by GBM, and is temporally pre-binned. The two types of PHAII data are CSPEC and CTIME, where CSPEC has 128 energy channels and CTIME has 8 energy channels. Learn how to plot the lightcurves and count spectra of gamma-ray bursts using CSPEC and CTIME data.

* **Fermi GBM TTE Data Tutorial:** TTE (Time-Tagged Event) data is one of the primary types of science data provided by GBM, and is temporally unbinned. Learn how to plot the lightcurves and count spectra of gamma-ray bursts using TTE data.

* **Fermi GBM Detector Responses Tutorial:** Detector response files allow you to compare a theoretical photon spectrum to an observed count spectrum. The two types of detector response files are rsp and rsp2 files, where rsp contain a single DRM (Detector Response Matrix), and rsp2 contain more than one DRM. Learn how to read, manipulate, and plot detector response 'rsp2' data.

* **Fermi GBM Position History Data Tutorial:** Position history (POSHIST) data contains the spacecraft location in orbit and pointing information for an entire day. Learn how to open, read, and plot position history data.

* **Fermi GBM Trigger Data Tutorial:** Trigger data (TRIGDAT) was designed to contain the minimum amount of data required for rapid on-ground characterization and localization of triggers. It contains pre-binned lightcurve data for each detector along with spacecraft position and attitude information. Learn how to plot the lightcurves of gamma-ray bursts and the positional data of the Fermi instrument using trigger data.

* **Fermi GBM Localizations and Skymaps Tutorial:** GBM produces localizations for GRBs that contain the best-modeled systematic uncertainty in the localization, individual detector pointings, and the geocenter location as observed by Fermi. Localization information is contained in HEALPix FITS files. Learn how to access localizations for gamma-ray bursts and to plot skymaps, both of the localizations and of other known points, using HEALPix files.

* **Finding GBM Data Tutorial:** GBM Data is hosted publicly on the HEASARC FTP server via the Fermi Science Support Center, and is stored in a consistent directory structure. Learn how to earch the HEASARC FTP server for trigger data and continuous data, as well as how to search the GBM catalogs.

* **Reduction and Export Tutorial:** Often, we would like to reduce GBM data and model the background so that we can examine the source spectrum. The following workflow will show you how you can do this and export the relevant data. Learn how to plot a lightcurve, complete a background fit, and write and export the results.

* **Spectral Analysis Tutorial:** Often, one would like to perform a spectral fit on GBM data. The following workflow will guide you through a simple example of this process. Learn how to fit GBM spectral data.

* **Fermi GBM Spectral Catalog Data Tutorial:** GBM provides standard spectral fits for each GRB triggered on-board, which are hosted as Spectral Catalog (SCat) files that contain data such as fit parameters, uncertainties, fluxes/fluences, fit statistic, and covariance. Learn how to download, access, and manipulate SCat data.

* **Data Primitives Tutorial:** Data primitives are the data classes that define the datatypes within the GDT at the most basic level. They handle the properties and manipulation of the fundamental types of data that GBM produces. Learn how to understand and manipulate the data primitives used in the Fermi Gamma-Ray Tools (GDT) toolkit.

******************
User Documentation
******************

Fermi Definitions
=================
.. toctree::
   :maxdepth: 1

   missions/fermi/time
   missions/fermi/frame
   missions/fermi/mcilwainl
   missions/fermi/plot

Fermi GBM
=========

Instrument Definitions
----------------------

.. toctree::
   :maxdepth: 1

   missions/fermi/gbm/detectors
   missions/fermi/gbm/collection
   missions/fermi/gbm/saa
   missions/fermi/gbm/headers

Data Types
----------

.. toctree::
   :maxdepth: 1

   missions/fermi/gbm/phaii
   missions/fermi/gbm/tte
   missions/fermi/gbm/response
   missions/fermi/gbm/localization/localization
   missions/fermi/gbm/localization/dol/dol
   missions/fermi/gbm/trigdat
   missions/fermi/gbm/poshist
   missions/fermi/gbm/scat
   missions/fermi/gbm/tcat



Data Finders and Catalogs
-------------------------

.. toctree::
   :maxdepth: 1

   missions/fermi/gbm/finders
   missions/fermi/gbm/catalogs

----

*******
License
*******
.. toctree::
   :maxdepth: 1
   
   license


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
