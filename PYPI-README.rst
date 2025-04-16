=========
GDT-Fermi
=========

The GDT-Fermi is an extension to Gamma-ray Data Tools that adds functions specific to the Fermi mission (GBM and LAT).

This software is not subject to EAR.

Normal Installation
-------------------

If you don't plan to contribute code to the project, the recommended install method is installing from PyPI using:

.. code-block:: sh

   pip install astro-gdt-fermi
   gdt-data init

The ``gdt-data init`` is required to initialize the library after installation of astro-gdt. You do not need to
perform the initialization again if astro-gdt was already installed and initialized.  There is no harm in running
it again "just in case".

Contributing Code or Documentation
----------------------------------

If you plan to help with the development or documentation of astro-gdt, then please visit our github site at
https://github.com/USRA-STI/gdt-fermi.