=========
GDT-Fermi
=========

The GDT-Fermi is an extension to Gamma-ray Data Tools that adds functions specific to the Fermi mission (GBM and LAT).

Normal Installation
-------------------

If you don't plan to contribute code to the project, the recommended install method is installing from PyPI using:

.. code-block:: sh

   pip install astro-gdt-fermi

Setting up a development environment
------------------------------------

If you do want to contribute code to this project (and astro-gdt), you can use the following commands to quickly setup a
development environment:

.. code-block:: sh

   mkdir gdt-devel
   cd gdt-devel
   python -m venv venv
   . venv/bin/activate
   pip install --upgrade pip setuptools wheel
   git clone git@github.com:USRA-STI/gdt-core.git
   git clone git@github.com:USRA-STI/gdt-fermi.git
   pip install -e gdt-core/
   pip install -r gdt-core/requirements.txt
   pip install -e gdt-fermi/
   pip install -r gdt-fermi/requirements.txt

This should result in git-devel having the following directory structure::

   .
   ├── venv
   ├── gdt-core
   └── gdt-fermi

and both gdt-core and gdt-fermi installed in the virtual environment named venv.

Namespace Packaging
-------------------

GDT-Fermi can be used as an example of how we expect other missions to contribute extensions to the Gamma-ray Data Tools.
GDT-Fermi uses namespace packaging which requires that src/gdt and src/gdt/missions to **NOT** contain a file named
``__init__.py``. You can learn more about Namespace packages by reading `PEP-420 <https://peps.python.org/pep-0420/>`_.