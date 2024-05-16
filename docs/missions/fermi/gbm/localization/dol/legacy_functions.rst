.. _gbm-dol-legacy-functions:

********************************************************************************************************************
Legacy DoL Functions
********************************************************************************************************************
(:mod:`gdt.missions.fermi.gbm.localization.dol.legacy_functions`)

The methods provided in this module handle the basic functionality needed
to operate the legacy DoL localization algorithm in a way that replicates
the original Fortran code. These include legacy angular operations
as well as coordinate conversions that would typically be handled by
astropy in other parts of the GDT. A number of these methods are performed
with 32-bit floating point accuracy to match the original Fortran code.
A full list of available methods is shown in the :ref:`Reference/API
section<gbm-dol-legacy-functions_api>` at the bottom of this page.

Here we will focus on describing the core localization fit procedure :meth:`~gdt.missions.fermi.gbm.localization.dol.legacy_functions.find_best_location`.
This method is a :math:`\chi^2` minimization routine with :math:`\chi^2` defined as

.. math::

    \chi^2 = \Sigma_i \frac{(M_i - E_i)^2}{E_i}

where :math:`M_i` are the measured counts and :math:`E_i`
are the expected counts in detector :math:`i`. This assumes the measured counts
follow a Gaussian distribution with symmetric variance :math:`E_i`, which is
true for GBM given that typical backgrounds alone provide a few hundred measured
counts across 50-300 keV in each detector during most GRB emission timescales.

The expected counts :math:`E_i` are determined from the sum of the source
:math:`S_i` and background :math:`B_i` expectations for each detector

.. math::

    E_i = S_i + B_i

:math:`B_i` are typically determined from a polynomial fit to background periods (:ref:`See Localization Example<dol-example-background>`).
:math:`S_i` consists of the response of each detector :math:`R_i` multiplied by the flux :math:`f` of the source

.. math::

    S_i = f \times R_i


where the source flux :math:`f` is found by solving for the value that minimizes :math:`\chi^2`. This
is done by analytically setting the derivative of :math:`\chi^2` equal
to zero and solving for :math:`f`

.. math::

    \frac{d}{df} \chi^{2} = 0

To make this process easier, the :math:`\chi^2` is modified to
exchange the model variance :math:`var(E_i) = E_i` with the variance
estimated from the measured counts :math:`var(M_i) = M_i`

.. math::

    \chi^{2}_{\, var(M_i)} = \Sigma_i \frac{(M_i - E_i)^2}{M_i}

In this case, only the term :math:`E_i` in the numerator depends on :math:`f`

.. math::

    \frac{d}{df} \chi^{2}_{\, var(M_i)} = \frac{d}{df} \Sigma_i \frac{(M_i - E_i)^2}{M_i} = -2 \Sigma_i \frac{(M_i - E_i)}{M_i} \frac{d E_i}{df} = 0

Substituting :math:`E_i = f \times R_i + B_i` then yields

.. math::

    \frac{d E_i}{df} = R_i

.. math::

    f = \frac{a}{b}, \quad a \equiv \Sigma_i \frac{(M_i - B_i)}{M_i} R_i, \quad b \equiv \Sigma_i \frac{R_i^2}{M_i}

.. _gbm-dol-legacy-functions_api:

Reference/API
=============

.. automodapi:: gdt.missions.fermi.gbm.localization.dol.legacy_functions
   :inherited-members:

