# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: William Cleveland
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# This software is not subject to EAR.
#
# Very closely based on the program GBM Trigger Classification.
# Written by:
#               Michael Briggs
#               University of Alabama in Huntsville (UAH)
#
# Based on algorithms developed by:
#               Chip Meegan, David Perrin and Eli Sidman.
#
# Published papers on the algorith used:
#
#      Perrin, D. J., Sidman, E. D., Meegan, C. A., Briggs, M. S., Connaughton, V.
#      "GLAST Burst Monitor Trigger Classification Algorithm",
#      American Astronomical Society, HEAD meeting #8, id.18.14; Bulletin of the American Astronomical Society,
#      Vol. 36, p.943 2004HEAD....8.1814P
#
#      Briggs, M. S., Connaughton, V., Paciesas, W., Preece, R., Meegan, C. A., Fishman, G., Kouveliotou, C.,
#      Wilson-Hodge, C., Diehl, R., Greiner, J., von Kienlin, A., Lichti, G., Steinle, H., Kippen, R. M.
#      "GLAST Burst Monitor On-Board Triggering, Locations and Event Classification",
#      THE FIRST GLAST SYMPOSIUM. AIP Conference Proceedings, Volume 921, pp. 450-451 (2007). 2007AIPC..921..450B
#
#      Briggs, M. S., Pendleton, G. N., Kippen, R. M., Brainerd, J. J., Hurley, K. Connaughton, V., Meegan, C. A.
#      "The Error Distribution of BATSE GRB Locations",
#      The Astrophysical Journal Supplement Series, 122, 503-518, 1999 [astro-ph/9901111].
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import numpy as np
import math
from enum import IntEnum

from scipy.stats import vonmises_fisher
from astropy.coordinates import get_sun, SkyCoord, GCRS, UnitSphericalRepresentation
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.trigdat import Trigdat
from gdt.missions.fermi.frame import FermiFrame
from gdt.missions.fermi.mcilwainl import calc_mcilwain_l
from gdt.missions.fermi.gbm.detectors import GbmDetectors

__all__ = ['Classification', 'TriggerClasses']

# Configuration values
ZENITH_TO_HORIZON_ARCMIN = 6780
TGF_ALGORITHM_NUMBERS = [43, 50, 116, 117, 118, 119]
GBM_ARCMIN_TO_RADIANS = 2.908882E-4
EGROUPS = [(3, 4), (1, 2), (4, 7), (2, 2), (5, 7)]
GBM_FLOAT_PI = 3.14159

# Used to determine local particles.
DET_SRC_COUNTS_RATIO_THRESHOLD_X100 = 100

# Probabilities based on hardness ratio
HARDNESS_RATIO_PROB_X10K = {
    'GRBS': [257, 393, 589, 830, 1073, 1261, 1341, 1286, 1114, 873, 625, 85, 61, 54, 47, 40, 34, 28],
    'PARTICLES': [10, 128, 448, 1025, 1602, 1666, 1474, 1089, 576, 641, 769, 384, 128, 64, 10, 10, 10, 10],
    'SOLARFLARES': [10, 10, 10, 10, 10, 10, 10, 10, 10, 200, 750, 1175, 1175, 1200, 1300, 1150, 850, 2110],
    'GENERIC_SGRS': [177, 194, 212, 231, 252, 275, 299, 325, 353, 382, 413, 588, 785, 973, 1118, 1188, 1167, 1059],
    'XRAY_TRANSIENTS': [10, 80, 403, 1048, 1854, 2258, 1854, 1129, 725, 483, 161, 10, 10, 10, 10, 10, 10, 10],
    'CYGNUS_X1': [634, 793, 952, 1587, 1428, 1111, 1269, 1111, 793, 317, 10, 10, 10, 10, 10, 10, 10, 10],
    'SGR_1806_20': [177, 194, 212, 231, 252, 275, 299, 325, 353, 382, 413, 588, 785, 973, 1118, 1188, 1167, 1059],
    'GROJ_0422_32': [10, 80, 403, 1048, 1854, 2258, 1854, 1129, 725, 483, 161, 10, 10, 10, 10, 10, 10, 10],
}

# Probabilities based on McIlwain L
MCILWAIN_L_PROB_X10K = {
    'PARTICLES': [10, 10, 20, 40, 60, 80, 110, 1900, 460, 720, 1050, 980, 1280, 2120, 190, 780, 130, 50],
    'NOT_PARTICLES': [1260, 1870, 1660, 1360, 840, 480, 480, 360, 280, 290, 320, 150, 220, 160, 140, 50, 50, 20],
}

# The initial probabilities used by GBM
INITIAL_PROBABILITIES = [
    0.00,  # ERROR
    0.00,  # UNRELIABLE_LOCATION
    0.00,  # LOCAL_PARTICLES
    0.00,  # BELOW_HORIZON
    0.86,  # GRB
    0.02,  # GENERIC_SGR
    0.02,  # GENERIC_TRANSIENT
    0.03,  # DISTANT_PARTICLES
    0.03,  # SOLAR_FLARE
    0.02,  # CYG_X1
    0.01,  # SGR_1806_20
    0.01,  # GROJ_0422_32
    0.00,  # RESERVED_1
    0.00,  # RESERVED_2
    0.00,  # RESERVED_3
    0.00,  # RESERVED_4
    0.00,  # RESERVED_5
    0.00,  # RESERVED_6
    0.00,  # RESERVED_7
    0.00,  # TGF
]

# RA and DECs of known sources that cause frequent triggers.
J2000 = GCRS(obstime="J2000")
SOURCES = {
    'CYG_X1': SkyCoord(5.22884, 0.614384, unit='rad', frame=J2000),
    'SGR_1806_20': SkyCoord(4.75016, -0.356240, unit='rad', frame=J2000),
    'GROJ_0422_32': SkyCoord(1.1419, 0.5743, unit='rad', frame=J2000),
}


class TriggerClasses(IntEnum):
    """Definition of numeric trigger classes"""
    # And names for the classes:
    # "single classes" where the probabilities of these will be either 0 or 1.
    # TGF is also a "single class" and listed as below as the highest value
    ERROR = 0
    UNRELIABLE_LOCATION = 1
    LOCAL_PARTICLES = 2
    BELOW_HORIZON = 3

    # Can be anywhere on the sky, above the Earth's horizon
    GRB = 4
    GENERIC_SGR = 5
    GENERIC_TRANSIENT = 6

    # distant particles are assumed to come from the Earth's horizon
    DISTANT_PARTICLES = 7

    # point sources of known location (moving in one case only)
    SOLAR_FLARE = 8
    CYG_X1 = 9
    SGR_1806_20 = 10
    GROJ_0422_32 = 11

    # Reserved trigger classes
    RESERVED_1 = 12
    RESERVED_2 = 13
    RESERVED_3 = 14
    RESERVED_4 = 15
    RESERVED_5 = 16
    RESERVED_6 = 17
    RESERVED_7 = 18

    # Largest trigger class ID allowed
    TGF = 19


class Classification:
    """Class for calculating source classification probabilities"""

    def __init__(self):
        self.probs = np.array(INITIAL_PROBABILITIES)
        self.dirty = False

    def sum_channels(self, arr: np.array, low: int, high: int):
        """Sums array values over a given energy channel range

        Args:
            arr (np.array): The array to sum
            low (int):      Lowest energy channel included in the sum
            high (int):     Highest energy channel included in the sum

        Returns:
            (float)
        """
        return arr[:, low:high + 1].sum(axis=1)

    def classify_trigdat(self, trigdat: Trigdat, use_loc: bool = True, 
                         verbose: bool = False):
        """Convenience method for quickly applying classification method
        to GBM triggered data.

        Args:
            trigdat (:class:`Trigdat`): The triggered data object
            use_loc (bool):    Use location info when True
            verbose (bool):    Display trigdat info to screen when True

        Returns:
            (list)
        """
        # gather trigdat info
        trigger_alg = trigdat.headers['PRIMARY']['TRIG_ALG']
        fsw_loc = trigdat.fsw_locations[-1]
        location = fsw_loc.location[:2]
        location_err = fsw_loc.location[2]

        if verbose:
            print("\nTrigger Algorithm = {}".format(trigger_alg))
            print("Location = (RA {:.3f}, Dec {:.3f}) deg".format(*location))
            print("Location Err = {:.3f} deg".format(location_err))

        # latest maxrates packet
        maxrates = trigdat.maxrates[-1]
        rates = maxrates.all_rates.reshape(maxrates.num_dets, maxrates.num_chans)
        maxrates_timescale = min(maxrates.timescale, 1024)
        time = Time(maxrates.time_range[-1], format='fermi')
        frame = trigdat.poshist.at(time)

        # background information scaled to match MAXRATES timescale
        bkg = trigdat.backrates
        bkg_timescale = min(bkg.time_range[1] - bkg.time_range[0], 1.024) * 1000
        bkg_rates = bkg.all_rates.reshape(bkg.num_dets, bkg.num_chans) * bkg_timescale / maxrates_timescale

        if verbose:
            print("\nRates\n", rates)
            print("\nBackground\n", bkg_rates)

        return self.classify(trigger_alg, frame, location, location_err, rates, bkg_rates, use_loc)

    def classify(self, trig_alg: int, frame: FermiFrame, loc: tuple,
                 loc_total_error: float, rates: np.array, bkg: np.array, use_loc: bool = False):
        """Method to calculate source classification probabilities using
        information from a trigger. Returns the two detectors closest to
        the source location when use_loc=True, otherwise returns top
        two detectors with highest rates in energy channels (3, 4)

        Args:
            trig_alg (int):          Trigger algorithm number
            frame (:class:`~..FermiFrame`): Spacecraft frame at the trigger time
            loc (tuple):             Best-fit source location given as
                                     (ra, dec) in degrees
            loc_total_error (float): Radius of 68% containment for
                                     the localization in degrees
            rates (np.array):        Detector rates from the latest
                                     maxrates data packet
            bkg (np.array):          Background rates
            use_loc (bool):          Uses location when True

        Returns:
            (list)
        """
        # Format location as SkyCoord for easier coordinate transformations later on
        loc = SkyCoord(*loc, unit='deg', frame=J2000)

        # Check for TGF
        if self.is_tgf(trig_alg):
            return []

        # Compute background subtracted counts
        counts = (rates - bkg).clip(min=0)

        # Check for local particles
        if self.is_local_particle(rates, bkg):
            return []

        # Get the counts for the energy group 0 and energy group 1
        e0 = self.sum_channels(counts, *EGROUPS[0])
        e1 = self.sum_channels(counts, *EGROUPS[1])

        if use_loc:
            # Choose the two closest detectors to the source location
            angles = [det.skycoord(frame).separation(loc, origin_mismatch="ignore").deg[0] for det in GbmDetectors if det.is_nai()]
            det_order = np.argsort(angles)
        else:
            # Chooses the two detector hardness by the number of counts in e0
            nai_e0 = [e0[det.number] for det in GbmDetectors if det.is_nai()]
            det_order = np.flip(np.argsort(nai_e0))

        # Hardness ratio of the highest detector and second highest detector
        first = e1[det_order[0]] / e0[det_order[0]]
        second = e1[det_order[1]] / e0[det_order[1]]

        if first == float("inf") or second == float("inf"):
            raise ZeroDivisionError("Detectors {} & {} caused a divide by zero error".format(det_order[0], det_order[1]))

        two_det_hardness_ratio = (first + second) / 2.0

        # Perform bayesian probability
        self.bayesian(frame, loc, loc_total_error, two_det_hardness_ratio)

        return det_order[:2]

    def bayesian(self, frame: FermiFrame, loc: tuple, loc_total_error: float,
                 two_det_hardness_ratio: float):
        """Computes bayesian probabilities for all source classes

        Args:
            frame (:class:`~..FermiFrame`): Spacecraft frame at trigger time
            loc (tuple):                    Best-fit source location given as
                                            (ra, dec) in degrees
            loc_total_error (float):        68% containment (stat + sys) radius
                                            for the localization in degrees
            two_det_hardness_ratio (float): Average hardness ratio in
                                            the closest two detectors
        """

        # We don't want to perform any bayesian calculations on a used array
        if self.dirty:
            raise ValueError("Probability array is dirty. Create a new one")

        # Check the trigger location against the S/C horizon
        dist_from_zenith = frame.sc_gcrs.separation(loc, origin_mismatch="ignore").rad
        if dist_from_zenith > ZENITH_TO_HORIZON_ARCMIN * GBM_ARCMIN_TO_RADIANS:
            self.single_possibility(TriggerClasses.BELOW_HORIZON)
            return

        # Calculate McIlwain L
        mcilwain_l = calc_mcilwain_l(frame.earth_location.lat.value, 
                                     frame.earth_location.lon.value)

        # calculate the likelihood for each model based on the localization
        # For sources with unknown locations, the probability is uniform across 
        # the visible sky.
        solid_angle_above_horizon = 2.0 * math.pi * (1.0 - math.cos(ZENITH_TO_HORIZON_ARCMIN * GBM_ARCMIN_TO_RADIANS))
        uniform_prob = 1.0 / solid_angle_above_horizon

        self.probs[TriggerClasses.GRB] *= uniform_prob
        self.probs[TriggerClasses.GENERIC_SGR] *= uniform_prob
        self.probs[TriggerClasses.GENERIC_TRANSIENT] *= uniform_prob

        # Relation between Von Mises-Fisher kapp to sigma: 
        # Eq. 8 in Briggs et al 1999 ApJS 122 503; good out ~20 degrees
        kappa = 1.0 / (0.66 * np.radians(loc_total_error))

        # Distant Particles are assumed to be on the Earth's limb.
        # (Calculated as the nearest point on the limb to the trigger location)
        ang_dist = ZENITH_TO_HORIZON_ARCMIN * GBM_ARCMIN_TO_RADIANS - dist_from_zenith
        self.probs[TriggerClasses.DISTANT_PARTICLES] *= vonmises_fisher.pdf(
            SkyCoord(0, 0.5 * np.pi - ang_dist, unit='rad').cartesian.xyz,
            SkyCoord(0, 0.5 * np.pi, unit='rad').cartesian.xyz, kappa)

        # Point sources of known Locations
        sun_coord = get_sun(frame.obstime).represent_as(UnitSphericalRepresentation)
        self.probs[TriggerClasses.SOLAR_FLARE] *= vonmises_fisher.pdf(
            sun_coord.to_cartesian().xyz,
            loc.cartesian.xyz, kappa)

        self.probs[TriggerClasses.CYG_X1] *= vonmises_fisher.pdf(
            SOURCES['CYG_X1'].cartesian.xyz,
            loc.cartesian.xyz, kappa)

        self.probs[TriggerClasses.SGR_1806_20] *= vonmises_fisher.pdf(
            SOURCES['SGR_1806_20'].cartesian.xyz,
            loc.cartesian.xyz, kappa)

        self.probs[TriggerClasses.GROJ_0422_32] *= vonmises_fisher.pdf(
            SOURCES['GROJ_0422_32'].cartesian.xyz,
            loc.cartesian.xyz, kappa)

        # calculate the likelihood for each model based on the spectral hardness
        # hardness ratio is source counts in 15-50 keV / source counts in 50-30
        # keV.
    
        # The table has two bin widths, with widths of 0.1 below HR = 1.0 and 
        #  0.5 widths above.
        
        # Bin   HR Range
        #  0   0.0 -- 0.1
        #  1   0.1 -- 0.2
        #  2   0.2 -- 0.3
        # ...
        #  9   0.9 -- 1.0
        # 10   1.0 -- 1.5
        # 11   1.5 -- 2.0
        # ...
        # 16   4.0 -- 4.5
        # 17   >= 4.5
        if two_det_hardness_ratio < 1.0:
            hr_bin = int(10 * two_det_hardness_ratio)
        else:
            hr_bin = int(2 * two_det_hardness_ratio + 8)

        if two_det_hardness_ratio <= 0.0:
            hr_bin = 0
        if hr_bin > 17:
            hr_bin = 17

        self.probs[TriggerClasses.GRB] *= HARDNESS_RATIO_PROB_X10K['GRBS'][hr_bin]
        self.probs[TriggerClasses.DISTANT_PARTICLES] *= HARDNESS_RATIO_PROB_X10K['PARTICLES'][hr_bin]
        self.probs[TriggerClasses.SOLAR_FLARE] *= HARDNESS_RATIO_PROB_X10K['SOLARFLARES'][hr_bin]
        self.probs[TriggerClasses.GENERIC_SGR] *= HARDNESS_RATIO_PROB_X10K['GENERIC_SGRS'][hr_bin]
        self.probs[TriggerClasses.GENERIC_TRANSIENT] *= HARDNESS_RATIO_PROB_X10K['XRAY_TRANSIENTS'][hr_bin]
        self.probs[TriggerClasses.CYG_X1] *= HARDNESS_RATIO_PROB_X10K['CYGNUS_X1'][hr_bin]
        self.probs[TriggerClasses.SGR_1806_20] *= HARDNESS_RATIO_PROB_X10K['SGR_1806_20'][hr_bin]
        self.probs[TriggerClasses.GROJ_0422_32] *= HARDNESS_RATIO_PROB_X10K['GROJ_0422_32'][hr_bin]

        # calculate the likelihood for each model based on the McIlwain L
        # This table has bins all with width 0.05, starting at McIlwain L value 1.0
        
        # Bin  McIlwain L range
        #  0   1.00 -- 1.05
        #  1   1.05 -- 1.10
        #  2   1.10 -- 1.15
        # ...
        # 16   1.80 -- 1.95
        # 17   >= 1.95
        ml_bin = int(20 * (mcilwain_l - 1.0))

        if mcilwain_l <= 0.0:
            ml_bin = 0
        if ml_bin > 17:
            ml_bin = 17

        self.probs[TriggerClasses.DISTANT_PARTICLES] *= MCILWAIN_L_PROB_X10K['PARTICLES'][ml_bin]
        self.probs[TriggerClasses.GRB] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.SOLAR_FLARE] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.GENERIC_SGR] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.GENERIC_TRANSIENT] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.CYG_X1] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.SGR_1806_20] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]
        self.probs[TriggerClasses.GROJ_0422_32] *= MCILWAIN_L_PROB_X10K['NOT_PARTICLES'][ml_bin]

        # Normalize the probabilities
        self.probs /= self.probs.sum()

        self.dirty = True

    def single_possibility(self, trig_class: int):
        """Set the probability array to a single possibility

        Args:
            trig_class (int): Index of the source class to assign
                              as the single possibility
        """
        self.probs = np.zeros(len(INITIAL_PROBABILITIES))
        self.probs[trig_class] = 1.0
        self.dirty = True

    def is_tgf(self, trig_alg: int):
        """Test if this is a TGF and set the probability array accordingly

        Args:
            trig_alg (int): Trigger algorithm number

        Returns:
            (bool)
        """
        if trig_alg in TGF_ALGORITHM_NUMBERS:
            self.single_possibility(TriggerClasses.TGF)
            return True
        return False

    def set_unreliable_location(self):
        """Set the probability array to unreliable location"""
        self.single_possibility(TriggerClasses.UNRELIABLE_LOCATION)

    def is_local_particle(self, trig_rates: np.array, bkg_rates: np.array):
        """Check if this is a local particle trigger

        Args:
            trig_rates (np.array): Detector rates in the trigger window
            bkg_rates (np.array): Background rates

        Returns:
            (bool)
        """
        # Use the ratio of the count rate of the detector with the min 
        # source count rate / the count rate of the detector with the 
        # max source count rate.

        # For point sources, some detectors aren't illuminated, so the minimum 
        # count rate will tend to be low, and the min_max ratio is low. T
        # The ratio tends approach one for local particles  (non-point source).

        egroup0_counts = self.sum_channels(trig_rates, *EGROUPS[0])
        egroup0_background = self.sum_channels(bkg_rates, *EGROUPS[0])

        min_source_counts = None
        max_source_counts = None
        for idx in range(trig_rates.shape[0]):
            src_counts = egroup0_counts[idx] - egroup0_background[idx]
            if min_source_counts is None or src_counts < min_source_counts:
                min_source_counts = src_counts
            if max_source_counts is None or src_counts > max_source_counts:
                max_source_counts = src_counts
        min_max_src_count_ratio = min_source_counts / max_source_counts

        if min_max_src_count_ratio > DET_SRC_COUNTS_RATIO_THRESHOLD_X100 / 100.0:
            self.single_possibility(TriggerClasses.LOCAL_PARTICLES)
            return True
        return False

    def rankings(self):
        """Returns the list of ranked source classes and probability
        
        Returns:
            (list)
        """
        rank = np.flip(np.argsort(self.probs))
        result = list()
        for r in rank:
            result.append((TriggerClasses(r), self.probs[r]))
        return result

    def __repr__(self):
        """(str): Representation of the object in a python console"""
        s = 'Classification('
        for tc in TriggerClasses:
            s += '\n               {:19s} {:.4f}'.format(tc.name, self.probs[tc])
        s += ')'
        return s

    def _repr_html_(self):
        """(str): Representation of the object in a python console"""
        s = '<p>Classification:</p><table>'
        s += '<tr><th>classification</th><th>probability</th></tr>'
        for tc in TriggerClasses:
            s += '<tr><td>{}</td><td>{:.4f}</td></tr>'.format(tc.name, self.probs[tc])
        s += '</table>'
        return s
