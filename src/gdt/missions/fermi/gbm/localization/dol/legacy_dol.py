# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Developed by: Joshua Wood
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Very closely based on the program DoL (Daughter of Locburst).
# Written by:
#               Valerie Connaughton
#               University of Alabama in Huntsville (UAH)
#
# Included in the Gamma-ray Data Tools with permission from UAH
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

import os
import sys
import argparse
import numpy as np

from . import legacy_spectral_models, legacy_functions
from ... import localization

legacy_dtorad = legacy_functions.legacy_dtorad

class legacy_DoL():
    r""" The legacy DoL localization code for Fermi GBM with float32 support
    
    Parameters:
        spec (list):  List of spectrum definitions given as (name, spectrum).
                      Example:
                      
                      .. code-block:: python
                        
                        spec = [("Hard", "comp,index=-0.25,epeak=1000.,amp=10.0"),
                                ("Normal", "comp,index=-1.15,epeak=350.0,amp=10.0"),
                                ("Soft", "comp,index=-1.95,epeak=50.0,amp=10.0")]
        
        locrates (list): List of paths to locrates file for each entry in ``spec``.
                         Example:

                         .. code-block:: python

                           locrates = ["locrates_1deg_50_300_hard_n.npy",
                                       "locrates_1deg_50_300_normal_n.npy",
                                       "locrates_1deg_50_300_soft_n.npy"]

        usedet  (np.array, optional): A boolean array containing True (1) if 
                                      the detector is to be included or False 
                                      (0) if a detector is not to be included.  
                                      If not set, then all detectors (including 
                                      BGOs) are used
        nen (int, optional): Number of energy bins. Default is 8.
        ndet (int, optional): Total number of detectors. Default is 14.
        ldet (int, optional): Max number of detectors for use in localization. 
                              Default is 12.
    """
    ver_no = "4.15a"
    """(str): Version number of the routine"""

    # default number of detectors and energy bins
    default_nen = 8
    """(int): Default number of energy bins for localization"""
    default_ndet = 14
    """(int): Default total number of detectors"""
    default_ldet = 12
    """(int):  Default number of detectors for localization"""

    def __init__(self, spec=None, locrates=None, usedet=None,
                 nen=default_nen, ndet=default_ndet, ldet=default_ldet,
                 verbose=False):
        # use default spectra when none are provided
        if spec is None:
            spec = [("Hard: ", legacy_spectral_models.band_hard),
                    ("Normal: ", legacy_spectral_models.band_norm),
                    ("Soft: ", legacy_spectral_models.band_soft)]

        # use default locrates files when none are provided
        if locrates is None:
            locrates = [legacy_spectral_models.band_hard_50_300,
                        legacy_spectral_models.band_norm_50_300,
                        legacy_spectral_models.band_soft_50_300]

        # set class properties
        self.nen = nen
        self.ndet = ndet
        self.ldet = ldet
        self.spec = spec
        """(list): List of spectrum definitions given as (name, spectrum)."""
        self.locrates = locrates
        """(list): List of locrates tables for each spectrum"""
        self.verbose = verbose
        """(bool): Print things to screen when True"""

        # detector selection
        if usedet is None:
            self.usedet = np.ones(14).astype('bool')
            self.usedet[-2:] = False  # default selection excludes BGO
        else:
            self.usedet = usedet.astype('bool')

        self.nai_az, self.nai_zen, self.nai_unit_vec, self.back_unit_vec = \
            legacy_functions.initialize_det_geometry(self.verbose)

        self.initialize_sky_grid(self.locrates[0])

        # check if user set a bap directory
        try:
            self.bapdir = os.environ["DOL_BAPDIR"]
        except BaseException:
            self.bapdir = ""

    # END __init__()

    def eval(self, crange, mrates, brates, sduration, bgduration,
             sc_pos, sc_quat, energies, fra, fdec, sc_time,
             scat_opt, fname="", odir=""):
        r""" Evaluate DoL localization and write results to files.

        Args:
            crange (np.ndarray(2, int32)): Array of energy channel IDs chosen 
                                           for localization
            mrates (np.ndarry(``ndet`` * ``nen``, int32)): Array of measured detector 
                                                   counts in each energy channel
            brates (np.ndarry(``ndet`` * ``nen``, int32)): Array of estimated background 
                                                   counts in each energy channel
            sduration (float32): Timescale of measured detector counts
            bgduration (float32): Timescale of estimated background counts
            sc_pos (np.ndarray(3, float32)): Spacecraft position x, y, z 
                                             coordinates
            sc_quat (np.ndarray(4, float64)): Quaternion of the space craft. 
                                              Last element is the scalar field.
            energies (np.ndarray(``nen`` + 1, float32)): Array with end points of 
                                                     energy channel bins
            fra (float32): Initial right ascension value for localization in 
                           degrees
            fdec (float32): Initial declination value for localization in degrees
            sc_time (int64): Spacecraft time in MET seconds
            scat_opt (int32): Scattering option. 1 means include scattering when 
                              computing the expected detector counts for a given 
                              spectrum. 0 means do not include scattering when 
                              computing expected detector counts.
            fname (str): Output file name string
            odir (str): Directory path to output file

        Returns:
            (dict): Dictionary object with information about the best localization
        """
        loc = self.localize(crange, mrates, brates, sduration, bgduration,
                            sc_pos, sc_quat, energies, sc_time, scat_opt)

        scx, scy, scz = loc["scx"], loc["scy"], loc["scz"]

        # information about the initial provided location
        ini = {"ra": fra / legacy_dtorad, "dec": fdec / legacy_dtorad}
        ini["pos"] = legacy_functions.ang_to_cart_dec(ini["ra"], ini["dec"])
        ini["az"], ini["zen"] = legacy_functions.j2000_to_sc(scx, scy, scz, ini["pos"])
        ini["xyz"] = legacy_functions.ang_to_cart_zen(ini["az"], ini["zen"])
        loc["initial"] = ini

        # information about the best estimated location
        best = loc["best"]
        best["pos"] = legacy_functions.ang_to_cart_zen(best["az"], best["zen"])
        best["ra"], best["dec"] = legacy_functions.sc_to_j2000(scx, scy, scz, best["pos"])
        best["xyz"] = legacy_functions.ang_to_cart_dec(best["ra"], best["dec"])
        best["lii"], best["bii"] = legacy_functions.eq2000_to_gal_r(best["ra"], best["dec"])

        # information about the sun
        sun = {k: v for k, v in zip(["ra", "dec"], legacy_functions.sun_loc(sc_time))}
        sun["xyz"] = legacy_functions.ang_to_cart_dec(sun["ra"], sun["dec"])
        sun["angle"] = legacy_functions.get_good_angle(sun["xyz"], best["xyz"])
        loc["sun"] = sun

        # information about the Earth center
        geo = loc["geo"]
        geo["angle"] = legacy_functions.get_good_angle(-loc["sc_pos"], best["xyz"])

        # compute angles between detectors and initial, best, geo locations
        loc["det_ang_initial"] = legacy_functions.get_det_geometry(
            ini["az"], ini["zen"], self.nai_az, self.nai_zen)
        loc["det_ang_best"] = legacy_functions.get_det_geometry(
            best["az"], best["zen"], self.nai_az, self.nai_zen)
        loc["det_ang_geo"] = legacy_functions.get_det_geometry(
            geo["az"], geo["zen"], self.nai_az, self.nai_zen)

        maxdet = loc["maxdet"]
        loc["err_chip"] = np.float32(10.) \
                          / np.sqrt((loc["c_mrates"][maxdet] - loc["c_brates"][maxdet]) \
                                    / np.cos(loc["det_ang_best"][maxdet]) \
                                    / np.float32(126.) / sduration)

        if self.verbose:
            self.stdout(loc)

        if len(fname):
            path = self.bapdir if len(self.bapdir) else odir
            path = os.path.join(path, "chi2grid_bn%s_v00.dat" % fname.strip())
            self.write_chi2grid(path, scx, scy, scz, sc_pos, best["chi2"])

            path = os.path.join(odir, "dol_%s.txt" % fname.strip())
            self.write_summary(path, loc)

        return loc

    # END eval()

    def localize(self, crange, mrates, brates, sduration, bgduration,
                 sc_pos, sc_quat, energies, sc_time, scat_opt):
        r""" Perform the localization routine over all spectra.

        Args:
            crange (np.ndarray(2, int32): Array of energy channel IDs chosen 
                                          for localization
            mrates (np.ndarry(``ndet`` * ``nen``, int32)):
                Array of measured detector counts in each energy channel
            brates (np.ndarry(``ndet`` * ``nen``, int32)):
                Array of estimated background counts in each energy channel
            sduration (float32): Timescale of measured detector counts
            bgduration (float32): Timescale of estimated background counts
            sc_pos (np.ndarray(float32)): Spacecraft position x, y, z 
                                          coordinates
            sc_quat (np.ndarray(4, float64)): Quaternion of the spacecraft. 
                                              Last element is the scalar field.
            energies (np.ndarray(``nen`` + 1, float32)): Array with end points 
                                                         of energy channel bins
            sc_time (int64): Spacecraft time in MET seconds
            scat_opt (int32): Scattering option. 1 means include scattering when 
                              computing the expected detector counts for a given 
                              spectrum. 0 means do not include scattering when 
                              computing expected detector counts.

        Returns:
            (dict): Dictionary object with information about best localization
        """
        # correct rates for detector deadtimes
        c_mrates, c_brates, cenergies, usedet, maxdet, signif, deadtime = \
            legacy_functions.deadtime_correct(crange, mrates, brates, sduration,
                                              bgduration, energies, self.verbose)

        # determine location of Earth center
        geodir, geo_az, geo_zen, scx, scy, scz = \
            legacy_functions.get_geocenter(sc_quat, sc_pos, self.verbose)

        if self.verbose:
            print(" geocenter az and zen")
            print("    %.14f    %.14f" %
                  (geo_az * legacy_dtorad, geo_zen * legacy_dtorad))

        # initial values for variables in the loop
        scattered_rates = None
        best = {'nchi2': np.float32(9999999.)}

        # compute localization for each spectrum
        for ispec, (name, spec) in enumerate(self.spec):

            # compute scattering matrix when requested
            if (scat_opt == 1) & (crange[1] > 2) and (scattered_rates is None):
                scattered_rates, _ = legacy_functions.compute_scat(
                    self.npoints, self.sky_grid, cenergies, geodir,
                    self.nai_az, self.nai_zen, self.nai_unit_vec,
                    self.back_unit_vec, front_only=True)

            if scattered_rates is None:
                # do not include scattering in predicted detector rates
                loctable_entries = self.locrates[ispec]
            else:
                # include scattering effect in predicted detector rates
                atm_scattered_rates = legacy_functions.add_scat(
                    self.npoints, self.sky_grid, scattered_rates, None, cenergies, spec)
                loctable_entries = np.int32(
                    np.float32(self.locrates[ispec]) + atm_scattered_rates)

            # find best guess for source location based on chi-square metric
            gindex, chi2, nchi2, rchi2, guess_az, guess_zen, guess_loc_err, \
                loc_reliable = legacy_functions.find_best_location(
                self.ndet, self.npoints, usedet.size, usedet,
                loctable_entries, sduration, c_mrates, c_brates,
                None, self.verbose)

            # store best answer as a dict object
            if nchi2 < best["nchi2"]:
                best["chi2"] = chi2
                best["nchi2"] = nchi2
                best["err"] = guess_loc_err
                best["ispec"] = ispec
                best["index"] = gindex

            if self.verbose:
                print(" {0} {1} {2:.9} {3:.9} {4:.9}".format(
                    name, gindex + 1, nchi2, guess_az, guess_zen))

        # END for (ispec)

        # NOTE: This is a logic error to do this here because it will always
        # grab the status of the final localization, not the best one
        best["loc_reliable"] = loc_reliable

        if len(self.spec) == 1:
            best["ispec"] = 0
        else:
            best["ispec"] += 1

        best["az"] = np.float32(loctable_entries[0][best["index"]]) / np.float32(60.) / np.float32(legacy_dtorad)
        best["zen"] = np.float32(loctable_entries[1][best["index"]]) / np.float32(60.) / np.float32(legacy_dtorad)

        if self.verbose:
            dtorad = np.float32(legacy_dtorad)
            print(" Chosen: {0} {1:.9} {2:.9} {3:.9}".format(
                best["index"] + 1, best["nchi2"],
                best["az"] * dtorad, best["zen"] * dtorad))

        # return a dictionary with relevant information       
        return {"best": best, "deadtime": deadtime,
                "geo": {"dir": geodir, "az": geo_az, "zen": geo_zen},
                "scx": scx, "scy": scy, "scz": scz, "sc_pos": sc_pos,
                "c_mrates": c_mrates, "c_brates": c_brates,
                "cenergies": cenergies, "signif": signif, "maxdet": maxdet}

    # END localize()

    def write_chi2grid(self, path, scx, scy, scz, sc_pos, chi2):
        r"""Write grid of chi2 values on the sky to a text file.

        Args:
            path (str): Path to output file
            scx (np.ndarray(3, float64)): Vector along +x axis of spacecraft
            scy (np.ndarray(3, float64)): Vector along +y axis of spacecraft
            scz (np.ndarray(3, float64)): Vector along +z axis of spacecraft
            sc_pos (np.ndarray(3, float64)): Spacecraft position
            chi2 (np.ndarray(``npoints``, float32)): Grid of chi2 values on the 
                                                     sky
        """
        # cartesian position vectors of every point on the sky
        pos = legacy_functions.ang_to_cart_zen(self.sky_grid[0], self.sky_grid[1])
        # convert position vectors to ra & dec grid
        j2000grid = np.array(legacy_functions.sc_to_j2000(scx, scy, scz, pos))
        # determine which points are visible
        visible = legacy_functions.get_occult(sc_pos, self.npoints, j2000grid)

        # open file and write values
        f = open(path, 'w')
        f.write("%12d\n" % self.npoints)
        for i in range(self.npoints):
            f.write("%6.1f%7.1f%12.2f%2d%7.1f%7.1f\n" %
                    (self.locrates[0][0][i] / np.float32(60.),
                     self.locrates[0][1][i] / np.float32(60.),
                     chi2[i], visible[i],
                     legacy_dtorad * j2000grid[0][i],
                     legacy_dtorad * j2000grid[1][i]))
        # END for (i)
        f.close()

    # END write_chi2grid()

    def write_summary(self, path, loc):
        r""" Write summary text file with important info.

        Args:
            path (str): Path to output file
            loc (dict): Dictionary with information about best localization
        """
        # gather the variables we'll need
        best, geo, sun = loc["best"], loc["geo"], loc["sun"]
        loc_pointer = 1 if best["loc_reliable"] else 0
        db_no = "%5d" % self._idb_no
        locver = self.ver_no.strip() + "db" + db_no.strip()

        # write variables to text file
        f = open(path, 'w')
        f.write("%8.2f %8.2f %8.2f %3d %9.1f %7.2f %8.2f %-10s\n" %
                (best["ra"] * legacy_dtorad, best["dec"] * legacy_dtorad,
                 best["err"], loc_pointer, loc["signif"],
                 best["az"] * legacy_dtorad, best["zen"] * legacy_dtorad, locver))
        f.write("%8.2f %8.2f %8.2f\n" %
                (geo["az"] * legacy_dtorad, geo["zen"] * legacy_dtorad,
                 geo["angle"]))
        f.write("%8.2f %8.2f\n" %
                (best["lii"] * legacy_dtorad, best["bii"] * legacy_dtorad))
        f.write("%8.2f %8.2f %8.2f\n" %
                (sun["ra"] * legacy_dtorad, sun["dec"] * legacy_dtorad,
                 sun["angle"]))
        if len("%9.2f" % best["nchi2"]) <= 9:
            f.write("%9.2f %3d\n" % (best["nchi2"], best["ispec"]))
        else:
            f.write("********* %3d\n" % best["ispec"])
        f.write(self.join_float(loc["det_ang_best"], legacy_dtorad) + "\n")
        f.write(self.join_float(loc["det_ang_geo"], legacy_dtorad) + "\n")
        f.write(self.join_float(loc["deadtime"][:self.ldet]) + "\n")
        f.close()

    # END write_summary()

    def to_GbmHealPix(self, loc, frame, grid_nearest=True, **kwargs):
        r""" Converts a DoL localization dictionary to a
        :class:`~gdt.missions.fermi.gbm.localization.GbmHealPix` probability map

        Args:
            loc (dict): Dictionary containing a chi2 array returned by the
                        :meth:`~gdt.missions.fermi.gbm.localization.dol.legacy_dol.legacy_DoL.eval` 
                        method
            frame (:class:`~gdt.missions.fermi.frame.FermiFrame`):
                Frame object with spacecraft position and rotation
            grid_nearest (bool): Perform approximate nearest pixel interpolation 
                                 between chi2 grid and GbmHealPix map when True.

        Returns:
            (:class:`~gdt.missions.fermi.gbm.localization.GbmHealPix`):
                GbmHealPix probability map describing the localization
        """
        # get chi-squared values from best localization fit and spacecraft coord
        chi2 = loc["best"]["chi2"]
        az, zen = self.sky_grid

        # wrap at 2 pi to facilitate simple grid interpolation onto healpix map
        mask = (az == 0.0)
        az = np.concatenate([az, np.full(mask.sum(), 2 * np.pi)])
        zen = np.concatenate([zen, zen[mask]])
        chi2 = np.concatenate([chi2, chi2[mask]])

        # calculate correponding ra, dec with legacy method
        pos = legacy_functions.ang_to_cart_zen(az, zen)
        ra, dec = legacy_functions.sc_to_j2000(loc["scx"], loc["scy"], loc["scz"], pos)

        c2g = localization.Chi2Grid.from_data(*np.degrees([az, zen, ra, dec]), chi2)
        c2g._quaternion = frame.quaternion
        c2g._scpos = frame.obsgeoloc.xyz.value
        c2g._trigtime = frame.obstime.fermi

        return localization.GbmHealPix.from_chi2grid(c2g, grid_nearest=grid_nearest, **kwargs)

    # END to_GbmHealPix()

    @property
    def locrates(self):
        r"""(list): List containing the numpy tables for each spectral shape"""
        return self._locrates

    # END locrates()

    @locrates.setter
    def locrates(self, val):
        r""" Function to set locrates member of DoL class

        Parameters
        ----------
        val : list
            List of string paths to locrates files
        """
        self._locrates = []
        self._idb_no = None
        for path in val:
            if not os.path.exists(path):
                raise ValueError("%s does not exist" % path)
            if path[-4:] != ".npy":
                raise ValueError("%s does not have .npy extension" % path)

            table, idb_no = legacy_functions.read_table(path, None, None)

            if self._idb_no is None:
                self._idb_no = idb_no
            elif idb_no != self._idb_no:
                raise ValueError("idb_no = %d of %s does not match previous "
                                 "idb_no = %d" % (idb_no, path, self._idb_no))

            self._locrates.append(table)

    # END locrates.setter()

    def initialize_sky_grid(self, loctable):
        r""" Build grid of sky positions using grid_res resolution.
        First find out how many points will fill grid = npoints.

        Args:
            loctable (np.ndarray): Table from locrates file 
        """
        self.npoints = loctable.shape[1]
        self.sky_grid = loctable[0:2] / legacy_functions.legacy_arcmin2rad
        self.sky_grid = self.sky_grid.astype(np.float32)

    # END initialize_sky_grid()

    def stdout(self, loc):
        r""" Routine to print standard output from DoL.

        Args:
            loc (dict): Dictionary object from eval() with localization info
        """
        # pick out some variables that we'll use frequently
        best, initial = loc["best"], loc["initial"]
        geo, sun = loc["geo"], loc["sun"]

        print(" Distance between them: ")
        print(" Initial position Ra, Dec: {0:.9} {1:.9}".format(
            initial["ra"] * legacy_dtorad, initial["dec"] * legacy_dtorad))
        print(" Initial position Az, Zen: {0:.17} {1:.17}".format(
            initial["az"] * legacy_dtorad, initial["zen"] * legacy_dtorad))
        print(" This position Ra, Dec: {0:.17} {1:.17}".format(
            best["ra"] * legacy_dtorad, best["dec"] * legacy_dtorad))
        print(" This position Az, Zen: {0:.17} {1:.17}".format(
            best["az"] * legacy_dtorad, best["zen"] * legacy_dtorad))

        good_ang = legacy_functions.get_good_angle(best["xyz"], initial["pos"])
        print(" J2000: {0:.9}".format(good_ang))
        good_ang = legacy_functions.get_good_angle(best["pos"], initial["xyz"])
        print(" SC: {0:.9}".format(good_ang))
        good_ang = legacy_functions.get_good_angle(-loc["sc_pos"], initial["pos"])
        print(" J2000: dist from source to geocenter (initial) {0:.9}".format(good_ang))
        geo_ang = loc["geo"]["angle"]
        print(" J2000: dist from source to geocenter (this prog) {0:.9}".format(geo_ang))
        norm_vec = initial["pos"] / np.sqrt(legacy_functions.dot(initial["pos"], initial["pos"]))
        print(" J2000: dist from source to z-axis (initial) {0:.17}".format(
            legacy_dtorad * np.arccos(legacy_functions.dot(loc["scz"], norm_vec))))
        norm_vec = best["xyz"] / np.sqrt(legacy_functions.dot(best["xyz"], best["xyz"]))
        print(" J2000: dist from source to z-axis (this prog) {0:.17}".format(
            legacy_dtorad * np.arccos(legacy_functions.dot(loc["scz"], norm_vec))))

        # angles of detectors to earth given geo az and zen
        print(" Detector Earth angles:")
        self.display_nai_angles(loc["det_ang_geo"])

        # angles of detectors to initial source position in sc az and zen
        print(" Source angles (initial pos):")
        self.display_nai_angles(loc["det_ang_initial"])

        # angles of detectors to best source position in sc az and zen
        print(" Source angles (this prog):")
        self.display_nai_angles(loc["det_ang_best"])

        # Begin standard output
        print(" Begin Standard output (PYTHON)")
        print(" Dol Version " + self.ver_no)
        print("GRB Ra, Dec, error:  %6.2f %6.2f %6.2f" %
              (best["ra"] * legacy_dtorad, best["dec"] * legacy_dtorad, best["err"]))
        print("GRB Az, Zen (Fermi):        %6.2f %6.2f" %
              (best["az"] * legacy_dtorad, best["zen"] * legacy_dtorad))
        print("Galactic coords:            %6.2f %6.2f" %
              (best["lii"] * legacy_dtorad, best["bii"] * legacy_dtorad))
        print("Sun RA, Dec:                %6.2f %6.2f" %
              (sun["ra"] * legacy_dtorad, sun["dec"] * legacy_dtorad))
        print("Sun Distance (deg)          %6.2f" % sun["angle"])
        # print off-zenith angle
        print("Angle to local zenith:      %6.2f" %
              (np.float32(180.) - geo_ang))
        # estimate error based on Chip's prescription
        print("Intensity (cts/cm^2/sec):   %6.2f" %
              (np.float32(100.) / loc["err_chip"] ** 2))
        maxdet = loc["maxdet"]
        diff = loc["c_mrates"][maxdet] - loc["c_brates"][maxdet]
        print("Fluence (cts/cm^2):         %6.2f" %
              (diff / np.cos(loc["det_ang_best"][maxdet]) / np.float32(126.)))
        print(" End Standard output (PYTHON)")

    # END stdout()

    def display_nai_angles(self, ang):
        r""" Function to show angles of NaI detectors relative to location.

        Args:
            ang (np.ndarray): Array with detector angles
        """
        for j in range(self.nai_az.size):
            print("{0:9} {1:0<18.17} {2:0<10.9} {3:0<10.9}".format(
                j, legacy_dtorad * ang[j], self.nai_az[j], self.nai_zen[j]))

    # END display_nai_angles()

    def join_float(self, val, unit=None):
        r""" Helper function to join float values into space separated str.
            
        Args:
            val (np.ndarray): Array of floating point values
            unit (float32/64): Unit used for display string

        Returns:
            (str): String of values separated by spaces
        """
        if unit is None:
            return " ".join(["%8.4f" % i for i in val])
        return " ".join(["%8.4f" % (i * unit) for i in val])

    # END join_float()


# END class legacy_DoL()

###############################################################################

def main():
    r""" Command line interface for the legacy DoL routine """

    nen = legacy_DoL.default_nen
    ndet = legacy_DoL.default_ndet

    #############
    # ARGUMENTS #
    #############

    # need to add spaces to ensure negative numbers aren't considered options
    for i, arg in enumerate(sys.argv):
        if (arg[0] == '-') and arg[1].isdigit():
            sys.argv[i] = ' ' + arg

    p = argparse.ArgumentParser(description="Calculates localization error.",
                                formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--crange", type=np.int32, nargs=2,
                   help="Energy bin range to use for localization")
    p.add_argument("--mrates", type=np.int32, nargs=ndet * nen,
                   help="Measured counts for all detectors and energy bins")
    p.add_argument("--brates", type=np.int32, nargs=ndet * nen,
                   help="Estimated background counts for all detectors and energy bins")
    p.add_argument("--sduration", type=np.float32,
                   help="Source exposure timescale in seconds")
    p.add_argument("--bgduration", type=np.float32,
                   help="Background exposure timescale in seconds")
    p.add_argument("--sc_pos", type=np.float32, nargs=3,
                   help="Spacecraft position [x, y, x] in km to Earth center")
    p.add_argument("--sc_quat", type=np.float64, nargs=4,
                   help="Spacecraft rotation in quaternion format with scalar last")
    p.add_argument("--energies", type=np.float32, nargs=nen + 1,
                   help="Energy bin boundaries")
    p.add_argument("--fra", type=np.float32,
                   help="Initial right ascension from flight software in degrees")
    p.add_argument("--fdec", type=np.float32,
                   help="Initial declination from flight software in degrees")
    p.add_argument("--sc_time", type=np.int64,
                   help="Time of the localization in Fermi mission elapsed seconds")
    p.add_argument("--scat_opt", type=np.int32, default=1,
                   help="Apply atmospheric scattering effects when set to 1")
    p.add_argument("--fname", type=str, default="",
                   help="Output file name")
    p.add_argument("--dir", type=str, default="",
                   help="Output directory")

    args = p.parse_args()

    #############
    # ARGUMENTS #
    #############

    ###########################################################################

    fra = args.fra
    fdec = args.fdec
    odir = args.dir
    fname = args.fname
    sc_time = args.sc_time
    scat_opt = args.scat_opt
    sduration = args.sduration
    bgduration = args.bgduration
    crange = np.array(args.crange)
    sc_pos = np.array(args.sc_pos)
    sc_quat = np.array(args.sc_quat)
    energies = np.array(args.energies)

    mrates = np.zeros((ndet, nen), np.int32)
    for i in range(ndet):
        for j in range(nen):
            mrates[i][j] = args.mrates[i * nen + j]

    brates = np.zeros((ndet, nen), np.int32)
    for i in range(ndet):
        for j in range(nen):
            brates[i][j] = args.brates[i * nen + j]

    if crange[0] in [0, 1] and crange[1] == 2:
        spec = [("Hard: ", legacy_spectral_models.band_soft)]
        locrates = [legacy_spectral_models.band_soft_5_50]
    else:
        spec = None
        locrates = None

    dol = legacy_DoL(spec=spec, locrates=locrates, verbose=True)

    dol.eval(crange, mrates, brates, sduration, bgduration,
             sc_pos, sc_quat, energies, fra, fdec, sc_time,
             scat_opt, fname, odir)


# END main()

###############################################################################

if __name__ == "__main__":
    main()
