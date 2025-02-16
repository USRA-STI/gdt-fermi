# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
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
import warnings
import astropy.io.fits as fits
import numpy as np
from scipy.stats import chi2
import scipy.interpolate
import healpy as hp
from astropy.coordinates import get_sun, SkyCoord, angular_separation
from astropy.coordinates.representation import CartesianRepresentation
from astropy.units import Quantity

from gdt.core.coords import Quaternion
from gdt.core.file import FitsFileContextManager
from gdt.core.healpix import HealPixLocalization
from gdt.missions.fermi.frame import *
from gdt.missions.fermi.time import Time
from gdt.missions.fermi.gbm.detectors import GbmDetectors
from gdt.missions.fermi.gbm.headers import HealpixHeaders

__all__ = ['GbmHealPix', 'Chi2Grid', 'ga_model', 'gbuts_o3_model', 'hitl_model',
           'robo_ba_model', 'untargeted_search_model']


class GbmHealPix(HealPixLocalization, FitsFileContextManager):
    """Class for GBM HEALPix localization files.
    """

    def __init__(self):
        HealPixLocalization.__init__(self)
        FitsFileContextManager.__init__(self)
        self._frame = None
        self._geo_loc = None
        self._geo_rad = None
        self._quat = None
        self._scpos = None
        self._sun_loc = None
        self._trigtime = None

    def convolve(self, model, *args, **kwargs):
        """Convolve the map with a model kernel.  The model can be a Gaussian
        kernel or any mixture of Gaussian kernels. Uses `healpy.smoothing 
        <https://healpy.readthedocs.io/en/latest/generated/healpy.sphtfunc.smoothing.html>`_.
        
        An example of a model kernel with a 50%/50% mixture of two Gaussians,
        one with a 1-deg width, and the other with a 3-deg width::
            
            def gauss_mix_example():
                sigma1 = np.deg2rad(1.0)
                sigma2 = np.deg2rad(3.0)
                frac1 = 0.50
                return ([sigma1, sigma2], [frac1])
        
        Args: 
            model (<function>): The function representing the model kernel
            *args: Arguments to be passed to the model kernel function
        
        Returns:
            (:class:`GbmHealPix`)
        """
        return super().convolve(model, *args, headers=self.headers, 
                                quaternion=self.quaternion, scpos=self.scpos, 
                                **kwargs)
        
    @property
    def frame(self):
        """(:class:`~gdt.core.coords.SpacecraftFrame`): The spacecraft frame at
        the time of the localization"""
        return self._frame

    @property
    def geo_location(self):
        """(astropy.coordinates.SkyCoord): The geocenter location at
        :attr:`trigtime`"""
        return self._geo_loc

    @property
    def geo_probability(self):
        """(float): The amount of localization probability on the Earth"""
        if self.geo_location is None:
            return None
        prob_mask, geo_mask = self._earth_mask()
        return self.prob[prob_mask][geo_mask].sum()

    @property
    def geo_radius(self):
        """(astropy.units.Quantity): The apparent angular radius of the Earth at
        :attr:`trigtime`.

        Note:
            If a :attr:`scpos` isn't set, then an average 67.5 deg is returned
        """
        # if the radius isn't known, use the average 67.5 deg radius
        if self._geo_rad is not None:
            return self._geo_rad
        else:
            return Quantity(67.5, unit='deg')

    @property
    def quaternion(self):
        """(:class:`~gdt.core.coords.Quaternion`): The spacecraft attitude
        quaternion"""
        return self._quat

    @property
    def scpos(self):
        """(astropy.coordinates.CartesianRepresentation):
           The spacecraft position in Earth inertial coordinates"""
        return self._scpos

    @property
    def sun_location(self):
        """(astropy.coordinates.SkyCoord): The Sun location at
        :attr:`trigtime`"""
        return self._sun_loc

    @classmethod
    def from_chi2grid(cls, chi2grid, nside=128, headers=None, filename=None, exact_nearest=False, grid_nearest=False):
        """Create a GbmHealPix object from a :class:`Chi2Grid` object.

        Args:
            chi2grid (:class:`Chi2Grid`): The chi2grid object containing the
                                         chi-squared/log-likelihood info.
            nside (int, optional): The nside resolution to use. Default is 128.
            headers (:class:`~gdt.core.headers.FileHeaders`, optional):
                The file headers
            filename (str, optional): The filename
            exact_nearest (bool): Use exact nearest pixel interpolation when True.
                                  This method is slow O(minute) due to
                                  angular difference calculation.
            grid_nearest (bool): Use approximate nearest pixel interpolation when True.
                                 This method is fast O(second) by using 2D grid.

        Returns:
            (:class:`GbmHealPix`)
        """
        if not isinstance(chi2grid, Chi2Grid):
            raise TypeError('chi2grid must be a Chi2Grid object')

        # convert chisq map to probability map assuming Wilk's theorem applies
        loglike = -chi2grid.chisq / 2.0
        probs = np.exp(loglike - np.max(loglike))

        # fill a low-resolution healpix map with probability values
        lores_nside = 64
        lores_npix = hp.nside2npix(lores_nside)
        lores_array = np.zeros(lores_npix)
        theta = cls._dec_to_theta(chi2grid.dec)
        phi = cls._ra_to_phi(chi2grid.ra)
        if exact_nearest:
            # exact nearest pixel method using angular difference (slow)
            lores_pix = np.arange(lores_npix)
            proj_theta, proj_phi = hp.pix2ang(lores_nside, lores_pix)
            idx = [angular_separation(proj_phi[i], 0.5 * np.pi - proj_theta[i],
                                      phi, 0.5 * np.pi - theta).argmin() for i in lores_pix]
            lores_array = probs[idx]
        elif grid_nearest:
            # approximate nearest pixel method using 2D grid (fast)
            lores_pix = np.arange(lores_npix)
            proj_theta, proj_phi = hp.pix2ang(lores_nside, lores_pix)
            lores_array = scipy.interpolate.griddata((phi, theta), probs, (proj_phi, proj_theta), method='nearest')
        else:
            # basic healpix conversion (can result in missing pixels)
            idx = hp.ang2pix(lores_nside, theta, phi)
            lores_array[idx] = probs

        # upscale to high-resolution
        hires_nside = nside
        hires_npix = hp.nside2npix(hires_nside)
        theta, phi = hp.pix2ang(hires_nside, np.arange(hires_npix))
        prob_array = hp.get_interp_val(lores_array, theta, phi)

        quat = Quaternion(chi2grid.quaternion)
        scpos = CartesianRepresentation(chi2grid.scpos, unit='m')

        obj = cls.from_data(prob_array, trigtime=chi2grid.trigtime,
                            headers=headers, filename=filename,
                            scpos=scpos, quaternion=quat)
        return obj

    @classmethod
    def from_data(cls, prob_arr, trigtime=None, headers=None, filename=None,
                  quaternion=None, scpos=None):
        """Create a GbmHealPix object from a healpix array.

        Args:
            prob_arr (np.array):
                The HEALPix array containing the probability/pixel
            trigtime (float, optional):
                The time corresponding to the localization
            headers (:class:`~gdt.core.headers.FileHeaders`, optional):
                The file headers
            filename (str, optional): The filename
            quaternion (:class:`~gdt.core.coords.Quaternion`, optional):
                The associated spacecraft quaternion used to determine the
                detector pointings in equatorial coordinates
            scpos (astropy.coordinates.representation.CartesianRepresentation, optional):
                The associated spacecraft position in Earth inertial coordinates
                used to determine the geocenter location in equatorial
                coordinates

        Returns:
            (:class:`GbmHealPix`)
        """
        obj = super().from_data(prob_arr, trigtime=trigtime, filename=filename)

        if headers is not None:
            if not isinstance(headers, HealpixHeaders):
                raise TypeError('headers must be a HealpixHeaders object')
        else:
            headers = cls._none_default_headers()
        obj._headers = headers

        if quaternion is not None:
            if not isinstance(quaternion, Quaternion):
                raise TypeError('quaternion must be a Quaternion object')
            obj._quat = quaternion

        if scpos is not None:
            if not isinstance(scpos, CartesianRepresentation):
                raise TypeError('scpos must be a CartesianRepresentation object')
            obj._scpos = scpos
            
        # if we have a trigtime, calculate sun position
        if trigtime is not None:
            trigtime = Time(trigtime, format='fermi')
            obj._sun_loc = get_sun(trigtime)
        elif obj._headers[1]['SUN_RA'] is not None:
            obj._sun_loc = SkyCoord(obj._headers[1]['SUN_RA'],
                                    obj._headers[1]['SUN_DEC'], unit='deg',
                                    frame='gcrs')
        else:
            obj._sun_loc = None
        
        # create the spacecraft frame with the info that we have at hand
        obj._frame = FermiFrame(obstime=trigtime,
                                quaternion=obj._quat, obsgeoloc=obj._scpos,
                                detectors=GbmDetectors)
        
        # if have scpos, create geocenter location and earth radius
        # if not, then try to pull from header
        if obj._scpos is not None:
            obj._geo_loc = obj._frame.geocenter
            obj._geo_rad = obj._frame.earth_angular_radius
        elif obj._headers[1]['GEO_RA'] is not None:
            obj._geo_loc = SkyCoord(obj._headers[1]['GEO_RA'],
                                    obj._headers[1]['GEO_DEC'], unit='deg',
                                    frame='gcrs')
            obj._geo_rad = Quantity(obj._headers[1]['GEO_RAD'], unit='deg')

        # if have trigtime and quaternion, create detector pointings
        # if not, then try to pull from header
        if (trigtime is not None) and (quaternion is not None):
            
            for det in obj._frame.detectors:
                pointing = (det.azimuth, det.elevation)
                det_coord = SkyCoord(*pointing, frame=obj._frame).gcrs[0]
                setattr(obj, det.name.lower() + '_pointing', det_coord)
        elif obj._headers[1]['N0_RA'] is not None:

            for det in GbmDetectors:
                ra_key = det.name.upper() + '_RA'
                dec_key = det.name.upper() + '_DEC'
                det_coord = SkyCoord(obj._headers[1][ra_key],
                                     obj._headers[1][dec_key], unit='deg',
                                     frame='gcrs')
                setattr(obj, det.name.lower() + '_pointing', det_coord)
        
        # build headers
        obj._headers = obj._build_headers(obj.trigtime, obj.nside)

        return obj

    @classmethod
    def multiply(cls, healpix1, healpix2, primary=0, output_nside=128):
        """Multiply two GbmHealPix maps and return a new map.

        Note:
            Either ``healpix1`` *or* ``healpix2`` can be a non-GbmHealPix
            object, however at least one of them must be a GbmHealPix object
            **and** the ``primary`` argument must be set to the appropriate
            GbmHealPix object otherwise a TypeError will be raised.

        Args:
            healpix1 (:class:`~.gdt.core.healpix.HealPix` or :class:`GbmHealPix`):
                One of the HEALPix maps to multiply
            healpix2 (:class:`~.gdt.core.healpix.HealPix` or :class:`GbmHealPix`):
                The other HEALPix map to multiply
            primary (int, optional): If 0, use the first map header information,
                                     or if 1, use the second map header
                                     information. Default is 1.
            output_nside (int, optional): The nside of the multiplied map.
                                          Default is 128.
        Returns
            :class:`GbmHealPix`: The multiplied map
        """

        if primary == 0:
            if not isinstance(healpix1, cls):
                raise TypeError('Primary HealPix (healpix1) is not of class {}. '
                                'Perhaps try setting healpix2 as the primary'.format(cls.__name__))
        else:
            if not isinstance(healpix2, cls):
                raise TypeError('Primary HealPix (healpix2) is not of class {}. '
                                'Perhaps try setting healpix1 as the primary'.format(cls.__name__))

        if primary == 0:
            headers = healpix1.headers
            quat = healpix1.quaternion
            scpos = healpix1.scpos
        else:
            headers = healpix2.headers
            quat = healpix2.quaternion
            scpos = healpix2.scpos

        obj = super().multiply(healpix1, healpix2, primary=primary,
                               output_nside=output_nside, quaternion=quat,
                               scpos=scpos)

        return obj

    def observable_fraction(self, healpix):
        """The observable fraction of a healpix probability region on the sky.
        Non-observable regions are ones that are behind the Earth.

        Args:
            healpix (:class:`HealPix`): The healpix region for which to
                                        calculate the observable fraction.
        Returns:
            (float)
        """
        if self.geo_location is None:
            raise RuntimeError('Location of geocenter is not known')

        # speed things up a bit by only considering pixels with non-zero prob
        prob_mask = (healpix.prob > 0.0)
        # get ra, dec coords for pixels and calculate angle from geocenter
        theta, phi = hp.pix2ang(healpix.nside, np.arange(healpix.npix))
        ra = self._phi_to_ra(phi)[prob_mask]
        dec = self._theta_to_dec(theta)[prob_mask]
        pts = SkyCoord(ra, dec, frame='icrs', unit='deg')
        # the mask of everything with prob > 0.0 and is visible
        ang = self.geo_location.separation(pts)
        geo_mask = (ang > self.geo_radius)

        # sum it up and divide by total prob (should be 1, but good to be sure)
        temp = np.copy(healpix.prob)
        temp = temp[prob_mask]
        frac = temp[geo_mask].sum() / healpix.prob.sum()
        return frac

    @classmethod
    def open(cls, file_path, **kwargs):
        """Open a GBM HEALPix FITS file and return the GbmHealPix object

        Args:
            file_path (str): The file path of the FITS file

        Returns:
            (:class:`GbmHealPix`)
        """
        # ignore comment length warnings
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=UserWarning)
            obj = super().open(file_path, **kwargs)

        # get the headers
        hdrs = [hdu.header for hdu in obj.hdulist]

        # some older files do not have these keywords
        if 'GEO_RAD' not in hdrs[1]:
            hdrs[1]['GEO_RAD'] = 67.5
        if 'COMMENT' not in hdrs[1]:
            hdrs[1]['COMMENT'] = ''
            hdrs[1]['COMMENT'] = ''
        headers = HealpixHeaders.from_headers(hdrs)

        trigtime = headers['PRIMARY']['TRIGTIME']

        # quaternion and scpos are stored as comments
        try:
            headers[1]['COMMENT'][0] = obj.hdulist[1].header['COMMENT'][0]
            headers[1]['COMMENT'][1] = obj.hdulist[1].header['COMMENT'][1]
        except:
            headers[1]['COMMENT'][0] = ''
            headers[1]['COMMENT'][1] = ''

        try:
            scpos_comment = headers[1]['COMMENT'][0]
            scpos = scpos_comment.split('[')[1].split(']')[0]
            scpos = np.array([float(el) for el in scpos.split()])
            scpos = CartesianRepresentation(scpos, unit='m')
        except:
            scpos = None
        try:
            quat_comment = headers[1]['COMMENT'][1]
            quat = quat_comment.split('[')[1].split(']')[0]
            quat = np.array([float(el) for el in quat.split()])
            quat = Quaternion(quat)
        except:
            quat = None

        # get the probability and significance arrays
        prob = obj.column(1, 'PROBABILITY').flatten()
        sig = obj.column(1, 'SIGNIFICANCE').flatten()
        if headers[1]['ORDERING'] == 'NESTED':
            npix = prob.size
            idx = hp.ring2nest(hp.npix2nside(npix), np.arange(npix))
            prob = prob[idx]
            sig = sig[idx]

        obj.close()

        obj = cls.from_data(prob, trigtime=trigtime, quaternion=quat,
                            scpos=scpos, filename=obj.filename,
                            headers=headers)
        obj._sig = sig

        return obj

    def region_probability(self, healpix, prior=0.5):
        r"""The probability that the localization is associated with
        the localization region from another map.  This is calculated
        against the null hypothesis that the two maps represent
        unassociated sources:

        :math:`P(A | \mathcal{I}) =
        \frac{P(\mathcal{I} | A) \ P(A)}
        {P(\mathcal{I} | A) \ P(A) + P(\mathcal{I} | \neg A) \ P(\neg A)}`

        where

        * :math:`P(\mathcal{I} | A)` is the integral over the overlap of the two
          maps once the Earth occultation has been removed for *this* map.
        * :math:`P(\mathcal{I} | \neg A)` is the integral over the overlap of
          *this* map with a uniform distribution on the sky (i.e. the probability
          the localization is associated with a random point on the sky)
        * :math:`P(A)` is the prior probability that *this* localization is
          associated with the *other* HEALPix map.

        Note:
            The localization region of *this* map overlapping the Earth will be
            removed and the remaining unocculted region is used for the
            calculation.  The *other* map is assumed to have no exclusionary
            region.

        Args:
            healpix (:class:`~gdt.core.healpix.HealPixLocalization`):
                The healpix map for which to calculate the spatial association.
            prior (float, optional): The prior probability that the localization
                                     is associated with the source. Default is
                                     0.5.

        Returns:
            (float)
        """
        if (prior < 0.0) or (prior > 1.0):
            raise ValueError('Prior probability must be within 0-1, inclusive')

        # convert uniform prob/sr to prob/pixel
        u = 1.0 / (4.0 * np.pi)

        # get the non-zero probability and earth masks
        try:
            prob_mask, geo_mask = self._earth_mask()
        except:
            prob_mask = np.ones_like(self.prob, dtype=bool)
            geo_mask = np.zeros_like(self.prob, dtype=bool)
        probmap1 = np.copy(self.prob)
        temp = probmap1[prob_mask]
        temp[geo_mask] = 0.0
        probmap1[prob_mask] = temp
        probmap1 = self._assert_prob(probmap1)

        # ensure maps are the same resolution and convert uniform prob/sr to
        # prob/pixel
        probmap2 = np.copy(healpix.prob)
        if self.nside > healpix.nside:
            probmap2 = hp.ud_grade(probmap2, nside_out=self.nside)
            probmap2 = self._assert_prob(probmap2)
            u *= hp.nside2resol(self.nside) ** 2
        elif self.nside < healpix.nside:
            probmap1 = hp.ud_grade(probmap1, nside_out=healpix.nside)
            probmap1 = self._assert_prob(probmap1)
            u *= hp.nside2resol(healpix.nside) ** 2
        else:
            u *= hp.nside2resol(self.nside) ** 2

        # alternative hypothesis: they are related
        alt_hyp = np.sum(probmap1 * probmap2)
        # null hypothesis: one of the maps is from an unassociated source
        # (uniform spatial probability)
        null_hyp = np.sum(probmap1 * u)

        # since we have an exhaustive and complete list of possibilities, we can
        # easily calculate the probability
        prob = (alt_hyp * prior) / ((alt_hyp * prior) + (null_hyp * (1.0 - prior)))
        return prob

    def remove_earth(self):
        """Return a new GbmHealPix with the probability on the Earth masked out.
        The remaining probability on the sky is renormalized.

        Note:
            The :attr:`geo_location` attribute must be available to use this
            function.


        Returns:
            (:class:`GbmHealPix`)
        """
        if self.geo_location is None:
            raise RuntimeError('Location of geocenter is not known')

        # get the non-zero probability and earth masks
        prob_mask, geo_mask = self._earth_mask()

        # zero out the probabilities behind the earth
        new_prob = np.copy(self.prob)
        temp = new_prob[prob_mask]
        temp[geo_mask] = 0.0
        new_prob[prob_mask] = temp
        # renormalize
        new_prob = self._assert_prob(new_prob)
        obj = type(self).from_data(new_prob, trigtime=self.trigtime,
                                   headers=self.headers, scpos=self.scpos,
                                   filename=self.filename,
                                   quaternion=self.quaternion)
        return obj

    def source_probability(self, ra, dec, prior=0.5):
        r"""The probability that the GbmHealPix localization is associated with
        a known point location.  This is calculated against the null hypothesis
        that the localization originates from an unassociated random source
        that has equal probability of origination anywhere in the sky:

        :math:`P(A | \mathcal{I}) =
        \frac{P(\mathcal{I} | A) \ P(A)}
        {P(\mathcal{I} | A) \ P(A) + P(\mathcal{I} | \neg A) \ P(\neg A)}`

        where

        * :math:`P(\mathcal{I} | A)` is the probability of the localization at
          the point source once the Earth occultation has been removed
        * :math:`P(\mathcal{I} | \neg A)` is the probability per pixel assuming
          a uniform distribution on the sky (i.e. the probability the
          localization is associated with a random point on the sky)
        * :math:`P(A)` is the prior probability that the localization is
          associated with the point source

        Note:
            If the point source is behind the Earth, then it is assumed that
            GBM could not observe it, therefore the probability will be zero.

        Args:
            ra (float): The RA of the known source location
            dec (float): The Dec of the known source location
            prior (float, optional): The prior probability that the localization
                                     is associated with the source. Default is
                                     0.5.

        Returns:
            (float)
        """
        if (prior < 0.0) or (prior > 1.0):
            raise ValueError('Prior probability must be within 0-1, inclusive')

        # convert uniform prob/sr to prob/pixel
        u = 1.0 / (4.0 * np.pi)
        u *= hp.nside2resol(self.nside) ** 2

        # the pixel probability of the skymap at the location of the point source
        try:
            p = self.remove_earth().probability(ra, dec, per_pixel=True)
        except:
            p = self.probability(ra, dec, per_pixel=True)

        # if we know the location of the earth and it's behind the earth,
        # then we obviously couldn't have seen it
        if self.geo_location is not None:
            pt = SkyCoord(ra, dec, frame='icrs', unit='deg')
            ang = self.geo_location.separation(pt)
            if ang < self.geo_radius:
                p = 0.0

        # null hypothesis is that they are not associated, therefore the sky map
        # is result of some source that has uniform probability on the sky
        prob = (p * prior) / ((p * prior) + (u * (1.0 - prior)))
        return prob

    def _build_hdulist(self):

        # create FITS and primary header
        hdulist = fits.HDUList()
        primary_hdu = fits.PrimaryHDU(header=self.headers['PRIMARY'])
        for key, val in self.headers['PRIMARY'].items():
            primary_hdu.header[key] = val
        hdulist.append(primary_hdu)

        # the healpix extension
        prob_arr = hp.reorder(self.prob, r2n=True).reshape(-1, 1024)
        prob_col = fits.Column(name='PROBABILITY', format='1024E',
                               array=prob_arr)
        sig_arr = hp.reorder(self.sig, r2n=True).reshape(-1, 1024)
        sig_col = fits.Column(name='SIGNIFICANCE', format='1024E',
                              array=sig_arr)
        hpx_hdu = fits.BinTableHDU.from_columns([prob_col, sig_col],
                                                header=self.headers['HEALPIX'])
        cflag = 0
        for key, val in self.headers['HEALPIX'].items():
            if key == 'COMMENT':
                hpx_hdu.header['COMMENT'][cflag] = val
                cflag += 1
                continue
            hpx_hdu.header[key] = val
        hdulist.append(hpx_hdu)
        hpx_hdu.header.comments['TTYPE1'] = 'Differential probability per pixel'
        hpx_hdu.header.comments['TTYPE2'] = 'Integrated probability'

        return hdulist

    def _build_headers(self, trigtime, nside):

        headers = self.headers.copy()
        headers['PRIMARY']['TRIGTIME'] = trigtime
        headers['HEALPIX']['NSIDE'] = nside
        headers['HEALPIX']['LASTPIX'] = hp.nside2npix(nside)
        headers['HEALPIX']['OBJECT'] = 'FULLSKY'

        if self.trigtime is not None:
            headers['HEALPIX']['SUN_RA'] = self.sun_location.ra.value
            headers['HEALPIX']['SUN_DEC'] = self.sun_location.dec.value

        if self.scpos is not None:
            headers['HEALPIX']['COMMENT'][0] = 'SCPOS: ' + \
                                               np.array2string(self.scpos.xyz.value)
            headers['HEALPIX']['GEO_RA'] = self.geo_location.ra.value
            headers['HEALPIX']['GEO_DEC'] = self.geo_location.dec.value
            headers['HEALPIX']['GEO_RAD'] = self.geo_radius.value

        if self.quaternion is not None:
            quat = np.append(self.quaternion.xyz, self.quaternion.w)
            headers['HEALPIX']['COMMENT'][1] = 'QUAT: ' + np.array2string(quat)
            for det in GbmDetectors:
                pointing = getattr(self, det.name.lower() + '_pointing')
                headers['HEALPIX'][det.name.upper() + '_RA'] = pointing.ra.value
                headers['HEALPIX'][det.name.upper() + '_DEC'] = pointing.dec.value

        return headers

    def _earth_mask(self):
        # speed things up a bit by only considering pixels with non-zero prob
        mask = (self.prob > 0.0)
        # get ra, dec coords for pixels and calculate angle from geocenter
        theta, phi = hp.pix2ang(self.nside, np.arange(self.npix))
        pts = SkyCoord(self._phi_to_ra(phi)[mask],
                       self._theta_to_dec(theta)[mask], frame='icrs',
                       unit='deg')
        ang = self.geo_location.separation(pts)

        # the mask of the non-zero probability pixels that are behind the earth
        geo_mask = (ang <= self.geo_radius)

        return mask, geo_mask

    @staticmethod
    def _none_default_headers():
        hdr = HealpixHeaders()
        hdr[0]['TRIGTIME'] = None
        hdr[1]['SUN_RA'] = None
        hdr[1]['SUN_DEC'] = None
        hdr[1]['GEO_RA'] = None
        hdr[1]['GEO_DEC'] = None
        hdr[1]['GEO_RAD'] = None
        hdr[1]['N0_RA'] = None
        hdr[1]['N0_DEC'] = None
        hdr[1]['N1_RA'] = None
        hdr[1]['N1_DEC'] = None
        hdr[1]['N2_RA'] = None
        hdr[1]['N2_DEC'] = None
        hdr[1]['N3_RA'] = None
        hdr[1]['N3_DEC'] = None
        hdr[1]['N4_RA'] = None
        hdr[1]['N4_DEC'] = None
        hdr[1]['N5_RA'] = None
        hdr[1]['N5_DEC'] = None
        hdr[1]['N6_RA'] = None
        hdr[1]['N6_DEC'] = None
        hdr[1]['N7_RA'] = None
        hdr[1]['N7_DEC'] = None
        hdr[1]['N8_RA'] = None
        hdr[1]['N8_DEC'] = None
        hdr[1]['N9_RA'] = None
        hdr[1]['N9_DEC'] = None
        hdr[1]['NA_RA'] = None
        hdr[1]['NA_DEC'] = None
        hdr[1]['NB_RA'] = None
        hdr[1]['NB_DEC'] = None
        hdr[1]['B0_RA'] = None
        hdr[1]['B0_DEC'] = None
        hdr[1]['B1_RA'] = None
        hdr[1]['B1_DEC'] = None

        return hdr

    def __repr__(self):
        s = '<{0}: {1}\n'.format(self.__class__.__name__, self.filename)
        s += ' NSIDE={0}; trigtime={1};\n'.format(self.nside, self.trigtime)
        s += ' centroid={}>'.format(self.centroid)
        return s


class Chi2Grid():
    """Class for the GBM legacy internal Chi2Grid localization files/objects.
    """

    def __init__(self):
        self._az = np.array([])
        self._zen = np.array([])
        self._ra = np.array([])
        self._dec = np.array([])
        self._chisq = np.array([])
        self._quaternion = None
        self._scpos = None
        self._trigtime = None

    @property
    def azimuth(self):
        """(np.array): The spacecraft azimuth grid points"""
        return self._az

    @property
    def chisq(self):
        """(np.array): The chi-squared value at each grid point"""
        return self._chisq

    @property
    def dec(self):
        """(np.array): The Dec grid points"""
        return self._dec

    @property
    def numpts(self):
        """(int): Number of sky points in the Chi2Grid"""
        return self._az.size

    @property
    def quaternion(self):
        """(np.array): The spacecraft attitude quaternion"""
        return self._quaternion

    @quaternion.setter
    def quaternion(self, val):
        if len(val) != 4:
            raise ValueError('quaternion must be a 4-element array')
        self._quaternion = np.asarray(val)

    @property
    def ra(self):
        """(np.array): The RA grid points"""
        return self._ra

    @property
    def scpos(self):
        """(np.array): The spacecraft position in Earth inertial coordinates"""
        return self._scpos

    @scpos.setter
    def scpos(self, val):
        if len(val) != 3:
            raise ValueError('scpos must be a 3-element array')
        self._scpos = np.asarray(val)

    @property
    def significance(self):
        """(np.array): The significance value at each point"""
        min_chisq = np.min(self.chisq)
        return 1.0 - chi2.cdf(self.chisq - min_chisq, 2)

    @property
    def trigtime(self):
        """(float): The trigger time"""
        return self._trigtime

    @trigtime.setter
    def trigtime(self, val):
        try:
            val = float(val)
        except:
            raise ValueError('trigtime must be a float')
        self._trigtime = val

    @property
    def zenith(self):
        """(np.array): The spacecraft zenith grid points"""
        return self._zen

    @classmethod
    def open(cls, filename):
        """Read a chi2grid file and create a Chi2Grid object

        Args:
            filename (str): The filename of the chi2grid file

        Returns:
           :class:`Chi2Grid`: The Chi2Grid object
        """
        with open(filename, 'r') as f:
            txt = list(f)

        obj = cls()

        numpts = int(txt[0].strip())
        txt = txt[1:]
        obj._az = np.empty(numpts)
        obj._zen = np.empty(numpts)
        obj._ra = np.empty(numpts)
        obj._dec = np.empty(numpts)
        obj._chisq = np.empty(numpts)
        for i in range(numpts):
            line = txt[i].split()
            obj._az[i] = float(line[0].strip())
            obj._zen[i] = float(line[1].strip())
            obj._chisq[i] = float(line[2].strip())
            obj._ra[i] = float(line[4].strip())
            obj._dec[i] = float(line[5].strip())

        return obj

    @classmethod
    def from_data(cls, az, zen, ra, dec, chisq):
        """Create a Chi2Grid object from arrays

        Args:
            az (np.array): The azimuth grid points
            zen (np.array): The zenith grid points
            ra (np.array): The RA grid points
            dec (np.array): The Dec grid points
            chisq (np.array): The chi-squared values at each grid point

        Returns:
            :class:`Chi2Grid`: The Chi2Grid object
        """
        az = np.asarray(az, dtype=float).flatten()
        zen = np.asarray(zen, dtype=float).flatten()
        ra = np.asarray(ra, dtype=float).flatten()
        dec = np.asarray(dec, dtype=float).flatten()
        chisq = np.asarray(chisq, dtype=float).flatten()

        numpts = az.size
        if (zen.size != numpts) or (ra.size != numpts) or (dec.size != numpts) \
                or (chisq.size != numpts):
            raise ValueError('All inputs must have same size')

        obj = cls()
        obj._az = az
        obj._zen = zen
        obj._ra = ra
        obj._dec = dec
        obj._chisq = chisq
        return obj


# Systematic Model definitions using healpy.smoothing
# --------------------------------------------------------
def ga_model():
    """The localization systematic model for the Ground-Automated localization:
    A mixture of a 3.72 deg Gaussian (80.4% weight) and a 13.7 deg Gaussian.

    References:
        `Connaughton, V. et al. 2015, ApJ, 216, 32
        <https://iopscience.iop.org/article/10.1088/0067-0049/216/2/32>`_
    """
    sigma1 = np.deg2rad(3.72)
    sigma2 = np.deg2rad(13.7)
    frac1 = 0.804
    return ([sigma1, sigma2], [frac1])


def gbuts_o3_model():
    """The localization systematic model for the targeted search during O3:
    a 2.7 deg Gaussian.

    References:
        `Goldstein, A. et al. 2019, arXiv: 1903.12597
        <https://arxiv.org/abs/1903.12597>`_
    """
    sigma = np.deg2rad(2.7)
    return ([sigma], [1.0])


def hitl_model(az):
    """The localization systematic model for the human-in-the loop localization:
    A mixture of a 4.17 deg Gaussian (91.8% weight) and a 15.3 deg Gaussian
    for a centroid between azimuth 292.5 - 67.5 or azimuth 112.5 - 247.5,
    otherwise a mixture of a 2.31 deg Gaussian (88.4% weight) and a
    13.2 deg Gaussian.

    References:
        `Connaughton, V. et al. 2015, ApJ, 216, 32
        <https://iopscience.iop.org/article/10.1088/0067-0049/216/2/32>`_

    Args:
        az (float): The localization centroid in spacecraft azimuth
    """
    if (az > 292.5) or (az <= 67.5) or ((az > 112.5) and (az < 247.5)):
        sigma1 = np.deg2rad(4.17)
        sigma2 = np.deg2rad(15.3)
        frac1 = 0.918
    else:
        sigma1 = np.deg2rad(2.31)
        sigma2 = np.deg2rad(13.2)
        frac1 = 0.884
    return ([sigma1, sigma2], [frac1])


def robo_ba_model(grb_type):
    """The localization systematic model for the RoboBA localization:
    A mixture of a 1.86 deg Gaussian (57.9% weight) and a 4.14 deg Gaussian
    for a "long" GRB, and a mixture of a 2.55 deg Gaussian (39.0% weight) and a
    4.43 deg Gaussian for a "short" GRB.

    References:
        `Goldstein, A. et al. 2020, ApJ, 895, 40
        <https://iopscience.iop.org/article/10.3847/1538-4357/ab8bdb>`_

    Args:
        grb_type (str): The type of GRB, either 'long' or 'short'
    """
    if grb_type == 'long':
        sigma1 = np.deg2rad(1.86)
        sigma2 = np.deg2rad(4.14)
        frac1 = 0.579
    elif grb_type == 'short':
        sigma1 = np.deg2rad(2.55)
        sigma2 = np.deg2rad(4.43)
        frac1 = 0.39
    else:
        raise ValueError("grb_type must either be 'long' or 'short'")
    return ([sigma1, sigma2], [frac1])


def untargeted_search_model():
    """The localization systematic model for the Untargeted Search:
    A 5.53 deg Gaussian
    """
    sigma = np.deg2rad(5.53)
    return ([sigma], [1.0])
