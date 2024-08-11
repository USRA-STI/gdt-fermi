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

import numpy as np

from . import __data_dir__
from . import legacy_spectral_models

#############
# CONSTANTS #
#############

legacy_rpi = np.float64(np.arccos(-1.))
legacy_twopi = np.float32(6.2831853)
legacy_pi = np.float64(np.float32(3.14159265358973))
legacy_dtorad = np.float64(180. / np.arccos(-1.))
legacy_deg_to_rad = np.float64(legacy_pi / np.float32(180.))
legacy_arcmin2rad = np.float32(60.) * legacy_dtorad

legacy_tenergies = np.array([
    9.88152, 21.9039, 30.6248, 39.6809, 53.7016,
    72.9625, 97.6607, 122.844, 163.847, 231.738,
    316.693, 424.295, 587.606, 741.422, 1096.22,
    1813.54, 2749.15]).astype(np.float32)


#############
# CONSTANTS #
#############

def arctan2(x, y):
    """ Compute arctan2. Account for special case of small x,y

    Parameters
    ----------
    x : float32, np.ndarray(float32)
        X coordinate(s)
    y : float32, np.ndarray(float32)
        Y coordinate(s)

    Returns
    -------
    a : float32, np.ndarray(float32)
        arctan2(s) of x and y
    """
    if isinstance(x, np.ndarray) or isinstance(y, np.ndarray):
        return_single = False
    else:
        return_single = True

    x = np.atleast_1d(x)
    y = np.atleast_1d(y)

    min_val = np.float32(1e-6)
    mask = (np.fabs(x) >= min_val) | (np.fabs(y) >= min_val)

    a = np.zeros(x.size, np.float32)
    a[mask] = np.float32(np.arctan2(y[mask], x[mask]))
    a[a < np.float32(0.0)] += legacy_twopi

    return a[0] if return_single else a


# END arctan2()

def choose_energy_data(ndet, nen, crange, mrates, brates, energies):
    """ Function to choose energy bins

    Parameters
    ----------
    ndet : int
        Number of detectors
    nen : int
        Number of energy bins
    crange : np.ndarray(2, int32)
        Array of energy channel IDs chosen for localization
    mrates : np.ndarray(ndet * nen, int32)
        Array of measured detector counts in each energy channel
    brates : np.ndarray(ndet * nen, int32)
        Array of estimated background counts in each energy channel
    energies : np.ndarray(nen + 1, float32)

    Returns
    -------
    c_mrates : np.ndarray(ndet, int32)
        Measured detector counts summed over chosen energy channels
    c_brates : np.ndarray(ndet, int32)
        Estimated background counts summed over chosen energy channels
    cenergies : np.ndarray(2, float32)
        Array with first and last index of chosen energy channel range
    """
    tmrates = np.zeros((ndet, nen), np.int32)
    tbrates = np.zeros((ndet, nen), np.int32)

    c_mrates = np.zeros(ndet, np.int32)
    c_brates = np.zeros(ndet, np.int32)
    cenergies = np.zeros(2, np.float32)

    # select only energy bins within crange
    for i in range(ndet):
        for j in range(nen):
            if (j < crange[0]) or (j > crange[1]):
                continue
            tmrates[i][j] = mrates[i][j]
            tbrates[i][j] = brates[i][j]

    for i in range(ndet):
        c_mrates[i] = tmrates[i].sum()
        c_brates[i] = tbrates[i].sum()

    cenergies[0] = energies[crange[0]]
    cenergies[1] = energies[crange[1] + 1]

    return c_mrates, c_brates, cenergies


# END choose_energy_data()

def get_geocenter(sc_quat, sc_pos, verbose=False):
    """ Get Earth center and spacecraft direction cosine matrix.

    Transform from "Compendium of Co-ordinate Transformations"
    B Shivakumar et al. National Aerospace Laboratories.
    Rewritten 05/23 using MSB's matrix (seems different)

    Tested sc cosines using quaternion & x,z cosines provided by
    T. Burnett (toby_quat.png) 05/07.  Input Toby's quaternion
    to fakerates_q.txt, return correct scx, scz (scy untested).

    Parameters
    ----------
    sc_quat : np.ndarray(4, float64)
        Quaternion of the space craft. Last element is the scalar field.
    sc_pos : np.ndarray(3, float64)
        Space craft position
    verbose : bool
        Print info to screen when True

    Returns
    -------
    geodir : np.ndarray(3, np.float64)
        Vector pointing from spacecraft to center of the Earth
    geo_az : float32
        Azimuth of geodir relative to spacecraft
    geo_zen : float32
        Zenith of geodir relative to spacecraft
    scx : np.ndarray(3, np.float64)
        Vector along +x axis of spacecraft
    scy : np.ndarray(3, np.float64)
        Vector along +y axis of spacecraft
    scz : np.ndarray(3, np.float64)
        Vector along +z axis of spacecraft
    """
    scx = np.zeros(3, np.float64)
    scx[0] = (sc_quat[0] ** 2 - sc_quat[1] ** 2 - sc_quat[2] ** 2 + sc_quat[3] ** 2)
    scx[1] = 2. * (sc_quat[0] * sc_quat[1] + sc_quat[3] * sc_quat[2])
    scx[2] = 2. * (sc_quat[0] * sc_quat[2] - sc_quat[3] * sc_quat[1])

    scy = np.zeros(3, np.float64)
    scy[0] = 2. * (sc_quat[0] * sc_quat[1] - sc_quat[3] * sc_quat[2])
    scy[1] = (-sc_quat[0] ** 2 + sc_quat[1] ** 2 - sc_quat[2] ** 2 + sc_quat[3] ** 2)
    scy[2] = 2. * (sc_quat[1] * sc_quat[2] + sc_quat[3] * sc_quat[0])

    scz = np.zeros(3, np.float64)
    scz[0] = 2. * (sc_quat[0] * sc_quat[2] + sc_quat[3] * sc_quat[1])
    scz[1] = 2. * (sc_quat[1] * sc_quat[2] - sc_quat[3] * sc_quat[0])
    scz[2] = (-sc_quat[0] ** 2 - sc_quat[1] ** 2 + sc_quat[2] ** 2 + sc_quat[3] ** 2)

    # Calculate geocenter relative to spacecraft pointing
    # using sc direction cosines and sc cartesian position.
    # Tested this going back and forth between coordinate
    # systems & frames, but no formal checks with other
    # data sets.  070815.

    geodir = np.array([np.dot(scx, -sc_pos),
                       np.dot(scy, -sc_pos),
                       np.dot(scz, -sc_pos)]).astype(np.float32)
    geodir = geodir / np.sqrt((geodir ** 2).sum())

    # Calculate geocenter azimuth and zenith angle from geocenter x,y,z

    geo_az = arctan2(geodir[0], geodir[1])
    geo_zen = np.float32(np.arctan2(
        np.linalg.norm(geodir[0:2]), geodir[2]))

    sdir2 = sc_pos / np.sqrt((sc_pos ** 2).sum())
    if sdir2[2] < np.float32(-1.):
        sdir2[2] = np.float64(-1.)
    if sdir2[2] > np.float32(1.):
        sdir2[2] = np.float64(1.)

    scra = arctan2(sdir2[0], sdir2[1])
    scdec = np.float32(np.arcsin(sdir2[2]))

    if verbose:
        print(" %.17f %.17f %.17f %.17f" % tuple(sc_quat))
        print("sc RA new %.15f" % (legacy_dtorad * scra))
        print("sc Dec %.16f" % (legacy_dtorad * scdec))
        print("rad %.5f" % np.sqrt((sc_pos ** 2).sum()))

    return geodir, geo_az, geo_zen, scx, scy, scz


# END get_geocenter()

def read_table(path, ndet, npoints):
    """ Read database table with detector response for each point on the sky.

    NOTE: This method can read from newer numpy files and older text files.

    Parameters
    ----------
    path : str
        Path to file containing the table
    ndet : int
        Number of detectors
    npoints : int32
        Number of points in the sky grid

    Returns
    -------
    t : np.ndarray((ndet, npoints), int32)
        Table with detector response for each sky location.
    idb_no : int32
        Database version number for this table
    """
    # quick load method using numpy file
    if "npy" in path:
        data = np.load(path, allow_pickle=True, encoding='bytes').item()
        return data[b"table"], data[b"idb_no"]

    # otherwise load from text file
    t = np.zeros((ndet, npoints), np.int32)
    f = open(path)

    for i in range(npoints):
        for j, val in enumerate(f.readline().split()):
            t[j][i] = np.int32(val)
        # END for (j, val)
    # END for (i)

    # read database version number
    try:
        idb_no = np.int32(f.readline().strip())
    except BaseException:
        idb_no = 0

    f.close()
    return t, idb_no


# END read_table()

def ang_to_cart_dec(ra, dec):
    """ Convert right ascension and declination to Cartesian coordinates.

    Parameters
    ----------
    ra : float32, np.ndarray(float32)
        Right ascension(s) in radians
    dec : float32, np.ndarray(float32)
        Declination(s) in radians

    Returns
    -------
    pos : np.ndarray(3, float32)
        X, Y, Z position vector(s)
    """
    if isinstance(ra, np.ndarray) or isinstance(dec, np.ndarray):
        return_single = False
    else:
        return_single = True

    ra = np.atleast_1d(ra)
    dec = np.atleast_1d(dec)

    pos = np.zeros((3, ra.size), np.float32)
    pos[0] = np.cos(dec) * np.cos(ra)
    pos[1] = np.cos(dec) * np.sin(ra)
    pos[2] = np.sin(dec)
    pos = pos.T

    return pos[0] if return_single else pos


# END ang_to_cart_dec()

def ang_to_cart_zen(az, zen):
    """ Convert azimuth and zenith to Cartesian coordinates.

    Parameters
    ----------
    az : float32, np.ndarray(float32)
        Azimuth angle(s) in radians
    zen : float32, np.ndarray(float32)
        Zenith angle(s) in radians

    Returns
    -------
    pos : np.ndarray(3, float32)
        X, Y, Z position vector(s)
    """
    if isinstance(az, np.ndarray) or isinstance(zen, np.ndarray):
        return_single = False
    else:
        return_single = True

    az = np.atleast_1d(az)
    zen = np.atleast_1d(zen)

    pos = np.zeros((3, az.size), np.float32)
    pos[0] = np.sin(zen) * np.cos(az)
    pos[1] = np.sin(zen) * np.sin(az)
    pos[2] = np.cos(zen)
    pos = pos.T

    return pos[0] if return_single else pos


# END ang_to_cart_zen()

def j2000_to_sc(scx, scy, scz, pos):
    """ Convert from j2000 position to space craft azimuth and zenith.

    NOTE: Based on sc direction cosines derived from quaternion,
    and a given  j2000 position(x,y,z), return azimuth and zenith angle 
    in sc frame.

    Parameters
    ----------
    scx : np.ndarray(3, float64)
        Space craft X direction cosines
    scy : np.ndarray(3, float64)
        Space craft Y direction cosines
    scz : np.ndarray(3, float64)
        Space craft Z direction cosines
    pos : np.ndarray(3, float32)
        J2000 position in cartesian (x,y,z)

    Returns
    -------
    az : float32
        Azimuth angle in radians
    zen : float32
        Zenith angle in radians
    """

    # Convert source pos to sc frame
    source_pos_sc = np.zeros(3, np.float64)
    source_pos_sc[0] = np.dot(scx.astype(np.float32), pos)
    source_pos_sc[1] = np.dot(scy.astype(np.float32), pos)
    source_pos_sc[2] = np.dot(scz.astype(np.float32), pos)

    # Transform source Cartesian source_pos to az and zenith in sc frame
    az = arctan2(source_pos_sc[0], source_pos_sc[1])
    zen = np.float32(np.arccos(source_pos_sc[2]))

    return az, zen


# END j2000_to_sc()

def sc_to_j2000(scx, scy, scz, pos):
    """ Convert Cartesian source position from space craft frame
    to right ascension and declination in Earth-centered frame.

    NOTE: Based on sc direction cosines derived from quaternion,
    and a given X,Y,Z in space craft frame

    Parameters
    ----------
    scx : np.ndarray(3, float64)
        Space craft X direction cosines
    scy : np.ndarray(3, float64)
        Space craft Y direction cosines
    scz : np.ndarray(3, float64)
        Space craft Z direction cosines
    pos : np.ndarray(3, float32)
        (az, zen) position in cartesian (x,y,z)

    Returns
    -------
    ra : float32
        Right ascension in radians
    dec : float32
        Declination in radians
    """

    if len(pos.shape) == 1:
        pos = pos[np.newaxis, :]
        return_single = True
    else:
        return_single = False

    # First transpose the direction cosines.
    iscx = np.array([scx[0], scy[0], scz[0]], np.float64)
    iscy = np.array([scx[1], scy[1], scz[1]], np.float64)
    iscz = np.array([scx[2], scy[2], scz[2]], np.float64)
    iscx = [scx[0], scy[0], scz[0]]
    iscy = [scx[1], scy[1], scz[1]]
    iscz = [scx[2], scy[2], scz[2]]

    # Convert source pos to J2000 frame
    source_pos_j2 = np.zeros(pos.T.shape, np.float32)
    source_pos_j2[0] = np.dot(pos, iscx)
    source_pos_j2[1] = np.dot(pos, iscy)
    source_pos_j2[2] = np.dot(pos, iscz)

    # Transform source Cartesian source_pos to ra and dec in j2000 frame

    # added vc -- rounding errors can make this less than -1. and
    # dec bombs.
    source_pos_j2[2][source_pos_j2[2] < -1.] = np.float64(-1.)
    source_pos_j2[2][source_pos_j2[2] > 1.] = np.float64(1.)

    dec = np.float32(np.arcsin(source_pos_j2[2]))
    ra = arctan2(source_pos_j2[0], source_pos_j2[1])

    if return_single:
        return ra[0], dec[0]
    return ra, dec


# END sc_to_j2000()

def get_occult(sc_pos, npoints, points_array):
    """ Get list of points occulted by the Earth.

    Parameters
    ----------
    sc_pos : np.ndarray(np.float32)
        Array with X,Y,Z of space craft
    npoints : int32
        Number of points in the sky grid
    points_array : np.ndarray((2, npoints), float32)
        Array with ra/dec locations

    Returns
    -------
    visible : np.ndarray((2, npoints), int32)
        Array with visibility status for each location
        (1 == visible, 0 == occulted)
    """
    sdir = ang_to_cart_dec(points_array[0], points_array[1])
    angl = get_good_angle(sdir, -sc_pos)

    min_vis = np.int32(68.5)
    visible = np.zeros(npoints, np.int32)
    visible[angl > min_vis] = 1

    return visible


# END get_occult()

def get_good_angle(input_xyz1, input_xyz2):
    """ Given 2 Cartesian (unit or not) vectors, return angle between them.

    Parameters
    ----------
    input_xyz1 : np.ndarray(float32)
        First positional vector
    input_xyz2 : np.ndarray(float32)
        Second positional vector

    Returns
    -------
    good_angle : float32
        Angle between the vectors in degrees
    """
    if (len(input_xyz1.shape) > 1) or (len(input_xyz2.shape) > 1):
        return_single = False
    else:
        return_single = True

    # ensure we have the correct dimensions for vector multiplication
    v1 = input_xyz1 if len(input_xyz1.shape) > 1 else input_xyz1[np.newaxis, :]
    v2 = input_xyz2 if len(input_xyz2.shape) > 1 else input_xyz2[np.newaxis, :]

    arg_to_ang = np.sum(
        (v1 / np.linalg.norm(v1, axis=1)[:, np.newaxis]) *
        (v2 / np.linalg.norm(v2, axis=1)[:, np.newaxis]), axis=1)

    arg_to_ang[arg_to_ang < np.float32(-1.)] = np.float64(-1.)
    arg_to_ang[arg_to_ang > np.float32(1.)] = np.float64(1.)

    if return_single:
        arg_to_ang = arg_to_ang[0]

    return np.float32(legacy_dtorad * np.arccos(arg_to_ang).astype(np.float64))


# END get_good_angle()

def compute_scat(npoints, rgrid, cenergies, geodir,
                 nai_az, nai_zen, nai_unit_vec, back_unit_vec,
                 front_only=False):
    """ Function to compute response of NaI detectors over a selected set
    of energy channels. NOTE: Output of this function must be multiplied by
    a spectral shape in order to provide expected detector rates.

    Parameters
    ----------
    npoints : int32
        Number of points in the sky grid
    rgrid : np.array((2, npoints), float32)
        Grid of (az, zen) positions on the sky
    cenergies : np.ndarray(2, float32)
        Array with first and last index of chosen energy channel range
    geodir : np.ndarray(3, np.float64)
        Vector pointing from spacecraft to center of the Earth
    nai_az : np.ndarray(12, float32)
        Azimuth position of each NaI detector on spacecraft in degrees
    nai_zen : np.ndarray(12, float32)
        Zenith position of each NaI detector on spacecraft in degrees
    nai_unit_vec : np.ndarray(12, float32)
        Unit vectors along the direction of each NaI detector on spacecraft
    back_unit_vec : np.ndarray(12, float32)
        Unit vector along the anti-direction of each NaI detector on spacecraft
    front_only : bool
        Only compute forward scattering geometry when True

    Returns
    -------
    front_scattered_rates : np.ndarray((12, npoints, 16), float32)
        Array with energy response to forward scattering computed
        separately for each location on the sky for each NaI detector
    back_scattered_rates : np.ndarray((12, npoints, 16), float32)
        Array with energy response to backward scattering computed
        separately for each location on the sky for each NaI detector
    """
    response_res = np.int32(20)
    front_scattered_rates = np.zeros((12, npoints, 16), np.float32)
    back_scattered_rates = np.zeros((12, npoints, 16), np.float32)

    # define grid for direct response
    grid_spacing = np.float32(190. / (180 / np.arccos(-1)) / (response_res - 1))
    grid_points = grid_spacing * np.arange(response_res, dtype=np.float32)

    # read atmospheric scattering data
    scatterdata = read_scatter_data()

    # choose CONT energy channels from CTIME energy range.
    erange = np.zeros(2, np.int32) - 1
    for i in range(15):
        if (legacy_tenergies[i] <= cenergies[0]) and \
                (legacy_tenergies[i + 1] > cenergies[0]):
            erange[0] = i
        if (legacy_tenergies[i] < cenergies[1]) and \
                (legacy_tenergies[i + 1] >= cenergies[1]):
            erange[1] = i
    # END for (i)

    if erange[0] == -1:
        erange[0] = 0
    if erange[1] == -1:
        erange[1] = 15

    idet = range(12)
    rcart = ang_to_cart_zen(rgrid[0], rgrid[1])

    # scattering for front detectors
    geom_fac_front = np.float32(126 * 50 / (2025. * 0.8))
    earthpoints, gperp, gpmag, elev = read_earthpoints(
        geodir, nai_unit_vec, back_unit_vec, 0)
    front_scattered_rates = geom_fac_front * get_scattered_rates(
        rcart, gperp, gpmag, elev, earthpoints,
        response_res, grid_points, scatterdata, 0,
        nai_az, nai_zen, nai_unit_vec, back_unit_vec)

    if front_only:
        return front_scattered_rates, back_scattered_rates

    # scattering for back detectors
    geom_fac_back = np.float32(126 * 50 / (2025. * 0.8))
    earthpoints, back_gperp, back_gpmag, back_elev = read_earthpoints(
        geodir, nai_unit_vec, back_unit_vec, 1)
    back_scattered_rates = geom_fac_back * get_scattered_rates(
        rcart, gperp, gpmag, back_elev, earthpoints,
        response_res, grid_points, scatterdata, 1,
        nai_az, nai_zen, nai_unit_vec, back_unit_vec)

    return front_scattered_rates, back_scattered_rates


# END compute_scat()

def read_earthpoints(geodir, nai_unit_vec, back_unit_vec, scat_geom,
                     path=__data_dir__ + 'earth_points.npy'):
    """ Read array of position vectors used to compute
    scattering of photons off Earth's atmosphere.

    NOTE: This method can read from newer numpy files and older text files.

    Parameters
    ----------
    geodir : np.ndarray(3, np.float64)
        Vector pointing from spacecraft to center of the Earth
    nai_unit_vec : np.ndarray(12, float32)
        Unit vectors along the direction of each NaI detector on spacecraft
    back_unit_vec : np.ndarray(12, float32)
        Unit vector along the anti-direction of each NaI detector on spacecraft
    scat_geom : int32
        Scattering geometry. 0=forward, 1=backward.

    Returns
    -------
    earthpoints : np.ndarray((3,236), float32)
        Array of position vectors relative to the Earth that were used to 
        build atmospheric scattering matrix obtained with read_scatter_data()
    gperp : np.ndarray(3, 12), float32)
        Direction perpendicular to line between
        Earth center and each NaI detector
    gpmag : np.ndarray(12, float32)
        Magnitude of each gperp vector
    elev : np.ndarray(12, float32)
        Elevation angle for each detector in radians
    """
    if 'npy' in str(path):
        # quick load method using numpy file
        earthpoints = np.load(path)
    else:
        # otherwise load from text file
        earthpoints = np.loadtxt(path, np.float32).transpose()

    gperp = np.zeros((3, 12), np.float32)
    gpmag = np.zeros(12, np.float32)
    elev = np.zeros(12, np.float32)

    # define earthpoints for detectors given geodir
    for idet in range(12):
        if scat_geom == 0:
            this_unit_vec = np.array([nai_unit_vec[i][idet] for i in range(3)])
        else:
            this_unit_vec = np.array([back_unit_vec[i][idet] for i in range(3)])

        this_angle = np.float32(
            np.float64(get_good_angle(geodir, this_unit_vec)) / legacy_dtorad)
        elev[idet] = legacy_rpi / np.float64(2.) - this_angle

        proj = dot(geodir, this_unit_vec)
        gperp[0][idet] = geodir[0] - proj * this_unit_vec[0]
        gperp[1][idet] = geodir[1] - proj * this_unit_vec[1]
        gperp[2][idet] = geodir[2] - proj * this_unit_vec[2]
        gpmag[idet] = np.sqrt(np.sum(gperp[:, idet] * gperp[:, idet], dtype=np.float32))

    return earthpoints, gperp, gpmag, elev


# END read_earthpoints()

def read_scatter_data(path=__data_dir__ + 'alocdat_comp.npy'):
    """ Read atmospheric scattering matrix.

    NOTE: This method can read from newer numpy files and older text files.

    Returns
    -------
    scatterdata : np.ndarray((16, 236, 19), float32)
        Array with scattering matrix
    """
    # quick load method using numpy file
    if 'npy' in path:
        return np.load(path)

    # otherwise load from text file
    loadtxt = np.loadtxt(path, np.float32)
    scatterdata = np.zeros((16, 236, 19), np.float32)

    j = np.arange(236)
    for k in range(19):
        scatterdata.T[k] = loadtxt[236 * k + j]

    return scatterdata


# END read_scatter_data()

def get_scattered_rates(rcart, gperp, gpmag, elev, earthpoints, response_res,
                        grid_points, scatterdata, scat_geom,
                        nai_az, nai_zen, nai_unit_vec, back_unit_vec):
    """Calculate amount of scattering given geometry.

    Parameters
    ----------
    rcart : np.array((2, npoints), float32)
        Grid of (az, zen) positions on the sky
    gperp : np.ndarray(3, 12), float32)
        Direction perpendicular to line between
        Earth center and each NaI detector
    gpmag : np.ndarray(12, float32)
        Magnitude of each gperp vector
    elev : np.ndarray(12, float32)
        Elevation angle for each detector in radians
    earthpoints : np.ndarray((3,236), float32)
        Array of position vectors relative to the Earth that were used to 
        build atmospheric scattering matrix obtained with read_scatter_data()
    response_res : int32
        Respolution of the response matrix as an integer bin number
    grid_points : np.ndarray(response_res, float32)
        Array of zenith angles for each response point in degrees
    scatterdata : np.ndarray((16, 236, 19), float32)
        Array with scattering matrix
    scat_geom : int32
        Scattering geometry. 0=forward, 1=backward.
    nai_az : np.ndarray(12, float32)
        Azimuth position of each NaI detector on spacecraft in degrees
    nai_zen : np.ndarray(12, float32)
        Zenith position of each NaI detector on spacecraft in degrees
    nai_unit_vec : np.ndarray(12, float32)
        Unit vectors along the direction of each NaI detector on spacecraft
    back_unit_vec : np.ndarray(12, float32)
        Unit vector along the anti-direction of each NaI detector on spacecraft
 
    Returns
    -------
    scattered_rates : np.ndarray
        Rates of scatter photons in each detector
    """
    if scat_geom == 0:
        scat_unit_vec = nai_unit_vec
        scat_geom_fac = np.float32(1.00)
    else:
        scat_unit_vec = back_unit_vec
        scat_geom_fac = np.float32(0.0)

    c = np.array([np.cross(rcart, scat_unit_vec.T[i]) for i in range(12)])
    a = np.array([np.cross(c[i], scat_unit_vec.T[i]) for i in range(12)])
    amag = np.linalg.norm(a, axis=2)

    min_val = np.float32(1e-5)
    mask = ((gpmag[:, np.newaxis] > min_val) & (amag > min_val))

    x = np.zeros(a.shape[:2], a.dtype)
    for i in range(3):
        x += a[:, :, i] * gperp[i][:, np.newaxis]
    x /= amag * gpmag[:, np.newaxis]

    x[x > np.float32(1.0)] = np.float32(1.0)
    x[x < np.float32(-1.0)] = np.float32(-1.0)

    az = np.zeros_like(mask, np.float32)
    az[mask] = np.arccos(x[mask])

    min_elev = np.float32(-1.48353)
    max_elev = np.float32(+1.48353)

    elev[elev < min_elev] = min_elev
    elev[elev > max_elev] = max_elev

    eindex = np.argmin(earthpoints.T[:, 2, np.newaxis] >
                       elev[np.newaxis, :], axis=0)
    eindex[elev == max_elev] = np.argmax(earthpoints[2, :] == max_elev)

    heindex = eindex - 1

    arrshape = (earthpoints.shape[1], earthpoints.shape[0], elev.shape[0])
    eearthpoints = np.zeros(arrshape, earthpoints.dtype) + np.float32(99.)
    heearthpoints = eearthpoints.copy()

    emask = (earthpoints.T[:, 2, np.newaxis] ==
             earthpoints.T[np.newaxis, eindex, 2])
    emask = np.repeat(emask[:, np.newaxis, :], 3, axis=1)

    hemask = (earthpoints.T[:, 2, np.newaxis] ==
              earthpoints.T[np.newaxis, heindex, 2])
    hemask.T[heindex == -1] = False
    hemask = np.repeat(hemask[:, np.newaxis, :], 3, axis=1)

    eearthpoints[emask] = np.repeat(earthpoints.T[:, :, np.newaxis], 12, axis=2)[emask]
    heearthpoints[hemask] = np.repeat(earthpoints.T[:, :, np.newaxis], 12, axis=2)[hemask]

    aindex = np.argmin(np.fabs(eearthpoints[:, 1, :, np.newaxis] - az[np.newaxis, :, :]), axis=0)
    haindex = np.argmin(np.fabs(heearthpoints[:, 1, :, np.newaxis] - az[np.newaxis, :, :]), axis=0)

    felev = np.zeros_like(aindex, dtype=np.float32)
    min_val = np.float32(0.001)
    mask = (earthpoints.T[haindex, 2] - earthpoints.T[aindex, 2]) > min_val
    felev[mask] = ((elev[:, np.newaxis] - earthpoints.T[aindex, 2])[mask] /
                   (earthpoints.T[haindex, 2] - earthpoints.T[aindex, 2])[mask])
    felev[mask] = ((elev[:, np.newaxis] - earthpoints[2, aindex])[mask] /
                   (earthpoints[2, haindex] - earthpoints[2, aindex])[mask])

    bdangle = np.float32(
        np.arccos(np.sum(scat_unit_vec[np.newaxis, :, :] * \
                         rcart[:, :, np.newaxis], axis=1)).T)

    scattered_rates = np.zeros((12, rcart.shape[0], 16), dtype=np.float32)
    for ien in range(16):
        r1 = interpolatex(response_res, grid_points,
                          scatterdata[ien, aindex, :], bdangle)
        r2 = interpolatex(response_res, grid_points,
                          scatterdata[ien, haindex, :], bdangle)
        scattered_rates[:, :, ien] = scat_geom_fac \
                                     * ((np.float32(1.0) - felev) * r1 + felev * r2)
    # END for (ien)

    return scattered_rates


# END get_scattered_rates()

def dot(x, y):
    """ Dot product of two vectors

    Parameters
    ----------
    x : np.ndarray(float32/64)
        First vector
    y : np.ndarray(float32/64)
        Second vector

    Returns
    -------
    dprod : float32/64
        Dot product of x and y
    """
    return x[0] * y[0] + x[1] * y[1] + x[2] * y[2]


def crossprod(x, y):
    """ Cross product of two vectors

    Parameters
    ----------
    x : np.ndarray(float32)
        First vector
    y : np.ndarray(float32)
        Second vector

    Returns
    -------
    cprod : np.ndarray(float32)
        Cross product of x cross y
    """
    cprod = np.zeros(3, np.float32)
    cprod[0] = x[1] * y[2] - x[2] * y[1]
    cprod[1] = x[2] * y[0] - x[0] * y[2]
    cprod[2] = x[0] * y[1] - x[1] * y[0]

    return cprod


# END crossprod()

def all_idx(idx, axis=0):
    """Make a multidimensional index mask into another array
    
    Parameters:
    -----------
    idx: np.array
        The index array
    axis: int, optional
        The axis over which the index array will be applied
        
    Returns:
    --------
    grid: tuple
        The multidimensional index mask
    """
    grid = list(np.ogrid[tuple(map(slice, idx.shape))])
    grid.insert(axis, idx)
    return tuple(grid)


# END all_idx()

def interpolatex(int_size, inter_array, apply_array, bit_along):
    """ Interpolate between 2 values by weighting by distance between 2 pts.

    Parameters
    ----------
    int_size : int32
        Length of inter_arr between which we're interpolating
    inter_array : np.ndarray(int_size, float32)
        Array of bin endpoints between which we'll interpolate
    apply_array : np.ndarray((ldet, npoints, int_size - 1), float32)
        Array with contents for each bin. 
    bit_along : np.ndarray((ldet, npoints), float32)
        Distance of bin along inter_array

    Returns
    -------
    value : float32, np.ndarray(float32)
        Interpolated value(s)
    """
    x = np.fabs(inter_array[:, np.newaxis, np.newaxis] -
                bit_along[np.newaxis, :, :]).argmin(axis=0)
    x[inter_array[x] > bit_along] -= 1

    f_scale = (bit_along - inter_array[x]) \
              / (inter_array[x + 1] - inter_array[x])
    aidx = all_idx(x)
    aidx = (aidx[1], aidx[2], aidx[0])
    aidx1 = all_idx(x + 1)
    aidx1 = (aidx1[1], aidx1[2], aidx1[0])
    return f_scale * apply_array[aidx1] + \
        (np.float32(1.) - f_scale) * apply_array[aidx]


# END interpolatex()

def initialize_det_geometry(verbose=False):
    """ Contains geometry from rmk of 07/07
    Given nai el & az, calculate unit vectors in SC frame of each NaI
    Tested 08/14/07 -- produces correct unit vector angles
    as per geometry from rmk of 07/07.

    Parameters
    ----------
    verbose : bool
        Print info to screen when True

    Returns
    -------
    nai_az : np.ndarray(12, float32)
        Azimuth position of each NaI detector on spacecraft in degrees
    nai_zen : np.ndarray(12, float32)
        Zenith position of each NaI detector on spacecraft in degrees
    nai_unit_vec : np.ndarray(12, float32)
        Unit vectors along the direction of each NaI detector on spacecraft
    back_unit_vec : np.ndarray(12, float32)
        Unit vector along the anti-direction of each NaI detector on spacecraft
    """
    nai_az = np.array(
        [45.89, 45.11, 58.44, 314.87, 303.15, 3.35,
         224.93, 224.62, 236.61, 135.19, 123.73, 183.74], np.float32)
    nai_zen = np.array(
        [20.58, 45.31, 90.21, 45.24, 90.27, 89.79,
         20.43, 46.18, 89.97, 45.55, 90.42, 90.32], np.float32)

    nai_az_rad = np.float32(np.float64(nai_az) / legacy_dtorad)
    nai_zen_rad = np.float32(np.float64(nai_zen) / legacy_dtorad)

    nai_unit_vec = np.zeros((3, 12), np.float32)
    back_unit_vec = np.zeros((3, 12), np.float32)
    for j in range(12):
        vec = ang_to_cart_zen(nai_az_rad[j], nai_zen_rad[j])
        for i in range(3):
            nai_unit_vec[i][j] = vec[i]
            back_unit_vec[i][j] = -nai_unit_vec[i][j]
        # END for (i)

        # Transform source Cartesian source_pos to az and zenith  in sc frame
        zen = np.float32(np.arccos(back_unit_vec[2][j]))
        az = arctan2(back_unit_vec[0][j], back_unit_vec[1][j])

        if verbose:
            print(" az and el of front and back of dets")
            print("    {0:.9} {1:.9} {2:.17} {3:.17}".format(nai_az[j],
                                                             nai_zen[j], legacy_dtorad * az, legacy_dtorad * zen))
            print("")
    # END for (j)

    return nai_az, nai_zen, nai_unit_vec, back_unit_vec


# END initialize_det_geometry()

def get_spec(spec, energies, mid_energies, erange):
    """ Get spectral shape in each energy bin

    Parameters
    ----------
    spec : str
        Spectrum definition string formatted like so:
        comp function - "comp,index=-1.15,epeak=350.0"
        band function - "band,alpha=-1.0,beta=-2.3,epeak=230.0"
        pl   function - "pl,index=-2"
    energies : np.ndarray(float32, nen + 1)
        End points of energy bins
    mid_energies : np.ndarray(float32, nen)
        Energy describing the geometric mean between the ends of an energy bin
    erange : list
        A 2 element list containing the first and last indices used for
        normalizing the spectral shape

    Returns
    -------
    spec_cts : np.ndarray(float32)
        Array with normalized photon counts in each energy bin for spec
    """
    # list of available spectral functions
    spec_func = {"pl": legacy_spectral_models.legacy_pl,
                 "comp": legacy_spectral_models.legacy_comp,
                 "band": legacy_spectral_models.legacy_band}

    # list of parameters required by each spectral definition
    required_param = {
        "pl": ["index"],
        "comp": ["index", "epeak"],
        "band": ["alpha", "beta", "epeak"],
    }

    # parse spectrum type from spec string
    spec_type = spec.split(',')[0]
    if spec_type not in list(spec_func.keys()):
        raise ValueError("Spectral type '%s' not recognized. Use -h to see a list"
                         " of spectral shapes allowed by --spec option")

    # parse spectral parameters from spec string
    spec_param = {v.split("=")[0]: np.float32(v.split("=")[1])
                  for v in spec.split(',')[1:]}
    fnorm = spec_param.pop("amp", np.float32(1.0))

    # check for unrecognized parameters
    for key in list(spec_param.keys()):
        if key not in required_param[spec_type]:
            raise ValueError("Unrecognized parameter '%s' for '%s' spectrum" %
                             (key, spec_type))

    # check for missing parameters
    for req in required_param[spec_type]:
        if req not in list(spec_param.keys()):
            raise ValueError("Missing parameter '%s' for '%s' spectrum" %
                             (req, spec_type))

    spec_cts = spec_func[spec_type](mid_energies, **spec_param)
    spec_cts *= (energies[1:] - energies[:-1])

    f = np.float32(0.)
    for i in range(erange[0], erange[1] + 1):
        f += spec_cts[i]

    spec_cts = spec_cts * fnorm / f
    return spec_cts


# END get_spec()

def add_scat(npoints, rgrid, front_scattered_rates, back_scattered_rates,
             cenergies, spec):
    """ Add scattering data for each spectrum 

    Parameters
    ----------
    npoints : int32
        Number of points in the sky grid
    rgrid : np.array((2, npoints), float32)
        Grid of (az, zen) positions on the sky
    front_scattered_rates : np.ndarray((12, npoints, 16), float32)
        Array with energy response to forward scattering computed
        separately for each location on the sky for each NaI detector
    back_scattered_rates : np.ndarray((12, npoints, 16), float32)
        Array with energy response to backward scattering computed
        separately for each location on the sky for each NaI detector
    cenergies : np.ndarray(2, float32)
        Array with first and last index of chosen energy channel range
    spec : str
        Spectrum definition string formatted like so:
        comp function - "comp,index=-1.15,epeak=350.0"
        band function - "band,alpha=-1.0,beta=-2.3,epeak=230.0"
        pl   function - "pl,index=-2"

    Returns
    -------
    atm_scattered_rates : np.ndarray((14, npoints), np.float32)
        Scattered rates in each detector for a source coming
        from each point on the sky
    """
    atm_scattered_rates = np.zeros((14, npoints), np.float32)

    # Choose CONT energy channels from CTIME energy range.
    erange = np.zeros(2, np.int32) - 1
    for i in range(15):
        if (legacy_tenergies[i] < cenergies[0]) and \
                (legacy_tenergies[i + 1] > cenergies[0]):
            erange[0] = i
        if (legacy_tenergies[i] < cenergies[1]) and \
                (legacy_tenergies[i + 1] > cenergies[1]):
            erange[1] = i
    # END for (i)

    if erange[0] == -1:
        erange[0] = 0
    if erange[1] == -1:
        erange[1] = 15
    # print("ERANGE: ", erange, cenergies, tenergies[erange[0]], tenergies[erange[1]+1])

    # Calculate spectra for atmospheric scattering
    mid_energies = np.sqrt(legacy_tenergies[:-1] * legacy_tenergies[1:])
    scat_spec = get_spec(spec, legacy_tenergies, mid_energies, erange)

    for i in range(16):
        for idet in range(12):
            atm_scattered_rates[idet + 2][:] += \
                front_scattered_rates[idet, :, i] * scat_spec[i]
    # END for (idet)
    # END for (i)

    return atm_scattered_rates


# END add_scat()

def find_best_location(ndet, vpoints, udet, usedet, loctable_entries,
                       sduration, c_mrates, c_brates, visible, verbose=False):
    """ Calculate chi2 for each visible location and find minimum.
    Algorithm for normalization from MSB memo of Jul 27 2005

    Parameters
    ----------
    ndet : int
        Number of detectors
    vpoints : int32
        Number of points in the sky grid
    udet : int32
        Length os usedet arg
    usedet: np.array()
        A boolean array containing True (1) if the detector is to be
        included or False (0) if a detector is not to be included.  If not
        set, then all detectors (including BGOs) are used
    loctable_entries : np.ndarray((ndet, npoints), int32)
        Table with detector response for each sky location.
    sduration : float32
        Timescale of measured detector counts
    c_mrates : np.ndarray(ndet, int32)
        Measured detector counts summed over chosen energy channels
    c_brates : np.ndarray(ndet, int32)
        Estimated background counts summed over chosen energy channels
    visible : np.ndarray((2, npoints), int32)
        Array with visibility status for each location
        (1 == visible, 0 == occulted)
    verbose : bool
        Print things to screen when True

    Returns
    -------
    gindex : int32
        Index of best localization position in the sky
    chi2 : np.ndarray(vpoints, np.float32)
        Array of chi-square values computed by taking the difference between
        the observed and expected excesses for each point on the sky
    nchi2 :
        Array of chi-square values computed by taking the difference between
        the observed and expected excesses for each point on the sky using
        an alternative normalization method for the expected excess
    rchi2 : float32
        Value of chi2 array for best location on the sky
    az : float32
        Azimuth of rchi2 in degrees
    zen : float32
        Zenith of rchi2 in degrees
    loc_err : float32
        Localization error in degrees
    loc_reliable : bool
        Status of best localization error.
        True=GOOD, False=BAD, cannot be localized
    """
    # compute normalization
    a = np.zeros(vpoints, np.float32)
    b = np.zeros(vpoints, np.int64)
    for d in usedet:
        loc_int64 = np.int64(loctable_entries[d + 2, :])
        a += np.float32(loc_int64) * np.float32(c_mrates[d] - c_brates[d]) / np.float32(c_mrates[d])
        b += np.int64(loc_int64 ** 2 / c_mrates[d])
    # END for (d)
    fnorm = np.float32(a / b.astype(np.float32))

    # compute chi2 relative to background + flux normalization
    chi2 = np.zeros(vpoints, np.float32)
    for d in usedet:
        loc_float32 = np.float32(loctable_entries[d + 2, :])
        chi2 += (np.float32(c_mrates[d] - c_brates[d]) - fnorm * loc_float32) ** 2 \
                / (np.float32(c_brates[d]) + fnorm * loc_float32)
    # END for (d)

    chi2[(chi2 < np.float32(0.)) |
         (chi2 > np.float32(9999999.))] = np.float32(9999999.)

    i = chi2.argmin()
    gindex = i
    rchi2 = chi2[gindex]
    az = loctable_entries[0][i] / np.float32(60.)
    zen = loctable_entries[1][i] / np.float32(60.)
    #    print("min chi2: {0} {1:.9} {2:.9} {3:.9}".format(i, chi2[i], az, zen))
    mask = (chi2 != chi2[i])
    jmask = chi2[mask].argmin()
    j = np.arange(vpoints)[mask][jmask]
    jaz = loctable_entries[0][j] / np.float32(60.)
    jzen = loctable_entries[1][j] / np.float32(60.)
    #    print("next min chi2: {0} {1:.9} {2:.9} {3:.9}".format(j, chi2[j], jaz, jzen))

    a1 = np.float32(az / legacy_dtorad)
    a2 = np.float32(jaz / legacy_dtorad)
    z1 = np.float32(zen / legacy_dtorad)
    z2 = np.float32(jzen / legacy_dtorad)
    t1 = ang_to_cart_zen(a1, z1)
    t2 = ang_to_cart_zen(a2, z2)
    laz = np.float32(loctable_entries[0, :] / np.float32(60.) / legacy_dtorad)
    lzen = np.float32(loctable_entries[1, :] / np.float32(60.) / legacy_dtorad)

    # test code for fiducial intensity check for chi2 reliability
    rates = c_mrates[:12] - c_brates[:12]
    x = rates.argmax()
    mask = rates != rates[x]
    y = np.arange(rates.size)[mask][rates[mask].argmax()]
    sc_mrates = np.float32(c_mrates[:12]) / sduration
    #    print("sc_mrates", ['{0:.9}'.format(sc_mrates[i]) for i in range(sc_mrates.size)])
    sc_brates = np.float32(c_brates[:12]) / sduration
    #    print(["{0:.9}".format(sc_mrates[i] - sc_brates[i]) for i in range(sc_brates.size)])
    norm = np.float32(1000.) / (sc_mrates[x] - sc_brates[x] + sc_mrates[y] - sc_brates[y])
    nc_mrates = np.int32((sc_mrates - sc_brates) * norm + sc_brates)

    a = np.float32(0.)
    b = np.float32(0.)
    for d in usedet:
        loc_int64 = np.int64(loctable_entries[d + 2, gindex])
        a += np.float32(loc_int64) * (np.float32(nc_mrates[d]) - sc_brates[d]) / np.float32(nc_mrates[d])
        b += np.float32(loc_int64 ** 2) / np.float32(nc_mrates[d])
    # END for (d)
    nfnorm = np.float32(a / b)

    nchi2 = np.float32(0.)
    for d in usedet:
        a = np.float32((np.float32(nc_mrates[d] - sc_brates[d])
                        - nfnorm * np.float32(loctable_entries[d + 2, gindex])) ** 2)
        b = (sc_brates[d] + nfnorm * np.float32(loctable_entries[d + 2, gindex]))
        nchi2 += a / b
    # END for (d)

    includethispoint = np.zeros(vpoints, np.int32)
    error_determined = False
    offset_chi2_delta = np.float32(0.)
    loc_reliable = True
    max_rel_chi2 = np.int32(500)

    loc_err_vsmall = False
    if nchi2 >= max_rel_chi2:
        loc_reliable = False
    loc_err = np.float32(50.)
    if chi2[j] - chi2[i] > np.float32(2.3):
        loc_err_vsmall = True
        loc_err = np.float32(1.)

    while not error_determined and loc_reliable and not loc_err_vsmall:

        loc_err = np.float32(0.)

        mask = ((chi2 - chi2[gindex] >= np.float32(2.28) - offset_chi2_delta) &
                (chi2 - chi2[gindex] <= np.float32(2.32) + offset_chi2_delta))

        t3 = ang_to_cart_zen(laz[mask], lzen[mask])
        includethispoint[mask] = 1
        for distt3 in get_good_angle(t1, t3):
            loc_err = loc_err + distt3
        # END for (distt3)

        if includethispoint.sum() > 1:
            error_determined = True
        if offset_chi2_delta >= np.float32(2.38):
            loc_err = np.float32(50.)
            loc_reliable = False
            error_determined = True
            if (verbose):
                print("location unreliable")

        offset_chi2_delta = offset_chi2_delta + np.float32(0.05)

    # END while (not error_determined)

    if loc_reliable and not loc_err_vsmall:
        loc_err = loc_err / np.float32(includethispoint.sum())
        # added for small error --- if it can't find delta chi2>2.3 it gets back to itself.
        if loc_err < np.float32(1.0):
            loc_err = np.float32(1.)

    return gindex, chi2, nchi2, rchi2, az, zen, loc_err, loc_reliable


# END find_best_location()

def get_det_geometry(input_az, input_zen, nai_az, nai_zen):
    """ Get detector geometry relative to source.

    Based on  geometry provided by T. Morse, calculate
    angles to detectors from calling az and zen.

    Parameters
    ----------
    input_az : float32
        Source azimuth in radians
    input_zen : float32
        Source zenith in radians
    nai_az : np.ndarray(12, float32)
        Azimuth position of each NaI detector on spacecraft in degrees
    nai_zen : np.ndarray(12, float32)
        Zenith position of each NaI detector on spacecraft in degrees

    Returns
    -------
    det_angs : np.ndarray(ldet, float32)
        Detector angles to source in radians
    """
    # Calculate angles from detectors using 
    # detector az and zen's from T.Morse.
    det_angs = np.float32(np.arccos(
        np.cos(np.float64(nai_zen) / legacy_dtorad) * np.cos(input_zen) +
        np.sin(np.float64(nai_zen) / legacy_dtorad) * np.sin(input_zen) *
        np.cos(np.float64(nai_az) / legacy_dtorad - input_az)))

    return det_angs


# END get_det_geometry()

def eq2000_to_gal_r(ra2000, dec2000):
    """ Function to convert from equatorial (ra, dec) coordinates in J2000
    epoch to Galactic coordinates (lii, bii)

    Parameters
    ----------
    ra2000 : float32
        J2000 right ascension in radians 
    dec2000 : float32
        J2000 declination in radians

    Returns
    -------
    lii : float32
        Galactic longitude in radians
    bii : float32
        Galactic latitude in radians

    Version V1.1

    Subroutine to convert equatorial coordinates, RA_2000 & DEC_2000
    (epoch 2000), to galactic coordinates, Lii and Bii.
    All arguments are in RADIANS.
    Michael S. Briggs, 16 July 1992, UAH / MSF! ES-62.

    Method:
    Step 1: Equatorial (epoch 2000) --> Equatorial (epoch 1950),
    Step 2: Equatorial (epoch 1950) --> Galactic.
    These two steps are required because GRO uses epoch 2000 equatorial
    coordinates and transformations from equatorial to galactic are
    for epoch 1950 equatorial coordinates.
    In essence, the galactic system is defined by the location of the
    galactic center and galactic North pole, and these locations are defined
    by the values they were assigned in 1950 coordinates.  See J. Meus.
    Step 1, precession from epoch 2000 to 1950, is somewhat simplified
    from the procedure given in The Astronomical Almanac.   I have left
    out corrections which are probably unnecessary for BATSE's position
    accuracy, e.g. I have left out the E-terms, which take care of the
    abberation effects due to the ellipticty of the earth's orbit.

    Step 1 is based upon p. B43 of the 1991 The Astronomical Almanac.
    Step 2 is baed upon a rotation matrix created based upon the definition
    of galactic coordinates.   See subroutine GET_ETOG_MATRIX.

    References:
    Astronomical Algorithms, J. Meeus, pp. 89 & 90, pp. 123-130 (esp.
    pp. 129 & 130).
    Astrophysical Formulae, 2nd. ed., K. R. Lang, pp. 498 & 504.
    Practical Astronomy with your !alculator, 3rd ed., P. Duffett-Smith,
    pp. 43 & 50.
    Third Reference !atalogue of Bright Galaxies, de Vancouleurs et al.,
    vol. 1, p. 11.
    Third Reference !atalogue of Bright Galaxies, de Vancouleurs et al.,
    vol. 2 & 3.
    The Astronomical Almanac for the year 1991, p. B43.

    Tested as follows:
       The Third Reference !atalogue of Bright Galaxies, volumes 2 & 3,
    list the coordinates of many galaxies.   The coordinates listed
    include RA & DEC epoch 2000, RA & DEC epoch 1950, and Lii and Bii,
    hence these volumes may serve as a Rosetta Stone for the testing
    of coordinate systems.  The volumes list RA to 0.1 sec of time and
    DEC to 1 sec of arc.  They list Lii and Bii both to 0.01 degrees.
       Running this subroutine with input the RA and DEC of randomly
    selected galaxies from these volumes, the output Lii and Bii were
    always found to exactly agree with the numbers in the volumes.
    Hence this subroutine is shown to be accurate to at least about
    0.01 degrees.

    Further testing: Mark Finger wrote a program for the same purpose,
    known as J2000_TO_GALII.   The programs were compared using 5000
    points randomly located on the sphere.   The maximum discrepancy
    was 8 milli-arc seconds!

    Input arguments:
    REAL*4 RA_2000                 ! right ascension, RADIANS
    REAL*4 DEC_2000                ! declination, RADIANS

    Output arguments:
    REAL*4 LII                     ! galactic longitude, RADIANS
    REAL*4 BII                     ! galactic latitude, RADIANS

    Matrix to convert equatorial (epoch 2000) to equatorial (epoch 1950).
    Remember that FORTRAN stores matrices backwards.   Taken from page B43
    of The Astronomical Almanac.
    """

    # NOTE: first set minv type as float32 and then caste to float64 because
    # this is how it was done in the original FORTRAN code. Need to replicate
    # legacy behavior.    
    minv = np.array(
        [[+0.9999256795, +0.0111814828, +0.0048590039],
         [-0.0111814828, +0.9999374849, -0.0000271771],
         [-0.0048590040, -0.0000271557, +0.9999881946]], np.float32).astype(np.float64)

    # spherical to cartesian:
    ra = np.float64(ra2000)
    dec = np.float64(dec2000)

    x1 = np.zeros(3, np.float64)
    x1[0] = np.cos(dec) * np.cos(ra)
    x1[1] = np.cos(dec) * np.sin(ra)
    x1[2] = np.sin(dec)

    # equatorial, epoch 2000 --> epoch 1950:
    # do the transformation via matrix multiplication:
    x2 = np.zeros(3, np.float64)
    for i in range(3):
        for j in range(3):
            x2[i] += minv[i][j] * x1[j]
        # END for (j)
    # END for (i)

    # obtain the matrix for equatorial
    etog = ETOG().value

    # equatorial, epoch 1950 --> galactic:
    # do the transformation via matrix multiplication:
    x3 = np.zeros(3, np.float64)
    for i in range(3):
        for j in range(3):
            x3[i] += etog[i][j] * x2[j]
        # END for (j)
    # END for (i)

    # cartesian to spherical:
    lii = np.float32(np.arctan2(x3[1], x3[0]))
    if lii < np.float32(0.0):
        lii += np.float32(np.float32(2.0) * legacy_pi)
    bii = np.float32(np.arcsin(x3[2]))

    return lii, bii


# END eq2000_to_gal_r()

def get_etog_matrix():
    """ Function to calculate equatorial (ra, dec) to Galactic (lii, bii)
    coordinate rotation matrix.

    Returns
    -------
    etog : np.ndarray((3,3), float64)
        Rotation matrix for transforming equatorial to galactic vectors

    Version V1.0

    This subroutine creates the 3x3 matrix that converts 1950 equatorial
    coordinates into galactic coordinates lii and bii.
    This subroutine is called by EQ1950_TO_GAL and by EQ2000_TO_GAL.
    Michael S. Briggs, 16 July 1992, UAH / MSFC ES-62.

    The transformation is based upon the following, which is the definition
    of the galatic coordinate system:
        Location                 RA 1950               DEC 1950
    North Galatic Pole    12h49m = 192.25 deg    +27d24' = +27.4 deg
      Galactic Origin     17h42.4m = 265.6 deg   -28d55' = -28.9166666..
    Also the ascending node of the galactic plane on the celestial equator
    is lii = 33.0 degrees.   All the above numbers are exact by definition
    in 1950 coordinates.
    See:
    1) Practical Astronomy with your Calculator, 3rd ed., 
       P. Duffett-Smith, p. 43 & 50.
    2) Astrophysical Formulae, 2nd ed., K. R. Lang, p. 498 & 504.
    3) Third Reference Catalogue of Bright Galaxies, vol. 1, p. 11.
    4) Astronomical Algorithms, J. Meeus, p. 89 & 90.

    Verification of this subroutine:
    The matrix calculated by this subroutine is a higher precision version
    of the matrix given on page 50 of Practical Astron. with your Calculator.
    The numbers in the matrix in that book agree with those calculated
    by this subroutine to +/-1 on the last digit.
    """

    etog = np.zeros((3, 3), np.float64)

    a = np.zeros((3, 3), np.float64)
    b = np.zeros((3, 3), np.float64)
    c = np.zeros((3, 3), np.float64)
    d = np.zeros((3, 3), np.float64)

    # Rotate about equatorial Z-axis so that projection of galactic Z-axis
    # onto equatorial XY-plane is perpendicular to Y-axis:
    # Angle comes from RA of galactic N pole:
    ang1 = np.float64(np.float32(-192.25) * legacy_deg_to_rad)

    a[0][0] = np.cos(ang1)
    a[0][1] = -np.sin(ang1)
    a[0][2] = 0.
    a[1][0] = np.sin(ang1)
    a[1][1] = np.cos(ang1)
    a[1][2] = 0.
    a[2][0] = 0.
    a[2][1] = 0.
    a[2][2] = 1.

    # Rotate about Y-axis so that current Z-axis is rotate to direction of
    # galactic Z-axis.   Angle comes from declination of galactic N pole:
    ang2 = np.float64(np.float32(90. - 27.4) * legacy_deg_to_rad)

    b[0][0] = np.cos(ang2)
    b[0][1] = 0.
    b[0][2] = -np.sin(ang2)
    b[1][0] = 0.
    b[1][1] = 1.
    b[1][2] = 0.
    b[2][0] = np.sin(ang2)
    b[2][1] = 0.
    b[2][2] = np.cos(ang2)

    # Above rotations establish galactic latitude b, now rotate about
    # galactic Z-axis to get the longitudes l correct:
    # Angle comes from the ascending node:
    ang3 = (np.float64(np.float32(33.) * legacy_deg_to_rad
                       - legacy_pi / np.float32(2.)))

    c[0][0] = np.cos(ang3)
    c[0][1] = -np.sin(ang3)
    c[0][2] = 0.
    c[1][0] = np.sin(ang3)
    c[1][1] = np.cos(ang3)
    c[1][2] = 0.
    c[2][0] = 0.
    c[2][1] = 0.
    c[2][2] = 1.

    # Multiply the 3 rotation matrices to obtain the transformation
    # matrix:
    for i in range(3):
        for j in range(3):
            for k in range(3):
                d[i][j] += b[i][k] * a[k][j]
            # END for (k)
        # END for (j)
    # END for (i)

    for i in range(3):
        for j in range(3):
            for k in range(3):
                etog[i][j] += c[i][k] * d[k][j]
            # END for (k)
        # END for (j)
    # END for (i)

    return etog


# END get_etog_matrix()

class ETOG(object):
    """ Definition of ETOG class which allows caching of ETOG value """

    # global etog matrix shared by all class instances.
    # This way we only need to make the matrix once.
    _etog = None

    def __init__(self):
        """ Construct etog matrix if it isn't already cached """
        if self._etog is None:
            self._etog = get_etog_matrix()

    @property
    def value(self):
        return self._etog


# Here follows translation of c-code described below by MSB. Translation
# by vc 081219. Translation by jw 190819.
#  A portion of the GLAST Burst Monitor Flight Software, managed by
#  NASA Marshall Space Flight Center.
#   Written by and Copyright by Michael S. Briggs, 2005, an employee of the
#   University of Alabama in Huntsville, and a member of the National Space
#   Science and Technology Center.
#   Filename: SunLoc.c

#   Michael S. Briggs, PC test version: 2005 August 8;
#   RTEMS / SPARC version: 2005 August 25.


#   Unless otherwise stated, all equations, page and chapter numbers
#   are from Astronomical Algorithms by Jean Meeus, copyright 1991
#   by Willmann-Bell, Inc.

#   Caution: the equations in that book are in DEGREES, unlike most
#   of the GBM FSW, and the mathematical functions of C and Newlib !
#   Most of the variables in this routine follow the equations of
#   Meeus and have units of degrees -- these variables have "_deg"
#   appended to their names.

#   Also consulted: Practical Astronomy with your Calculator, 3rd ed.,
#   by Peter Duffett-Smith.

#   Most of the variables are double precision because we don't want to
#   loose precision over the long time span of the GLAST mission.
#   However, we use single precision trig functions so that we don't
#   have to include the double precision versions in the executable --
#   we don't want the executable to become larger by also including
#   double precision trig functions.   It makes sense that single
#   precision angles are good enough and testing with the PC
#   version confirmed this.

#   Testing (with the PC version) was done by comparing to the
#   RA and Dec values listed n The Astronomical Almanac for the year 2000.
#   A special test version was used in which the input was days
#   (of 86400 s) (double type) rather than seconds.
#   The input values ranged from
#   0.0  for 2000 Dec 32 = Julian Date 2451910.5 (= GLAST Epoch !)
#   to -367.0 for 2000 Jan 0 = Julian Date 2451543.5
#   The worst case disagreement was approx
#   several seconds of RA (i.e., few x E-5 deg) and
#   25 arc seconds of dec (i.e., 7E-3 deg).
#   This test did NOT test UT vs Dynamical Time since it implicitly
#   assumed that the input times were Dynamical Times, matching
#   the times of the table in the Astronomical Almanac.
#   See comments below about why the difference doesn't matter for
#   the accuracy needed for GBM (UT vs GPS vs Dynamical are all close
#   enough).

def sun_loc(sc_time_sec):
    """ Compute the location of the sun in equatorial (ra, dec) coordinates

    NOTE: We calculate JulianDay as an offset from the Julian Day of the
    GLAST Epoch: JD_OF_GLAST_EPOCH is defined in gbm_central.h -- comments
    there describe the calculation of its value.

    Parameters
    ----------
    sc_time_sec : int64
        Spacecraft time in MET seconds

    Returns
    -------
    sun_ra : float32
        Right ascension of the sun in radians
    sun_dec : float32
        Declination of the sun in radians
    """
    jd_of_glast_epoch = np.float32(2451910.50)
    sec_per_day = np.float32(86400.00)
    jd = np.float64(jd_of_glast_epoch + np.float32(sc_time_sec) / sec_per_day)

    # We neglect the correction from S/C Time (GPS Time ??) to
    # Dynamical Time because it doesn't matter for the accuracy
    # that we need.  The Sun (really the earth) only moves 360 degrees
    # in one year, and the difference between UT or GPS times
    # and Dynamical or Ephemeris Times will be roughly 60 s
    # during the GLAST era.  The difference arises because the
    # rotation of the Earth is slowing down -- UT is adjusted as
    # the earth slows down, while dynamical (aka ephemeris)
    # advances continuously and is more suitable for astronomical
    # calculations.

    # Equation 24.1, Julian centuries of 36525 ephemeris days from
    # J2000.0 (2000 January 1.5 TD):
    t_cen = np.float64((jd - np.float32(2451545.0)) / np.float32(36525.0))

    # Equation 24.2: Geometric Mean [ecliptic] Longitude of the Sun,
    # in degrees.  "mean" means pretending a circular orbit.
    long_mean_deg = np.float32(280.46645) + np.float32(36000.76983) * t_cen + np.float32(0.0003032) * t_cen * t_cen

    # Equation 24.3: Mean anomaly of Sun in degress.
    # "anomaly" is the angle the sun has moved.
    # [mean ==> viewed from the center of the putatively circularly orbit.]
    # ["anomaly" is referenced to perigee.]
    mean_anomaly_deg = np.float32(357.52910) + np.float32(35999.05030) * t_cen - np.float32(0.0001559) * t_cen ** 2 - \
                       np.float32(4.8E-7) * t_cen ** 3

    # Equation 24.4: eccentricity of the Earth's orbit:
    eccen = np.float32(0.016708617) - np.float32(0.000042037) * t_cen - np.float32(1.236E-7) * t_cen * t_cen

    # The equation of the center -- approx solution of ellipical orbit:
    equation_center_deg = (
            (np.float32(1.914600) - np.float32(0.004817) * t_cen - np.float32(0.000014) * t_cen ** 2) * np.sin(
        mean_anomaly_deg / legacy_dtorad)
            + (np.float32(0.019993) - np.float32(0.000101) * t_cen) * np.sin(
        np.float32(2.0) * mean_anomaly_deg / legacy_dtorad)
            + np.float32(0.000290) * np.sin(np.float32(3.0) * mean_anomaly_deg / legacy_dtorad))

    # The Sun's True Longitude:
    true_long_deg = long_mean_deg + equation_center_deg

    # Since we want J2000 coordinates, we do not apply the corrections
    # for nutation and aberration to obtain apparent longitude.

    # The ecliptic latitude of the Sun, referenced to the Ecliptic of the date !,
    # never exceeds 1.2 arc secs, which we can surely neglect.

    # To convert from Ecliptic to Equatorial Coordinates, we need the
    # obliquity of the ecliptic, eq. 21.2 (but converted to decimal degrees):
    obliquity_ecliptic_deg = np.float32(23.4392911) - np.float32(46.8150) / np.float32(3600.0) * t_cen \
                             - np.float32(0.00059) / np.float32(3600.0) * t_cen ** 2 \
                             + np.float32(0.001813) / np.float32(3600.0) * t_cen ** 3

    # The one bad case for atan2 should be impossible -- the sun
    # cannot be at either pole of the Ecliptic coordinate system:
    # Inverse tangent of ( cos (obliquity_ecliptic_deg) * sin (TrueLong_deg) ) /  cos (TrueLong_deg) :
    sun_ra = np.float32(np.arctan2(np.cos(obliquity_ecliptic_deg / legacy_dtorad) *
                                   np.sin(true_long_deg / legacy_dtorad), np.cos(true_long_deg / legacy_dtorad)))

    # Change range of RA from -PI to +PI to 0 to 2PI:
    if sun_ra < np.float64(0.):
        sun_ra = np.float32(sun_ra + np.float64(2.) * legacy_rpi)

    sin_Dec_2000 = np.float32(np.sin(obliquity_ecliptic_deg / legacy_dtorad) * np.sin(true_long_deg / legacy_dtorad))

    if sin_Dec_2000 > np.float64(1.000) and sin_Dec_2000 < np.float64(1.001):
        sin_Dec_2000 = np.float32(np.float64(-1.0))
    if sin_Dec_2000 < np.float64(-1.000) and sin_Dec_2000 > np.float64(-1.001):
        sin_Dec_2000 = np.float32(np.float64(-1.0))

    sun_dec = np.float32(np.arcsin(sin_Dec_2000))

    return sun_ra, sun_dec


# END sun_loc()

def deadtime_correct(crange, mrates, brates, sduration, bgduration,
                     energies, verbose=False):
    """ Compute deadtime corrected rates

    Parameters
    ----------
    crange : np.ndarray(2, int32)
        Array of energy channel IDs chosen for localization
    mrates : np.ndarry(ndet * nen, int32)
        Array of measured detector counts in each energy channel
    brates : np.ndarry(ndet * nen, int32)
        Array of estimated background counts in each energy channel
    sduration : float32
        Timescale of measured detector counts
    bgduration : float32
        Timescale of estimated background counts
    energies : np.ndarray(nen + 1, float32)
        Array with end points of energy channel bins
    verbose : bool
        Print things to screen when True

    Returns
    -------
    c_mrates : np.ndarray(ndet, int32)
        Measured detector counts summed over chosen energy channels
    c_brates : np.ndarray(ndet, int32)
        Estimated background counts summed over chosen energy channels
    cenergies : np.ndarray(2, float32)
        Array with first and last index of chosen energy channel range
    usedet : np.ndarray(ldet, int32)
        Array determining whether we use detector in localization.
        0=Not Used, 1=Used
    maxdet : int
        Array index for detector with the largest count excess
        relative to background during the signal duration 
    data_signif : float32
        Gaussian significance estimate of this signal detection
    deadtime : np.ndarray(ndet, float32)
        Array of cumulative deadtimes for each detector
    """
    nen = energies.size - 1
    ndet = mrates.size // nen

    if nen * ndet != mrates.size:
        raise ValueError("Rates not provided for all energy bins")

    if verbose:
        print(" DEADTIME: ")

    deadtime = np.zeros(ndet, np.float32)
    bdeadtime = np.zeros(ndet, np.float32)
    for i in range(ndet):
        for j in range(nen):
            if j == (nen - 1):
                factor = 1.0e-5
            else:
                factor = 2.6e-6
            deadtime[i] += mrates[i][j] * factor
            bdeadtime[i] += brates[i][j] * factor
        # END for (j)
        deadtime[i] = 1. - deadtime[i]
        bdeadtime[i] = 1. - bdeadtime[i]

        if verbose:
            print("%12d %12.9f %16.9f" % (i + 1, deadtime[i], bdeadtime[i]))

        mrates[i] = mrates[i] / deadtime[i]
        brates[i] = brates[i] / bdeadtime[i]

    # END for (i)

    mrates = np.int32(sduration * mrates.astype(np.float32))
    brates = np.int32(sduration * brates.astype(np.float32))

    c_mrates, c_brates, cenergies = choose_energy_data(
        ndet, nen, crange, mrates, brates, energies)
    c_diff_rates = c_mrates - c_brates

    detmask = np.float32(c_brates[0:ndet - 2]) / bgduration > 100
    usedet = np.where(detmask == True)[0].astype(np.int32)

    i = c_diff_rates[:12].argmax()
    data_signif = np.float32(c_diff_rates[i]) / np.sqrt(np.float32(c_brates[i]))

    diff_rates = mrates - brates
    sratio = np.float32(diff_rates[i][1] + diff_rates[i][2]) \
             / np.float32(diff_rates[i][3] + diff_rates[i][4])

    if verbose:
        print(" Rates: ")
        print("%12d" * c_mrates.size % tuple(c_mrates))
        print("  ")
        print("%12d" * c_brates.size % tuple(c_brates))
        print("  ")
        print("%12d" * c_diff_rates.size % tuple(c_diff_rates))
        print("  ")
        print(" Energies: %13.7f%18.7f" % (cenergies[0], cenergies[1]))
        print("  ")
        print(" duration   {0:.9}".format(sduration))
        print(" data signif {0:11} {1:.9}".format(i + 1, data_signif))
        print(" Softness ratio %13.9f" % sratio)

    return c_mrates, c_brates, cenergies, usedet, i, data_signif, deadtime

# END deadtime_correct()
