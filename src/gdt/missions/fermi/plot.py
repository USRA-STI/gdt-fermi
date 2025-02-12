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
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path

from gdt.core.plot.plot import EarthPoints, PlotElement, GdtCmap
from gdt.core.plot.earthplot import EarthPlot
from .mcilwainl import calc_mcilwain_l

__all__ = ['FermiEarthPlot', 'FermiIcon', 'McIlwainL', 'mcilwain_map']

class FermiIcon(EarthPoints):
    """Plot a Fermi icon on the Earth.

    Parameters:
        lat (np.array): The latitude value
        lon (np.array): The longitude value
        proj (GeoAxesSubplot): The Cartopy projection
        alpha (float, optional): The alpha opacity
        **kwargs: Other plotting options
    """
    def __init__(self, lat, lon, proj, alpha=1.0, **kwargs):
        self._norm_width = 31.4 * 2.0
        self._norm_height = 7.0 * 2.0
        super().__init__(lat, lon, proj, color=None, alpha=alpha, **kwargs)

    @property
    def lat(self):
        """(np.array): Normalized plot coordinates for the LAT"""
        return self._normalize(self._lat())

    @property
    def gbm_side(self):
        """(np.array): Normalized plot coordinates for the GBM side panel"""
        return self._normalize(self._gbm_side())

    @property
    def left_panel(self):
        """(np.array): Normalized plot coordinates for left panel + solar array"""
        return self._normalize(self._left_panel())

    @property
    def right_panel(self):
        """(np.array): Normalized plot coordinates for right panel + solar array"""
        return self._normalize(self._right_panel())

    @property
    def antenna(self):
        """(np.array): Normalized plot coordinates for antenna"""
        return self._normalize(self._antenna())

    def _create(self, lat, lon, proj):
    
        lon = np.asarray(lon)
        lat = np.asarray(lat)
        
        lon[(lon > 180.0)] -= 360.0
        x, y = (lon, lat)
        z = 10
        factor = 50.
        fermilat = self.lat * factor
        fermilat[:, 0] += x
        fermilat[:, 1] += y
        path1 = Path(fermilat, closed=True)
        patch1 = patches.PathPatch(path1, facecolor='#DCDCDC', zorder=z)
        proj.add_patch(patch1)

        gbm = self.gbm_side * factor
        gbm[:, 0] += x
        gbm[:, 1] += y
        path2 = Path(gbm, closed=True)
        patch2 = patches.PathPatch(path2, facecolor='#B79241', zorder=z)
        proj.add_patch(patch2)

        panel = self.left_panel * factor
        panel[:, 0] += x
        panel[:, 1] += y
        path3 = Path(panel, closed=True)
        patch3 = patches.PathPatch(path3, facecolor='#45597C', zorder=z)
        proj.add_patch(patch3)

        panel = self.right_panel * factor
        panel[:, 0] += x
        panel[:, 1] += y
        path4 = Path(panel, closed=True)
        patch4 = patches.PathPatch(path4, facecolor='#45597C', zorder=z)
        proj.add_patch(patch4)

        antenna = self.antenna * factor
        antenna[:, 0] += x
        antenna[:, 1] += y
        path5 = Path(antenna, closed=True)
        patch5 = patches.PathPatch(path5, facecolor='#546165', zorder=z)
        proj.add_patch(patch5)

        return [patch1, patch2, patch3, patch4, patch5]

    def _normalize(self, pts):
        return (pts / self._norm_width)

    def _lat(self):
        pts = [[-2.5, 3.5], [-2.5, 1.2], [2.5, 1.2], [2.5, 3.5], [-2.5, 3.5]]
        pts = np.array(pts)
        return pts

    def _gbm_side(self):
        pts = [[-2.5, 1.2], [-2.5, -2.1], [2.5, -2.1], [2.5, 1.2], [-2.5, 1.2]]
        pts = np.array(pts)
        return pts

    def _left_panel(self):
        pts = [[-2.5, -1.0], [-4.5, -2.5], [-15.7, -2.5], [-15.7, 0.5],
               [-4.5, 0.5], [-2.5, -1.0]]
        pts = np.array(pts)
        return pts

    def _right_panel(self):
        pts = [[2.5, -1.0], [4.5, -2.5], [15.7, -2.5], [15.7, 0.5],
               [4.5, 0.5], [2.5, -1.0]]
        pts = np.array(pts)
        return pts

    def _antenna(self):
        pts = [[0.5, -2.1], [0.5, -3.5], [1.5, -3.5], [1.5, -2.1], [0.5, -2.1]]
        pts = np.array(pts)
        return pts


class McIlwainL(PlotElement):
    """Plot class for the McIlwain L heatmap.

    Parameters:
        lat_range (float, float): The latitude range
        lon_range (float, float): The longitude range
        proj (Cartopy Projection): The Cartopy projection
        ax (:class:`matplotlib.axes`): The axis on which to plot
        colorbar (bool, optional): If True, create a colorbar for the heatmap. 
                                   Default is True
        color (:class:`~gdt.plot.plot.GdtCmap`): The colormap of the heatmap. 
                                  Default is 'viridis'
        alpha (float, optional): The alpha opacity of the heatmap
        norm (:class:`matplotlib.colors.Normalize` or similar, optional):
            The normalization used to scale the colormapping to the heatmap 
            values. This can be initialized by the defined matplotlib 
            normalizations or a custom normalization.
        levels (list, optional): The number of plot levels.  If not set,
                                 the default is to plot from 0.9 to 1.7 in 
                                 increments of 0.1.
        **kwargs: Other plotting options
    """
    def __init__(self, lat_range, lon_range, proj, colorbar=True, 
                 color=GdtCmap('viridis'), alpha=None, norm=None, 
                 levels=None, **kwargs):
        
        self._colorbar = colorbar
        self._norm = norm
        self._levels = levels
        if self._levels is None:
            self._levels = [0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7]
        
        super().__init__(color=color, alpha=alpha)
        self._kwargs = kwargs
        artists = self._create(lat_range, lon_range, proj)
        self._artists = self._sanitize_artists(artists)

        # set the colormap
        self.color = color
        self.color.set_callback(self._artists[0].changed)
        
    @property
    def alpha(self):
        """(float): The alpha opacity value"""
        return self._alpha
    @alpha.setter
    def alpha(self, alpha):
        self._artists[0].set_alpha(alpha)
        if len(self._artists) == 2:
            self._artists[1].set_alpha(alpha)
            self._artists[1]._draw_all()
        self._alpha = alpha

    @property
    def color(self):
        """(:class:`~gdt.plot.plot.GdtCmap`): The colormap"""
        return self._color
    @color.setter
    def color(self, color):
        if not isinstance(color, GdtCmap):
            raise TypeError('color must be of type GdtCmap')
        self._color = color
        self._artists[0].set_cmap(color.cmap)

    @property
    def colorbar(self):
        """(matplotlib.colorbar.Colorbar): The colorbar object"""
        if self._colorbar:
            return self._artists[-1]

    @property
    def levels(self):
        """(list): The contour levels"""
        # mark TODO: Add a setter, which will have to do a replot
        return self._levels
    
    @property
    def norm(self):
        """(matplotlib.colors.Normalize or similar): The colormap normalization"""
        return self._norm
    @norm.setter
    def norm(self, norm):
        self._artists[0].set_norm(norm)
        self._norm = norm

    def _create(self, lat_range, lon_range, proj):
        artists = mcilwain_map(lat_range, lon_range, proj, 
                               cmap=self._color.cmap, alpha=self._alpha, 
                               norm=self._norm, levels=self._levels, 
                               **self._kwargs)
        artists = [artists]
        if self._colorbar:
            artists.append(self._make_colorbar(proj, artists[0]))

        return artists

    def _make_colorbar(self, ax, artist):
        cb = plt.colorbar(artist, label='McIlwain L', ax=ax, shrink=0.6,
                          pad=0.2, orientation='horizontal')
        cb._draw_all()
        return cb

    def __repr__(self):
        spaces = ' ' * 12
        s = "<McIlwainL: color='{0}';\n{1}".format(self.color.name, spaces) 
        s += 'alpha={0};\n{1}'.format(self.alpha, spaces)
        s += 'num_contours={0};\n{1}'.format(len(self.levels), spaces)
        s += 'colorbar={0}>'.format(self._colorbar)
        return s


def mcilwain_map(lat_range, lon_range, proj, saa_mask=None, cmap=None,
                 alpha=0.5, num_lat_points=108, num_lon_points=720,
                 levels=None, **kwargs):
    """Plot a McIlwain L heatmap on the Earth.
    
    Args:
        lat_range (float, float): The latitude range
        lon_range (float, float): The longitude range
        proj (Cartopy Projection): The Cartopy projection
        ax (matplotlib.axes): The plot axes references
        saa_mask (:class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`, optional)
            Mask out the SAA from the heatmap
        color (str, optional): The color of the heatmap
        alpha (float, optional): The alpha opacity of the heatmap
        kwargs (optional): Other plotting keywords
    
    Returns:
        (matplotlib.contour.QuadContourSet)
    """
    # do create an array on the earth
    lat_array = np.linspace(*lat_range, num_lat_points)
    lon_array = np.linspace(*lon_range, num_lon_points)
    LAT, LON = np.meshgrid(lat_array, lon_array)

    # mcilwain l over the grid
    mcl = calc_mcilwain_l(LAT, LON)

    # if we want to mask out the SAA
    if saa_mask is not None:
        xy = list(zip(saa_mask.longitude, saa_mask.latitude))
        saa_path = patches.Polygon(xy).get_path()
        mask = saa_path.contains_points(np.array((LON.ravel(), LAT.ravel())).T)
        mcl[mask.reshape(mcl.shape)] = 0.0

    # do the plot  
    image = proj.contourf(LON, LAT, mcl, levels=levels, alpha=alpha, cmap=cmap, 
                          **kwargs)
    return image


class FermiEarthPlot(EarthPlot):
    """Class for plotting Fermi's orbit, including the McIlwain L heatmap.
    
    Note:
        This class requires installation of Cartopy.
    
    Parameters:
        saa (:class:`~gdt.core.geomagnetic.SouthAtlanticAnomaly`, optional): 
            If set, displays the SAA polygon.
        mcilwain (optional, bool): Set to True if plotting the McIlwain L 
                                   heatmap, otherwise False.
        **kwargs: Options to pass to :class:`~gdt.plot.plot.GdtPlot`        
    """
    def __init__(self, saa=None, mcilwain=True, **kwargs):
        lat_range = (-30.00, 30.00)
        lon_range = (-180.0, 180.0)
        super().__init__(lat_range=lat_range, lon_range=lon_range, saa=saa, 
                         **kwargs)
        self._trig_mcilwain = None
        self._mcilwain = None

        if mcilwain:
            self._mcilwain = McIlwainL(lat_range, lon_range, self._m, 
                                       alpha=0.5, saa_mask=saa)

    @property
    def mcilwainl(self):
        """(:class:`McIlwainL`): The McIlwain L heatmap"""
        return self._mcilwain

    def add_spacecraft_frame(self, *args, **kwargs):
        super().add_spacecraft_frame(*args, icon=FermiIcon, **kwargs)
    
    def standard_title(self):
        """Add a standard plot title containing orbital position and McIlwain L
        """
        if self.spacecraft is not None:
            coord = self.spacecraft.coordinates
            title = 'Latitude, East Longitude: ({0}, {1})\n'.format(*coord)
            lat = float(coord[0][:-1]) * (-1 if "S" in coord[0] else 1)
            lon = float(coord[1][:-1]) * (-1 if "W" in coord[1] else 1)
            mcl = calc_mcilwain_l(lat, lon)
            title += 'McIlwain L: {:.2f}'.format(mcl)
            self._m.set_title(title)
