#
# Copyright 2008-2011 Sergio Pascual
# 
# This file is part of PyEmir
# 
# PyEmir is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# PyEmir is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with PyEmir.  If not, see <http://www.gnu.org/licenses/>.
# 

'''Recipe for the reduction of multiobject spectroscopy.

Recipe to reduce observations obtained in multiobject spectroscopy,
considering different possibilities depending on the observing mode. In
particular, the following observing modes are considered: stare spectra,
dithered/nodded spectra, and offset spectra.

The special flexibility of the Cold Slit Unit (CSU) allows to obtain slit
configurations in which some of the individual slits can be perfectly aligned
with some other slits. When these aligned slits are contiguous, they constitute
what are going to be considered as **pseudo-long-slits**. In the particular
case in which a given individual slit is not aligned with any other, the
pseudo-long-slit will have the shortest length *Lslit*, the length of a single
individual slit. In the rest of the cases, where the number of individual slits
that are aligned consecutively is *Nslit*>1, the total length of the
pseudo-long-slit will be *Nslit* x *Lslit*.

In order to unambiguously proceed with a proper sky subtraction for different
observing modes, a table specifying which pseudo-long-slits with *sky*
spectra are associated to the pseudo-long-slits with *science* spectra must be
provided. In addition, this table must specify the characteristics of the data
obtained through each pseudo-long-slit. More precisely, in the case of
*science* spectra, it must indicate whether the target is puntual or extended.
The *sky* spectra must indicate the corresponding position and range (pixels
along the slit, referred to the central pseudo-long-slit coordinates) of the
valid sky regions. In this way, both *sky* and *science* spectra can be
extracted from the same pseudo-long-slit, if this is required.

The following possibilities must be considered at the time of creating the
association table among pseudo-long-slits:

 * A single *sky* pseudo-long-slit can be associated to several *science*
   pseudo-long-slits.
 * The association between pseudo-long-slits is not restricted to
   pseudo-long-slits of the same spectroscopic image. For example, in the case
   of separate exposures for *science* and *sky* spectra, a given
   pseudo-long-slit can contain the *science* spectra in one exposure and the
   *sky* spectra in another exposure.
 * The configuration of the CSU must be identical when pseudo-long-slits of
   different images are associated.

**Observing modes:**

 * Stare spectra
 * Dithered/Nodded spectra
 * Offset Spectra 

**Inputs:**

 * Science frames + [Sky Frames]
 * An indication of the observing mode: **stare spectra**, 
   **dithered/nodded spectra**, or **offset spectra**. In the case of
   dithered/nodded spectra, the corresponding offsets must be indicated (these
   offsets must be integer and can be zero)
 * The table relating the associations among the different pseudo-long-slits
   and the *science* and *sky* intervals that will be extracted at the end of
   the data reduction (TBD if this table must be in the FITS header and/or in
   other format)
 * Master Dark 
 * Bad pixel mask (BPM) 
 * Non-linearity correction polynomials 
 * Spectroscopic Master flat (twilight/dome flats)
 * Master background (thermal background, only in K band)
 * Exposure Time (must be the same in all the frames, TBD)
 * Airmass for each frame
 * Detector model (gain, RN, lecture mode)
 * Average extinction curve the spectral range
 * Response curve for each pseudo-long-slit.

**Outputs:**

 * 2D image with two extensions: final 2D image and associated 2D variance
   image, before sky-subtaction (TBD)
 * Final sky-subtracted and extracted spectra. TBD the way in which these
   spectra (and associated variance spectra) will be arranged (single 2D image,
   individual 1D spectra in different FITS extensions, individual 1D FITS files
   for each pseudo-long-slit)

**Procedure:**

The reduction will be carried out by extracting the different pseudo-long-slits
and following a *traditional* long-slit reduction, making the appropriate
treatment for the sky subtraction. The basic steps must include:

 * Data modeling (if appropriate/possible) and variance frame creation (from
   first principles)

 * Correction of non-linearity

 * Dark correction: this step can be avoided in most cases, since the
   extracted sky spectra may contain that information. However, if observing
   modes with different exposure times for *science* and *sky* frames are
   allowed, this correction may be necessary.

 * Flatfielding: distinguish between high frequency (pixel-to-pixel) and 
   low-frequency (overall response and slit illumination) corrections. 
   Lamp flats are adequate for the former and twilight flats for the second.

 * Detection and extraction of slits: apply Border_Detection algorithm, from
   own frames or from flatfields.

 * Cleaning

   * Single spectroscopic image: sigma-clipping algorithm removing 
     local background in pre-defined direction(s).
   * Multiple spectroscopic images: sigma-clipping from comparison between 
     frames.
 
 * Wavelength calibration and C-distortion correction of each slit. 
   Double-check with available sky lines.

 * Sky-subtraction (number of sources/slit will be allowed to be > 1?).

   * Subtraction using sky signal at the borders of the same slit.
   * Subtraction using sky signal from other(s) slit(s), not necessarily 
     adjacent.
 
 * Spectrophotometric calibration of each slit, using the extinction correction 
   curve and the master spectrophotometric calibration curve.

 * Spectra extraction: define optimal, average, peak, FWHM.

'''

