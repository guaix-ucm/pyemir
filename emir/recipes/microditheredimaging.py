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

'''Recipe for the reduction of microdithering imaging.

Recipe to reduce observations obtained in imaging mode with microdithering.
A critical piece of information
here is a table that clearly specifies which images can be labelled as
*science*, and which ones as *sky*. Note that some images are used both as
*science* and *sky* (when the size of the targets are small compared to the
offsets).

**Observing modes:**
 * Micro-dithered images 

**Inputs:**

 * Science frames + [Sky Frames]
 * A table relating each science image with its sky image(s)  (TBD if it's in the FITS header and/or in other format)
 * Offsets between them
 * Master Dark 
 * Bad pixel mask (BPM) 
 * Non-linearity correction polynomials 
 * Master flat (twilight/dome flats)
 * Master background (thermal background, only in K band)
 * Exposure Time (must be the same in all the frames)
 * Airmass for each frame
 * Detector model (gain, RN, lecture mode)
 * Average extinction in the filter
 * Astrometric calibration (TBD)

**Outputs:**

 * Image with three extensions: final image scaled to the individual exposure
   time, variance  and exposure time map OR number of images combined (TBD).

**Procedure:**

Images are regridded to a integer subdivision of the pixel and then they are
corrected from dark, non-linearity and flat. It should be desirable that the
microdithering follows a pattern that can be easily translated to a subdivision
of the pixel size (by an integer *n* = 2, 3, 4,...) that does not requires a
too high *n* value. An iterative process starts:

 * Sky is computed from each frame, using the list of sky images of each
   science frame. The objects are avoided using a mask (from the second
   iteration on).

 * The relatiev offsets are the nominal from the telescope. From the second
   iteration on, we refine them using bright objects.

 * We combine the sky-subtracted images, output is: a new image, a variance
   image and a exposure map/number of images used map.

 * An object mask is generated.

 * We recompute the sky map, using the object mask as an additional input. From
   here we iterate (typically 4 times).

 * Finally, the images are corrected from atmospheric extinction and flux
   calibrated.

 * A preliminary astrometric calibration can always be used (using the central
   coordinates of the pointing and the plate scale in the detector). A better
   calibration might be computed using available stars (TBD).

'''
