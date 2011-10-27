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

'''Spectral Flatfield Recipe.

Recipe to process spectral flat-fields. The flat-on and flat-off images are
combined (method?) separately and the subtracted to obtain a thermal subtracted
flat-field.

**Observing modes:**

    * Multislit mask Flat-Field
     
**Inputs:**

 * A list of lamp-on flats
 * A model of the detector 

**Outputs:**

 * A combined spectral flat field with with variance extension and quality flag.

**Procedure:**

 * TBD

'''

