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

'''Detector Nonlinearity Recipe.

Recipe to calibrate the detector nonlinearity.

**Observing modes:**

 * Not considered

**Inputs:**

 * A list of dome flat fields with the adequate exposure times.

**Outputs:**

 * The polynomial in a given format (TBD)

**Procedure:**

The median of the reduced flatfields is computed. For each exposure time, the
average and standard deviation of the ADUs is computed. A linear fit is
computed for the lowest ADU counts, with standard deviations as weights. A new
fit is performed, this time of ADU_linear/ADU_obs, versus ADU_obs. The
polynomial terms of the fit are computed.

'''

