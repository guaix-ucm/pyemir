#
# Copyright 2011 Sergio Pascual
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

def redimensioned_images(label, iteration, ext='.fits'):
    dn = '%s_r_i%01d%s' % (label, iteration, ext)
    mn = '%s_mr_i%01d%s' % (label, iteration, ext)
    return dn, mn

def object_mask(label, iteration, ext='.fits'):
    return '%s_mro_i%01d%s' % (label, iteration, ext)

def skyflat_proc(label, iteration, ext='.fits'):
    dn = '%s_rf_i%01d%s' % (label, iteration, ext)
    return dn

def skybackground(label, iteration, ext='.fits'):
    dn = '%s_sky_i%01d%s' % (label, iteration, ext)
    return dn

def skybackgroundmask(label, iteration, ext='.fits'):
    dn = '%s_skymask_i%01d%s' % (label, iteration, ext)
    return dn


def skysub_proc(label, iteration, ext='.fits'):
    dn = '%s_rfs_i%01d%s' % (label, iteration, ext)
    return dn

def skyflat(label, iteration, ext='.fits'):
    dn = 'superflat_%s_i%01d%s' % (label, iteration, ext)
    return dn

def segmask(iteration, ext='.fits'):
    return "check_i%01d%s" % (iteration, ext)

