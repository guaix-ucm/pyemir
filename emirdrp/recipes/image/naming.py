#
# Copyright 2011-2014 Universidad Complutense de Madrid
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

'''

Naming intermediate images

'''


def name_redimensioned_frames(label, step, ext='.fits'):
    dn = '%s_r%s' % (label, ext)
    mn = '%s_mr%s' % (label, ext)
    return dn, mn


def name_object_mask(label, step, ext='.fits'):
    return '%s_mro_i%01d%s' % (label, step, ext)


def name_skybackground(label, step, ext='.fits'):
    dn = '%s_sky_i%01d%s' % (label, step, ext)
    return dn


def name_skybackgroundmask(label, step, ext='.fits'):
    dn = '%s_skymask_i%01d%s' % (label, step, ext)
    return dn


def name_skysub_proc(label, step, ext='.fits'):
    dn = '%s_rfs_i%01d%s' % (label, step, ext)
    return dn


def name_skyflat(label, step, ext='.fits'):
    dn = 'superflat_%s_i%01d%s' % (label, step, ext)
    return dn


def name_skyflat_proc(label, step, ext='.fits'):
    dn = '%s_rf_i%01d%s' % (label, step, ext)
    return dn


def name_segmask(step, ext='.fits'):
    return "check_i%01d%s" % (step, ext)
