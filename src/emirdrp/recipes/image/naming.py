#
# Copyright 2011-2023 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0-or-later
# License-Filename: LICENSE.txt
#

"""

Naming intermediate images

"""


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
