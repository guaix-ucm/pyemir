#
# Copyright 2008-2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division

import enum
import numpy

from numina.array.bbox import BoundingBox

EMIR_NBARS = 55

class CSUConf(object):
    """Information about the configuration of slits in the CSU"""
    def __init__(self, barmodel):
        self.name = 'CSU'
        self.conf_id = 'v1'
        self.conf_f = 'UNKNOWN'
        self.slits = {}
        # Indices
        self.LBARS = list(range(1, EMIR_NBARS + 1))
        self.RBARS = [(i + EMIR_NBARS) for i in self.LBARS]
        self.NBARS = EMIR_NBARS
        self.bars = barmodel
        # CSU is open if
        self.L_LIMIT = 6
        self.R_LIMIT = 2057
        # Generate default slits, without header information
        self._update_slits({})

    def is_open(self):
        lopen = all(self.bars[barid].x2 <= self.L_LIMIT for barid in self.LBARS)
        ropen = all(self.bars[barid].x1 >= self.R_LIMIT for barid in self.RBARS)

        return lopen and ropen

    def is_closed(self):
        return not self.is_open()

    def set_state_from_img(self, img):
        """Set state using a FITS image"""

        # FIXME: CSUCONFF is in primary
        hdr0 = img[0].header
        self.conf_f = hdr0['CSUCONFF']

        if 'MECS' in img:
            # Get header from extension
            hdr_csu = img['MECS'].header
        else:
            # Get header from main header
            hdr_csu = img[0].header
        self.set_state(hdr_csu)

    def set_state(self, hdr):
        """Read CSU information and slits from header"""

        # FIXME: this keyword in only in primary
        # it should be duplicated in MECS
        try:
            self.conf_f = hdr['CSUCONFF']
        except KeyError:
            pass
        # Read CSUPOS and set position in model
        # The bars
        for idx in self.bars:
            key = "CSUP{}".format(idx)
            # UNIT is mm

            # set CSUPOS for BAR
            # Using the model, this sets the X,Y coordinate
            self.bars[idx].csupos = hdr[key]

        self._update_slits(hdr)

    def _update_slits(self, hdr):
        """Recreate slits"""

        # Slits.
        # For the moment we will fuse only reference slits
        OUTPOS = -100
        # clean existing slits
        self.slits = {}
        mm = []
        for idx in self.LBARS:

            # References from header
            try:
                slit_t = hdr["SLIFL%d" % idx]
                target_type = TargetType(slit_t)
            except KeyError:
                target_type = TargetType.UNKNOWN

            xref = hdr.get("XRSLI%d" % idx, OUTPOS) - 1
            yref = hdr.get("YRSLI%d" % idx, OUTPOS) - 1
            target_coordinates = (xref, yref)

            xref = hdr.get("XVSLI%d" % idx, OUTPOS) - 1
            yref = hdr.get("YVSLI%d" % idx, OUTPOS) - 1
            target_coordinates_v = (xref, yref)

            mm.append((idx, target_type, target_coordinates, target_coordinates_v))

        bag = merge_slits(mm, max_slits=3, tol=1e-2)

        for idx, l_ids in bag.items():
            r_ids = [lid + EMIR_NBARS for lid in l_ids]

            lbars = {lid: self.bars[lid] for lid in l_ids}
            rbars = {rid: self.bars[rid] for rid in r_ids}

            this = LogicalSlit(idx, lbars=lbars, rbars=rbars)
            # References from header
            ref = mm[l_ids[0] - 1]
            this.target_type = ref[1]
            this.target_coordinates = ref[2]
            this.target_coordinates_v = ref[3]
            self.slits[idx] = this


class TargetType(enum.Enum):
    """Possible targets in a slit"""
    UNASSIGNED = 0
    SOURCE = 1
    REFERENCE = 2
    SKY = 3
    UNKNOWN = 4


class CSUBarModel(object):
    def __init__(self, idx, xs, xdelt, yc, active=True):
        self.idx = idx
        self.xs = xs
        self.xdelt = xdelt
        # FIXME: pixels?
        self.slit_h_virt = 16.242
        self.y0 = yc
        self.y1 = yc - self.slit_h_virt - 3
        self.y2 = yc + self.slit_h_virt + 3
        # Height in virt pixels
        self._x1 = 0
        self._x2 = 2048
        self.active = active
        self.csupos = 0.0

    def position(self, csupos):
        # start in 1 position
        return self.xs + self.xdelt * csupos

    @property
    def current_pos(self):
        # start in 1 position
        return self.position(self.csupos)

    @property
    def xpos(self):
        # start in 1 position
        return self.position(self.csupos)

    @property
    def x2(self):
        return self._x2

    @property
    def x1(self):
        return self._x1

    def bbox(self):
        # origin is 1
        x1 = max(self.x1 - 1, 0)
        y1 = max(self.y1 - 1, 0)
        x2 = max(self.x2 - 1, 0, x1)
        y2 = max(self.y2 - 1, 0, y1)
        return BoundingBox.from_coordinates(x1, x2, y1, y2)


class CSUBarModelL(CSUBarModel):

    @property
    def x2(self):
        return self.position(self.csupos)


class CSUBarModelR(CSUBarModel):

    @property
    def x1(self):
        return self.position(self.csupos)


class PhysicalBar(object):
    def __init__(self, idx, xpos, y1, y2, active=True):
        self.idx = idx
        self.xpos = xpos
        self.y1 = y1
        self.y2 = y2
        self.active = active


class PhysicalBarL(PhysicalBar):
    def __init__(self, idx, xpos, y1, y2, active=True):
        super(PhysicalBarL, self).__init__(idx, xpos, y1, y2, active)


class PhysicalBarR(PhysicalBar):
    def __init__(self, idx, xpos, y1, y2, active=True):
        super(PhysicalBarR, self).__init__(idx, xpos, y1, y2, active)


class LogicalSlit(object):
    """Slit formed from combination of PhysicalBarL and PhysicalBarR"""
    def __init__(self, idx, lbars, rbars, target_type=TargetType.UNKNOWN):
        UNDEFPOS = -100
        self.target_type = target_type
        self.target_coordinates = (UNDEFPOS, UNDEFPOS)
        self.target_coordinates_v = (UNDEFPOS, UNDEFPOS)
        self.idx = idx
        self.lbars = lbars
        self.rbars = rbars
        # Bar ids
        lbars_ids = list(lbars.keys())
        rbars_ids = list(rbars.keys())
        lbars_ids.sort()
        rbars_ids.sort()

        # Must have the same size
        if len(lbars) != len(rbars):
            raise ValueError('len() are different')

        # bars must be paired
        for l, r in zip(lbars_ids, rbars_ids):
            if r != l + EMIR_NBARS:
                raise ValueError('not paired {} and {}'.format(l, r))

        for a, b in zip(lbars_ids, lbars_ids[1:]):
            if b != a + 1:
                raise ValueError('not contiguous {} and {}'.format(a, b))

        self.lbars_ids = lbars_ids
        self.rbars_ids = rbars_ids
        # BBOX
        id_first = self.lbars_ids[0]
        id_last = self.lbars_ids[-1]

        bbox_y1 = self.lbars[id_first].y1
        bbox_y2 = self.lbars[id_last].y2

        bbox_x1 = min(self.lbars.values(), key=lambda obk: obk.xpos).xpos
        bbox_x2 = max(self.rbars.values(), key=lambda obk: obk.xpos).xpos
        # Origin is 1
        self.bbox_int = BoundingBox.from_coordinates(bbox_x1-1, bbox_x2-1, bbox_y1-1, bbox_y2-1)

    def bbox(self):
        return self.bbox_int

    def width(self):
        bar_x1 = min(self.lbars.values(), key=lambda obk: obk.xpos)
        bar_x2 = max(self.rbars.values(), key=lambda obk: obk.xpos)
        return bar_x2.xpos - bar_x1.xpos


def create_bar_models(barstab):
    """Create a dictionary of models of bars using a table"""
    bars = {}
    for params in barstab:
        barid = int(params[0])
        x0 = params[3]
        xdelt = params[2]
        y0 = params[1]
        if barid > EMIR_NBARS:
            bar = CSUBarModelR(barid, x0, xdelt, y0)
        else:
            bar = CSUBarModelL(barid, x0, xdelt, y0)
        bars[barid] = bar
    return bars


def read_csu_from_header(barmodel, hdr):
    """Read CSU information and slits from header"""
    # FIXME: CSUCONFF is in primary
    csu_conf = CSUConf(barmodel)
    csu_conf.set_state(hdr)
    return csu_conf


def read_csu_from_image(barmodel, img):
    """Read CSU information and slits from FITS"""
    csu_conf = CSUConf(barmodel)
    csu_conf.set_state_from_img(img)
    return csu_conf


def merge_slits(mm, max_slits=3, tol=1e-2):
    """"Merge contiguous REFERENCE slits"""
    cx1 = 0
    bag = {}
    while True:
        t_id, t_type, t_coor, t_coor_v = mm[cx1]
        rest = mm[cx1 + 1: cx1 + max_slits]
        bag[t_id] = [t_id]
        cx2 = cx1 + 1
        if t_type == TargetType.REFERENCE:
            # For REFERENCE
            for r_id, r_type, r_coor, r_coor_v in rest:
                dis = numpy.hypot(t_coor[0] - r_coor[0], t_coor[1] - r_coor[1])
                if r_type == TargetType.REFERENCE and (dis < tol) :
                    bag[t_id].append(r_id)
                    cx2 += 1
                else:
                    break
        cx1 = cx2
        if cx1 >= len(mm):
            break
    return bag
