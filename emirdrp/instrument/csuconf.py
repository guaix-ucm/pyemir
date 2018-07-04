#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division

import enum

from numina.array.bbox import BoundingBox

from emirdrp.core import EMIR_NBARS


class CSUConf(object):
    """Information about the configuration of slits in the CSU"""
    def __init__(self):
        self.name = 'CSU'
        self.conf_id = 'v1'
        self.conf_f = 'UNKNOWN'
        self.slits = {}
        # Indices
        self.lbars = list(range(1, 55 + 1))
        self.rbars = [(i + 55) for i in self.lbars]
        self.bars = {}
        # CSU is open if
        self.L_LIMIT = 6
        self.R_LIMIT = 2057

    def is_open(self):
        lopen = all(self.bars[barid].x2 <= self.L_LIMIT for barid in self.lbars)
        ropen = all(self.bars[barid].x1 >= self.R_LIMIT for barid in self.rbars)

        return lopen and ropen

    def is_closed(self):
        return not self.is_open()


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
        self.slit_h_virt = 16.242
        self.y0 = yc
        self.y1 = yc - self.slit_h_virt-3
        self.y2 = yc + self.slit_h_virt+3
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
    """Slit formed from combination of PysicalBarL  and PhysicalBarR"""
    def __init__(self, idx, lbars, rbars, target_type=TargetType.UNKNOWN):
        self.target_type = target_type
        self.target_coordinates = (-100, -100)
        self.target_coordinates_v = (-100, -100)
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


def create_bar_models(barstab):
    bars = {}
    for params in barstab:
        barid = int(params[0])
        x0 = params[3]
        xdelt = params[2]
        y0 = params[1]
        if barid > 55:
            bar = CSUBarModelR(barid, x0, xdelt, y0)
        else:
            bar = CSUBarModelL(barid, x0, xdelt, y0)
        bars[barid] = bar
    return bars


def read_csu_2(hdr, barmodel):
    """Read CSU information and slits from header"""
    conf = CSUConf()
    conf.name = 'NAME'
    conf.conf_id = 'v1'
    conf.conf_f = hdr.get('CSUCONFF', 'UNKNOWN')
    conf.bars = barmodel
    # Read CSUPOS and set position in model
    # The bars
    for idx in conf.bars:
        key = "CSUP{}".format(idx)
        # UNIT of CSUPOS?
        # who knows

        # set CSUPOS for BAR
        # Using the model, this sets the X,Y coordinate
        conf.bars[idx].csupos = hdr[key]
        # print('read_csu {}, position is {}'.format(idx, conf.bars[idx].current_pos))

    # Slits. We dont have keywords to fuse slitlets into larger logical slits
    # so there is one slitslet for each pair of bars
    for idx in range(1, 55 + 1):
        l_ids = [idx]
        r_ids = [idx + 55]

        lbars = {lid: conf.bars[lid] for lid in l_ids}
        rbars = {rid: conf.bars[rid] for rid in r_ids}

        this = LogicalSlit(idx, lbars=lbars, rbars=rbars)
        # References from header
        try:
            slit_t = hdr["SLIFL%d" % idx]
            this.target_type = TargetType(slit_t)
        except KeyError as err:
            print('warning', err)
            this.target_type = TargetType.UNKNOWN

        xref = hdr.get("XRSLI%d" % this.idx, -100) - 1
        yref = hdr.get("YRSLI%d" % this.idx, -100) - 1
        this.target_coordinates = (xref, yref)

        xref = hdr.get("XVSLI%d" % this.idx, -100) - 1
        yref = hdr.get("YVSLI%d" % this.idx, -100) - 1
        this.target_coordinates_v = (xref, yref)
        conf.slits[idx] = this

    return conf