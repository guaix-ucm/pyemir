#
# Copyright 2008-2018 Universidad Complutense de Madrid
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


from __future__ import division


import sys
import enum
import math
import logging
import re

import numpy as np
import matplotlib.pyplot as plt
import sep
from numina.array.utils import coor_to_pix_1d
from numina.array.offrot import fit_offset_and_rotation
from numina.core import Requirement, Product, Parameter, RecipeError
from numina.core.products import ArrayType
from numina.core.requirements import ObservationResultRequirement
from numina.types.qc import QC

import emirdrp.datamodel as datamodel
from emirdrp.core import EmirRecipe, EMIR_NBARS, EMIR_PLATESCALE
from emirdrp.processing.bars import slits_to_ds9_reg, find_bars
from emirdrp.processing.combine import basic_processing_with_combination
from emirdrp.products import DataFrameType, NominalPositions

from emirdrp.requirements import MasterBadPixelMaskRequirement
from emirdrp.requirements import MasterBiasRequirement
from emirdrp.requirements import MasterDarkRequirement
from emirdrp.requirements import MasterIntensityFlatFieldRequirement
from emirdrp.requirements import MasterSkyRequirement


class MaskCheckRecipe(EmirRecipe):

    """
    Acquire a target.

    Recipe for the processing of multi-slit/long-slit check images.

    **Observing modes:**

        * MSM and LSM check

    """

    # Recipe Requirements
    #
    obresult = ObservationResultRequirement()
    master_bpm = MasterBadPixelMaskRequirement()
    master_bias = MasterBiasRequirement()
    master_dark = MasterDarkRequirement()
    master_flat = MasterIntensityFlatFieldRequirement()
    master_sky = MasterSkyRequirement()

    bars_nominal_positions = Requirement(NominalPositions,
                                         'Nominal positions of the bars'
                                         )
    median_filter_size = Parameter(5, 'Size of the median box')
    average_box_row_size = Parameter(7, 'Number of rows to average for fine centering (odd)')
    average_box_col_size = Parameter(21, 'Number of columns to extract for fine centering (odd)')
    fit_peak_npoints = Parameter(3, 'Number of points to use for fitting the peak (odd)')

    # Recipe Products
    frame = Product(DataFrameType)
    # derivative = Product(DataFrameType)
    slits = Product(ArrayType)
    positions3 = Product(ArrayType)
    positions5 = Product(ArrayType)
    positions7 = Product(ArrayType)
    positions9 = Product(ArrayType)
    DTU = Product(ArrayType)
    ROTANG = Product(float)
    TSUTC1 = Product(float)
    csupos = Product(ArrayType)
    csusens = Product(ArrayType)
    offset = Product(list, description='In arcseconds')
    rotang = Product(float, description='Counterclockwise, in degrees')

    def run(self, rinput):
        self.logger.info('starting processing for bars detection')

        flow = self.init_filters(rinput)

        hdulist = basic_processing_with_combination(rinput, flow=flow)

        hdr = hdulist[0].header
        self.set_base_headers(hdr)

        self.save_intermediate_img(hdulist, 'reduced_image.fits')

        try:
            rotang = hdr['ROTANG']
            tsutc1 = hdr['TSUTC1']
            dtub, dtur = datamodel.get_dtur_from_header(hdr)
            csupos = datamodel.get_csup_from_header(hdr)
            if len(csupos) != 2 * EMIR_NBARS:
                raise RecipeError('Number of CSUPOS != 2 * NBARS')
            csusens = datamodel.get_cs_from_header(hdr)

        except KeyError as error:
            self.logger.error(error)
            raise RecipeError(error)

        self.logger.debug('start finding bars')

        allpos, slits = find_bars(hdulist,
                                  rinput.bars_nominal_positions,
                                  csupos,
                                  dtur,
                                  average_box_row_size=rinput.average_box_row_size,
                                  average_box_col_size=rinput.average_box_col_size,
                                  fit_peak_npoints=rinput.fit_peak_npoints,
                                  median_filter_size=rinput.median_filter_size,
                                  logger=self.logger
                                  )

        self.logger.debug('end finding bars')

        if self.intermediate_results:
            import numpy
            with open('ds9.reg', 'w') as ds9reg:
                slits_to_ds9_reg(ds9reg, slits)
            numpy.savetxt('slits.txt', slits)

        off, angle, qc = compute_off_rotation(hdulist, slits, logger=self.logger)
        off_arc = [o * EMIR_PLATESCALE for o in off]
        self.logger.info('offset %s (pix), rotang %6.3f (deg)', off, angle)
        self.logger.info('offset %s (arcsec), rotang %6.3f (deg)', off_arc, angle)

        result = self.create_result(frame=hdulist,
                                    slits=slits,
                                    positions9=allpos[9],
                                    positions7=allpos[7],
                                    positions5=allpos[5],
                                    positions3=allpos[3],
                                    DTU=dtub,
                                    ROTANG=rotang,
                                    TSUTC1=tsutc1,
                                    csupos=csupos,
                                    csusens=csusens,
                                    qc=qc,
                                    offset=off_arc,
                                    rotang=angle
                                    )
        return result


class CSUConf(object):
    """Information about the configuration of slits in the CSU"""
    def __init__(self):
        self.name = 'CSU'
        self.conf_id = 'v1'
        self.slits = {}
        self.lbars = {}
        self.pbars = {}


def read_csu(hdr, slit_table):
    """Read CSU information and slits from header"""
    conf = CSUConf()
    conf.name = hdr.get('INSMODE', 'NULL')
    conf.conf_id = hdr.get('CONFID', 'v1')

    slits_lg = conf.slits
    lbars = conf.lbars
    pbars = conf.pbars

    for idx, slit in enumerate(slit_table, 1):
        xpos1, y2, xpos2, y2, xpos2, y1, xpos1, y1 = slit
        lbars[idx] = PhysicalBarL(idx, xpos1, y1, y2)
        pbars[idx + EMIR_NBARS] = PhysicalBarR(idx + EMIR_NBARS, xpos2, y1, y2)

    # loop over everything, count SLITFL
    # pattern1 = re.compile(r"SLITFL(\d+)")
    pattern2 = re.compile(r"BARLSL(\d+)")

    # Find bars that belong to a slit
    available_bars = range(1, EMIR_NBARS + 1)
    used_bars = [0] * len(available_bars)

    slit_lbar_ids = {}
    for key in hdr:
        bel_re = pattern2.match(key)
        if bel_re:
            bun_idx = int(bel_re.group(1))
            value = hdr[key]
            slit_lbar_ids.setdefault(value, []).append(bun_idx)
            # print(bun_idx, value)

    if len(slit_lbar_ids) == 0:
        for idx in range(1, EMIR_NBARS + 1):
            slit_lbar_ids[idx] = [idx]

    # Check
    for slit_id, bars_ids in slit_lbar_ids.items():
        for bar_id in bars_ids:
            used = used_bars[bar_id - 1]
            if used == 0:
                used_bars[bar_id - 1] = bar_id
            else:
                raise ValueError('bar used in slits {} and {}'.format(bar_id, used))

    for idx, ls in slit_lbar_ids.items():
        bars_l = {bidx: lbars[bidx] for bidx in ls}
        bars_r = {bidx + EMIR_NBARS: pbars[bidx + EMIR_NBARS] for bidx in ls}

        slit_t = hdr.get("SLITFL%02d" % idx, 0)
        target_coordinates = None
        target_type = TargetType.UNASSIGNED
        if slit_t == 2:
            xref = hdr["XREASL%02d" % idx] - 1
            yref = hdr["YREASL%02d" % idx] - 1
            target_type = TargetType.REFERENCE
            target_coordinates = (xref, yref)
        elif slit_t == 1:
            target_type = TargetType.SOURCE
        elif slit_t == 0:
            target_type = TargetType.UNASSIGNED
        else:
            raise ValueError("invalid value '{}' for SLITFLXX".format(slit_t))

        slits_lg[idx] = LogicalSlit(
            idx, bars_l, bars_r,
            target_type=target_type,
            target_coordinates=target_coordinates
        )

    return conf


class TargetType(enum.Enum):
    """Possible targets in a slit"""
    SOURCE = 1
    UNKNOWN = 2
    UNASSIGNED = 3
    SKY = 4
    REFERENCE = 5
    # aliases for the other fields
    STAR = 5
    BLANK = 4


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
    def __init__(self, idx, lbars, rbars, target_type=TargetType.UNKNOWN, target_coordinates=None):
        self.target_type = target_type
        self.target_coordinates = target_coordinates
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
            if (r != l + EMIR_NBARS):
                raise ValueError('not paired {} and {}'.format(l, r))

        for a, b in zip(lbars_ids, lbars_ids[1:]):
            if (b != a + 1):
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

        self.bbox_int = [coor_to_pix_1d(w) for w in [bbox_x1, bbox_x2, bbox_y1, bbox_y2]]

    def bbox(self):
        return self.bbox_int


def compute_off_rotation(hdulist, slits, logger=None):

    if logger is None:
        logger = logging.getLogger(__name__)

    data = hdulist[0].data
    hdr = hdulist[0].header
    swapped_code = (sys.byteorder == 'little') and '>' or '<'
    if data.dtype.byteorder == swapped_code:
        data = data.byteswap().newbyteorder()

    # add_to_header(hdr)
    conf = read_csu(hdr, slits)

    logger.info('we have %s slits', len(conf.slits))
    ref_slits = [slit for slit in conf.slits.values() if slit.target_type is TargetType.REFERENCE]
    logger.info('we have %s reference slits', len(ref_slits))

    p1 = [] # Reference positions
    q1 = [] # Measured positions
    for slit_lg in ref_slits:
        bbox = slit_lg.bbox()
        region = np.s_[bbox[2]:bbox[3] + 1, bbox[0]:bbox[1] + 1]
        res = comp_centroid(data, region, debug_plot=False, plot_reference=slit_lg.target_coordinates)
        if res is None:
            logger.warning('no object found in slit %s, skipping', slit_lg.idx)
            continue
        m_x = res[0] + region[1].start
        m_y = res[1] + region[0].start
        logger.debug('in slit %s, reference is %s', slit_lg.idx, slit_lg.target_coordinates)
        logger.debug('in slit %s, object is (%s, %s)', slit_lg.idx, m_x, m_y)
        p1.append(slit_lg.target_coordinates)
        q1.append((m_x, m_y))

    logger.info('compute offset and rotation with %d points', len(p1))
    if len(p1) == 0:
        logger.warning('cant compute offset and rotation with 0 points')
        offset = [0.0, 0.0]
        rot = [[1, 0], [0, 1]]
        angle = 0.0
        qc = QC.BAD
    else:
        offset, rot = fit_offset_and_rotation(np.array(p1), np.array(q1))
        angle = math.atan2(rot[1, 0], rot[0, 0])
        angle = np.rad2deg(angle)
        qc = QC.GOOD
    logger.info('offset is %s', offset)
    logger.info('rot matrix is %s', rot)
    logger.info('rot angle %5.2f deg', angle)
    return offset, angle, qc


def comp_centroid(data, region, debug_plot=True, plot_reference=None, logger=None):
    from matplotlib.patches import Ellipse

    if logger is None:
        logger = logging.getLogger(__name__)

    ref_x = region[1].start
    ref_y = region[0].start
    logger.debug('region ofset is %s, %s', ref_x, ref_y)
    subimage = data[region].copy()
    bkg = sep.Background(subimage)
    data_sub = subimage - bkg
    objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
    # Select brightest object
    logger.debug('%d object found', len(objects))

    if len(objects) == 0:
        # print('No objects')
        return None

    iadx = objects['flux'].argmax()
    # plot background-subtracted image

    if debug_plot:
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       vmin=m - s, vmax=m + s, origin='lower')

        print(plot_reference)
        if plot_reference:
            e = Ellipse(xy=(plot_reference[0] - ref_x, plot_reference[1] - ref_y),
                        width=6,
                        height=6,
                        angle=0)
            e.set_facecolor('none')
            e.set_edgecolor('green')
            ax.add_artist(e)

        # plot an ellipse for each object
        for idx, obj in enumerate(objects):
            e = Ellipse(xy=(obj['x'], obj['y']),
                        width=6 * obj['a'],
                        height=6 * obj['b'],
                        angle=obj['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            if idx == iadx:
                e.set_edgecolor('blue')
            else:
                e.set_edgecolor('red')
            ax.add_artist(e)

        plt.show()
    maxflux = objects[iadx]
    return maxflux['x'], maxflux['y']


def add_to_header(hdr):

    logical_slits = [(5, 6, 7, 8, 9), (10, 11, 12, 13), (15, 16, 17),
                     (25, 26, 27), (35, 36, 37, 38), (42, 43, 44, 45),
                     (47, 48, 49, 50)]

    ref = [(679.72959501, 247.518684593),
           (987.645411484, 394.558690093),
           (880.403465621, 572.913094434),
           (1156.4863229, 933.602503354),
           (1399.31597101, 1333.42710666),
           (1412.65781051, 1600.04941192),
           (686.282485999, 1782.41253144)
           ]

    for idx in range(1, 55 + 1):
        hdr["BARSLI%02d" % idx] = 0

    for idx, ls in enumerate(logical_slits, 1):
        hdr["SLITFL%02d" % idx] = 2
        hdr["XVIRSL%02d" % idx] = ref[idx - 1][0]
        hdr["YVIRSL%02d" % idx] = ref[idx - 1][1]
        hdr["XREASL%02d" % idx] = ref[idx - 1][0]
        hdr["YREASL%02d" % idx] = ref[idx - 1][1]

        for bidx in ls:
            hdr["BARLSL%02d" % bidx] = idx


def add_to_header1(hdr):

    for idx in range(1, 55 + 1):
        hdr["SLITFL%02d" % idx] = 1
        hdr["XVIRSL%02d" % idx] = 0
        hdr["YVIRSL%02d" % idx] = 0
        hdr["XREASL%02d" % idx] = 0
        hdr["YREASL%02d" % idx] = 0


def add_to_header2(hdr):

    for idx in range(1, 55 + 1):
        hdr["SLITFL%02d" % idx] = 1
        hdr["XVIRSL%02d" % idx] = 0
        hdr["YVIRSL%02d" % idx] = 0
        hdr["XREASL%02d" % idx] = 0
        hdr["YREASL%02d" % idx] = 0

    hdr["SLITFL33"] = 2