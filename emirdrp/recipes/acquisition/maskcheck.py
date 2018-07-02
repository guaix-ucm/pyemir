#
# Copyright 2008-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

from __future__ import division


import sys
import enum
import math
import logging
import itertools

import numpy as np
import matplotlib.pyplot as plt
import sep
import skimage.filters as filt

import skimage.feature as feat
from scipy import ndimage as ndi
from numina.array.utils import coor_to_pix_1d
from numina.array.bbox import BoundingBox
from numina.array.offrot import fit_offset_and_rotation
from numina.array.peaks.peakdet import refine_peaks
from numina.core import Requirement, Result, Parameter
from numina.types.qc import QC
from numina.core.query import ResultOf
import numina.types.array as tarray

from emirdrp.core import EMIR_NBARS, EMIR_PLATESCALE_PIX, EMIR_PLATESCALE
from emirdrp.core import EMIR_PIXSCALE
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import basic_processing_, combine_images
from emirdrp.processing.combine import process_ab, process_abba
import emirdrp.requirements as reqs
import emirdrp.products as prods
import emirdrp.instrument.distortions as dist


_NOMINALPOS = \
"""
   1  11.21484759	6.050300152	1.006931003
   2  49.30894249	6.048943068	1.305389218
   3  86.41217695	6.049308669	1.158653968
   4  123.7548247	6.049885147	1.111796273
   5  161.6850211	6.047962374	2.132388709
   6  199.1181558	6.048789102	2.361405842
   7  236.7785942	6.050007163	1.834087278
   8  274.6029416	6.047577921	1.861126798
   9  311.9311216	6.049781377	1.263138212
  10  349.2233764	6.049463532	0.599254159
  11  386.7347265	6.050508173	0.982413022
  12  424.365245	6.050510187	0.918456545
  13  461.9852561	6.049841237	1.514800167
  14  499.4136008	6.050851183	1.519270219
  15  536.5579487	6.051029481	1.034473702
  16  573.4715124	6.050678957	1.275667546
  17  610.9907427	6.052483829	0.248068988
  18  648.5853825	6.050819015	0.784814234
  19  685.9746342	6.050694356	1.20143528
  20  723.7704432	6.052307258	0.990134405
  21  761.1686225	6.052449288	0.629286656
  22  798.6618566	6.051362842	0.637959623
  23  836.2400178	6.051799848	1.211803226
  24  873.8721331	6.051498233	0.644248767
  25  911.3460353	6.054077871	0.011146096
  26  948.5305385	6.052221584	0.963347496
  27  986.0982998	6.053368433	0.87872476
  28  1023.731596	6.052427027	0.801232259
  29  1061.196977	6.052319943	0.055943108
  30  1098.571986	6.052740583	0.731693911
  31  1136.074888	6.05409902	1.276135495
  32  1173.758293	6.05130894	0.506021149
  33  1211.303227	6.052692448	1.072851013
  34  1248.523557	6.051314344	0.472204979
  35  1285.859445	6.051979653	0.050803575
  36  1323.604486	6.051208944	0.806198309
  37  1361.754394	6.051232903	-0.080860162
  38  1399.180986	6.050721818	0.656130308
  39  1436.370347	6.051125813	0.990730245
  40  1474.099407	6.04943513	1.325242595
  41  1511.057571	6.051369181	0.401007343
  42  1548.268634	6.04911982	0.204480737
  43  1585.711148	6.049386564	0.475008698
  44  1623.06518	6.048524412	0.980650651
  45  1660.729215	6.050010384	0.867948914
  46  1698.173454	6.049537497	0.974817179
  47  1735.751871	6.05016223	0.662702411
  48  1773.172551	6.04877701	0.872045426
  49  1810.665282	6.049070212	1.146726653
  50  1848.385975	6.048937228	0.69930393
  51  1886.140714	6.047995261	0.964154069
  52  1923.900424	6.048294844	1.02590626
  53  1961.181748	6.052027118	1.011679953
  54  1998.66557	6.049786384	0.997423159
  55  2037.108886	6.051019236	-0.018442041
  56  11.21484759	-6.050949702	2062.952966
  57  49.30894249	-6.049057793	2061.96369
  58  86.41217695	-6.049509347	2062.632646
  59  123.7548247	-6.050344768	2062.755618
  60  161.6850211	-6.050647553	2063.635573
  61  199.1181558	-6.05032146	2063.088107
  62  236.7785942	-6.050630119	2063.617342
  63  274.6029416	-6.049497138	2062.404227
  64  311.9311216	-6.050546071	2063.159375
  65  349.2233764	-6.050469039	2063.200361
  66  386.7347265	-6.052716316	2063.336789
  67  424.365245	-6.049246635	2062.7177
  68  461.9852561	-6.050867087	2063.514774
  69  499.4136008	-6.0503855	2062.576461
  70  536.5579487	-6.050332102	2063.164178
  71  573.4715124	-6.050531515	2062.804584
  72  610.9907427	-6.052611793	2064.386082
  73  648.5853825	-6.052401462	2063.360315
  74  685.9746342	-6.051752845	2061.927832
  75  723.7704432	-6.051841443	2062.896491
  76  761.1686225	-6.052234416	2063.100489
  77  798.6618566	-6.054154105	2063.507102
  78  836.2400178	-6.053753227	2063.34106
  79  873.8721331	-6.05139629	2062.943969
  80  911.3460353	-6.053373651	2063.323767
  81  948.5305385	-6.052794361	2063.214531
  82  986.0982998	-6.053463769	2062.839414
  83  1023.731596	-6.05172025	2062.749323
  84  1061.196977	-6.053160768	2062.943505
  85  1098.571986	-6.052837995	2063.288001
  86  1136.074888	-6.051994274	2063.316395
  87  1173.758293	-6.050987612	2062.682207
  88  1211.303227	-6.051955824	2063.046156
  89  1248.523557	-6.049989693	2063.177237
  90  1285.859445	-6.052119441	2062.920172
  91  1323.604486	-6.05073856	2062.956267
  92  1361.754394	-6.051192485	2062.853586
  93  1399.180986	-6.049418457	2063.128934
  94  1436.370347	-6.050694762	2062.976157
  95  1474.099407	-6.049564403	2062.128429
  96  1511.057571	-6.050042144	2062.17244
  97  1548.268634	-6.0511598	2062.868856
  98  1585.711148	-6.048381044	2062.005524
  99  1623.06518	-6.04875246	2062.378611
 100  1660.729215	-6.05004154	2062.383435
 101  1698.173454	-6.047694054	2062.168699
 102  1735.751871	-6.048362116	2062.125675
 103  1773.172551	-6.046601223	2061.752033
 104  1810.665282	-6.048933664	2062.238642
 105  1848.385975	-6.047375916	2062.206302
 106  1886.140714	-6.048968933	2062.148436
 107  1923.900424	-6.048072168	2062.682802
 108  1961.181748	-6.049388185	2062.875553
 109  1998.66557	-6.051123287	2061.869519
 110  2037.108886	-6.050515352	2062.176973
"""


class CSUConf(object):
    """Information about the configuration of slits in the CSU"""
    def __init__(self):
        self.name = 'CSU'
        self.conf_id = 'v1'
        self.slits = {}
        # Indices
        self.lbars = list(range(1, 55 + 1))
        self.rbars = [(i + 55) for i in self.lbars]
        self.bars = {}


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
        return BoundingBox.from_coordinates(
            x1=max(self.x1-1, 0), x2=max(self.x2-1, 0),
            y1=max(self.y1-1, 0), y2=max(self.y2-1, 0)
        )


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


def comp_centroid(data, bounding_box, debug_plot=False, plot_reference=None, logger=None):
    from matplotlib.patches import Ellipse

    if logger is None:
        logger = logging.getLogger(__name__)

    region = bounding_box.slice
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
    maxflux = objects[iadx]

    if debug_plot:
        fig, ax = plt.subplots()
        m, s = np.mean(data_sub), np.std(data_sub)
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       vmin=m - s, vmax=m + s, origin='lower',
                       extent=bounding_box.extent)
        if plot_reference:
            e = Ellipse(xy=(plot_reference[0], plot_reference[1]),
                        width=6,
                        height=6,
                        angle=0)
            e.set_facecolor('none')
            e.set_edgecolor('green')
            ax.add_artist(e)

        # plot an ellipse for each object
        for idx, obj in enumerate(objects):
            e = Ellipse(xy=(obj['x'] + ref_x, obj['y'] + ref_y),
                        width=6 * obj['a'],
                        height=6 * obj['b'],
                        angle=obj['theta'] * 180. / np.pi)
            e.set_facecolor('none')
            if idx == iadx:
                e.set_edgecolor('blue')
            else:
                e.set_edgecolor('red')
            ax.add_artist(e)
        return maxflux['x'], maxflux['y'], ax
    else:
        return maxflux['x'], maxflux['y']


class MaskCheckRecipe(EmirRecipe):

    """
    Acquire a target.

    Recipe for the processing of multi-slit/long-slit check images.

    **Observing modes:**

        * MSM and LSM check

    """

    # Recipe Requirements
    #
    obresult = reqs.ObservationResultRequirement()
    master_bpm = reqs.MasterBadPixelMaskRequirement()

    #bars_nominal_positions = Requirement(
    #    prods.NominalPositions,
    #    'Nominal positions of the bars'
    #)

    # Recipe Products
    slit_image = Result(prods.ProcessedImage)
    object_image = Result(prods.ProcessedImage)
    offset = Result(tarray.ArrayType)
    angle = Result(float)

    def run(self, rinput):
        self.logger.debug('workaround for bars_nominal_positions')
        try:
            import StringIO as S
        except ImportError:
            import io as S
        import numpy
        ss = S.StringIO(_NOMINALPOS)
        rinput.bars_nominal_positions = numpy.loadtxt(ss)

        self.logger.info('starting processing for bars detection')
        # Combine and masking
        flow = self.init_filters(rinput)

        # count frames
        frames = rinput.obresult.frames

        nframes = len(frames)
        if nframes == 1:
            slit_frames = frames
            object_frames = frames
            background_subs = False
        elif nframes == 2:
            slit_frames = frames[0:]
            object_frames = frames
            background_subs = True
        elif nframes == 4:
            slit_frames = frames[0::3]
            object_frames = frames
            background_subs = True
        else:
            raise ValueError("expected 1, 2 or 4 frames, got {}".format(nframes))

        interm = basic_processing_(frames, flow, self.datamodel)

        if nframes == 1:
            hdulist_slit = combine_images(interm[:], self.datamodel)
            hdulist_object = hdulist_slit
            background_subs = False
        elif nframes == 2:
            hdulist_slit = combine_images(interm[0:], self.datamodel)
            hdulist_object = process_ab(interm, self.datamodel)
            background_subs = True
        elif nframes == 4:
            hdulist_slit = combine_images(interm[0::3], self.datamodel)
            hdulist_object = process_abba(interm, self.datamodel)
            background_subs = True
        else:
            raise ValueError("expected 1, 2 or 4 frames, got {}".format(nframes))

        self.set_base_headers(hdulist_slit[0].header)
        self.set_base_headers(hdulist_object[0].header)

        self.save_intermediate_img(hdulist_slit, 'slit_image.fits')

        self.save_intermediate_img(hdulist_object, 'object_image.fits')

        # Get slits
        # Rotation around (0,0)
        # For other axis, offset is changed
        # (Off - raxis) = Rot * (Offnew - raxis)
        crpix1 = hdulist_slit[0].header['CRPIX1']
        crpix2 = hdulist_slit[0].header['CRPIX2']

        rotaxis = np.array((crpix1 - 1, crpix2 - 1))

        self.logger.debug('center of rotation (from CRPIX) is %s', rotaxis)

        csu_conf, slits_bb = self.compute_slits(hdulist_slit, rinput.bars_nominal_positions)

        image_sep = hdulist_object[0].data.astype('float32')

        self.logger.debug('center of rotation (from CRPIX) is %s', rotaxis)

        offset, angle, qc = compute_off_rotation(image_sep, csu_conf, slits_bb,
                                                 rotaxis=rotaxis, logger=self.logger,
                                                 debug_plot=True, intermediate_results=True
                                                 )

        # Convert mm to m
        offset_out = np.array(offset) / 1000.0
        # Convert DEG to RAD
        angle_out = np.deg2rad(angle)
        result = self.create_result(
            slit_image=hdulist_slit,
            object_image=hdulist_object,
            offset=offset_out,
            angle=angle_out,
            qc=qc
        )

        return result

    def compute_slits(self, hdulist, bars_nominal_positions):
        # Get slits
        hdr = hdulist[0].header
        # Extract DTU and CSU information from headers

        dtuconf = self.datamodel.get_dtur_from_header(hdr)

        # coordinates transformation from DTU coordinates
        # to image coordinates
        # Y inverted
        # XY switched
        # trans1 = [[1, 0, 0], [0,-1, 0], [0,0,1]]
        # trans2 = [[0,1,0], [1,0,0], [0,0,1]]
        trans3 = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]  # T3 = T2 * T1

        vec = np.dot(trans3, dtuconf.coor_r) / EMIR_PIXSCALE
        self.logger.debug('DTU shift is %s', vec)

        self.logger.debug('create bar model')
        barmodel = create_bar_models(bars_nominal_positions)
        csu_conf = read_csu_2(hdr, barmodel)

        if self.intermediate_results:
            self.logger.debug('create bar mask from predictions')
            mask1 = np.ones_like(hdulist[0].data)
            for i in itertools.chain(csu_conf.lbars, csu_conf.rbars):
                bar = csu_conf.bars[i]
                mask1[bar.bbox().slice] = 0

            self.save_intermediate_array(mask1, 'mask1.fits')

            self.logger.debug('create slit mask from predictions')
            mask1 = np.zeros_like(hdulist[0].data)
            for slit in csu_conf.slits.values():
                mask1[slit.bbox().slice] = slit.idx
            self.save_intermediate_array(mask1, 'mask2.fits')

            self.logger.debug('create slit reference mask from predictions')
            mask1 = np.zeros_like(hdulist[0].data)
            for slit in csu_conf.slits.values():
                if slit.target_type == TargetType.REFERENCE:
                    mask1[slit.bbox().slice] = slit.idx
            self.save_intermediate_array(mask1, 'mask3.fits')

        self.logger.debug('finding borders of slits')
        self.logger.debug('not strictly necessary...')
        data = hdulist[0].data
        self.logger.debug('dtype of data %s', data.dtype)

        self.logger.debug('median filter (3x3)')
        image_base = ndi.filters.median_filter(data, size=3)

        # Cast as original type for skimage
        self.logger.debug('casting image to unit16 (for skimage)')
        iuint16 = np.iinfo(np.uint16)
        image = np.clip(image_base, iuint16.min, iuint16.max).astype(np.uint16)

        self.logger.debug('compute Sobel filter')
        # FIXME: compute sob and sob_v is redundant
        sob = filt.sobel(image)
        self.save_intermediate_array(sob, 'sobel_image.fits')
        sob_v = filt.sobel_v(image)
        self.save_intermediate_array(sob_v, 'sobel_v_image.fits')

        # Compute detector coordinates of bars
        all_coords_virt = np.empty((110, 2))
        all_coords_real = np.empty((110, 2))

        # Origin of coordinates is 1
        for bar in csu_conf.bars.values():
            all_coords_virt[bar.idx - 1] = bar.xpos, bar.y0

        # Origin of coordinates is 1 for this function
        _x, _y = dist.exvp(all_coords_virt[:, 0], all_coords_virt[:, 1])
        all_coords_real[:, 0] = _x
        all_coords_real[:, 1] = _y

        # FIXME: hardcoded value
        h = 16
        slit_h_virt = 16.242
        slit_h_tol = 3
        slits_bb = {}

        mask1 = np.zeros_like(hdulist[0].data)

        for idx in range(EMIR_NBARS):
            lbarid = idx + 1
            rbarid = lbarid + EMIR_NBARS
            ref_x_l_v, ref_y_l_v = all_coords_virt[lbarid - 1]
            ref_x_r_v, ref_y_r_v = all_coords_virt[rbarid - 1]

            ref_x_l_d, ref_y_l_d = all_coords_real[lbarid - 1]
            ref_x_r_d, ref_y_r_d = all_coords_real[rbarid - 1]

            width_v = ref_x_r_v - ref_x_l_v
            # width_d = ref_x_r_d - ref_x_l_d

            if (ref_y_l_d >= 2047 + h) or (ref_y_l_d <= 1 - h):
                # print('reference y position is outlimits, skipping')
                continue

            if width_v < 5:
                # print('width is less than 5 pixels, skipping')
                continue

            plot = False
            regionw = 12
            px1 = coor_to_pix_1d(ref_x_l_d) - 1
            px2 = coor_to_pix_1d(ref_x_r_d) - 1
            prow = coor_to_pix_1d(ref_y_l_d) - 1

            comp_l, comp_r = calc0(image, sob_v, prow, px1, px2, regionw, h=h,
                                   plot=plot, lbarid=lbarid, rbarid=rbarid,
                                   plot2=False)

            region2 = 5
            px21 = coor_to_pix_1d(comp_l)
            px22 = coor_to_pix_1d(comp_r)

            comp2_l, comp2_r = calc0(image, sob_v, prow, px21, px22, region2,
                                     refine=True,
                                     plot=plot, lbarid=lbarid, rbarid=rbarid,
                                     plot2=False)

            if np.any(np.isnan([comp2_l, comp2_r])):
                self.logger.warning("converting NaN value, border of=%d", idx + 1)
                comp2_l, comp2_r = comp_l, comp_r
            # print('slit', lbarid, '-', rbarid, comp_l, comp_r)
            # print('pos1', comp_l, comp_r)
            # print('pos2', comp2_l, comp2_r)

            xpos1_virt, _ = dist.pvex(comp2_l + 1, ref_y_l_d)
            xpos2_virt, _ = dist.pvex(comp2_r + 1, ref_y_r_d)

            y1_virt = ref_y_l_v - slit_h_virt - slit_h_tol
            y2_virt = ref_y_r_v + slit_h_virt + slit_h_tol
            _, y1 = dist.exvp(xpos1_virt + 1, y1_virt)
            _, y2 = dist.exvp(xpos2_virt + 1, y2_virt)
            # print(comp2_l, comp2_r, y1 - 1, y2 - 1)
            cbb = BoundingBox.from_coordinates(comp2_l, comp2_r, y1 - 1, y2 - 1)
            slits_bb[lbarid] = cbb
            mask1[cbb.slice] = lbarid

        self.save_intermediate_array(mask1, 'mask3.fits')
        return csu_conf, slits_bb


def pix2virt(pos, origin=1):
    ddef_o = 1
    off = ddef_o - origin
    pos = np.atleast_2d(pos) + off
    nx, ny = dist.pvex(pos[:, 0], pos[:, 1])
    res = np.stack((nx, ny), axis=1)
    return res - off


def compute_off_rotation(data, csu_conf, slits_bb, rotaxis=(0, 0),
                         logger=None, debug_plot=False,
                         intermediate_results=True
                         ):

    if logger is None:
        logger = logging.getLogger(__name__)

    swapped_code = (sys.byteorder == 'little') and '>' or '<'
    if data.dtype.byteorder == swapped_code:
        data = data.byteswap().newbyteorder()

    logger.info('we have %s slits', len(csu_conf.slits))
    refslits = [slit for slit in csu_conf.slits.values() if slit.target_type is TargetType.REFERENCE]
    logger.info('we have %s reference slits', len(refslits))

    p2 = []  # Pos REF VIRT
    p1 = []  # Pos REF DET
    q1 = []  # Pos Measured DET
    q2 = []  # Pos Measured VIRT
    EMIR_REF_IPA = 90.0552

    for this in refslits:
        bb = slits_bb[this.idx]
        region = bb.slice
        target_coordinates = this.target_coordinates
        res = comp_centroid(data, bb, debug_plot=debug_plot, plot_reference=target_coordinates)

        logger.debug('in slit %s, reference is %s', this.idx, target_coordinates)

        if res is None:
            logger.warning('no object found in slit %s, skipping', this.idx)
            continue

        if debug_plot and intermediate_results:
            ax = res[2]
            ax.set_title('slit %s' % this.idx)
            plt.savefig('centroid_slit_%s.png' % this.idx)
            # plt.show()

        m_x = res[0] + region[1].start
        m_y = res[1] + region[0].start

        logger.debug('in slit %s, object is (%s, %s)', this.idx, m_x, m_y)
        p1.append(this.target_coordinates)
        p2.append(this.target_coordinates_v)
        q1.append((m_x, m_y))

    logger.info('compute offset and rotation with %d points', len(p1))
    qc = QC.BAD

    if len(p2) == 0:
        logger.warning("can't compute offset and rotation with 0 points")
        offset = np.array([0.0, 0.0])
        rot = np.array([[1, 0], [0, 1]])
        angle = 0.0
    else:
        logger.debug('convert coordinates to virtual, ie, focal plane')
        q2 = pix2virt(q1, origin=0)
        # Move from objects to reference
        logger.debug('compute transform from measured objects to reference coordinates')
        offset, rot = fit_offset_and_rotation(q2, p2)
        logger.debug('rotation matrix')
        logger.debug('%s', rot)
        logger.debug('translation (with rotation around 0)')
        logger.debug('%s', offset)

        logger.debug('center of rotation (from CRPIX) is %s', rotaxis)
        logger.debug('translation (with rotation around rotaxis)')
        newoff = np.dot(rot, offset - rotaxis) + rotaxis
        logger.debug('%s', newoff)

        offset = newoff
        angle = math.atan2(rot[1, 0], rot[0, 0])
        angle = np.rad2deg(angle)
        qc = QC.GOOD
    logger.info('offset is %s', offset)
    logger.info('rot matrix is %s', rot)
    logger.info('rot angle %5.2f deg', angle)

    o_mm = offset * EMIR_PLATESCALE_PIX / EMIR_PLATESCALE
    angle = np.deg2rad(EMIR_REF_IPA)
    ipa_rot = create_rot2d(angle)
    logger.info('OFF (mm) %s', o_mm)
    logger.info('Default IPA is %s', EMIR_REF_IPA)
    o_mm_ipa = np.dot(ipa_rot, o_mm)

    logger.info('=========================================')
    logger.info('Offset Target in Focal Plane Frame %s mm', o_mm_ipa)
    logger.info('=========================================')

    logger.debug('MEAN of REF-MEASURED (ON DETECTOR) %s', np.subtract(p1, q1).mean(axis=0))
    logger.debug('MEAN pf REF-MEASURED (VIRT) %s', np.subtract(p2, q2).mean(axis=0))

    return offset, angle, qc


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


def calc0(image, sob, prow, px1, px2, regionw, h=16, refine=False,
          plot=True, lbarid=0, rbarid=0, plot2=False):
    borders = [px1 - regionw, px1 + regionw, px2 - regionw, px2 + regionw,
               prow - h, prow + h]
    borders = np.clip(borders, 0, 2047)
    cc1l = borders[0]
    cc1r = borders[1]
    cc2l = borders[2]
    cc2r = borders[3]
    py1 = borders[4]
    py2 = borders[5]

    bb = BoundingBox.from_coordinates(cc1l, cc2r, py1, py2)

    lres = []

    # print(cc1l, cc2r, py1, py2)

    steps = range(-3, 3 + 1)
    scale = 3
    for s in steps:
        pr1 = prow + scale * s
        if pr1 < 0 or pr1 >= 2048:
            continue
        val = calc1b(image, sob, pr1, cc1l, cc1r, refine=refine,
                     plot=plot, barid=lbarid, px=px1)
        if val is not None:
            lres.append(val)

    # Average left bar position
    if lres:
        row_l, col_l, colr_l = zip(*lres)
    else:
        row_l, col_l, colr_l = [], [], []
    if refine:
        comp_l = np.mean(colr_l)
    else:
        comp_l = np.median(col_l)

    rres = []
    scale = 3
    for s in steps:
        pr1 = prow + scale * s
        if pr1 < 0 or pr1 >= 2048:
            continue

        val = calc1b(image, sob, pr1, cc2l, cc2r, sign=-1, refine=refine,
                     plot=plot, barid=rbarid, px=px2)
        if val is not None:
            rres.append(val)
    # Average left bar position
    if rres:
        row_r, col_r, colr_r = zip(*rres)
    else:
        row_r, col_r, colr_r = [], [], []

    if refine:
        comp_r = np.mean(colr_r)
    else:
        comp_r = np.median(col_r)

    if plot2:
        plt.imshow(image[bb.slice], extent=bb.extent)
        plt.scatter(colr_l, row_l, marker='+', color='black')
        plt.scatter(colr_r, row_r, marker='+', color='black')
        plt.axvline(px1, color='blue')
        plt.axvline(px2, color='blue')
        plt.axvline(comp_l, color='red')
        plt.axvline(comp_r, color='red')
        plt.show()

    return comp_l, comp_r


def calc1b(image, sob, pr, ccl, ccr, sign=1, refine=False, plot=True, barid=0, px=0):
    cut1 = sob[pr, ccl:ccr + 1]
    mcut1 = image[pr, ccl:ccr + 1]

    res = feat.peak_local_max(sign * cut1, exclude_border=4, num_peaks=2)

    if len(res) != 0:
        res1 = res[:, 0]
    else:
        res1 = []

    if sign > 0:
        bar_s = 'l'
    else:
        bar_s = 'r'

    if plot:
        xdummy = np.arange(ccl, ccr + 1)
        plt.title("{}bar {}, prow={} px={}".format(bar_s, barid, pr, px))
        plt.plot(xdummy, sign * cut1)
        plt.axvline(px, color='g')
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        ax2.plot(xdummy, mcut1, '--')
        for r in res1:
            plt.axvline(ccl + r, color='red')
        for r in res[:, 0]:
            plt.axvline(ccl + r, color='red', linestyle='--')
    xpeak = None
    if len(res1) != 0:
        idx = (sign * cut1)[res1].argmax()
        xpeak = ccl + res1[idx]

    if plot:
        if xpeak:
            plt.axvline(xpeak, color='black')
        plt.show()

    if xpeak:
        if refine:
            x_t, y_t = refine_peaks(sign * sob[pr], np.array([xpeak]), window_width=3)
            xpeak_ref = x_t
        else:
            xpeak_ref = xpeak
        return pr, xpeak, xpeak_ref
    else:
        return None


def create_rot2d(angle):
    ca = math.cos(angle)
    sa = math.sin(angle)
    return np.array([[ca, -sa], [sa, ca]])
