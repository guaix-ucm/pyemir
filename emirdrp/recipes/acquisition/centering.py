#
# Copyright 2020 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""Fine centering of the MOS"""

from collections import defaultdict

from numina.core import Parameter, Requirement, Result
from numina.core.requirements import ObservationResultRequirement

from emirdrp.core.recipe import EmirRecipe
from emirdrp.instrument.dtuconf import DtuConf
import emirdrp.products as prods


class FineCenteringRecipe(EmirRecipe):
    obresult = ObservationResultRequirement()
    bars_nominal_positions = Requirement(
        prods.NominalPositions, 'Nominal positions of the bars',
        optional=True, default=None
    )
    box_row_size = Parameter(7, 'Number of rows to sum for fine centering (odd)')
    box_col_size = Parameter(7, 'Number of columns to sum for fine centering (odd)')

    img_max = Result(int)
    img_max_uuid = Result(str)

    def run(self, rinput):
        self.logger.info('starting processing for fine centering')
        obresult = rinput.obresult
        if rinput.bars_nominal_positions is None:
            bars_nominal_positions = prods.default_nominal_positions()
        else:
            bars_nominal_positions = rinput.bars_nominal_positions

        uuids = dict()
        for idx, frame in enumerate(obresult.frames, 1):
            with frame.open() as hdulist:
                img_uuid = hdulist[0].header['UUID']
                uuids[idx] = img_uuid
            self.logger.info('image {} has UUID {}'.format(idx, img_uuid))

        images = dict()
        slitids_col = set()

        hrows = rinput.box_row_size // 2
        hcols = rinput.box_col_size // 2
        for idx, frame in enumerate(obresult.frames, 1):
            with frame.open() as hdulist:
                dtuconf = DtuConf.from_img(hdulist)
                vec = dtuconf.vector_shift()
                self.logger.debug('DTU shift is %s', vec)
                csu_conf = self.load_csu_conf(bars_nominal_positions, hdulist)
                self.logger.info('image {} has CSU conf from file: {}'.format(idx, csu_conf.conf_f))
                if not csu_conf.is_open():
                    slitdic = compute_flux(hdulist[0].data, csu_conf, hcols=hcols, hrows=hrows)
                    slitids_col.update(slitdic.keys())
                    images[idx] = slitdic
                else:
                    self.logger.info('CSU is open, not detecting slits')

        # Counting maxima
        cd = defaultdict(int)
        for key in sorted(slitids_col):
            fluxes = []
            for img, slits in images.items():
                flux = slits[key]
                fluxes.append((flux, img))
            res = max(fluxes)
            msg = "For slit {0}, max is in image {1[1]} with value {1[0]}".format(key, res)
            cd[res[1]] += 1
            self.logger.info(msg)
        # Find image with more maxima
        img_max = 0
        for k, v in cd.items():
            if v > img_max:
                img_max = k
            msg = "**Image {} has {} maximum values**".format(k, v)
            self.logger.info(msg)

        result = self.create_result(
            img_max=img_max,
            img_max_uuid=uuids[img_max]
        )
        self.logger.info('end processing for fine centering')
        return result

    def run_qc(self, recipe_input, recipe_result):
        from numina.types.qc import QC
        recipe_result.qc = QC.GOOD
        return recipe_result

    def load_csu_conf(self, bars_nominal_positions, hdulist):
        import emirdrp.instrument.csuconf as csuconf

        self.logger.debug('create bar model')
        barmodel = csuconf.create_bar_models(bars_nominal_positions)
        return csuconf.read_csu_from_image(barmodel, hdulist)


def compute_flux(data, csu_conf, hcols=2, hrows=4, logger=None):
    # We can use a dictionary of slits_bb created by measuring the slits
    # in the image or we can use the slits given by the model of the CSU
    # which is correct around ~ 1.5 pix
    import logging

    from numina.array.utils import image_box2d

    from emirdrp.instrument.csuconf import TargetType

    result = dict()
    if logger is None:
        logger = logging.getLogger(__name__)

    logger.debug('we have %s slits', len(csu_conf.slits))
    accept_type = [TargetType.REFERENCE, TargetType.SOURCE]
    refslits = [slit for slit in csu_conf.slits.values() if slit.target_type in accept_type]
    logger.debug('we have %s reference and source slits', len(refslits))

    for slit in refslits:
        target_coordinates = slit.target_coordinates
        logger.debug('slit %s is formed by bars %s %s', slit.idx, slit.lbars_ids, slit.rbars_ids)
        logger.debug('slit %s has reference %s', slit.idx, target_coordinates)
        # There is a function for this:

        x0 = target_coordinates[0]
        y0 = target_coordinates[1]
        # box is inverted!!!!
        box = (hrows, hcols)
        sf = image_box2d(x0, y0, data.shape, box)
        result[slit.idx] = data[sf].sum()

    return result
