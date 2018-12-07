#
# Copyright 2016-2018 Universidad Complutense de Madrid
#
# This file is part of PyEmir
#
# SPDX-License-Identifier: GPL-3.0+
# License-Filename: LICENSE.txt
#

"""
Spectroscopy mode, ABBA
"""

import astropy.io.fits as fits
import numina.core
import numina.exceptions
import numpy
import numina.core.query as qmod
import numina.ext.gtc

from numina.array import combine
from numina.core import Result, Requirement
from numina.exceptions import RecipeError
from numina.core.requirements import ObservationResultRequirement

import emirdrp.datamodel
import emirdrp.decorators
import emirdrp.products as prods
from emirdrp.core.recipe import EmirRecipe
from emirdrp.processing.combine import basic_processing


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


class BaseABBARecipe(EmirRecipe):
    """Process images in ABBA mode"""

    obresult = ObservationResultRequirement(
        query_opts=qmod.ResultOf(
            'STARE_SPECTRA.stare',
            node='children',
            id_field="stareSpectraIds"
        )
    )
    accum_in = Requirement(
        prods.DataFrameType,
        description='Accumulated result',
        optional=True,
        destination='accum',
        query_opts=qmod.ResultOf(
            'LS_ABBA.accum',
            node='prev'
        )
    )

    spec_abba = Result(prods.ProcessedMOS)
    # Accumulate 'spec_abba' results
    accum = Result(prods.ProcessedMOS, optional=True)

    def build_recipe_input(self, obsres, dal, pipeline='default'):
        if numina.ext.gtc.check_gtc():
            self.logger.debug('running in GTC environment')
            return self.build_recipe_input_gtc(obsres, dal)
        else:
            self.logger.debug('running outside of GTC environment')
            return super(BaseABBARecipe, self).build_recipe_input(
                obsres, dal
            )

    def build_recipe_input_gtc(self, obsres, dal):
        self.logger.debug('start recipe input builder')
        stareImagesIds = obsres.stareSpectraIds
        self.logger.debug('Stare Spectra images IDS: %s', stareImagesIds)
        stareImages = []
        for subresId in stareImagesIds:
            subres = dal.getRecipeResult(subresId)
            stareImages.append(subres['elements']['stare'])

        naccum = obsres.naccum
        self.logger.info('naccum: %d', naccum)
        if naccum != 1:  # if it is not the first dithering loop
            self.logger.info("SEARCHING LATEST RESULT LS_ABBA TO ACCUMULATE")
            latest_result = dal.getLastRecipeResult("EMIR", "EMIR", "LS_ABBA")
            accum_dither = latest_result['elements']['accum']
            self.logger.info("FOUND")
        else:
            self.logger.info("NO ACCUMULATION LS_ABBA")
            accum_dither = stareImages[0]

        newOR = numina.core.ObservationResult()
        newOR.frames = stareImages
        newOR.naccum = naccum
        newOR.accum = accum_dither
        newRI = self.create_input(obresult=newOR)
        self.logger.debug('end recipe input builder')
        return newRI

    def run(self, rinput):
        partial_result = self.run_single(rinput)
        new_result = self.aggregate_result(partial_result, rinput)
        return new_result

    def run_single(self, rinput):
        self.logger.info('starting spectroscopy ABBA reduction')

        flow = self.init_filters(rinput)
        nimages = len(rinput.obresult.frames)
        self.logger.info('we receive %d images', nimages)
        if nimages != 4:
            msg = 'Recipe expects 4 images, received %d' % nimages
            raise numina.exceptions.RecipeError(msg)

        procesed_hdulists = basic_processing(rinput, flow)

        # INPUTS are ABBA, so
        #
        hdulist = self.process_abba(procesed_hdulists)
        grism = hdulist[0].header.get('GRISM', 'unknown')
        if grism.lower() == 'open':
            # perform TEST10 in addition
            import emirdrp.recipes.acquisition.maskcheck as mk
            from numina.core import ObservationResult, DataFrame
            import numpy
            try:
                import StringIO as S
            except ImportError:
                import io as S

            self.logger.info('GRISM is OPEN, doing a RECIPE10')
            sub = mk.MaskCheckRecipe()
            sub.configure(instrument='EMIR', mode='TEST10')
            o = ObservationResult()
            o.__dict__ = rinput.obresult.__dict__
            o.frames = [DataFrame(frame=hdulist)]
            subd = {}
            subd['obresult'] = o
            ss = S.StringIO(_NOMINALPOS)
            subd['bars_nominal_positions'] = numpy.loadtxt(ss)

            subinput = mk.MaskCheckRecipe.RecipeInput(**subd)

            sub.run(subinput)

        result = self.create_result(spec_abba=hdulist)
        self.logger.info('end spectroscopy ABBA reduction')
        return result

    def process_abba(self, images):
        # Process four images in ABBA mode
        dataA0 = images[0][0].data.astype('float32')
        dataB0 = images[1][0].data.astype('float32')

        dataB1 = images[2][0].data.astype('float32')
        dataA1 = images[3][0].data.astype('float32')

        dataAB0 = dataA0 - dataB0
        dataAB1 = dataA1 - dataB1

        dataABBA = dataAB0 + dataAB1

        hdulist = self.create_proc_hdulist(images, dataABBA)
        self.logger.debug('update result header')
        hdu = hdulist[0]
        hdu.header['history'] = "Processed ABBA"
        hdu.header['NUM-NCOM'] = (2, 'Number of combined frames')
        dm = emirdrp.datamodel.EmirDataModel()
        for img, key in zip(images, ['A', 'B', 'B', 'A']):
            imgid = dm.get_imgid(img)
            hdu.header['history'] = "Image '{}' is '{}'".format(imgid, key)

        return hdulist

    def create_proc_hdulist(self, cdata, data_array):
        import astropy.io.fits as fits
        import uuid
        # Copy header of first image
        base_header = cdata[0][0].header.copy()

        hdu = fits.PrimaryHDU(data_array, header=base_header)
        self.set_base_headers(hdu.header)
        hdu.header['EMIRUUID'] = str(uuid.uuid1())
        # Update obsmode in header
        hdu.header['OBSMODE'] = 'LS_ABBA'
        # Headers of last image
        hdu.header['TSUTC2'] = cdata[-1][0].header['TSUTC2']
        result = fits.HDUList([hdu])
        return result

    def create_accum_hdulist(self, cdata, data_array_n,
                             method_name='unkwnow', use_errors=False):
        import uuid

        base_header = cdata[0][0].header.copy()
        hdu = fits.PrimaryHDU(data_array_n[0], header=base_header)
        hdr = hdu.header
        self.set_base_headers(hdr)
        hdu.header['EMIRUUID'] = str(uuid.uuid1())
        hdr['IMGOBBL'] = 0
        hdr['TSUTC2'] = cdata[-1][0].header['TSUTC2']

        hdu.header['history'] = "Combined %d images using '%s'" % (
            len(cdata),
            method_name
        )
        #hdu.header['history'] = 'Combination time {}'.format(
        #    datetime.datetime.utcnow().isoformat()
        #)
        # Update NUM-NCOM, sum of individual frames
        ncom = 0
        for hdul in cdata:
            ncom += hdul[0].header['NUM-NCOM']
        hdr['NUM-NCOM'] = ncom

        #
        if use_errors:
            varhdu = fits.ImageHDU(data_array_n[1], name='VARIANCE')
            num = fits.ImageHDU(data_array_n[2], name='MAP')
            hdulist = fits.HDUList([hdu, varhdu, num])
        else:
            hdulist = fits.HDUList([hdu])

        return hdulist

    def aggregate_result(self, partial_result, rinput):
        obresult = rinput.obresult
        # Check if this is our first run
        naccum = getattr(obresult, 'naccum', 0)
        accum = getattr(obresult, 'accum', None)
        # result to accumulate
        result_key = 'spec_abba'
        field_to_accum = getattr(partial_result, result_key)

        if naccum == 0:
            self.logger.debug('naccum is not set, do not accumulate')
            return partial_result
        elif naccum == 1:
            self.logger.debug('round %d initialize accumulator', naccum)
            newaccum = field_to_accum
        elif naccum > 1:
            self.logger.debug('round %d of accumulation', naccum)
            newaccum = self.aggregate_frames(accum, field_to_accum, naccum)
        else:
            msg = 'naccum set to %d, invalid' % (naccum,)
            self.logger.error(msg)
            raise RecipeError(msg)

        # Update partial result
        partial_result.accum = newaccum

        return partial_result

    def aggregate_frames(self, accum, frame, naccum):
        return self.aggregate2(accum, frame, naccum)

    def aggregate2(self, img1, img2, naccum):

        frames = [img1, img2]
        use_errors = True
        # Initial checks
        fframe = frames[0]
        # Ref image
        img = fframe.open()
        has_num_ext = 'NUM' in img
        has_bpm_ext = 'BPM' in img
        base_header = img[0].header
        baseshape = img[0].shape

        data_hdul = []
        for f in frames:
            img = f.open()
            data_hdul.append(img)

        if has_num_ext:
            self.logger.debug('Using NUM extension')
            masks = [numpy.where(m['NUM'].data, 0, 1).astype('uint8') for m in data_hdul]
        elif has_bpm_ext:
            self.logger.debug('Using BPM extension')
            masks = [m['BPM'].data for m in data_hdul]
        else:
            self.logger.warning('BPM missing, use zeros instead')
            false_mask = numpy.zeros(baseshape, dtype='int16')
            masks = [false_mask for _ in data_hdul]

        self.logger.info('Combine target images (final, aggregate)')

        weight_accum = 2 * (1 - 1.0 / naccum)
        weight_frame = 2.0 / naccum
        self.logger.debug("weights for 'accum' and 'frame', %s", [weight_accum, weight_frame])
        scales = [1.0 / weight_accum, 1.0 / weight_frame]
        method = combine.mean
        data_arr = [hdul[0].data for hdul in data_hdul]
        out = method(data_arr, masks=masks, scales=scales, dtype='float32')

        self.logger.debug('create result image')

        return self.create_accum_hdulist(
            data_hdul,
            out,
            method_name=method.__name__,
            use_errors=False
        )
