
try:
    import DF

    from numina.ext.gtc import register

    _more = {'emirdrp.dataproducts.MasterBias': DF.TYPE_FRAME,
             'emirdrp.dataproducts.MasterDark': DF.TYPE_FRAME,
             'emirdrp.dataproducts.MasterIntensityFlat': DF.TYPE_FRAME,
             'emirdrp.dataproducts.CoordinateList2DType': DF.TYPE_DOUBLE_ARRAY2D,
             }
    register(_more)
except ImportError:
    pass
