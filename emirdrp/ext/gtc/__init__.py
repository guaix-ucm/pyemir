
try:
    import DF

    from numina.ext.gtc import register

    _more = {'emirdrp.products.MasterBias': DF.TYPE_FRAME,
             'emirdrp.products.MasterDark': DF.TYPE_FRAME,
             'emirdrp.products.MasterIntensityFlat': DF.TYPE_FRAME,
             'emirdrp.products.CoordinateList2DType': DF.TYPE_DOUBLE_ARRAY2D,
             }
    register(_more)
except ImportError:
    pass
