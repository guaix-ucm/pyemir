
try:
    import DF


    from numina.ext.gtc import register

    _more = {
             'emir.dataproducts.MasterBias': DF.TYPE_FRAME,
             'emir.dataproducts.MasterDark': DF.TYPE_FRAME,
             'emir.dataproducts.MasterIntensityFlat': DF.TYPE_FRAME, 
             'emir.dataproducts.CoordinateList2DType': DF.TYPE_DOUBLE_ARRAY2D,
             }
    register(_more)
except ImportError:
    pass
