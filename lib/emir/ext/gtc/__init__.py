
import Whatever


from numina.ext.gtc import register

_more = {
'emir.dataproducts.MasterBias': 'DPK.DPKFrame'
'emir.dataproducts.MasterDark': 'DPK.FRame',
'emir.dataproducts.MasterIntensityFlat': 'DPK.FRame', 
'emir.dataproducts.CoordinateList2DType': 'DPK.Array2D',
}

register(_more)
