# -*- coding: utf-8 -*-
"""
GeoTIFF I/O and valid-pixel extraction used by the classification stage.
"""
import os
from osgeo import gdal
import numpy as np


def Get_Filename(path, filetype):
    """Return the names of every file under `path` carrying the given extension."""
    name = []
    for root, dirs, files in os.walk(path):
        for i in files:
            if filetype in i:
                name.append(i.replace(filetype, ''))
    return [item + filetype for item in name]


def Read_Geotiff(path):
    """Read band 1 of a GeoTIFF, returning the array, the geotransform and the projection."""
    ds = gdal.Open(path)
    geotrans = ds.GetGeoTransform()
    proj = ds.GetProjection()
    arr = ds.GetRasterBand(1).ReadAsArray()
    return arr, geotrans, proj


def Write_Geotiff(arr, rows, cols, geotrans, proj, outFileName):
    """Write a single-band UInt16 GeoTIFF carrying the supplied grid definition."""
    driver = gdal.GetDriverByName("GTiff")
    outdata = driver.Create(outFileName, cols, rows, 1, gdal.GDT_UInt16)
    outdata.SetGeoTransform(geotrans)
    outdata.SetProjection(proj)
    outdata.GetRasterBand(1).WriteArray(arr)
    outdata.FlushCache()


def Extract_ValidPixels(pm18V, pm23V, pm36V, bsvv, maskArray, validMYIMask):
    """
    Reduce the daily rasters to the pixels that take part in the classification.

    A pixel is valid when all three brightness temperatures are present, the land
    mask flags it as ocean, and the backscatter is neither a land flag nor the
    sensor fill value. `maskArray` is recoded in place to 0 = land, -1 = no data,
    1 = ocean outside the ice-covered area, 2 = valid pixel for classification.
    """
    FYGR37_18PRE = (pm36V - pm18V) / (pm36V + pm18V)
    FYGR23_18PRE = (pm23V - pm18V) / (pm23V + pm18V)
    iceIdxPRE = np.where((FYGR37_18PRE < 0.045) & (FYGR23_18PRE < 0.04))
    fullIceMask = maskArray.copy()
    fullIceMask[:] = 0
    fullIceMask[iceIdxPRE] = 1

    if -33 in pm18V:
        validIdx = np.where((pm18V != -33) & (pm23V != -33) & (pm36V != -33)
                            & (maskArray == 1)
                            & (bsvv != 1) & (bsvv != -50))
        oceanIdx = np.where((pm18V != -33) & (pm23V != -33) & (pm36V != -33) & (maskArray == 1))
    elif -50 in pm18V:
        validIdx = np.where((pm18V != -50) & (pm23V != -50) & (pm36V != -50)
                            & (maskArray == 1)
                            & (bsvv != 1) & (bsvv != -50))
        oceanIdx = np.where((pm18V != -50) & (pm23V != -50) & (pm36V != -50) & (maskArray == 1))
    else:
        print('Undefined Nodata value')
        validIdx = np.nan
        oceanIdx = np.nan

    pm18V = pm18V[validIdx]
    pm23V = pm23V[validIdx]
    pm36V = pm36V[validIdx]
    bsvv = bsvv[validIdx]
    validMYIMaskPartial = validMYIMask[validIdx]

    landIdx = np.where(maskArray == 0)
    maskArray[:, :] = -1
    maskArray[landIdx] = 0
    maskArray[oceanIdx] = 1
    maskArray[validIdx] = 2

    pm18V = pm18V.reshape(pm18V.shape[0], 1)
    pm23V = pm23V.reshape(pm18V.shape[0], 1)
    pm36V = pm36V.reshape(pm36V.shape[0], 1)
    bsvv = bsvv.reshape(bsvv.shape[0], 1)
    validMYIMaskPartial = validMYIMaskPartial.reshape(validMYIMaskPartial.shape[0], 1)
    return pm18V, pm23V, pm36V, bsvv, maskArray, fullIceMask, validMYIMaskPartial


def Normalize_byColumn(data):
    """Scale every column of `data` to the 0-1 range."""
    mins = data.min(0)
    maxs = data.max(0)
    ranges = maxs - mins
    row = data.shape[0]
    normData = data - np.tile(mins, (row, 1))
    normData = normData / np.tile(ranges, (row, 1))
    return normData
