# -*- coding: utf-8 -*-
"""
    @File       : post_NSIDC_SIE.py
    @Description: Fuse the daily GEOARC-SIT classification with the NSIDC sea ice extent
                  product. The extent product is read on the same 608x896 polar
                  stereographic grid, so pixels the classifier left empty because the
                  radiometer or the scatterometer did not cover them can still be given
                  an ice class, and ice that the classifier invented where the extent
                  product says open water is reverted.
    Paths below point at one particular run and have to be re-pointed for the next one.
"""
import os
import datetime
from osgeo import gdal, osr
import numpy as np
from scipy import ndimage

# Final maps from SIT_main.py, one folder per winter half-year
SITPath = r'D:\SITY\result\ResultMask\YYYYMMDD\finalResult'
# NSIDC extent resampled onto the product grid
SIEPath = r'D:\SITY\Pre_processed_Data\DataSeaiceClassification\SIE'
# Annual mean multi-year ice fraction, source of the climatological edge zone. One
# MYI_Probability_<episode>.tif per winter half-year, as written by the yearly comparison.
yearMeanSITPath = r'D:\SITY\analysis_own\AnalysisResult_YperY\YYYYMMDD'
run_date = datetime.datetime.now().strftime('%Y%m%d')
outPath = os.path.join(r'D:\SITY\result\Result_CorrectedBySIE', run_date)
if not os.path.exists(outPath):
    os.makedirs(outPath)

landMaskPath = r"D:\SITY\Mask\LandMask.tif"
landMaskNewPath = r"D:\SITY\Mask\landMaskNew.tif"

# Where the two land masks disagree, the SIT has to fall back to no data.
# In the difference raster the disagreeing pixels are 1.
geoData = gdal.Open(landMaskPath)
geotrans = geoData.GetGeoTransform()
proj = geoData.GetProjection()
landMask = geoData.GetRasterBand(1).ReadAsArray()
landMask = np.array(landMask, dtype=int)
landMaskNewFile = gdal.Open(landMaskNewPath)
landMaskNew = landMaskNewFile.GetRasterBand(1).ReadAsArray()
landMaskDiff = np.zeros((896, 608), dtype=int)
landMaskDiff[np.where((landMask == 1) & (landMaskNew == 0))] = 1


def lonlat_obtain(gdalVar, raster):
    """Per-pixel longitude and latitude of the grid covered by `raster`."""
    old_cs = osr.SpatialReference()
    old_cs.ImportFromWkt(gdalVar.GetProjectionRef())
    # NSIDC grid is defined on the Hughes 1980 ellipsoid, not on WGS 84
    new_cs_wkt = """
    GEOGCS["Hughes 80",
        DATUM["Hughes_1980",
            SPHEROID["Hughes_1980",6378273,298.279411123,
                AUTHORITY["EPSG","7058"]],
            AUTHORITY["EPSG","6054"]],
        PRIMEM["Greenwich",0,
            AUTHORITY["EPSG","8901"]],
        UNIT["degree",0.01745329251994328,
            AUTHORITY["EPSG","9122"]],
        AUTHORITY["EPSG","4054"]]"""
    new_cs = osr.SpatialReference()
    new_cs.ImportFromWkt(new_cs_wkt)

    transform = osr.CoordinateTransformation(old_cs, new_cs)
    nrows, ncols = raster.shape

    lons = np.zeros((nrows, ncols))
    lats = np.zeros((nrows, ncols))

    for i in range(0, nrows):
        for j in range(0, ncols):
            px = geotrans[0] + j * geotrans[1] + i * geotrans[2]
            py = geotrans[3] + j * geotrans[4] + i * geotrans[5]
            coords = transform.TransformPoint(px, py)
            lons[i, j] = coords[1]
            lats[i, j] = coords[0]
    return lons, lats


lons, lats = lonlat_obtain(geoData, landMask)

for _, yearEpisodes, _ in os.walk(SITPath):
    for yearEpisode in yearEpisodes:
        # Keep only folders named like 2023_2024
        if not (len(yearEpisode) == 9 and yearEpisode[:4].isdigit()
               and yearEpisode[5:].isdigit() and yearEpisode[4] == '_'):
            continue

        # Climatological multi-year ice edge zone of this episode: pixels around the
        # band where multi-year ice is found in 85% of the yearly maximum fraction maps.
        breakFlag = 0
        MYIPFile = ''
        for _, _, MYIPFiles in os.walk(yearMeanSITPath):
            if breakFlag == 1:
                break
            for MYIPFile in MYIPFiles:
                if yearEpisode in MYIPFile and 'MYI' in MYIPFile:
                    breakFlag = 1
                    break
        if 'MYI' not in MYIPFile:
            continue
        MYIPData = gdal.Open(os.path.join(yearMeanSITPath, MYIPFile))
        MYIP = MYIPData.GetRasterBand(1).ReadAsArray()
        MYIP = np.array(MYIP, dtype=int)
        MYIPMax = np.max(MYIP)
        MYIBinary = MYIP.copy()
        MYIBinary[:] = 0
        MYIEdgeZone = MYIBinary.copy()
        MYIBinary[np.where(MYIP >= int(MYIPMax * 0.85))] = 1
        MYIDilationKernel = np.ones((5, 5), dtype=int)
        dilatedMYI = ndimage.binary_dilation(MYIBinary, structure=MYIDilationKernel).astype(MYIDilationKernel.dtype)
        MYIEdgeZone[np.where((MYIBinary == 0) & (dilatedMYI == 1))] = 1

        for _, DOYs, _ in os.walk(os.path.join(SITPath, yearEpisode)):
            for DOY in DOYs:
                doyDirPath = os.path.join(SITPath, yearEpisode, DOY)
                for _, _, SITFiles in os.walk(doyDirPath):
                    for SITFile in SITFiles:
                        SITFilePath = os.path.join(doyDirPath, SITFile)
                        SITData = gdal.Open(SITFilePath)
                        if SITData is None:
                            print(f"File not found: {SITFilePath}")
                            continue
                        SIT = SITData.GetRasterBand(1).ReadAsArray()
                        geotrans = SITData.GetGeoTransform()
                        proj = SITData.GetProjection()
                        SIT = np.array(SIT, dtype=int)
                        date = SITFile[-12:-4]
                        # 1031 is the artificial end-of-season map, not a real day
                        if date[4:8] == '1031':
                            continue

                        SIEYEPath = os.path.join(SIEPath, yearEpisode)
                        breakFlag = 0
                        SIEFile = ''

                        for _, _, SIEFiles in os.walk(SIEYEPath):
                            if breakFlag == 1:
                                breakFlag = 0
                                break
                            for SIEFile in SIEFiles:
                                if date in SIEFile:
                                    breakFlag = 1
                                    break

                        SIE = []
                        if date in SIEFile:
                            SIEData = gdal.Open(os.path.join(SIEYEPath, SIEFile))
                            SIE = SIEData.GetRasterBand(1).ReadAsArray()
                            SIE = np.array(SIE, dtype=int)
                        if date not in SIEFile:
                            continue

                        # Ice edge band around today's classified ice, and, one erosion
                        # step further in, the area where an overestimation can be excused
                        binaryMask = SIT.copy()
                        binaryMask[np.where((SIT == 0) | (SIT == 1))] = 1
                        binaryMask[np.where(binaryMask != 1)] = 0
                        Kernel = np.ones((5, 5), dtype=int)
                        Kernel2 = np.ones((10, 10), dtype=int)
                        iceEdgeEroded = ndimage.binary_erosion(binaryMask, structure=Kernel).astype(Kernel.dtype)
                        iceEdgeDilated = ndimage.binary_dilation(binaryMask, structure=Kernel).astype(Kernel.dtype)
                        iceEdgeEroded2 = ndimage.binary_erosion(binaryMask, structure=Kernel).astype(Kernel2.dtype)
                        arcExcludedIce = ndimage.binary_dilation(iceEdgeEroded2, structure=Kernel).astype(Kernel2.dtype)
                        iceEdgeZone = SIT.copy()
                        iceEdgeZone[:] = 0
                        iceEdgeZone[np.where((iceEdgeEroded == 0) & (iceEdgeDilated == 1) & (lats < 86))] = 1

                        # Land mask disagreement becomes no data
                        SIT[np.where(landMaskDiff == 1)] = 5

                        SITCorr = SIT.copy()
                        SITCorr[:] = 3
                        SITCorr[np.where(SIT == 1)] = 1
                        SITCorr[np.where(SIT == 0)] = 0

                        # Extent says open water well inside the ice cover: the class is
                        # dropped. Extent says ice right at the classified edge: open
                        # water is turned into first-year ice.
                        SITCorr[np.where(((SIT == 0) | (SIT == 1)) & (SIE == 0) & (arcExcludedIce == 0))] = 3
                        SITCorr[np.where((SIT == 3) & (SIE == 1) & (iceEdgeZone == 1))] = 0

                        # Empty pixels inside the climatological multi-year ice edge zone
                        # take multi-year ice, the rest below 86N takes first-year ice
                        SITCorr[np.where((MYIEdgeZone == 1) & (SIT == 5) & (lats < 86) & (SIE == 1))] = 1
                        SITCorr[np.where((SIT == 5) & (SIE == 1) & (lats < 86) & (SITCorr != 1))] = 0

                        # The extent product keeps a ring of no data along its own land
                        # mask, which would leave the corrected ice detached from the
                        # coast. Ice dilated into that ring is given back a class.
                        ENBinary = SIE.copy()
                        ENBinary[:] = 0
                        ENBinary[np.where(SIE == 253)] = 1
                        TIBinary = SITCorr.copy()
                        TIBinary[:] = 0
                        TMBinary = TIBinary.copy()
                        TIBinary[np.where((SITCorr == 0) | (SITCorr == 1))] = 1
                        TMBinary[np.where(SITCorr == 1)] = 1
                        Kernel3 = np.ones((5, 5))
                        TIDilated = ndimage.binary_dilation(TIBinary, structure=Kernel3).astype(Kernel3.dtype)
                        TMDilated = ndimage.binary_dilation(TMBinary, structure=Kernel3).astype(Kernel3.dtype)
                        SITCorr[np.where((TMDilated == 1) & (ENBinary == 1))] = 1
                        SITCorr[np.where((TMDilated == 0) & (TIDilated == 1) & (ENBinary == 1))] = 0

                        SITCorr[np.where((SIE == 254) | (SIT == 4))] = 4
                        # Central Arctic gap stays empty
                        SITCorr[np.where((SIT == 5) & (lats >= 86))] = 5

                        driver = gdal.GetDriverByName("GTiff")
                        if not os.path.exists(os.path.join(outPath, yearEpisode)):
                            os.makedirs(os.path.join(outPath, yearEpisode))
                        outFilePath = os.path.join(outPath, yearEpisode, 'SIT_CorrectedBySIE_' + date + '.tif')
                        outdata = driver.Create(outFilePath, 608, 896, 1, gdal.GDT_UInt16)
                        outdata.SetGeoTransform(geotrans)
                        outdata.SetProjection(proj)
                        outdata.GetRasterBand(1).WriteArray(SITCorr)
                        outdata.FlushCache()

                        print(date + ' is done.')
