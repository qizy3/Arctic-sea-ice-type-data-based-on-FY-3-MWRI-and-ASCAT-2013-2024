# -*- coding: utf-8 -*-
"""
    @File       : SIT_main.py
    @Description: Main program of the GEOARC-SIT sea ice type product. Runs the
                  classification, the spatial initialisation, the air temperature
                  optimisation, the multi-year ice drift-range optimisation and the
                  marginal ice zone correction, one winter half-year at a time.
"""
import os
import datetime

import config

import pre_process_fun as fc_pre
import SIT_Classifi_fun as fc_classi
import post_process_fun as fc_post


def IsLeapYear(year):
    if str(year)[2:4] == '00':
        return year % 400 == 0
    return year % 4 == 0


def DOYToString(doy):
    """Day of year as the three digit folder name used by the input data."""
    return '%03d' % doy


def DOYToDate(year, doy):
    """Day of year as a YYYYMMDD string."""
    accumDays = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    accumDaysLeap = [0, 31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335]
    AccumArray = accumDaysLeap if IsLeapYear(year) else accumDays
    month = 0
    for accumday in AccumArray:
        if doy <= accumday:
            break
        month += 1
    day = doy - AccumArray[month - 1]
    return '%d%02d%02d' % (year, month, day)


def Process(year, DOYStart, DOYEnd, scatFolder, pmFolder, resultFolder, areaMask_path,
            validMYIMaskPath, CentralNodataMaskPath, t2mPath):
    # Winter half-year the episode belongs to, named by its two calendar years
    if DOYStart > 183:
        yearEpisode = str(year) + '_' + str(year + 1)
    else:
        yearEpisode = str(year - 1) + '_' + str(year)

    # One folder per processing stage: raw classification, initialisation, after the
    # temperature step, and the final map that feeds the NSIDC extent fusion.
    originalResultFolder = os.path.join(resultFolder, 'originalResult')
    initializedResultFolder = os.path.join(resultFolder, 'initializedResult')
    middleResultFolder = os.path.join(resultFolder, 'middleResult')
    finalResultFolder = os.path.join(resultFolder, 'finalResult')

    scatYearFolder = os.path.join(scatFolder, str(year))
    pmYearFolder = os.path.join(pmFolder, str(year))
    t2mYearFolder = os.path.join(t2mPath, str(year))

    ogResultYearFolder = os.path.join(originalResultFolder, yearEpisode)
    itResultYearFolder = os.path.join(initializedResultFolder, yearEpisode)
    mdResultYearFolder = os.path.join(middleResultFolder, yearEpisode)
    fnResultYearFolder = os.path.join(finalResultFolder, yearEpisode)

    # Part 1: classification
    # ----------------------------------------------------------------------------------
    for doy in range(DOYStart, DOYEnd + 1):
        doyStr = DOYToString(doy)
        dateStr = DOYToDate(year, doy)

        scatDOYFolder = os.path.join(scatYearFolder, doyStr)
        pmDOYFolder = os.path.join(pmYearFolder, doyStr)
        ogResultDOYFolder = os.path.join(ogResultYearFolder, doyStr)

        if not os.path.exists(scatDOYFolder) or not os.path.exists(pmDOYFolder):
            print('backScatter or brightTemperature data on ' + str(year) + ' ' + doyStr + 'th is/are missing')
            continue

        fileType = '.tif'
        bsvvList = fc_pre.Get_Filename(scatDOYFolder, fileType)
        tbList = fc_pre.Get_Filename(pmDOYFolder, fileType)

        if 'ASCAT' in scatFolder:
            bsvv_path = os.path.join(scatDOYFolder, bsvvList[0])
        else:
            print('Error, unavailable scatFolder')
            bsvv_path = ''

        if 'AMSR' in pmFolder:
            tb18v_path = os.path.join(pmDOYFolder, tbList[1])
            tb23v_path = os.path.join(pmDOYFolder, tbList[3])
            tb36v_path = os.path.join(pmDOYFolder, tbList[5])
        elif 'FY' in pmFolder:
            tb18v_path = os.path.join(pmDOYFolder, tbList[0])
            tb23v_path = os.path.join(pmDOYFolder, tbList[1])
            tb36v_path = os.path.join(pmDOYFolder, tbList[2])
        else:
            print('Error, unavailable pmFolder')
            tb18v_path = ''
            tb23v_path = ''
            tb36v_path = ''

        [bsvv_array, _, _] = fc_pre.Read_Geotiff(bsvv_path)
        [tb18v_array, _, _] = fc_pre.Read_Geotiff(tb18v_path)
        [tb23v_array, _, _] = fc_pre.Read_Geotiff(tb23v_path)
        [tb36v_array, _, _] = fc_pre.Read_Geotiff(tb36v_path)
        [maskArray, geotrans, proj] = fc_pre.Read_Geotiff(areaMask_path)
        [validMYIMask, _, _] = fc_pre.Read_Geotiff(validMYIMaskPath)
        [CentralNodataMask, _, _] = fc_pre.Read_Geotiff(CentralNodataMaskPath)

        # Keep only the ocean pixels where all four rasters carry a value
        tb18v_ocean, tb23v_ocean, tb36v_ocean, bsvv_ocean, maskArray, fullIceMask, validMYIMaskPartial = \
            fc_pre.Extract_ValidPixels(tb18v_array, tb23v_array, tb36v_array, bsvv_array, maskArray, validMYIMask)

        if len(tb18v_ocean) == 0 or len(tb23v_ocean) == 0 or len(tb36v_ocean) == 0 or len(bsvv_ocean) == 0:
            continue

        [sit_arr, rows, cols] = fc_classi.Type_Classification(bsvv_ocean, tb18v_ocean, tb23v_ocean, tb36v_ocean,
                                                             maskArray, fullIceMask, validMYIMaskPartial,
                                                             CentralNodataMask)

        if not os.path.exists(ogResultDOYFolder):
            os.makedirs(ogResultDOYFolder)
        ogResultPath = os.path.join(ogResultDOYFolder, 'original_SeaIceType_' + dateStr + '.tif')
        fc_pre.Write_Geotiff(sit_arr, rows, cols, geotrans, proj, ogResultPath)
        print(dateStr + ' classified')

    # Part 2: optimisation of the classified maps
    # ----------------------------------------------------------------------------------
    for doy in range(DOYStart, DOYEnd + 1):
        doyStr = DOYToString(doy)
        dateStr = DOYToDate(year, doy)

        ogResultDOYFolder = os.path.join(ogResultYearFolder, doyStr)
        itResultDOYFolder = os.path.join(itResultYearFolder, doyStr)
        mdResultDOYFolder = os.path.join(mdResultYearFolder, doyStr)
        fnResultDOYFolder = os.path.join(fnResultYearFolder, doyStr)
        t2mDOYFolder = os.path.join(t2mYearFolder, doyStr)

        # Look back up to ten days for the previous map that is actually usable
        formerExist = 0
        daysDiff = 0
        doyFormer = -1
        for i in range(1, 11):
            doyFormer = DOYToString(doy - i)

            # In January the first candidate is the last day of the previous calendar year,
            # except when the episode already starts in November, where the folder exists.
            if doy == DOYStart and DOYStart != 304 and DOYStart != 305:
                doyFormer = '366' if IsLeapYear(year - 1) else '365'

            if os.path.exists(os.path.join(ogResultYearFolder, doyFormer)):
                formerExist = 1
                daysDiff = i
                break

        fileType = '.tif'

        if formerExist == 0:
            # First day of the episode: no previous map exists, so the drift-range step
            # has nothing to compare against and only the marginal ice zone runs.
            if not os.path.exists(ogResultDOYFolder):
                continue

            ogList = fc_pre.Get_Filename(ogResultDOYFolder, fileType)
            ogResultPath = os.path.join(ogResultDOYFolder, ogList[0])
            [ogMask, geotrans, proj] = fc_pre.Read_Geotiff(ogResultPath)
            itMask = fc_classi.Initialize_MYI(ogMask)
            itMask = fc_post.Filter_atMIZ(itMask)

            for folder in (itResultDOYFolder, mdResultDOYFolder, fnResultDOYFolder):
                if not os.path.exists(folder):
                    os.makedirs(folder)

            rows = itMask.shape[0]
            cols = itMask.shape[1]
            itResultPath = os.path.join(itResultDOYFolder, 'initialized_SeaIceType_' + dateStr + '.tif')
            mdFirstDayPath = os.path.join(mdResultDOYFolder, 't2mOptimized_SeaIceType_' + dateStr + '.tif')
            fnFirstDayPath = os.path.join(fnResultDOYFolder, 'SIT_' + dateStr + '.tif')
            fc_pre.Write_Geotiff(itMask, rows, cols, geotrans, proj, itResultPath)
            fc_pre.Write_Geotiff(itMask, rows, cols, geotrans, proj, mdFirstDayPath)
            fc_pre.Write_Geotiff(itMask, rows, cols, geotrans, proj, fnFirstDayPath)
            print(dateStr + ' initialised')

        elif formerExist == 1:
            # ----- initialisation ----- #
            if not os.path.exists(ogResultDOYFolder):
                continue

            ogList = fc_pre.Get_Filename(ogResultDOYFolder, fileType)
            ogResultPath = os.path.join(ogResultDOYFolder, ogList[0])
            [ogMask, geotrans, proj] = fc_pre.Read_Geotiff(ogResultPath)
            itMask = fc_classi.Initialize_MYI(ogMask)

            if not os.path.exists(itResultDOYFolder):
                os.makedirs(itResultDOYFolder)
            itResultPath = os.path.join(itResultDOYFolder, 'initialized_SeaIceType_' + dateStr + '.tif')
            rows = itMask.shape[0]
            cols = itMask.shape[1]
            fc_pre.Write_Geotiff(itMask, rows, cols, geotrans, proj, itResultPath)
            print(dateStr + ' initialised')

            # ----- air temperature optimisation ----- #
            mdResultFormerDOYFolder = os.path.join(mdResultYearFolder, doyFormer)
            mdFormerList = fc_pre.Get_Filename(mdResultFormerDOYFolder, fileType)
            mdFormerPath = os.path.join(mdResultFormerDOYFolder, mdFormerList[0])
            [mdFormerMask, _, _] = fc_pre.Read_Geotiff(mdFormerPath)

            t2mList = fc_pre.Get_Filename(t2mDOYFolder, fileType)
            t2mPath = os.path.join(t2mDOYFolder, t2mList[0])
            [t2mArray, _, _] = fc_pre.Read_Geotiff(t2mPath)

            itList = fc_pre.Get_Filename(itResultDOYFolder, fileType)
            itPath = os.path.join(itResultDOYFolder, itList[0])
            [itMask, geotrans, proj] = fc_pre.Read_Geotiff(itPath)

            mdMask = fc_post.CorrectionOnT2m(mdFormerMask, itMask, t2mArray, daysDiff)

            if not os.path.exists(mdResultDOYFolder):
                os.makedirs(mdResultDOYFolder)
            mdResultPath = os.path.join(mdResultDOYFolder, 't2mOptimized_SeaIceType_' + dateStr + '.tif')
            rows = mdMask.shape[0]
            cols = mdMask.shape[1]
            fc_pre.Write_Geotiff(mdMask, rows, cols, geotrans, proj, mdResultPath)
            print(dateStr + ' temperature optimisation done')

            # ----- final optimisation ----- #
            # Drift range first, then the marginal ice zone
            fnResultFormerDOYFolder = os.path.join(fnResultYearFolder, doyFormer)
            fnFormerList = fc_pre.Get_Filename(fnResultFormerDOYFolder, fileType)
            fnFormerPath = os.path.join(fnResultFormerDOYFolder, fnFormerList[0])
            [fnFormerMask, _, _] = fc_pre.Read_Geotiff(fnFormerPath)

            mdList = fc_pre.Get_Filename(mdResultDOYFolder, fileType)
            mdPath = os.path.join(mdResultDOYFolder, mdList[0])
            [mdMask, geotrans, proj] = fc_pre.Read_Geotiff(mdPath)

            fnMask = fc_post.Filter_onMovement(fnFormerMask, mdMask)
            fnMask = fc_post.Filter_atMIZ(fnMask)

            if not os.path.exists(fnResultDOYFolder):
                os.makedirs(fnResultDOYFolder)
            fnResultPath = os.path.join(fnResultDOYFolder, 'SIT_' + dateStr + '.tif')
            rows = fnMask.shape[0]
            cols = fnMask.shape[1]
            fc_pre.Write_Geotiff(fnMask, rows, cols, geotrans, proj, fnResultPath)
            print(dateStr + ' finally optimised')


# Start from here
if __name__ == '__main__':
    run_date = datetime.datetime.now().strftime('%Y%m%d')
    resultFolder = os.path.join(config.RESULT_MASK_DIR, run_date)
    if not os.path.exists(resultFolder):
        os.makedirs(resultFolder)

    scatSourceFolder = config.SCAT_SOURCE
    pmSourceFolder = config.PM_SOURCE

    areaMask_path = config.LAND_MASK
    validMYIMaskPath = config.VALID_MYI_MASK
    CentralNodataMaskPath = config.CENTRAL_NODATA_MASK
    t2mPath = config.T2M_SOURCE

    yearBegin = 2023
    yearEnd = 2024

    # Each episode is processed as its second half first, then the following first half,
    # because the optimisation always looks back at the day before.
    for year in range(yearBegin, yearEnd + 1):
        for half in ['yearLower', 'yearUpper']:
            if half == 'yearUpper':
                year += 1

            leapYear = 1 if IsLeapYear(year) else 0

            # November to April of the winter half-year, as day of year numbers
            if half == 'yearLower':
                DOYStart = 305 if leapYear else 304
                DOYEnd = 366 if leapYear else 365
            else:
                DOYStart = 1
                DOYEnd = 121 if leapYear else 120

            if year < 2013:
                print('No input data configured for ' + str(year))
                continue

            scatFolder = scatSourceFolder
            # FY-3C MWRI covers 2013-2018, FY-3D MWRI from 2019 onwards
            satellite = '3C' if year < 2019 else '3D'
            pmFolder = os.path.join(pmSourceFolder, 'FYData', satellite)

            Process(year, DOYStart, DOYEnd, scatFolder, pmFolder, resultFolder, areaMask_path,
                    validMYIMaskPath, CentralNodataMaskPath, t2mPath)

    print('All processes were done.')
