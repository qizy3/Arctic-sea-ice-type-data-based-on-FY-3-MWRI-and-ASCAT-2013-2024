# -*- coding: utf-8 -*-
"""
Sea ice type classification and the spatial initialisation of day one.

Class codes carried through the whole workflow:
  0 = first-year ice, 1 = multi-year ice, 3 = open water, 4 = land, 5 = no data.
"""
import numpy as np
from sklearn.cluster import KMeans
import scipy.signal as signal
from scipy import ndimage

from pre_process_fun import Normalize_byColumn


def Type_Classification(bsvv_ocean, pm18v_ocean, pm23v_ocean, pm36v_ocean, maskArray, fullIceMask,
                        validMYIMaskPartial, CentralNodataMask):
    """
    Split the ice-covered pixels into first-year and multi-year ice.

    Ice is first separated from open water with two brightness temperature gradients,
    then the backscatter together with the 36/18 GHz gradient is split into two clusters
    by k-means inside the climatological multi-year ice area. Which of the two clusters
    is multi-year ice is decided afterwards from three fixed sample points near the pole.
    """
    # MWRI gradient ratios: 37V/18V below 0.045 and 23V/18V below 0.04 means ice
    FYGR37_18 = (pm36v_ocean - pm18v_ocean) / (pm36v_ocean + pm18v_ocean)
    FYGR23_18 = (pm23v_ocean - pm18v_ocean) / (pm23v_ocean + pm18v_ocean)
    fullIceIdx = np.where((FYGR37_18 < 0.045) & (FYGR23_18 < 0.04))
    iceIdx = np.where((FYGR37_18 < 0.045) & (FYGR23_18 < 0.04) & (validMYIMaskPartial == 1))

    pmGR36_18 = (pm36v_ocean - pm18v_ocean) / (pm36v_ocean + pm18v_ocean)
    bsvv_ice = bsvv_ocean[iceIdx]
    bsvv_ice = bsvv_ice.reshape(bsvv_ice.shape[0], 1)
    pmGR_ice = pmGR36_18[iceIdx]
    pmGR_ice = pmGR_ice.reshape(pmGR_ice.shape[0], 1)

    # Two-column sample matrix, one column per feature
    bsvv_n_ice = Normalize_byColumn(bsvv_ice)
    pmGR_n_ice = Normalize_byColumn(pmGR_ice)
    sample_n_matrix = np.hstack((bsvv_n_ice, pmGR_n_ice))
    kmeans = KMeans(n_clusters=2, random_state=0).fit(sample_n_matrix)

    # Flat index of the four background classes
    [rows, cols] = maskArray.shape
    std_vector = maskArray.reshape(rows * cols, 1)  # land 0, no data -1, open ocean 1, valid 2
    land_idx = np.where(std_vector == 0)
    nodata_idx = np.where(std_vector == -1)
    trueOcean = np.where(std_vector == 1)
    ocean_idx = np.where(std_vector == 2)

    sit_data = std_vector.copy()
    sit_data[:, 0] = 0

    # Start from open water, mark every ice pixel 10, then overwrite the clustered pixels
    sit_ocean = pmGR36_18.copy()
    sit_ocean[:, :] = 3
    sit_ocean[fullIceIdx] = 10
    sit_ocean[iceIdx] = kmeans.labels_

    sit_data[nodata_idx] = 5
    sit_data[land_idx] = 4
    sit_data[trueOcean] = 3
    sit_data[ocean_idx[0]] = sit_ocean

    # Three sample points near the pole are multi-year ice by definition. If the
    # clustering assigned them to 0, the two labels have to be swapped.
    hitPosition1 = int((896 / 2) * 608 + (608 / 16 * 7) + 608 * 50)
    hitPosition2 = int((896 / 2) * 608 + (608 / 32 * 12) + 608 * 20)
    hitPosition3 = int((896 / 2) * 608 + (608 / 32 * 16) + 608 * 60)
    status = \
        sit_data[hitPosition1, 0] * sit_data[hitPosition2, 0] + sit_data[hitPosition1, 0] * sit_data[hitPosition3, 0]\
        + sit_data[hitPosition2, 0] * sit_data[hitPosition3, 0]
    if status == 0:
        fyi_idx = np.where(sit_data == 1)
        myi_idx = np.where(sit_data == 0)
        sit_data[fyi_idx] = 0
        sit_data[myi_idx] = 1
    sit_data[np.where(sit_data == 10)] = 0

    sit_data = sit_data.reshape(rows, cols)
    test_data = sit_data.copy()
    # Ice-classified pixels that are not inside the land mask become no data
    sit_data[np.where((test_data == 3) & (fullIceMask == 1))] = 5
    sit_data[np.where(CentralNodataMask == 1)] = 5

    return sit_data, rows, cols


def Initialize_MYI(mask):
    """
    Build the multi-year ice extent of the first day of the episode.

    A 13x13 median filter on the binary multi-year ice layer gives a spatially coherent
    initial extent, and a 5x5 dilation of the land layer keeps a buffer where nothing is
    allowed to grow. Remaining isolated no data speckles are then absorbed into the
    neighbourhood they sit in.
    """
    land_idx = np.where(mask == 4)
    maskSmooth = signal.medfilt(mask, (3, 3))
    myi_idx = np.where(maskSmooth == 1)
    fyi_idx = np.where(maskSmooth == 0)
    nodata_idx = np.where(maskSmooth == 5)

    maskMedian = maskSmooth.copy()
    maskMedian[:, :] = 0
    maskMedian[myi_idx] = 1
    maskMedian = signal.medfilt(maskMedian, (13, 13))
    myi_init_idx = np.where(maskMedian == 1)

    kernelLand = np.ones((5, 5), dtype=int)
    LandMask = maskSmooth.copy()
    LandMask[:, :] = 0
    LandMask[land_idx] = 1
    dilatedLand = ndimage.binary_dilation(LandMask, structure=kernelLand).astype(kernelLand.dtype)
    preserveZone_idx = np.where((dilatedLand - LandMask) == 1)

    maskDoubleSmooth = signal.medfilt(mask, (13, 13))
    maskDoubleSmooth[fyi_idx] = 0
    maskDoubleSmooth[myi_idx] = 0
    maskDoubleSmooth[myi_init_idx] = 1
    maskDoubleSmooth[land_idx] = 4
    maskDoubleSmooth[nodata_idx] = 5
    maskDoubleSmooth[preserveZone_idx] = 5

    # Comparing the two layouts at each pixel tells whether a no data pixel is surrounded
    # by open water or by ice, and can therefore be relabelled instead of staying empty.
    maskToRemove = mask.copy()
    maskToRemove[:, :] = 500
    maskToRemove[np.where(mask == 5)] = 100
    maskRemoved = maskDoubleSmooth.copy()
    maskRemoved[:, :] = 500
    maskRemoved[np.where(maskDoubleSmooth == 3)] = 50
    maskRemoved[np.where((maskDoubleSmooth == 0) | (maskDoubleSmooth == 1))] = 15
    Difference = maskToRemove - maskRemoved
    NodataToOWIdx = np.where(Difference == 50)
    NodataToFMYIIdx = np.where(Difference == 85)
    mask[NodataToOWIdx] = maskDoubleSmooth[NodataToOWIdx]
    mask[NodataToFMYIIdx] = maskDoubleSmooth[NodataToFMYIIdx]
    mask[np.where(maskDoubleSmooth == 5)] = 5
    return mask
