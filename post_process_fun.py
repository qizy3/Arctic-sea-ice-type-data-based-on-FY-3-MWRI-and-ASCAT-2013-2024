# -*- coding: utf-8 -*-
"""
The three post-classification optimisations applied to every daily map.
"""
import numpy as np
import scipy.signal as signal
from scipy import ndimage


def CorrectionOnT2m(formerMask, mask, t2m, daysDiff):
    """
    Undo multi-year ice that was lost because the air temperature was above freezing.

    `daysDiff` is how many days separate the reference map from the one being corrected.
    The larger the gap, the coarser the erosion kernel used to describe the band the
    marginal ice zone may have advanced through, so that only pixels close to the interior
    are trusted as still being multi-year ice.
    """
    if 3 < daysDiff <= 6:
        kernelSize = 3
    elif 6 < daysDiff < 9:
        kernelSize = 5
    else:
        kernelSize = 7

    myiBinary = formerMask.copy()
    myiBinary[np.where(myiBinary != 1)] = 0

    # Outer band of yesterday's multi-year ice, treated as uncertain drift range
    erosionKernel = np.ones((kernelSize, kernelSize), dtype=int)
    erodedMYI = ndimage.binary_erosion(myiBinary, structure=erosionKernel).astype(erosionKernel.dtype)

    # Pixels that survive the erosion keep their multi-year ice label, the band loses it
    tempMask = formerMask.copy()
    if daysDiff > 3:
        tempMask[np.where((myiBinary == 1) & (erodedMYI == 0))] = -1

    # Today's first-year ice inside yesterday's trusted multi-year ice is put back
    mask[np.where((tempMask == 1) & (mask == 0) & (t2m >= 270.15))] = 1
    return mask


def Filter_onMovement(optimizedMask, rawMask):
    """
    Reject multi-year ice that appeared between two days without room to have drifted in.

    Yesterday's final map is dilated by 5x5; a pixel that falls outside that envelope but
    is nevertheless labelled multi-year ice by today's raw classification can only be a
    misclassification, so it is demoted to first-year ice.
    """
    LandAndNodataIdx = np.where((rawMask == 4) | (rawMask == 5))

    myiSmoothIdx = np.where(optimizedMask == 1)
    binaryMaskOptimized = optimizedMask.copy()
    binaryMaskOptimized[:, :] = 0
    binaryMaskOptimized[myiSmoothIdx] = 1
    kernel = np.ones((5, 5), dtype=int)
    dilatedOptimized = ndimage.binary_dilation(binaryMaskOptimized, structure=kernel).astype(kernel.dtype)

    myiIdxRaw = np.where(rawMask == 1)
    binaryMaskRaw = rawMask.copy()
    binaryMaskRaw[:, :] = 0
    binaryMaskRaw[myiIdxRaw] = 1

    # Where the predicted envelope falls short of what the raw classification already claims
    dilatedZone = dilatedOptimized - binaryMaskRaw
    wrongMyiIdx = np.where(dilatedZone == -1)

    movementOptimized = rawMask.copy()
    movementOptimized[wrongMyiIdx] = 0
    movementOptimized[LandAndNodataIdx] = rawMask[LandAndNodataIdx]
    return movementOptimized


def Filter_atMIZ(maskSmooth):
    """
    Remove multi-year ice wrongly sitting in the marginal ice zone.

    Two erosion bands are built around the ice edge: the 11x11 one marks the marginal ice
    zone itself, the 9x9 one the outer fringe where errors are judged. Within the outer
    fringe, multi-year ice pixels that a 5x5 median filter does not keep are isolated
    speckles and are reverted to first-year ice.
    """
    iceIdx = np.where((maskSmooth == 4) | (maskSmooth == 0) | (maskSmooth == 1))
    LandAndNodataIdx = np.where((maskSmooth == 4) | (maskSmooth == 5))

    binaryMask = maskSmooth.copy()
    binaryMask[:, :] = 0
    binaryMask[iceIdx] = 1

    largeKernel = np.ones((11, 11), dtype=int)
    largeEroded = ndimage.binary_erosion(binaryMask, structure=largeKernel).astype(largeKernel.dtype)
    largeMarginZone = binaryMask - largeEroded

    smallKernel = np.ones((9, 9), dtype=int)
    smallEroded = ndimage.binary_erosion(binaryMask, structure=smallKernel).astype(smallKernel.dtype)
    smallMarginZone = binaryMask - smallEroded

    myiMZIdx = np.where((maskSmooth == 1) & (largeMarginZone == 1))
    marginalZone = maskSmooth.copy()
    marginalZone[:, :] = 0
    marginalZone[myiMZIdx] = 1

    MZ_medfilted = signal.medfilt(marginalZone, (5, 5))
    noiseMZ = marginalZone - MZ_medfilted
    mismyi_idx = np.where((noiseMZ == 1) & (smallMarginZone == 1))

    MZ_Optimized = maskSmooth.copy()
    MZ_Optimized[mismyi_idx] = 0
    MZ_Optimized[LandAndNodataIdx] = maskSmooth[LandAndNodataIdx]
    return MZ_Optimized
