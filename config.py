# -*- coding: utf-8 -*-
"""
Global paths and shared parameters for the SIType-Arctic product generation.

SIT_main.py takes its inputs from here. The two later stages, post_NSIDC_SIE.py and
Export_SIT_NetCDF.py, carry their own run specific paths at the top of the file instead.
"""
import os

ROOT_DIR = r"D:\SITY"

# ---------------------------------------------------------------------------
# Input data
# ---------------------------------------------------------------------------
PRE_PROCESSED_DATA = os.path.join(ROOT_DIR, r"Pre_processed_Data\DataSeaiceClassification")
SCAT_SOURCE = os.path.join(PRE_PROCESSED_DATA, "ASCAT")   # backscatter, MetOp/ASCAT
PM_SOURCE = PRE_PROCESSED_DATA                            # brightness temperature, FY-3 MWRI
T2M_SOURCE = os.path.join(PRE_PROCESSED_DATA, "ERA5")     # ERA5 2 m air temperature

# ---------------------------------------------------------------------------
# Masks
# ---------------------------------------------------------------------------
MASK_DIR = os.path.join(ROOT_DIR, "Mask")
LAND_MASK = os.path.join(MASK_DIR, "LandMask.tif")
VALID_MYI_MASK = os.path.join(MASK_DIR, "ValidMYIMask.tif")
CENTRAL_NODATA_MASK = os.path.join(MASK_DIR, "CentralNodataMask.tif")

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
RESULT_MASK_DIR = os.path.join(ROOT_DIR, "result", "ResultMask")
