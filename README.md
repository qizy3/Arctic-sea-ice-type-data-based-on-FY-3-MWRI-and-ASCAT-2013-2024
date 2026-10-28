# SIType-Arctic: production code

Source code that produces the daily Arctic sea ice type dataset **SIType-Arctic**
(first-year ice / multi-year ice) for the November–April seasons of 2013–2024, on the
12.5 km NSIDC North Polar Stereographic grid (EPSG:3413).

The data archive itself, with the variable definitions and reading examples, is
documented separately in the dataset README delivered with the NetCDF files.
This file only describes the code.

## Contents

| File | Role |
|---|---|
| `config.py` | Input, mask and output locations read by the classification stage |
| `SIT_main.py` | Stage 1: classification, initialisation, warm-air correction, MYI drift correction, ice-edge correction |
| `pre_process_fun.py` | GeoTIFF read/write and valid-pixel extraction |
| `SIT_Classifi_fun.py` | K-means FYI/MYI classification and the spatial initialisation of day one |
| `post_process_fun.py` | The three optimisation steps, one function each |
| `post_NSIDC_SIE.py` | Stage 2: fusion with the NSIDC daily sea ice extent product |
| `Export_SIT_NetCDF.py` | Stage 3: conversion of the fused daily maps to CF-1.11 NetCDF-4 |

Run them in that order: `SIT_main.py`, then `post_NSIDC_SIE.py`, then `Export_SIT_NetCDF.py`.

## How the product is made

The code follows the processing scheme: preprocess, sea ice type classification, then
post-processing, which is split into ice type correction and data gap filling.

### Preprocess

ASCAT grid data and FY-3 MWRI swath data are brought onto the shared 608 × 896 EPSG:3413
grid, which yields a backscatter dataset and a brightness temperature dataset. The
regridding itself is outside this repository.

`pre_process_fun.Extract_ValidPixels` then reduces the daily pair to the pixels that take
part in the classification: all three brightness temperatures present, the land mask
flagging the pixel as ocean, and the backscatter carrying neither the land flag nor the
sensor fill value. The same pass recodes the mask into land, no data, open ocean and valid
pixel, and records which pixels are ice at all, so that a pixel classified as open water
inside a solid ice cover can later be recognised as an omission.

### Sea ice type classification

**Ice-water classification.** A pixel is ice when the MWRI gradient ratio 37V/18V is below
0.045 and 23V/18V below 0.04. Everything else is open water.

**Sea ice type clustering by K-means.** Only the ice pixels inside the climatological MYI
area are split further. Each is described by two features, ASCAT σ⁰(VV) and the MWRI
gradient ratio 37V/18V, both normalised column-wise to 0–1, and k-means separates them into
two groups. Which group is multi-year ice is settled afterwards from three fixed sample
points near the pole, where the ice is multi-year by definition.

The result is the initial sea ice type map, written to `originalResult`.

### Post-processing

#### Ice type correction

Four steps, in `post_process_fun.py`, applied to every day in the order below — the two
middle ones need the previous day's map, which is why they come after the initialisation
and before the edge is cleaned up.

1. **Spatial initialization.** On the first day of an episode there is no previous map to
   work from, so the MYI extent is built from the classification itself: a 13×13 median
   filter over the MYI layer, a 5×5 dilation of the land mask to reserve a buffer nothing
   may grow into, and relabelling of isolated no-data speckles to whatever surrounds them.
   On later days the same initialisation still provides the starting map the next three
   steps refine. Written to `initializedResult`.
2. **Warm-air intrusion correction.** MYI that the classifier lost on a day whose ERA5
   maximum 2 m air temperature reached 0 °C is put back wherever the pixel is still MYI in
   the previous day's corrected map. Once the two maps are more than three days apart, the
   outer band of yesterday's MYI is no longer trusted as a reference: the band is set by an
   erosion of 3, 5 or 7 pixels as the gap in days grows, and only the eroded interior keeps
   the MYI label. Written to `middleResult`.
3. **MYI drift correction.** MYI that appeared since the previous day outside a 5×5
   dilation of the previous day's MYI cannot have drifted there and is demoted to FYI.
4. **Ice edge correction.** Two erosion bands are built around the ice edge, 11×11 for the
   marginal ice zone proper and 9×9 for its outer fringe. MYI in the zone that a 5×5 median
   filter does not keep is an isolated speckle in the fringe and is reverted to FYI.

Together, steps 3 and 4 give the corrected sea ice type map with its FYI, MYI, open water,
land and no-data classes, written to `finalResult`. On the first day of an episode only the
initialisation and the ice edge correction run, there being no previous map to compare with.

#### Data gap filling

The corrected maps are fused with the NSIDC daily sea ice extent product, resampled onto the
same grid beforehand. This fills the gaps the sensor swaths left and removes the classes the
classification had invented:

- ice sitting well inside the cover where extent says open water is dropped;
- open water right at the classified edge where extent says ice becomes first-year ice;
- an empty pixel below 86°N becomes multi-year ice if the annual mean MYI fraction places it
  in the climatological MYI edge zone, and first-year ice otherwise; above 86°N the central
  Arctic gap stays empty;
- the ring of no data that extent keeps along its own land mask is given back an ice class,
  so the corrected ice stays attached to the coast.

The output, `SIT_CorrectedBySIE_YYYYMMDD.tif`, is the final sea ice type product.
`Export_SIT_NetCDF.py` turns it into one CF-1.11 NetCDF-4 file per day, in a folder per
month, carrying the class codes across unchanged.


## Inputs

The code expects these, already resampled onto the product grid:

| Input | `config.py` entry | Notes |
|---|---|---|
| ASCAT sigma-0, V-pol | `SCAT_SOURCE` | one GeoTIFF per day, under `ASCAT/YYYY/YYYYDDD/` |
| FY-3 MWRI Tb 18.7 / 23.8 / 36.5 GHz, V-pol | `PM_SOURCE` | FY-3C for 2013–2018, FY-3D from 2019 |
| ERA5 daily maximum 2 m air temperature | `T2M_SOURCE` | Kelvin |
| NSIDC daily sea ice extent | `SIEPath` in `post_NSIDC_SIE.py` | resampled; extent, not concentration |
| Annual mean MYI fraction | `yearMeanSITPath` in `post_NSIDC_SIE.py` | one `MYI_Probability_<episode>.tif` per half-year |
| Land mask, climatological MYI area, central no-data mask | `LAND_MASK`, `VALID_MYI_MASK`, `CENTRAL_NODATA_MASK` | |

Regridding the native sensor and product files onto this grid, and the yearly averaging
behind the MYI fraction, are separate steps and are not part of this repository.

## Configuration

`config.py` holds the paths for stage 1. Stages 2 and 3 keep their paths at the top of
their own file, since they point at one particular run; the `YYYYMMDD` there has to be
replaced with the folder written by the preceding stage. The year range and the two
satellite assignments are set at the bottom of `SIT_main.py`.

Class codes are fixed at 0 = first-year ice, 1 = multi-year ice, 3 = open water,
4 = land, 5 = no data (also the NetCDF `_FillValue`). Code 2 is unused.
`Export_SIT_NetCDF.py` sets the dataset version through its `VERSION` constant.

## Requirements

Python 3 with GDAL (`osgeo`), NumPy, SciPy, scikit-learn, rasterio, pyproj and netCDF4.
Stage 1 reads and writes GeoTIFF, so the GDAL Python bindings have to be installed
separately from pip/conda packages on most systems.

## Missing and degraded days

24 of the November–April calendar days have no file, because the scatterometer or the
radiometer did not cover the Arctic on those days; the list is in the dataset README.
No placeholder is generated, so a time series should be built from the dates that are
actually present.

On a further 148 of the 1,970 delivered maps the drift correction of step 3 was omitted
during production, on days whose input imagery carried a defect and on the day directly
after it, since the previous-day reference would have carried the same defect. This code
applies step 3 to every day, so re-running it reproduces the archive exactly except on
those 148 days.

## Citation

Qi, Z., Ye, Y., Hui, F., & Cheng, X. (2026). *Arctic sea ice type data based on FY-3 MWRI
and ASCAT, 2013–2024* [Data set]. Zenodo. https://doi.org/10.5281/zenodo.23489985

Ziyun Qi¹, Yufang Ye¹²³ ✉, Fengming Hui¹²³, and Xiao Cheng¹²³

1. School of Geospatial Engineering and Science, Sun Yat-sen University, Zhuhai, 519000, China
2. Southern Marine Science and Engineering Guangdong Laboratory (Zhuhai), Zhuhai, 519000, China
3. Key Laboratory of Comprehensive Observation of Polar Environment (Sun Yat-sen University), Ministry of Education, Zhuhai, 519000, China

✉ Corresponding author: yeyf8@mail.sysu.edu.cn

## License

Creative Commons Attribution 4.0 International (CC BY 4.0). The product is intended for
scientific research and is not valid for navigation.
