# -*- coding: utf-8 -*-
"""
    @File       : Export_SIT_NetCDF.py
    @Description: Convert the daily GeoTIFF maps that came out of the SIE fusion step
                  (SIType-Arctic, produced by the GEOARC-SIT workflow) into CF-1.11
                  NetCDF. One file per day, sorted into one folder per month.
"""
import os
import re
import sys
import glob
import datetime

import numpy as np
import rasterio
from pyproj import Transformer
import netCDF4 as nc

# Input: the daily maps written by post_NSIDC_SIE.py, one sub-folder per winter half-year.
# Replace YYYYMMDD with the folder of the run being converted.
SRC_DIR = r'D:\SITY\result\Result_CorrectedBySIE\YYYYMMDD'
OUT_DIR = r'D:\SITY\result\SIType-Arctic'
FILE_STEM = 'SIType-Arctic'


def nc_rel_path(yyyymmdd):
    """SIType-Arctic/YYYYMM/SIType-Arctic-YYYYMMDD.nc"""
    return os.path.join(yyyymmdd[:6], f'{FILE_STEM}-{yyyymmdd}.nc')

# Class codes, kept in step with SIT_Classifi_fun.Type_Classification and post_NSIDC_SIE.py
FLAG_VALUES = [0, 1, 3, 4]
FLAG_MEANINGS = 'first_year_ice multi_year_ice open_water land'
FILL_VALUE = 5

VERSION = '20261022'                  # NetCDF dataset version
FORCE = '--force' in sys.argv

DATE_RE = re.compile(r'SIT_CorrectedBySIE_(\d{8})\.tif$')


def collect_dates():
    months = {}
    for path in glob.glob(os.path.join(SRC_DIR, '*', 'SIT_CorrectedBySIE_*.tif')):
        m = DATE_RE.search(os.path.basename(path))
        if not m:
            continue
        date = m.group(1)
        months.setdefault(date[:6], []).append((date, path))
    return {k: sorted(v) for k, v in sorted(months.items())}


def grid_geometry(ref_path):
    """Read the grid definition from a reference GeoTIFF: x/y coordinates, EPSG code and 2D lat/lon."""
    with rasterio.open(ref_path) as ds:
        transform = ds.transform
        epsg = ds.crs.to_epsg()
        rows, cols = ds.height, ds.width

    x = np.array([transform.c + transform.a * (i + 0.5) for i in range(cols)], dtype='f8')
    y = np.array([transform.f + transform.e * (j + 0.5) for j in range(rows)], dtype='f8')

    tr = Transformer.from_crs(f'EPSG:{epsg}', 'EPSG:4326', always_xy=True)
    xv, yv = np.meshgrid(x, y)
    lon, lat = tr.transform(xv, yv)
    return x, y, epsg, rows, cols, lat.astype('f4'), lon.astype('f4'), transform


def projection_attrs(epsg):
    return {
        'grid_mapping_name': 'polar_stereographic',
        'latitude_of_projection_origin': 90.0,
        'straight_vertical_longitude_from_pole': -45.0,
        'standard_parallel': 70.0,
        'false_easting': 0.0,
        'false_northing': 0.0,
        'scale_factor_at_projection_origin': 1.0,
        'long_name': 'CF grid mapping - NSIDC Sea Ice Polar Stereographic North',
        'epsg_code': f'EPSG:{epsg}',
        'crs_wkt': 'PROJCS["WGS 84 / NSIDC Sea Ice Polar Stereographic North",'
                   'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563]],'
                   'PRIMEM["Greenwich",0],UNIT["degree",0.0174532925199433]],'
                   'PROJECTION["Polar_Stereographic"],PARAMETER["latitude_of_origin",70],'
                   'PARAMETER["central_meridian",-45],PARAMETER["false_easting",0],'
                   'PARAMETER["false_northing",0],UNIT["metre",1],AUTHORITY["EPSG","3413"]]',
        'proj4_params': '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 '
                        '+datum=WGS84 +units=m +no_defs',
    }


GLOBAL_ATTRS = {
    'title': 'SIType-Arctic daily Arctic sea ice type (FYI/MYI) product corrected with NSIDC sea ice extent',
    'product_name': 'SIType-Arctic',
    'version': f'v{VERSION}',
    'Conventions': 'CF-1.11, ACDD-1.3',
    'source': 'Daily sea ice type GeoTIFF maps produced by the FY-3 MWRI and ASCAT '
              'classification workflow',
    'institution': 'School of Geospatial Engineering and Science, Sun Yat-sen University, China',
    'creator_name': 'Ziyun Qi',
    'creator_email': 'qizy3@mail2.sysu.edu.cn',
    'creator_institution': 'School of Geospatial Engineering and Science, Sun Yat-sen University, China',
    'contact': 'Ziyun Qi (qizy3@mail2.sysu.edu.cn)',
    'license': 'Creative Commons Attribution 4.0 International (CC BY 4.0), '
               'https://creativecommons.org/licenses/by/4.0/',
    'data_repository': 'Zenodo',
    'doi': '10.5281/zenodo.23489985',
    'landing_page': 'https://doi.org/10.5281/zenodo.23489985',
    'disclaimer': 'Not valid for navigation.',
    'platform': 'MetOp ASCAT scatterometer; Fengyun-3D MWRI radiometer',
    'sensor': 'ASCAT sigma-0 (V-pol); MWRI 18.7/23.8/36.5 GHz brightness temperature (V-pol)',
    'ancillary_data': 'NSIDC daily sea ice extent (resampled to the product grid); '
                      'ERA5 2 m air temperature; climatological multi-year ice extent and '
                      'minimal MYI range mask; annual mean MYI fraction (AnalysisResult_YperY)',
    'processing_level': 'L3 - daily gridded sea ice class map',
    'spatial_resolution': '12500 m',
    'grid': 'NSIDC Sea Ice Polar Stereographic North (EPSG:3413), 608 columns x 896 rows',
    'crs': 'EPSG:3413',
    'geographic_scope': 'Arctic Ocean and adjacent sub-Arctic seas',
    'geospatial_lat_min': '37.36',
    'geospatial_lat_max': '90.0',
    'geospatial_lon_min': '-180.0',
    'geospatial_lon_max': '180.0',
    'time_resolution_period': '1 day',
    'temporal_coverage': '2013-11-01 to 2024-04-30, cold season only (November to April)',
    'date_created': datetime.datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ'),
    'summary': 'Daily sea ice type classification (first-year ice / multi-year ice) over the Arctic, '
               'derived from ASCAT scatterometry and FY-3D MWRI brightness temperatures, spatially '
               'filtered, temperature-optimised, and finally fused with the NSIDC daily sea ice extent '
               'product to fill underestimated and overestimated ice classes. Pixel values are integer '
               'class codes: 0 = first-year ice, 1 = multi-year ice, 3 = open water, 4 = land, '
               '5 = no data.',
    'keywords': 'Arctic, sea ice type, first-year ice, multi-year ice, ASCAT, FY-3D MWRI, NSIDC SIE, '
                'polar stereographic',
    'references': 'Dataset DOI: https://doi.org/10.5281/zenodo.23489985 (Zenodo). '
                  'See README.md delivered with this dataset.',
    'comment': 'This file was produced by converting the operational GeoTIFF outputs; class codes and '
               'grid are bit-identical to the source rasters.',
}


def write_day(date, tif_path, x, y, lat, lon, epsg, rows, cols, ref_transform):
    out_path = os.path.join(OUT_DIR, nc_rel_path(date))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    ds = nc.Dataset(out_path, 'w', format='NETCDF4')
    ds.createDimension('time', 1)
    ds.createDimension('nv', 2)
    ds.createDimension('y', rows)
    ds.createDimension('x', cols)

    v = ds.createVariable('time', 'f8', ('time',))
    v.units = 'days since 1970-01-01 00:00:00 UTC'
    v.calendar = 'standard'
    v.axis = 'T'
    v.standard_name = 'time'
    v.long_name = 'nominal time of the daily sea ice type composite'
    v.bounds = 'time_bounds'
    vb = ds.createVariable('time_bounds', 'f8', ('time', 'nv'))
    vb.units = 'days since 1970-01-01 00:00:00 UTC'
    vb.calendar = 'standard'
    vb.long_name = 'start and end time of each daily composite'

    day = datetime.datetime.strptime(date, '%Y%m%d')
    noon = day + datetime.timedelta(hours=12)
    units = 'days since 1970-01-01 00:00:00 UTC'
    v[:] = np.array([nc.date2num(noon, units, calendar='standard')], dtype='f8')
    vb[:] = np.array([[nc.date2num(day, units, calendar='standard'),
                       nc.date2num(day + datetime.timedelta(days=1), units, calendar='standard')]],
                     dtype='f8')

    vx = ds.createVariable('x', 'f8', ('x',))
    vx.units = 'm'
    vx.standard_name = 'projection_x_coordinate'
    vx.long_name = 'x coordinate of projection (EPSG:3413)'
    vx.axis = 'X'
    vx[:] = x
    vy = ds.createVariable('y', 'f8', ('y',))
    vy.units = 'm'
    vy.standard_name = 'projection_y_coordinate'
    vy.long_name = 'y coordinate of projection (EPSG:3413)'
    vy.axis = 'Y'
    vy[:] = y

    vg = ds.createVariable('polar_stereographic', 'i2', ())
    for k, val in projection_attrs(epsg).items():
        setattr(vg, k, val)

    vlat = ds.createVariable('latitude', 'f4', ('y', 'x'), zlib=True, complevel=4, shuffle=True)
    vlat.units = 'degrees_north'
    vlat.standard_name = 'latitude'
    vlat.long_name = 'latitude of each grid cell centre'
    vlat[:] = lat
    vlon = ds.createVariable('longitude', 'f4', ('y', 'x'), zlib=True, complevel=4, shuffle=True)
    vlon.units = 'degrees_east'
    vlon.standard_name = 'longitude'
    vlon.long_name = 'longitude of each grid cell centre'
    vlon[:] = lon

    vsit = ds.createVariable('sea_ice_type', 'u1', ('time', 'y', 'x'),
                             zlib=True, complevel=4, shuffle=True, chunksizes=(1, rows, cols),
                             fill_value=FILL_VALUE)
    vsit.long_name = 'sea ice type class code'
    vsit.units = '1'
    vsit.grid_mapping = 'polar_stereographic'
    vsit.coordinates = 'latitude longitude'
    vsit.comment = 'value 5 marks no-data cells (missing sensor swath, cloud/masked areas, or the ' \
                   'central Arctic gap); it is stored as the CF _FillValue and is auto-masked on read.'
    vsit.flag_values = np.array(FLAG_VALUES, dtype='u1')
    vsit.flag_meanings = FLAG_MEANINGS
    vsit.flag_descriptions = 'first-year ice | multi-year ice | open water (no ice cover) ' \
                             '| land mask'
    vsit.valid_range = np.array([0, 5], dtype='u1')
    vsit.resolution = '12.5 km, daily'

    with rasterio.open(tif_path) as src:
        if src.transform != ref_transform or (src.width, src.height) != (cols, rows):
            raise SystemExit(f'grid mismatch: {tif_path}')
        vsit[0, :, :] = src.read(1).astype('u1')

    iso = f'{date[:4]}-{date[4:6]}-{date[6:8]}'
    version_iso = f'{VERSION[:4]}-{VERSION[4:6]}-{VERSION[6:8]}'
    for k, val in GLOBAL_ATTRS.items():
        if k == 'temporal_coverage':
            val = iso
        elif k == 'geospatial_lat_min':
            val = f'{lat.min():.4f}'
        elif k == 'geospatial_lat_max':
            val = f'{lat.max():.4f}'
        elif k == 'geospatial_lon_min':
            val = f'{lon.min():.4f}'
        elif k == 'geospatial_lon_max':
            val = f'{lon.max():.4f}'
        setattr(ds, k, val)
    ds.temporal_coverage_start = iso
    ds.temporal_coverage_end = iso
    ds.production_date = version_iso
    ds.history = f'{datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")} ' \
                 f'Export_SIT_NetCDF.py converted the source daily GeoTIFF sea ice type map ' \
                 f'into a CF-compliant NetCDF-4 file'
    ds.close()

    return out_path


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    months = collect_dates()
    if not months:
        raise SystemExit(f'no source GeoTIFF found under {SRC_DIR}')

    ref = months[sorted(months)[0]][0][1]
    x, y, epsg, rows, cols, lat, lon, ref_transform = grid_geometry(ref)
    print(f'grid: {cols}x{rows}, EPSG:{epsg}, x {x[0]:.0f}..{x[-1]:.0f}, y {y[0]:.0f}..{y[-1]:.0f}')
    print(f'lat range {lat.min():.3f}..{lat.max():.3f}, lon range {lon.min():.3f}..{lon.max():.3f}')

    inventory = []
    skipped = 0
    for yyyymm, items in months.items():
        for date, tif_path in items:
            out_path = os.path.join(OUT_DIR, nc_rel_path(date))
            if os.path.exists(out_path) and not FORCE:
                skipped += 1
                continue
            path = write_day(date, tif_path, x, y, lat, lon, epsg, rows, cols, ref_transform)
            size_kb = os.path.getsize(path) / 1024.0
            inventory.append({'date': date, 'file': nc_rel_path(date).replace('\\', '/'),
                              'size_kb': round(size_kb, 1)})
        print(f'{yyyymm}: {len(items)} days written')

    total_mb = sum(it['size_kb'] for it in inventory) / 1024.0
    print('\nINVENTORY_SUMMARY')
    print(f'daily files written: {len(inventory)}, skipped (exists): {skipped}, '
          f'new size: {total_mb:.1f} MB')


if __name__ == '__main__':
    main()
