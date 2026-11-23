"""Numerical and metadata tests for the historical ozone preparation workflow."""
import importlib.util
import os
import subprocess
from pathlib import Path
import tempfile
import unittest

import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('prepare', ROOT / 'experiments/ozone_historical/prepare_ozone.py')
prepare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prepare)


class Preprocessing(unittest.TestCase):
    def test_log_pressure_and_latitude(self):
        p = np.array([100000., 10000., 1000.])  # deliberately descending
        lat = np.array([90., 0., -90.])
        lon = np.array([180., 0., 90., 270.])
        v = np.broadcast_to((1e-6 + 1e-7*np.log(p)[:, None, None] +
                              1e-9*lat[None, :, None]), (3, 3, 4)).copy()
        ps = np.array([[80000., 100000.], [90000., 95000.]])
        sigma = np.array([.1, .8])
        target_lat = np.array([30., -60.])
        actual = prepare.regrid(v, p, lat, lon, sigma, target_lat, [359., 1.], ps)
        expected = 1e-6 + 1e-7*np.log(sigma[:, None, None]*ps) + 1e-9*target_lat[None, :, None]
        np.testing.assert_allclose(actual, expected, rtol=1e-14)

    def test_periodic_longitude_and_pressure_edges(self):
        v = np.tile(np.array([0., 1., 0., 1.])*1e-6, (2, 2, 1))
        result = prepare.regrid(v, [1000., 100000.], [-80., 80.], [0., 90., 180., 270.],
                                [.001, 1.], [-89., 89.], [359., 1.], np.full((2, 2), 101325.))
        np.testing.assert_allclose(result, 1e-6/90, rtol=1e-13)
        with self.assertRaisesRegex(ValueError, 'unique'):
            prepare.regrid(v, [1000., 100000.], [-80., 80.], [0., 90., 180., 360.],
                           [.5], [0.], [0.], np.array([[1e5]]))

    def test_published_month_units_and_missing(self):
        with tempfile.TemporaryDirectory() as tmp:
            with nc.Dataset(Path(tmp)/'source.nc', 'w') as ds:
                ds.createDimension('time', 3)
                time = ds.createVariable('time', 'f8', ('time',))
                time.units = 'months since 1850-01-01 00:00'
                time.calendar = 'standard'
                time[:] = [0.5, 1.5, 12.5]
                np.testing.assert_array_equal(prepare.month_indices(time), [22200, 22201, 22212])
                time[0] = np.ma.masked
                with self.assertRaisesRegex(ValueError, 'missing'):
                    prepare.month_indices(time)

    def test_expected_files_required(self):
        with self.assertRaisesRegex(ValueError, 'four requested'):
            prepare.prepare([], Path('missing'), Path('missing'))


class FullConversion(unittest.TestCase):
    def test_four_file_conversion_and_fixed_seasonal_cycle(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            sigma, lat, lon = prepare.sigma_grid()
            pressure_file = root/'ps.nc'
            with nc.Dataset(pressure_file, 'w') as ds:
                for name, size in [('time', 12), ('lat', 32), ('lon', 64)]:
                    ds.createDimension(name, size)
                ds.createVariable('lat', 'f8', ('lat',))[:] = lat
                ds.createVariable('lon', 'f8', ('lon',))[:] = lon
                t = ds.createVariable('time', 'f8', ('time',))
                t.units, t.calendar = 'months since 1850-01-01 00:00', 'standard'
                t[:] = np.arange(12) + .5
                ps = ds.createVariable('ps', 'f8', ('time', 'lat', 'lon'))
                ps.units = 'Pa'
                ps[:] = 100000.
            files = []
            start = 0
            for period, count in zip(prepare.PERIODS, [600, 600, 600, 180]):
                path = root/(prepare.PREFIX + period + '.nc')
                files.append(path)
                with nc.Dataset(path, 'w') as ds:
                    ds.source_id = 'UReading-CCMI-1-0'
                    for name, size in [('time', count), ('plev', 2), ('lat', 2), ('lon', 2)]:
                        ds.createDimension(name, size)
                    for name, values in [('plev', [1., 1000.]), ('lat', [-90., 90.]), ('lon', [0., 180.])]:
                        ds.createVariable(name, 'f8', (name,))[:] = values
                    ds['plev'].units = 'hPa'
                    t = ds.createVariable('time', 'f8', ('time',))
                    t.units, t.calendar = 'months since 1850-01-01 00:00', 'standard'
                    t[:] = start + np.arange(count) + .5
                    o3 = ds.createVariable('vmro3', 'f4', ('time', 'plev', 'lat', 'lon'))
                    o3.units = 'mole mole-1'
                    o3[:] = (1e-6 + (start + np.arange(count)) * 1e-10)[:, None, None, None]
                start += count
            for fixed in (False, True):
                output = root/f'converted-{fixed}.nc'
                prepare.prepare(files, pressure_file, output, fixed)
                with nc.Dataset(output) as ds:
                    self.assertEqual(ds['o3'].shape, (1980, 10, 32, 64))
                    np.testing.assert_allclose(ds['o3'][1979], 1e-6 + (11 if fixed else 1979)*1e-10, rtol=1e-6)
                    self.assertEqual(int(ds['year'][-1]), 2014)
                    self.assertEqual(int(ds['month'][-1]), 12)
                    self.assertIn('sha256', ds.source_files)
                    self.assertEqual(ds.reference_ps_sha256, prepare.sha256(pressure_file))
                model = os.environ.get('PLASIC_TEST_MODEL')
                if model:
                    subprocess.run([str(Path(model).resolve()), '--co2-ppm', '284.316853',
                                    '--steps', '1', '--start-year', '1850',
                                    '--data-dir', str(ROOT/'src/data/surface'),
                                    '--ozone-forcing', str(output), '--output', str(root/'restart')],
                                   check=True, capture_output=True, timeout=120)
                    self.assertTrue((root/'restart').is_file())

if __name__ == '__main__':
    unittest.main()
