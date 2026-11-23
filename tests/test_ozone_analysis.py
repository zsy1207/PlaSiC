"""Test area/time weighting and the scientific contrast of experiment commands."""
import importlib.util
from datetime import datetime
from pathlib import Path
import tempfile
import unittest
import netCDF4 as nc
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


analysis = load('analysis', 'experiments/ozone_historical/fig1/compute_temperature.py')
runner = load('runner', 'experiments/ozone_historical/run_experiments.py')


class HistoricalAnalysis(unittest.TestCase):
    def test_gaussian_and_calendar_weighting(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'tas.nc'
            nodes, _ = np.polynomial.legendre.leggauss(32)
            bounds_dates = [datetime(1850, m, 1) for m in range(1, 13)] + [datetime(1851, 1, 1)]
            days = np.array([(b-a).days for a, b in zip(bounds_dates[:-1], bounds_dates[1:])])
            with nc.Dataset(path, 'w') as ds:
                for name, size in [('time', 12), ('lat', 32), ('lon', 64), ('bnds', 2)]:
                    ds.createDimension(name, size)
                ds.createVariable('lat', 'f8', ('lat',))[:] = np.rad2deg(np.arcsin(nodes[::-1]))
                t = ds.createVariable('time', 'f8', ('time',))
                t.units, t.calendar, t.bounds = 'days since 1850-01-01', 'proleptic_gregorian', 'time_bnds'
                edges = nc.date2num(bounds_dates, t.units, t.calendar)
                t[:] = (edges[:-1]+edges[1:])/2
                ds.createVariable('time_bnds', 'f8', ('time', 'bnds'))[:] = np.stack([edges[:-1],edges[1:]],axis=1)
                v = ds.createVariable('tas', 'f8', ('time', 'lat', 'lon'))
                v.units = 'K'
                v[:] = 280 + np.arange(12)[:,None,None] + nodes[::-1][None,:,None]**2
            years, values = analysis.annual_temperature(path, 1850, 1850)
            np.testing.assert_array_equal(years, [1850])
            np.testing.assert_allclose(values, 280 + np.average(np.arange(12),weights=days) + 1/3, rtol=1e-14)
            with nc.Dataset(path, 'a') as ds:
                ds['time_bnds'][-1, 1] -= 1
            with self.assertRaisesRegex(ValueError, 'complete calendar month'):
                analysis.annual_temperature(path, 1850, 1850)

    def test_forcing_contrasts_and_leap_days(self):
        for case in runner.CASES:
            cmd = runner.command(Path('/model'), Path('/ic'), Path('/o3'), Path('/output'), case)
            self.assertEqual(cmd[cmd.index('--steps')+1], '1928480')
            self.assertEqual('--co2-ppm' in cmd, case == 'ozone_only')
            self.assertEqual('--co2-forcing' in cmd, case != 'ozone_only')
            self.assertIn('--ozone-forcing', cmd)


if __name__ == '__main__':
    unittest.main()
