"""Run the three 1850–2014 forcing experiments from one common spin-up state."""
from __future__ import annotations
import argparse
import hashlib
import netCDF4 as nc
import numpy as np
from datetime import date
import json
import os
from pathlib import Path
import shlex
import subprocess
import time

CASES = ('co2_ozone', 'ozone_only', 'co2_only')
ROOT = Path(__file__).resolve().parents[2]
CO2 = ROOT / 'src/data/forcing/CO2/mole-fraction-of-carbon-dioxide-in-air_input4MIPs_GHGConcentrations_CMIP_UoM-CMIP-1-2-0_gn-15x360deg_000001-201412.csv'


def command(model, restart, ozone, directory, case, end_year=2014, mpirun='mpirun'):
    steps = (date(end_year + 1, 1, 1) - date(1850, 1, 1)).days * 32
    result = [*shlex.split(mpirun), '--bind-to', 'none', '-np', '2', str(model),
              '--restart', str(restart), '--reset-nstep', '21622368',
              '--steps', str(steps), '--ozone-forcing', str(ozone),
              '--monthly-output', str(directory / 'monthly'),
              '--progress', str(directory / 'progress.jsonl'), '--progress-interval', '3200',
              '--output', str(directory / 'final.restart')]
    return result + (['--co2-ppm', '284.316853'] if case == 'ozone_only'
                     else ['--co2-forcing', str(CO2)])


def run(args):
    for p in (args.model, args.restart, args.historical_ozone, args.fixed_ozone, CO2):
        if not p.is_file():
            raise FileNotFoundError(p)
    for path, expected in ((args.historical_ozone, 'historical 1850-2014'),
                           (args.fixed_ozone, 'repeated 1850 seasonal cycle')):
        with nc.Dataset(path) as ds:
            if getattr(ds, 'ozone_experiment', None) != expected:
                raise ValueError(f'{path}: expected {expected}')
            np.testing.assert_array_equal(ds['year'][:], np.repeat(np.arange(1850, 2015), 12))
            np.testing.assert_array_equal(ds['month'][:], np.tile(np.arange(1, 13), 165))
    inputs = {}
    for path in (args.model, args.restart, args.historical_ozone, args.fixed_ozone, CO2):
        with path.open('rb') as stream:
            inputs[str(path)] = hashlib.file_digest(stream, 'sha256').hexdigest()
    for case in CASES:
        if (args.output / case).exists():
            raise FileExistsError(args.output / case)
    # A run is never skipped on the mere presence of a restart or partial output.
    for case in CASES:
        directory = args.output / case
        directory.mkdir(parents=True, exist_ok=False)
        (directory / 'monthly').mkdir()
        ozone = args.fixed_ozone if case == 'co2_only' else args.historical_ozone
        cmd = command(args.model, args.restart, ozone, directory, case, args.end_year, args.mpirun)
        manifest = dict(command=cmd, cwd=str(ROOT / 'src'), end_year=args.end_year,
                        state='running', start_unix=time.time(), input_sha256=inputs)
        (directory / 'run.json').write_text(json.dumps(manifest, indent=2))
        with (directory / 'run.log').open('w') as log:
            result = subprocess.run(cmd, cwd=ROOT / 'src', stdout=log, stderr=subprocess.STDOUT)
        manifest.update(returncode=result.returncode, finish_unix=time.time(),
                        state='exited' if result.returncode == 0 else 'failed')
        (directory / 'run.json').write_text(json.dumps(manifest, indent=2))
        result.check_returncode()
        # The compute step independently checks full monthly coverage before plotting.
        if not (directory / 'final.restart').is_file():
            raise RuntimeError(f'{case}: no final restart')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', type=Path, default=ROOT / 'build/T21_L10_MPI_P2_O16/plasic.x')
    parser.add_argument('--restart', type=Path, default=ROOT / 'experiments/historical_3d_500y/initial_conditions/icA.restart')
    parser.add_argument('--historical-ozone', required=True, type=Path)
    parser.add_argument('--fixed-ozone', required=True, type=Path)
    parser.add_argument('--output', type=Path, default=Path(os.environ.get('PLASIC_OZONE_RUNS', '/workspace/plasic-runs/historical')))
    parser.add_argument('--end-year', type=int, default=2014)
    parser.add_argument('--mpirun', default='mpirun')
    args = parser.parse_args()
    if not 1850 <= args.end_year <= 2014:
        parser.error('--end-year must be between 1850 and 2014')
    for key in ('model', 'restart', 'historical_ozone', 'fixed_ozone', 'output'):
        setattr(args, key, getattr(args, key).resolve())
    run(args)
