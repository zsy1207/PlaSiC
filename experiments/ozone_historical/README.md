# Historical CO₂ and ozone forcing, 1850–2014

Three experiments compare annual global near-surface air temperature under changing CO₂ and ozone. They follow the [maintainer's five requests](https://github.com/sunmoumou1/PlaSiC/pull/9#issuecomment-7079546681) and use the grid, ocean and warm initial condition from the [reference historical experiment](../historical_3d_500y/).

| Experiment | CO₂ | Ozone |
|---|---|---|
| CO₂ + ozone | Historical UoM monthly concentrations | Historical UReading monthly fields |
| Ozone only | Fixed at 284.316853 ppmv | Historical UReading monthly fields |
| CO₂ only | Historical UoM monthly concentrations | Repeating external 1850 seasonal cycle |

All three start from `icA.restart`, spun up with idealized ozone. The common switch to external ozone can cause adjustment. Each configuration has one member. Own-series baselines remove a mean offset, but do not identify forced trends separately from initial adjustment, model drift or internal variability. The three trajectories describe the model under these configurations.

## Forcing data

The four original `vmro3` files are from CMIP6 input4MIPs, UReading-CCMI-1-0, native grid `gn`, version `v20160711`. They cover 1850–1899, 1900–1949, 1950–1999 and 2000–2014. Together they contain 1,980 monthly fields on 66 pressure levels and a 96 × 144 horizontal grid, totalling 5,927,876,651 bytes. SHA256 values in [SHA256SUMS](SHA256SUMS) match the ORNL and Argonne catalogue replicas. [Dataset identifier](https://hdl.handle.net/21.14100/ab27892c-eaa5-3430-a3c5-f9b0549a556f) · [ESGF search requested by the maintainer](https://github.com/sunmoumou1/PlaSiC/pull/9#issuecomment-7089574072)

Run the following commands from the repository root. The defaults keep downloads and model output in ignored experiment subdirectories. `PLASIC_OZONE_RUNS` can point to another disk. The model requires the repository's C, SHTns, NetCDF and MPI dependencies. The scripts use Python 3.10 or newer, NumPy, SciPy, netCDF4 and Matplotlib.

```bash
exp="$PWD/experiments/ozone_historical"
data="$exp/data"
export PLASIC_OZONE_RUNS="${PLASIC_OZONE_RUNS:-$exp/runs}"
mkdir -p "$data" "$PLASIC_OZONE_RUNS"
base=https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/input4MIPs/CMIP6/CMIP/UReading/UReading-CCMI-1-0/atmos/mon/vmro3/gn/v20160711
(
  cd "$data" || exit 1
  while read -r sha filename; do
    curl --fail --location --retry 3 --continue-at - --output "$filename" "$base/$filename" || exit 1
  done < "$exp/SHA256SUMS"
  sha256sum --check "$exp/SHA256SUMS"
)
```

## Prepare and run

The study uses T21, 10 stretched atmospheric sigma levels, the 16-level 3-D ocean and two MPI ranks. All runs use the Gregorian calendar and 32 steps per day. First obtain the twelve monthly surface-pressure fields used for remapping. This one-year reference uses idealized ozone and supplies the pressure climatology.

```bash
make -C src -j4 all MPI=1 NPRO=2 NLAT=32 NLEV=10
mkdir -p "$data/reference/monthly"
mpirun --bind-to none -np 2 build/T21_L10_MPI_P2_O16/plasic.x \
  --restart experiments/historical_3d_500y/initial_conditions/icA.restart \
  --reset-nstep 21622368 --steps 11680 --co2-ppm 284.316853 \
  --monthly-output "$data/reference/monthly" --output "$data/reference/final.restart"

python "$exp/prepare_ozone.py" "$data"/vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_*.nc \
  --reference-ps "$data/reference/monthly/ps.nc" --output "$data/ozone_historical.nc"
python "$exp/prepare_ozone.py" "$data"/vmro3_input4MIPs_ozone_CMIP_UReading-CCMI-1-0_gn_*.nc \
  --reference-ps "$data/reference/monthly/ps.nc" --fixed-1850 --output "$data/ozone_1850_cycle.nc"

python "$exp/run_experiments.py" --historical-ozone "$data/ozone_historical.nc" \
  --fixed-ozone "$data/ozone_1850_cycle.nc" --jobs 2
python "$exp/fig1/compute_temperature.py"
python "$exp/fig1/plot_temperature.py" --final
```

The launcher can be replaced in full, for example with `--mpirun 'srun -n 2'` under an appropriate Slurm allocation. Each case takes 1,928,480 steps. To exercise a short run, set `--end-year 1850`. To continue it, rerun with the same input files and output directory. Each completed year is atomically published at `<case>/<year>/` with its restart, monthly NetCDF files and provenance manifest. The current format checks the SHA256 of the restart and every monthly output on resume. Older manifests verify the restart and temperature file. An incomplete or changed checkpoint stops the run.

## Remapping and interpretation

- Source values are ozone mole fractions. The model converts mole to mass fraction using 47.9982/28.9644 for ozone/dry air. Latitude/longitude interpolation is bilinear and periodic in longitude. Vertically each pressure cell contributes in proportion to its overlap with the target layer. Use the deepest source value between the bottom source edge and local ground where necessary.
- Model interfaces use the same monthly reference surface pressure in every year. This conserves the vertical column on that reference pressure grid. The radiation calculation uses the model's instantaneous pressure. The converter reads the Gaussian horizontal grid from the pressure file. Set `--levels` to match atmospheric builds other than L10.
- Both external files use identical time interpolation. January 1850 holds its first value before mid-month, and December 2014 holds its last value after mid-month. In other months the model interpolates between adjacent monthly centres. The interpolated field generally has a different calendar-month average from its input anchor. The fixed file repeats the twelve 1850 monthly fields throughout.
- In 60 sampled months, changing from bilinear to area-conservative horizontal remapping changed the global column by at most 0.018 DU and individual cells by up to 3.70 DU. Substituting the 1850–1859 monthly pressure climatology from the CO₂-only experiment changed global columns by up to 0.080 DU and individual cells by up to 6.30 DU. These are offline sensitivities, with both methods evaluated at a common pressure in each comparison.
- The L10 top layer contains 74.6–81.0% of the global ozone column over the 1,980 months. In the sampled months, using layer-midpoint pressure in the longwave calculation overestimates the top-layer pressure-weighted ozone path by 43–48% and the full-column path by 11–15%. These percentages describe the optical-path approximation. They do not quantify the resulting flux or temperature errors.
- The 1851–1860 mean top-of-atmosphere net downward radiation and full-depth ocean heat storage have opposite signs in all three runs. Their differences are approximately +3.99, +4.72 and +4.03 W m⁻² for CO₂ + ozone, ozone only and CO₂ only. A matched one-year reference with idealized ozone also showed this sign difference. This partial budget excludes other reservoirs.

## Temperature analysis and validation

Use the model’s near-surface air temperature (`tas`), Gaussian quadrature in latitude, equal longitude weights and actual calendar-month lengths. The CSV contains each annual mean in kelvin and its anomaly relative to that case's inclusive 1850–1900 mean. The figure shows annual values and a separate centred nine-year mean for each case. Means span 1854–2010 with no endpoint padding. Model diagnostics, the CSV and PNG/PDF are written under `$PLASIC_OZONE_RUNS`. Final results require 165 complete years from each case.

```bash
make -C src test-data test-ozone test-ozone-mpi
```

For each model invocation with prescribed ozone, the first radiation step reports the area-weighted global ozone column in DU, calculated from the actual layer mass fractions and surface pressure. A quantitative test prescribes a field with a known 300 DU global mean under PlaSiC's STP density convention. In a matched two-rank replay, the first radiation step of 1850 reported 287.631 DU. An independent calculation from the same instantaneous pressure gave 287.630594 DU. The ozone-specific tests cover remapping and native-grid units, calendar edges, startup validation, two-rank error propagation and radiation responses. The general upstream `test-runtime` target in `f29ee8e` references a missing test script. GCC 14 also issues a pre-existing `-Wmaybe-uninitialized` warning in the ocean halo code under strict builds. The default model and ozone tests build without promoting that upstream warning to an error.
