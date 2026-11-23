"""Plot verified annual temperatures; --final exports PNG and editable PDF."""
import argparse
import csv

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

from compute_temperature import CASES, OUTPUT


def main(final=False):
    available = {font.name for font in font_manager.fontManager.ttflist}
    family = next(name for name in ('Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans')
                  if name in available)
    with (OUTPUT / 'annual_temperature.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    years = np.array([int(r['year']) for r in rows])
    np.testing.assert_array_equal(years, np.arange(1850, 2015))
    plt.rcParams.update({'font.family': family, 'font.size': 8, 'axes.linewidth': .6,
                         'xtick.major.width': .6, 'ytick.major.width': .6,
                         'legend.frameon': False, 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(170 / 25.4, 85 / 25.4), dpi=600, layout='constrained')
    for case, label, color in zip(CASES, ('CO₂ + ozone', 'Ozone only', 'CO₂ only'),
                                  ('#222222', '#0072B2', '#D55E00')):
        values = np.array([float(r[case + '_anomaly_K']) for r in rows])
        if not np.isfinite(values).all():
            raise ValueError('nonfinite temperature anomaly')
        ax.plot(years, values, label=label, color=color, linewidth=.9)
    ax.axhline(0, color='.6', linewidth=.5, zorder=0)
    ax.set(xlim=(1850, 2014), xlabel='Year', ylabel='Surface air temperature anomaly (K)')
    ax.legend(loc='best')
    for extension in (('png', 'pdf') if final else ('png',)):
        fig.savefig(OUTPUT / f'fig1_surface_temperature.{extension}', dpi=600)
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--final', action='store_true')
    main(parser.parse_args().final)
