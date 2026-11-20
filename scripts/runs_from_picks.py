"""Print one scripts/make_run.sh command per row of a picks CSV written by
selection/find_dust_analogues.py. Review the commands, then run them.

    python scripts/runs_from_picks.py dust_analogues.csv            # numbered 1..N
    python scripts/runs_from_picks.py dust_analogues.csv --start 6  # numbered 6..
    python scripts/runs_from_picks.py dust_analogues.csv | bash     # create them all
"""
import argparse
import os
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument('csv')
ap.add_argument('--start', type=int, default=1, help='first system number')
a = ap.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
d = pd.read_csv(a.csv)
for k, (_, r) in enumerate(d.iterrows()):
    ids = [int(round(r[c])) for c in ('BH_id', 'Neighbour1BH_id', 'Neighbour2BH_id')]
    pos = ', '.join(f'{r[c]:.1f}' for c in ('BH_x', 'BH_y', 'BH_z'))
    lbol = ' '.join(f'{r[c]:.4e}' for c in ('BH_Lbol', 'Neighbour1BH_Lbol', 'Neighbour2BH_Lbol'))
    print(f'{here}/make_run.sh {a.start + k} {int(r.snapshot):06d} {r.z:.3f} {int(r.pynbody_haloid)} '
          f'{ids[0]} {ids[1]} {ids[2]} "{pos}" {lbol}')
