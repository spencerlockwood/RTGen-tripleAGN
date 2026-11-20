"""Check, before submitting, that halo_number really contains the host BH.

    python scripts/check_halo.py SNAP HALO "X, Y, Z"
    python scripts/check_halo.py 003072 118727 "-3864.7, -1807.4, 4519.8"

Reads row HALO of the snapshot's AHF_halos file (line 0 is the header, so the
1-based row = tangos halo_number = pynbody halo_numbers='v1') and prints its
Mvir, Rvir and the distance of the host BH from its centre. Uses SIM_DIR,
SIM_PREFIX, SIM_H and SIM_BOX_KPC from scripts/config.sh (source it first, or
run through: bash -c 'source scripts/config.sh && python scripts/check_halo.py ...').
"""
import glob
import os
import sys

snap, row, pos = sys.argv[1], int(sys.argv[2]), [float(v) for v in sys.argv[3].split(',')]
h = float(os.environ.get('SIM_H', 0.6777))
L = float(os.environ.get('SIM_BOX_KPC', 25000.))
sim_dir, prefix = os.environ['SIM_DIR'], os.environ['SIM_PREFIX']

f = glob.glob(os.path.join(sim_dir, f'{prefix}.{snap}.z*.AHF_halos'))[0]
z = float(f.split('.z')[-1].split('.AHF')[0])
a = 1 / (1 + z)
Lp = L * a                                             # box size, proper kpc
v = open(f).readlines()[row].split()
c = [float(v[k]) / h * a - L / 2 * a for k in (5, 6, 7)]   # AHF kpc/h comoving -> proper kpc, box centred on 0
d = [((ci - pi + Lp / 2) % Lp) - Lp / 2 for ci, pi in zip(c, pos)]
r = sum(x * x for x in d) ** .5
rvir = float(v[11]) / h * a
print(f'{os.path.basename(f)} row {row}: Mvir = {float(v[3]) / h:.2e} Msun, Rvir = {rvir:.0f} kpc, '
      f'host BH offset = {r:.1f} kpc  ->  {"OK" if r < 0.5 * rvir else "CHECK: far from halo centre"}')
