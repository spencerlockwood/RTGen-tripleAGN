import numpy as np
import pynbody as pn
import pynbody.plot.sph as sph
from astropy.cosmology import FlatLambdaCDM
import input_params as pars
import constants as const

# INPUT
input_galaxy_dir = pars.input_galaxy_dir
filename = pars.input_galaxy_dir + pars.input_galaxy_filename

# ****************************************         
# ****************************************                                                             
# **************************************** 

# Define some important variables
zstart = pars.zstart

cosmo = FlatLambdaCDM(H0=pars.H0, Om0=pars.Om0, Tcmb0=pars.Tcmb)

tstart = cosmo.lookback_time(zstart).to_value()
age_universe = cosmo.age(0).to_value()
tstart = age_universe - tstart

zs = const.zs

# Read simulation snapshot
#
# IMPORTANT: We load lazily and immediately filter down to a single halo via
# snap.halos(). This keeps memory usage bounded to that halo's particles
# instead of materializing the whole 25 Mpc box (which OOMs, since Romulus25
# contains far more dark matter + gas + star particles than fit in 64GB).
snap = pn.load(filename, lazy=True)
snap.physical_units()

# Romulus25 snapshots carry an on-disk 'grp' array, so pynbody 2 would pick
# HaloNumberCatalogue by default. Force the AHF catalogue instead, and number
# halos the way tangos' halo_number does: 1-based file order, which is
# pynbody's halo_numbers='v1'. (The default 'ahf' numbering uses the MPI-AHF
# hash IDs, and 'file-order' is 0-based, i.e. off by one, giving a different
# halo hundreds of kpc away.)
halos = snap.halos(priority=['AHFCatalogue'], halo_numbers=pars.halo_numbering)
gal = halos[pars.halo_number]

# Compute time and redshift
#
# Romulus25 is a cosmological run, so the age of the universe at this snapshot
# and its redshift are read directly from the snapshot. tsnap replaces the
# hard-coded value in input_params.py (13.8 Gyr was only right for z=0) and is
# what stellar ages are measured from below.
ctx = snap.conversion_context()
tsnap = float(snap.properties['time'].in_units('Gyr', **ctx))
z_real = float(snap.properties['z'])
time_real = tsnap

print('Snapshot age of universe (Gyr): ', tsnap)
print('Snapshot redshift:              ', z_real)

# Write time and redshift to file

ff = open('time_and_z.txt', 'w')

ff.write('# This file contains information about time and redshift of the snapshot (read from the snapshot itself) \n')
ff.write('# sim time [Gyr]      real time [Gyr]      z\n')
ff.write(str(tsnap) + '      ' + str(time_real) + '     ' + str(z_real) + '\n')
ff.close()

# ****************************************
# Recenter the halo on its own coordinate origin
# ****************************************
#
# gal.gas['pos'] / gal.star['pos'] are still expressed in the *box's*
# coordinate frame at this point (halos() only filters particle membership,
# it does not recenter). Since Romulus25 is a full cosmological box, the
# halo can sit anywhere within it, potentially tens of thousands of kpc from
# (0,0,0) -- and sph_interpolation.py downstream assumes the incoming
# positions are already galaxy-centered (it only applies a fixed +lbox/2
# shift and a +/-lbox/2 threshold cut). Without recentering here, that
# threshold cut discards essentially all particles.
#
# IMPORTANT: this centering is done entirely with plain numpy on the small,
# already halo-filtered arrays (gal.gas['pos'], gal.star['pos']). Earlier
# attempts used pynbody.analysis.halo.center(gal, ...) and gal.ancestor.wrap(),
# both of which operate on (or fall back to) the *full snapshot's* particle
# arrays -- that is what caused the repeated OOM kills. boxsize is read from
# a snapshot-level property, which is cheap (a single scalar, no particle
# array is touched to get it).
boxsize = float(snap.properties['boxsize'].in_units('kpc', **ctx))
print('Box size (proper kpc): ', boxsize)

pos_gas_raw   = np.array(gal.gas['pos'])
pos_star_raw  = np.array(gal.star['pos'])
star_iord     = np.array(gal.star['iord'])
star_tform    = np.array(gal.star['tform'])

def periodic_delta(x, center, box):
    """Minimum-image displacement of x from center under box periodicity."""
    d = x - center
    d -= box * np.round(d / box)
    return d

# Centre the grid on the HOST black hole of the triple rather than on the
# halo's stellar centre. For the triple-AGN systems the BHs need not sit at
# the halo centre (several hosts are ~5e12 Msun groups), and the RT grid is
# only sizex = 30 kpc across, so the BH position is the reliable anchor.
# BHs are star-family particles in ChaNGa; they are found by iord.
bh_iords = np.asarray(pars.bh_iords, dtype=star_iord.dtype)
missing = bh_iords[~np.isin(bh_iords, star_iord)]
if len(missing) > 0:
    raise RuntimeError('BH iord(s) %s not found in halo %d -- wrong halo_number / halo_numbering?'
                       % (missing, pars.halo_number))
bh_idx = np.array([np.where(star_iord == i)[0][0] for i in bh_iords])
cen = pos_star_raw[bh_idx[0]].copy()

print('Host BH position, box frame (proper kpc): ', cen)
print('Catalogue host BH position (proper kpc):  ', np.array(pars.host_pos_kpc))
bh_offsets = periodic_delta(pos_star_raw[bh_idx], cen, boxsize)
for lab, iord, off in zip(['Host', 'N1  ', 'N2  '], bh_iords, bh_offsets):
    print('  %s BH %d: offset from grid centre (kpc) = %s  |r| = %.3f'
          % (lab, iord, np.round(off, 3), np.linalg.norm(off)))

# BH positions in the grid frame (simulation axes, proper kpc, host at 0),
# read by postproc/band_luminosity.py to place an aperture on each member.
bh_mass = np.array(gal.star['mass'].in_units('Msol'))[bh_idx]
with open('bh_info.txt', 'w') as fbh:
    fbh.write('# role  iord  dx_kpc  dy_kpc  dz_kpc  mass_Msun   (sim axes, origin = host BH)\n')
    for lab, iord, off, mbh in zip(['Host', 'N1', 'N2'], bh_iords, bh_offsets, bh_mass):
        fbh.write('%s %d %.5f %.5f %.5f %.5e\n' % (lab, iord, off[0], off[1], off[2], mbh))

# Unwrap star/gas positions relative to the host BH using the minimum-image
# convention, so a halo straddling the box boundary stays contiguous.
pos_star_unwrapped = periodic_delta(pos_star_raw, cen, boxsize)
pos_gas_unwrapped  = periodic_delta(pos_gas_raw, cen, boxsize)

# Black holes are star particles with tform < 0. Exclude them from the stellar
# arrays: otherwise each BH (up to ~5e8 Msun here) is treated as a single old
# stellar population and injects spurious starlight at the nucleus.
is_star = star_tform >= 0.
print('Excluding %d BH particles (tform<0) from the stellar arrays' % np.sum(~is_star))

# Read gas properties

met = gal.gas['metals'] / zs
met[np.where(met<0)] = 0

pos = pos_gas_unwrapped   # centred on the host BH, physical kpc
temp = gal.gas['temp']
vel = gal.gas['vel']
dens = gal.gas['rho']#.in_units('g cm**-3')
#mass = 2.3262e+05*gal.gas['mass']#.in_units('g')
mass = gal.gas['mass'].in_units('Msol')

# Read star properties
posstar = pos_star_unwrapped[is_star]   # centred on the host BH, physical kpc
massstar = np.array(gal.star['mass'].in_units('Msol'))[is_star]
metstar = np.array(gal.star['metals'])[is_star] / zs
#print('Diagnostics stellar metallicity: ', np.amin(metstar), np.amax(metstar), np.mean(metstar))
metstar[np.where(metstar<0.)] = 0.

print('Total M* (1e9 M_sum): ', np.sum(massstar)/1e9)
print('Total Mgas (1e9 M_sum): ', np.sum(mass)/1e9)
print('')

# Stellar age [Gyr] = age of the universe at the snapshot - formation time
tformstar = tsnap - star_tform[is_star]

xx = pos[:,0]
yy = pos[:,1]
zz = pos[:,2]

vx = vel[:,0]
vy = vel[:,1]
vz = vel[:,2]

# Sanity check: confirm most of the halo's mass now falls well within the
# +/-sizex/2 (etc.) grid used by sph_interpolation.py. If these numbers are
# not comparable to ~sizex/2 in kpc, halo_number points at the wrong object.
print('Gas |x|,|y|,|z| max (kpc):   ', np.amax(np.abs(xx)), np.amax(np.abs(yy)), np.amax(np.abs(zz)))

np.save('posx_gas.npy', xx)
np.save('posy_gas.npy', yy)
np.save('posz_gas.npy', zz)
np.save('vx_gas.npy', vx)
np.save('vy_gas.npy', vy)
np.save('vz_gas.npy', vz)
np.save('mass_gas.npy', mass)
np.save('temperature_gas.npy', temp)
np.save('metallicity.npy', met)

xx = posstar[:,0]
yy = posstar[:,1]
zz = posstar[:,2]

print('Star |x|,|y|,|z| max (kpc):  ', np.amax(np.abs(xx)), np.amax(np.abs(yy)), np.amax(np.abs(zz)))

np.save('posx_star.npy', xx)
np.save('posy_star.npy', yy)
np.save('posz_star.npy', zz)
np.save('mass_star.npy', massstar)
np.save('metallicity_star.npy', metstar)
np.save('tform_star.npy', tformstar)