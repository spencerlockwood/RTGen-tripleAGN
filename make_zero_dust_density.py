import numpy as np
 
# A numerically negligible floor rather than exact 0.0 -- consistent with
# how this pipeline already handles "should be zero" quantities elsewhere
# (see numsmall = 1e-8 in constants.py, used to avoid divide-by-zero /
# undefined-log issues). RADMC-3D's radiative transfer math (optical depth,
# scattering source function normalization) can behave badly with exact
# zero opacity in some code paths -- a tiny floor avoids that class of bug
# entirely while still being physically indistinguishable from "no dust."
#
# 1e-60 rather than 1e-30: with 2 species, 1e-30 each is 2e-30 g/cm^3, which is
# NOT negligible -- the median real dust density is ~1e-28 at z=0 and ~1e-29 at
# z~3 (where ~20% of cells hold less dust than that). Since this run reuses the
# dusty dust_temperature.dat, a 1e-30 floor glows in the far-IR and fakes an
# 'intrinsic' FIR bump. 1e-60 is still non-zero but truly negligible.
FLOOR_VALUE = 1e-60
 
SRC = 'dust_density.inp'          # your existing, real dust density file
DST = 'dust_density_zerodust.inp' # new file for the dust-free comparison run
 
with open(SRC, 'r') as f:
    iformat_line = f.readline()
    ncells_line  = f.readline()
    nspec_line   = f.readline()
 
ncells = int(ncells_line)
nspec  = int(nspec_line)
 
print(f'Read header from {SRC}: ncells={ncells}, nspec={nspec}')
print(f'Writing zero-dust (floor={FLOOR_VALUE:.1e}) density file to {DST} ...')
 
zeros = np.full(ncells, FLOOR_VALUE)
 
with open(DST, 'w') as f:
    f.write(iformat_line)
    f.write(ncells_line)
    f.write(nspec_line)
    for ispec in range(nspec):
        zeros.tofile(f, sep='\n', format="%13.6e")
        f.write('\n')
 
print('Done.')