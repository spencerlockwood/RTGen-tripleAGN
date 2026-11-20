# Band luminosities of the AGN and their host galaxies

This page shows how to get the luminosity of each triple-AGN member in a
specific band, for both the host galaxy and the AGN, from a finished run. Use
`postproc/band_luminosity.py`. The sections after the quick start explain what
it computes, so you can check it or do it by hand.

> **Read this first: the AGN are not in the radiative transfer.** RADMC-3D is
> given only stellar sources (Starburst99 spectra of the star particles). The
> BHs only set where the grid is centred. So the images contain
> **starlight, plus dust emission heated by starlight, and nothing else.**
>
> * The **host-galaxy** band luminosity comes straight from the images (§2).
> * The **AGN** band luminosity has to be estimated separately (§3), from:
>   * the AGN's bolometric luminosity in the catalogue,
>   * a bolometric correction,
>   * the dust column between the BH and the camera, measured through the same RT grid.

---

## 1. Quick start

```bash
cd /path/to/RTGen-tripleAGN
python postproc/band_luminosity.py $RUNS_DIR/sys2_snap003072 --band-rest 0.48 0.54
```

Each run of the script uses exactly one band:

| Option | Band |
|---|---|
| `--band-rest LO HI` | Top-hat in **rest-frame** wavelength, µm (e.g. `0.5 0.6` ≈ rest V) |
| `--band-obs LO HI` | Top-hat in **observed-frame** wavelength, µm (shifted to the rest frame by 1/(1+z)) |
| `--filter FILE` | A filter curve: two columns, wavelength in Å and throughput (the format of `filters/HST`, `filters/JWST`, `filters/Euclid`). Observed frame unless you add `--filter-frame rest` |
| `--xray EMIN EMAX` | X-ray band, keV. **AGN only**, because the RT has no X-ray emission |

Other options:

| Option | Default | Meaning |
|---|---|---|
| `--aperture-kpc R` | 1 | Radius of the circle around each BH. A pixel inside two circles goes to the nearer BH |
| `--lbol H N1 N2` | – | AGN bolometric luminosities, erg/s |
| `--catalog CSV` | – | Read L_bol from a picks CSV (`selection/find_dust_analogues.py` output), matching the row by host BH id |
| `--agn-kbol K` | built-in | Bolometric correction K = L_bol / νL_ν to use for this band (§3.2) |
| `--no-agn-dust` | off | Skip the line-of-sight dust column (faster; the AGN is then unattenuated) |
| `--json FILE` | – | Also write every number to a JSON file |

Where the AGN L_bol comes from, first match wins:

1. `--lbol`;
2. `bh_lbol_erg_s` in the run's `input/input_params.py`. `scripts/make_run.sh` fills this in when you give it the three L_bol values, and `scripts/runs_from_picks.py` always does;
3. `--catalog`.

Stage 4 (`submit_dust_obscuration.sh`) already runs the script once for rest-frame
0.5–0.6 µm and writes `band_luminosity_V_rest.txt` in the run directory. Run it
by hand for any other band.

### Example: sys2 (z = 1.185), rest-frame 4800–5400 Å, r = 1 kpc

```text
Stars (+dust emission) from the RT images, aperture r = 1 kpc around each BH
        L_band dusty L_band nodust  A_band  nuLnu dusty M_AB dusty
       [erg/s]       [erg/s]       [mag]   [erg/s]      [rest-band]
Host      1.1044e+43    1.5262e+43    0.35   9.3501e+43     -21.41
N1        1.3912e+42    1.4209e+42    0.02   1.1773e+43     -19.16
N2        1.3914e+42    1.4126e+42    0.02   1.1773e+43     -19.16
Total     3.6910e+43    4.0736e+43    0.11   3.1235e+44     -22.72

AGN estimate (NOT from the RT): L_bol from dust_analogues.csv; nuL_nu = L_bol / K, K = 9.26 (Shen+08 BC at 5100 A);
            L_bol  nuLnu,int L_band,int tau_band L_band,obs   Sigma_d   N_H(los) AGN frac
Host    8.268e+43  8.928e+42  1.052e+42     0.50  6.402e+41  2.60e-05  2.152e+21    0.055
N1      8.857e+40  9.565e+39  1.127e+39     0.77  5.191e+38  4.07e-05  2.536e+21    0.000
N2      1.651e+40  1.783e+39  2.100e+38     0.16  1.789e+38  8.40e-06  8.341e+20    0.000
```

How to read it. In this band at 1 kpc, the host galaxy's starlight is about
1.1e43 erg/s after dust. The host AGN adds 6.4e41 erg/s after a dust optical
depth of 0.5 along the camera's sightline, i.e. 5.5% of the total. The two
companions' AGN contribute nothing measurable in the optical. (The N1 and N2
starlight values agreeing to 0.01% is a coincidence: N1 loses about 13 pixels
to the host aperture.)

The same run, other bands:

```bash
# Observed HST/WFC3 F160W (prints apparent AB magnitudes too)
python postproc/band_luminosity.py RUN --filter filters/HST/HST_WFC3_IR.F160W.dat
# Hard X-ray, intrinsic AGN luminosity and line-of-sight N_H
python postproc/band_luminosity.py RUN --xray 2 10
# A band with no built-in AGN correction: supply your own K = L_bol / nuL_nu
python postproc/band_luminosity.py RUN --band-rest 2.0 2.4 --agn-kbol 12
```

---

## 2. Host galaxy: band luminosity from the RT images

### 2.1 What is in the image files

| File | Content |
|---|---|
| `image_cont.out` | Dusty continuum image (stage 3) |
| `dust_obscuration_analysis/image_cont_nodust.out` | The same camera with the dust removed (stage 4), i.e. intrinsic starlight |

Both are RADMC-3D ASCII images:

* **Header:** `iformat`, `nx ny`, `nλ`, pixel size in cm (x and y), then the nλ
  wavelengths in µm.
* **Body:** for each wavelength, `nx × ny` intensities. x runs fastest. Read it as
  `image[λ, iy, ix]`; `read_radmc3d_image()` in `postproc/plot_dust_obscuration.py` does this.
* **Units:** specific intensity I_ν, in erg s⁻¹ cm⁻² Hz⁻¹ sr⁻¹.
* **Frame:** **rest-frame**. RADMC-3D knows nothing about redshift.
* **Camera:** set by `incl 90 sizepc 20000 npixx 256 npixy 256`, i.e. 20 kpc
  across and 78 pc per pixel. It has 150 wavelengths from 0.01 to 10⁴ µm,
  logarithmically spaced (0.032 dex ≈ 7.6% apart).

### 2.2 From intensity to luminosity

The luminosity density seen along this sightline (isotropic-equivalent) is

```
L_ν(λ) = 4π · A_pix · Σ_pixels I_ν(λ)          [erg s⁻¹ Hz⁻¹],   A_pix = Δx · Δy  [cm²]
```

The distance drops out: flux at distance d is `Σ I_ν A_pix / d²`, and `L = 4π d² F`.
This is the same factor as in the original `postproc/compute_sed.py`. Sum over
only the pixels of an aperture to get that aperture's L_ν.

### 2.3 Band integrals

Take a throughput T(λ_rest). For a top-hat, T = 1 inside the band. The script reports:

| Quantity | Definition | Use |
|---|---|---|
| `L_band` | ∫ L_ν T dν  [erg/s] | Energy emitted in the band |
| ⟨L_ν⟩ | ∫ L_ν T dν/ν ⁄ ∫ T dν/ν  [erg/s/Hz] | Photon-counting mean, the AB-magnitude convention |
| `nuLnu` | ν_pivot ⟨L_ν⟩ [erg/s] | "νL_ν" at the band's pivot wavelength |
| `M_AB` | −2.5 log ⟨L_ν⟩/(4π (10 pc)²) − 48.60 | Absolute AB magnitude in the rest-frame band |
| `A_band` | −2.5 log (⟨L_ν⟩_dusty / ⟨L_ν⟩_dust-free) | Dust attenuation of the starlight in that band |

L_ν is interpolated in log–log from the 150 camera wavelengths onto a fine grid
inside the band before integrating.

**Observed-frame bands** (`--band-obs`, `--filter`):

* Light seen in observed band [λ₁, λ₂] left the galaxy in rest band
  [λ₁, λ₂]/(1+z). The script shifts the filter by 1/(1+z) and integrates the rest-frame spectrum.
* The result `L_band` is the luminosity that ends up in the observed band, so the
  observed band flux is F = L_band / (4π D_L²).
* The apparent magnitude uses F_ν,obs = (1+z) ⟨L_ν⟩ / (4π D_L²):
  m_AB = −2.5 log F_ν,obs − 48.60.
* D_L uses the Romulus25 cosmology (flat, H₀ = 67.77, Ω_m = 0.3086).

### 2.4 Which pixels belong to which galaxy

**Image orientation.** `sph_interpolation.py` writes the grid with axes
(sim x, sim z, sim y) → (RADMC x, y, z). With `incl 90`, the camera therefore
looks along the **simulation z axis**:

* image x = simulation x; image y = simulation y;
* the image is centred on the host BH;
* the observer is on the −z side.

This was checked on the outputs:

* The projected star-particle mass in (x, y) matches the dust-free 1.6 µm
  image with correlation 0.9998, better than any of the other 23 axis/sign
  combinations.
* An attenuation model with dust in front from −z reproduces the observed
  A_V map with correlation 0.98, against 0.55 for +z (sys1).

**Each BH's projected position** is its (Δx, Δy) offset from the host, taken from:

* `sph_interpolation/bh_info.txt`, written by `rewrite_sim.py`;
* or, for runs made before that file existed, the "offset from grid centre"
  lines of `sph_interpolation/job_sph_interpolation.out`.

Pixel (iy, ix) is centred at

```
x = (ix + 0.5 − nx/2) · Δx,   y = (iy + 0.5 − ny/2) · Δy
```

**Apertures.** A pixel belongs to member k if it lies within `--aperture-kpc`
of BH k *and* BH k is the nearest of the three. That way pixels aren't counted
twice when the circles overlap, which they do: the separations are about
1–2 kpc. `Total` is the whole 20 × 20 kpc image.

> In a late-stage merger, "the host galaxy of N1" isn't well defined. Report the
> aperture radius with any number, and check how much it changes between, e.g.,
> 0.5 and 1 kpc. Projected separations here are along sim z, so they differ
> from the catalogue's `Separation_*_2D_kpc` if the catalogue projected along
> another axis.

### 2.5 Doing it by hand

```python
import sys, numpy as np
sys.path.insert(0, 'postproc')
from plot_dust_obscuration import read_radmc3d_image

x0, y0 = 1.51, 0.89     # projected BH offset from the host [kpc] (bh_info.txt: dx, dy)
l1, l2 = 0.5, 0.6       # rest-frame band [µm]

im = read_radmc3d_image('RUN/image_cont.out')        # or .../image_cont_nodust.out
lam = im['lam']                                         # rest-frame µm
Lnu = 4*np.pi * im['sizepix_x']*im['sizepix_y'] * im['image'].sum(axis=(1, 2))   # erg/s/Hz, whole image

# aperture of radius 1 kpc around a BH at projected (x0, y0) kpc
kpc = 3.0857e21; n = im['nx']; px = im['sizepix_x']/kpc
c = (np.arange(n) + 0.5 - n/2) * px
X, Y = np.meshgrid(c, c)
m = (X - x0)**2 + (Y - y0)**2 <= 1.0**2
Lnu_ap = 4*np.pi * im['sizepix_x']*im['sizepix_y'] * im['image'][:, m].sum(axis=1)

# energy in a rest-frame top-hat band [l1, l2] µm: interpolate log L_nu onto a
# fine grid first (the camera has only ~3 samples per 10% in wavelength)
lf = np.geomspace(l1, l2, 1000)
Lf = np.exp(np.interp(np.log(lf), np.log(lam), np.log(np.maximum(Lnu_ap, 1e-300))))
L_band = np.trapz(Lf * 2.998e14 / lf, np.log(lf))          # = ∫ L_nu dnu   [erg/s]
```

---

## 3. AGN: band luminosity estimate

### 3.1 Bolometric luminosity

The catalogue's `BH_Lbol`, `Neighbour1BH_Lbol` and `Neighbour2BH_Lbol` are
L_bol = ε_r Ṁ c² with ε_r = 0.1. That is the Romulus value; the CSV has
L_bol / Ṁ = 0.1 c² exactly. This is the **intrinsic** luminosity. Ṁ is the
accretion rate at that one snapshot, so it can vary a lot from step to step.

### 3.2 Bolometric corrections (intrinsic band luminosity)

The intrinsic band luminosity is νL_ν = L_bol / K for optical/UV bands, and
L_X = L_bol / K_X for X-ray bands. Built-in corrections:

| Band | K | Source | When the script uses it |
|---|---|---|---|
| 1350 Å | 3.81 | Shen et al. 2008 (Richards et al. 2006 mean quasar SED) | rest-frame pivot within ±15% of 1350 Å |
| 3000 Å | 5.15 | same | within ±15% of 3000 Å |
| 5100 Å | 9.26 | same | within ±15% of 5100 Å |
| 2–10 keV | K_X = 10.96 [1 + (log(L_bol/L☉)/11.93)^17.79] | Duras et al. 2020 (general AGN sample) | `--xray` |
| other X-ray bands | 2–10 keV value × band fraction for photon index Γ = 1.8 | power law | `--xray E1 E2` (e.g. 2–8 keV is 0.84 × L(2–10)) |

For any other band (near-/mid-IR, other UV/optical wavelengths, radio), pass
your own `--agn-kbol K`, taken from the SED template or relation you want to
use, e.g. the Hopkins et al. 2007 or Shen et al. 2020 bolometric corrections.
For a top-hat band, the in-band intrinsic luminosity assumes a flat νL_ν across
the band: L_band = νL_ν · ln(λ₂/λ₁).

These are mean relations for unobscured, luminous AGN. Each has a scatter of
roughly 0.2–0.3 dex. They are least reliable at the low L_bol of the
companions (≲ 10⁴¹ erg/s).

### 3.3 Dust attenuation along the camera sightline

To be consistent with the images, the AGN is dimmed by the dust that is
actually in the RT grid between the BH and the camera:

```
τ_λ = Σ_species κ_ext,s(λ) · Σ_s,      Σ_s = ∫ ρ_dust,s dl   (from the BH towards −z_sim)
L_band,obs = L_band,int · ⟨e^(−τ_λ)⟩_band
```

* ρ_dust,s is read from `dust_density.inp`, the dust actually used in the RT,
  including the metallicity-dependent dust-to-gas ratio.
* κ_ext = κ_abs + κ_scat, per gram of dust, from the run's `dustkappa_*.inp`.
  For a point source, light scattered out of the sightline is lost.
* The column starts halfway through the BH's own cell. Cells are 117 pc (30 kpc / 256).
* `N_H(los)` is the same column in gas: X ρ_gas / m_p with X = 0.76, from `input/mass_gas.DAT`.
  The X-ray output reports it, but the script does **not** apply X-ray
  absorption. Use your own absorption model (e.g. photoelectric absorption for
  N_H) to go from intrinsic to observed L_X.

In the example above, the dust column towards the host AGN of sys2 gives
τ ≈ 0.5 at 5100 Å, and N_H ≈ 2×10²¹ cm⁻².

**What this column misses:**

* **Sub-resolution obscuration.** There is no torus and no nuclear gas below
  ~100 pc, and the SPH smoothing spreads out dense gas. This column is a
  **lower limit** on the true obscuration. For comparison, UGC 2369S's N core
  has N_H ~ 10^25.5 cm⁻².
* **Other viewing angles.** The column is for this one camera direction (sim −z);
  other directions need a new column.
* **Reprocessed AGN light.** Light the dust absorbs is not re-emitted in this
  estimate: no AGN-heated far-IR. That needs AGN sources in the RT (§5).

### 3.4 Totals

For each member, `L_band_total = L_band(stars, dusty, aperture) + L_band,obs(AGN)`,
and `AGN frac` is the AGN's share of that total. The `--json` output holds
every term separately.

---

## 4. Caveats (all bands)

* **Starlight only, no AGN emission, in the images.** See §3 and §5.
* **No nebular emission.** There are no emission lines (Hα, [O III], …) and no
  nebular continuum, so narrow-band and line luminosities can't come from these images.
* **Wavelength sampling.** The sampling is 0.032 dex. Bands narrower than about
  two samples (narrow filters, `N` filters) are only interpolated; the script
  warns. For those, re-render with more wavelengths: edit
  `camera_wavelength_micron.inp` and re-run the `radmc3d image loadlambda …`
  lines of stages 3 and 4.
* **Above 160 µm (rest).** Starburst99 has no data there, and RTGen
  extrapolates it. The dust-free image (and A_band) there is not physical; the
  dusty image (dust emission) is fine.
* **One sightline.** Everything is for the fixed `incl 90` camera
  (sim +z direction). Luminosities are isotropic-equivalent for that sightline.
* **Monte Carlo noise.** Faint outer pixels carry photon noise
  (nphot = 10⁷). Totals and kpc apertures are robust; single pixels are not.

---

## 5. Putting the AGN into the RT instead (not done here)

To get self-consistent AGN light, including dust heated by the AGN, its far-IR
emission, and scattered AGN light, the AGN would have to be a RADMC-3D source:

* Write a `stars.inp` with one point source per BH, at the grid positions from
  `bh_info.txt` in the RADMC axis order (x, z, y).
* Give each source an AGN SED normalised to its L_bol.
* Re-run `mctherm` and the images.

The dusty and dust-free images would then contain the AGN, and
`band_luminosity.py` would read the AGN directly from the images. A small
aperture around each BH would hold the AGN plus nuclear starlight.
