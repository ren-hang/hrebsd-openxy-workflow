"""Two-grain mixture map from real patterns, to measure what a mixed pattern does to OpenXY (truth: zero strain, rotation and GND).
Copies a prepared work dir. In each of the two grains every usable pattern becomes ONE interior pattern of that grain plus the measured
pixel noise (as make_null_map.py), mixed with a fraction f of ONE interior pattern of the other grain:
  --levels 0 0.01 0.02 ...         f in contiguous bands of columns (bias vs f), or
  --measured work\gbmix_Pt6\mixing.npz   f = the fraction measured at each pixel by gb_mixing.py (a strain-free copy of the real boundary)
Points with f > 0 get IQ = CI = 0, Fit = 99 in the work copies so OpenXY picks a pure reference. All other points keep their real patterns.
usage: make_mix_map.py --src work\hrebsd_Pt6gnd --stem Pt6gnd --out work\hrebsd_Pt6mix --name Pt6mix --grains 1 3 (--levels ... | --measured ...)
Writes mix_truth.npz (f per pixel, base pattern positions) next to the work files; run it with hrebsd.py --stage run, then --stage post."""
import argparse, os, shutil, numpy as np, tifffile
from scipy import ndimage as ndi

ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True); ap.add_argument("--stem", required=True); ap.add_argument("--out", required=True); ap.add_argument("--name", required=True)
ap.add_argument("--grains", nargs=2, type=int, required=True); ap.add_argument("--levels", nargs="+", type=float); ap.add_argument("--measured")
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args(); rng = np.random.default_rng(a.seed)
if (a.levels is None) == (a.measured is None): ap.error("give exactly one of --levels / --measured")
d = np.load(os.path.join(a.src, "maps.npz")); nr, nc = int(d["nrows"]), int(d["ncols"]); grain = d["grain"]; use = d["usable"]
os.makedirs(os.path.join(a.out, "Corrected_EBSPs"), exist_ok=True)
tif = lambda r, c, s=a.stem: os.path.join(a.src, "Corrected_EBSPs", f"{s}_r{r}c{c}.tif")
first = open(os.path.join(a.src, "first_image.txt")).read().strip().replace(a.stem, a.name); open(os.path.join(a.out, "first_image.txt"), "w").write(first)
base, sigma, pos = {}, {}, {}
for g in a.grains:                                       # deepest usable interior point of each grain, and its pixel noise
    m = use & (grain == g); inner = ndi.distance_transform_edt(grain == g) * m
    r0, c0 = np.unravel_index(np.argmax(inner), inner.shape); base[g] = tifffile.imread(tif(r0, c0)).astype(float); pos[g] = (int(r0), int(c0))
    rr, cc = np.nonzero(m); sig = []
    for r, c in zip(rr[::7], cc[::7]):
        if c + 1 < nc and m[r, c + 1]: sig.append(np.std(tifffile.imread(tif(r, c)).astype(float) - tifffile.imread(tif(r, c + 1)).astype(float)) / np.sqrt(2))
    sigma[g] = float(np.median(sig)); print(f"grain {g}: {m.sum()} usable pts, base r{r0}c{c0}, pixel noise sigma {sigma[g]:.2f}")
f = np.zeros((nr, nc))
if a.levels:
    for g in a.grains:                                   # contiguous column bands of equal point count, one level each
        rr, cc = np.nonzero(use & (grain == g)); k = np.argsort(cc * nr + rr, kind="stable")
        for i, chunk in enumerate(np.array_split(k, len(a.levels))): f[rr[chunk], cc[chunk]] = a.levels[i]
else:
    fm = np.load(a.measured)["f_other"]; f = np.where(np.isfinite(fm) & np.isin(grain, a.grains) & use, np.clip(fm, 0, 1), 0.0)
done = np.zeros((nr, nc), bool)
for g, h in (a.grains, a.grains[::-1]):
    for r, c in zip(*np.nonzero(use & (grain == g))):
        im = (1 - f[r, c]) * base[g] + f[r, c] * base[h] + rng.normal(0, sigma[g], base[g].shape)
        tifffile.imwrite(os.path.join(a.out, "Corrected_EBSPs", f"{a.name}_r{r}c{c}.tif"), np.clip(im, 0, 65535).astype(np.uint16)); done[r, c] = True
for r in range(nr):
    for c in range(nc):
        if not done[r, c]: shutil.copy(tif(r, c), os.path.join(a.out, "Corrected_EBSPs", f"{a.name}_r{r}c{c}.tif"))
mixed = (f > 0).ravel()
for ext, cols in ((".ang", (5, 6, -1)), (".txt", (5, 6, 7))):             # IQ, CI, Fit columns: keep mixed points from being the reference
    out, i = [], 0
    for L in open(os.path.join(a.src, a.stem + ext)).read().splitlines():
        if L.startswith("#") or not L.strip(): out.append(L); continue
        if mixed[i]:
            t = L.split(); t[cols[0]], t[cols[1]], t[cols[2]] = "0", "0", "99"; L = "  ".join(t)
        out.append(L); i += 1
    open(os.path.join(a.out, a.name + ext), "w").write("\n".join(out) + "\n")
idx = os.path.join(a.src, a.stem + "_index.txt")
if os.path.exists(idx): shutil.copy(idx, os.path.join(a.out, a.name + "_index.txt"))
np.savez_compressed(os.path.join(a.out, "mix_truth.npz"), f=f, replaced=done, grains=np.array(a.grains), base_pos=np.array([pos[g] for g in a.grains]),
                    levels=np.array(a.levels if a.levels else []), measured=a.measured or "")
print(f"wrote {a.out}: {int((f > 0).sum())} mixed points, f max {f.max():.3f}")
