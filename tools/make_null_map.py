"""Strain-free null map for noise-floor tests: copy a prepared work dir, replace every usable pattern of each large grain by ONE pattern of that grain
plus independent Gaussian noise with the pixel noise measured between real neighbours. True strain, rotation and GND are then zero.
usage: make_null_map.py --src work\hrebsd_Pt6gnd --stem Pt6gnd --out work\hrebsd_Pt6null --name Pt6null [--grains 1 3] [--noise-scale 1.0]"""
import argparse, os, shutil, numpy as np, tifffile
ap = argparse.ArgumentParser()
ap.add_argument("--src", required=True); ap.add_argument("--stem", required=True); ap.add_argument("--out", required=True); ap.add_argument("--name", required=True)
ap.add_argument("--grains", nargs="+", type=int, default=[1, 3]); ap.add_argument("--noise-scale", type=float, default=1.0); ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args(); rng = np.random.default_rng(a.seed)
d = np.load(os.path.join(a.src, "maps.npz")); nr, nc = int(d["nrows"]), int(d["ncols"]); grain = d["grain"]; use = d["usable"]
os.makedirs(os.path.join(a.out, "Corrected_EBSPs"), exist_ok=True)
for f in os.listdir(a.src):
    p = os.path.join(a.src, f)
    if os.path.isfile(p) and f.startswith(a.stem) and f.endswith((".ang", ".txt")) or f == "first_image.txt":
        dst = f.replace(a.stem, a.name) if f.startswith(a.stem) else f; shutil.copy(p, os.path.join(a.out, dst))
tif = lambda r, c, s=a.stem: os.path.join(a.src, "Corrected_EBSPs", f"{s}_r{r}c{c}.tif")
first = open(os.path.join(a.src, "first_image.txt")).read().strip().replace(a.stem, a.name); open(os.path.join(a.out, "first_image.txt"), "w").write(first)
done = np.zeros((nr, nc), bool)
for g in a.grains:
    m = use & (grain == g); rr, cc = np.nonzero(m)
    k = len(rr) // 2; r0, c0 = rr[k], cc[k]; base = tifffile.imread(tif(r0, c0)).astype(float)
    # pixel noise from horizontal neighbour pairs inside the grain (includes the tiny real rotation, so slightly conservative)
    sig = []
    for r, c in zip(rr[::7], cc[::7]):
        if c + 1 < nc and m[r, c + 1]: sig.append(np.std(tifffile.imread(tif(r, c)).astype(float) - tifffile.imread(tif(r, c + 1)).astype(float)) / np.sqrt(2))
    s = np.median(sig) * a.noise_scale; print(f"grain {g}: {len(rr)} pts, base r{r0}c{c0}, pixel noise sigma {s:.2f} (signal std {base.std():.1f}, mean {base.mean():.0f})")
    for r, c in zip(rr, cc):
        tifffile.imwrite(os.path.join(a.out, "Corrected_EBSPs", f"{a.name}_r{r}c{c}.tif"), np.clip(base + rng.normal(0, s, base.shape), 0, 65535).astype(np.uint16)); done[r, c] = True
for r in range(nr):
    for c in range(nc):
        if not done[r, c]: shutil.copy(tif(r, c), os.path.join(a.out, "Corrected_EBSPs", f"{a.name}_r{r}c{c}.tif"))
print("wrote", a.out)
