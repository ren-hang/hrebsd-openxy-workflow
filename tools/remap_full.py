"""Full-map remapping pass around OpenXY. Pass 1 = an existing hrebsd.py run (work\hrebsd_<NAME>). Every non-reference pattern is un-rotated by the rotation R1 of its
pass-1 F (rays in the crystal frame of its grain's reference pattern, OpenXY/gen_synthetic geometry), pass 2 is run through hrebsd.py, and the final F = R1 . F2.
Writes work\hrebsd_<NAME>_remap\maps.npz (same keys as pass 1, F / strain / rotation replaced) + png maps, and prints the pass-1 vs final difference.
usage: remap_full.py <MAP.XTAL folder> --name Pt6full --pc X Y L [--elev 11.04] [--tilt 70] [--stage prep|run|post|all]  (other args as hrebsd.py)"""
import argparse, glob, os, shutil, subprocess, sys, numpy as np, tifffile
import scipy.io as sio
from scipy.linalg import polar
from scipy.ndimage import map_coordinates
from scipy.spatial.transform import Rotation as SR
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_synthetic import euler2gmat
from hrebsd import brady, PY

ap = argparse.ArgumentParser(); ap.add_argument("map"); ap.add_argument("--name", required=True); ap.add_argument("--pc", nargs=3, type=float, required=True)
ap.add_argument("--elev", type=float, default=11.04); ap.add_argument("--tilt", type=float, default=70.0); ap.add_argument("--nproc", type=int, default=3)
ap.add_argument("--material", default="platinum"); ap.add_argument("--min-rot", type=float, default=0.5, help="only remap points whose pass-1 rotation is at least this many degrees (below it remapping only adds noise, NOTES 13.5)"); ap.add_argument("--stage", default="all", choices=["prep", "run", "post", "all"])
a = ap.parse_args()
w1 = os.path.join("work", "hrebsd_" + a.name); n2 = a.name + "_remap"; w2 = os.path.join("work", "hrebsd_" + n2)
mat = sorted(glob.glob(os.path.join(w1, "out", "AnalysisParams_*.mat")), key=os.path.getmtime)[-1]
S = sio.loadmat(mat, squeeze_me=True, struct_as_record=False)["Settings"]; F1 = np.transpose(S.data.F, (2, 0, 1)); n = len(F1); ref = np.asarray(S.RefInd).astype(int) - 1
ang = np.array([l.split()[:3] for l in open(os.path.join(w1, a.name + ".ang")) if not l.startswith("#") and l.strip()], float)

rot1 = np.array([np.degrees(np.linalg.norm(SR.from_matrix(polar(f)[0]).as_rotvec())) if np.isfinite(f).all() else 0.0 for f in F1])

def prep():
    xs, ys, zs = brady(*a.pc); al = np.pi / 2 - np.radians(a.tilt) + np.radians(a.elev)
    Qps = np.array([[0, -np.cos(al), -np.sin(al)], [-1, 0, 0], [0, np.sin(al), -np.cos(al)]])
    os.makedirs(os.path.join(w2, "Corrected_EBSPs"), exist_ok=True)
    for f in os.listdir(w1):
        if f.startswith(a.name) and f.endswith((".ang", ".txt")): shutil.copy(os.path.join(w1, f), os.path.join(w2, f.replace(a.name, n2, 1)))
    open(os.path.join(w2, "first_image.txt"), "w").write(open(os.path.join(w1, "first_image.txt")).read().replace(a.name, n2))
    src = lambda i: os.path.join(w1, "Corrected_EBSPs", f"{a.name}_r{i // nc}c{i % nc}.tif")
    t0 = tifffile.imread(src(0)); P = min(t0.shape); X0 = (t0.shape[1] - P) // 2
    rr, cc = np.mgrid[1:P + 1, 1:P + 1].astype(float); pts = np.stack([(xs * P - cc).ravel(), ((1 - ys) * P - rr).ravel(), np.full(cc.size, -zs * P)])
    Mcache = {}; done = 0
    for i in range(n):
        im = tifffile.imread(src(i))[:, X0:X0 + P]
        if i != ref[i] and np.isfinite(F1[i]).all() and np.abs(F1[i] - np.eye(3)).max() < 0.2 and rot1[i] >= a.min_rot:
            if ref[i] not in Mcache: Mcache[ref[i]] = euler2gmat(*ang[ref[i]]) @ Qps
            M = Mcache[ref[i]]; R1 = polar(F1[i])[0]; vp = M.T @ (R1 @ (M @ pts)); s = -zs * P / vp[2]
            c = xs * P - s * vp[0]; r = (1 - ys) * P - s * vp[1]
            im = np.clip(map_coordinates(im.astype(float), [(r - 1).reshape(P, P), (c - 1).reshape(P, P)], order=3, mode="nearest"), 0, 65535).astype(np.uint16); done += 1
        tifffile.imwrite(os.path.join(w2, "Corrected_EBSPs", f"{n2}_r{i // nc}c{i % nc}.tif"), im)
    print("remapped %d of %d patterns -> %s (square crop %d px)" % (done, n, w2, P))

def hr(stage):
    cmd = [PY, "-I", "tools/hrebsd.py", a.map, "--pc", *map(str, a.pc), "--name", n2, "--elev", str(a.elev), "--nproc", str(a.nproc), "--stage", stage]
    r = subprocess.run(cmd, capture_output=True, text=True); print(r.stdout[-600:], r.stderr[-600:])
    if r.returncode: sys.exit("hrebsd.py %s failed" % stage)

def post():
    from replot_maps import replot
    d1 = np.load(os.path.join(w1, "maps.npz")); d2 = dict(np.load(os.path.join(w2, "maps.npz"))); nr, nc_ = int(d2["nrows"]), int(d2["ncols"])
    F2 = d2["F"].reshape(-1, 3, 3); F = F1.copy()
    for i in range(n):
        R1 = polar(F1[i])[0] if (i != ref[i] and np.isfinite(F1[i]).all() and np.abs(F1[i] - np.eye(3)).max() < 0.2 and rot1[i] >= a.min_rot) else np.eye(3)
        F[i] = R1 @ F2[i]
    U = np.array([polar(f)[1] - np.eye(3) for f in F]); dev = U - np.trace(U, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
    rv = np.array([SR.from_matrix(polar(f)[0]).as_rotvec() for f in F]) * 180 / np.pi; g = lambda v: np.asarray(v).reshape((nr, nc_) + np.shape(v)[1:])
    d2.update(F=g(F), dev_strain=g(dev), rot_deg=g(np.linalg.norm(rv, axis=1)), rotx=g(rv[:, 0]), roty=g(rv[:, 1]), rotz=g(rv[:, 2]), strain_norm=g(np.linalg.norm(dev.reshape(n, 9), axis=1)))
    for k, (p, q) in dict(e11=(0, 0), e12=(0, 1), e13=(0, 2), e22=(1, 1), e23=(1, 2), e33=(2, 2)).items(): d2[k] = g(dev[:, p, q])
    d2["usable"] = d2["usable"] & d1["usable"]; d2["F_pass1"] = d1["F"]; d2["dev_strain_pass1"] = d1["dev_strain"]
    for k in [k for k in d2 if k.startswith("gnd") or k in ("nye_alpha",)]: del d2[k]
    np.savez_compressed(os.path.join(w2, "maps.npz"), **d2); replot(os.path.join(w2, "maps.npz"))
    u = d2["usable"]; dd = np.linalg.norm((d2["dev_strain"] - d1["dev_strain"]).reshape(-1, 9), axis=1).reshape(nr, nc_)
    print("usable in both %d; |final - pass1| deviatoric strain (x1e-3): median %.3f  p90 %.3f  max %.3f" % (u.sum(), np.median(dd[u]) * 1e3, np.percentile(dd[u], 90) * 1e3, dd[u].max() * 1e3))
    for gi in np.unique(d2["grain"][u]):
        m = u & (d2["grain"] == gi)
        if m.sum() > 50: print("  grain %d n=%d: median strain norm pass1 %.3f -> final %.3f (x1e-3); median rotation %.4f deg" % (gi, m.sum(), np.median(d1["strain_norm"][m]) * 1e3, np.median(d2["strain_norm"][m]) * 1e3, np.median(d2["rot_deg"][m])))

nc = int(np.load(os.path.join(w1, "maps.npz"))["ncols"])
if a.stage in ("prep", "all"): prep()
if a.stage in ("run", "all"): hr("run")
if a.stage in ("post", "all"): hr("post"); post()
