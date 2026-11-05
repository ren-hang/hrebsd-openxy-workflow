"""One-command OpenXY hrEBSD workflow for Thermo Pathfinder maps (backscatter, 70 deg tilt, 20 kV).
  python tools\\hrebsd.py MAP_FOLDER --pc X Y L [--name NAME] [--block 10] [--phase Platinum] [--material platinum] [--stage check|prepare|run|post|all]
MAP_FOLDER holds <name>.pat and <name>.ANG (the *.MAP.XTAL folder, a local copy; nothing is written there).
Stages: check   -> preflight report (pattern size, scan size, indexing, phase, PC sanity); exits 1 on problems that would give wrong numbers
        prepare -> OpenXY inputs in work\hrebsd_<NAME>: tifs, orientation-corrected .ang (Rz(az) Rcorr), grain file, reference kept away from boundaries
        run     -> OpenXY headless (MATLAB), log in logs\hrebsd_<NAME>.log (--gnd adds the dislocation-density calculation)
        post    -> maps.npz (all values, replottable), results.csv, summary.md, strain / rotation / overview PNGs, GND PNGs with --gnd
Strain values are produced, not interpreted. Run from the repository root."""
import argparse, csv, glob, os, re, struct, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

PY = os.path.join("work", "kpenv", "Scripts", "python.exe")
RMAT = os.path.join("inputs", "Deba_calibration", "EBSD_Rcorr_calibration_20260505.mat")
CODE = os.path.join("work", "OpenXY_BYU", "Code")
PHOSPHOR_UM = 93542.4                       # phosphor height in um: mperpix = PHOSPHOR_UM / pattern height (161.28 for 580 px, 80.64 for 1160 px)
DEFAULT_PC = (15.8, 7.9, 16.0)              # instrument default = "not entered"; never correct


def brady(X, Y, L):                         # Pathfinder X,Y,L (mm) -> TSL/OpenXY x*,y*,z*  (Brady E.1-E.3)
    return 0.5 + 0.5 * ((X - 15.8) / 20), 1 - Y / 20, L / 20


def find(folder, pat):
    f = [p for p in glob.glob(os.path.join(folder, "*")) if re.fullmatch(pat, os.path.basename(p), re.I)]
    return f[0] if f else None


def rows_of(path, ncol=10):
    return np.array([l.split()[:ncol] for l in open(path, errors="replace") if not l.startswith("#") and l.strip()], float)


def check(a):
    problems, notes = [], []
    pat, ang = find(a.map, r".*\.pat"), find(a.map, r".*\.ang")
    if not pat: problems.append("no .pat file in the folder")
    if not ang: problems.append("no .ANG file in the folder")
    if problems: return None, problems, notes
    with open(pat, "rb") as f: _, nc, nr, pw, ph = struct.unpack("<5I", f.read(20))
    size, exp = os.path.getsize(pat), 25 + nc * nr * (pw * ph * 2 + 8) - 8
    notes.append(f".pat: {nc} x {nr} points, patterns {pw} x {ph} px (OpenXY uses the central {ph} x {ph} square); file {size/1e9:.2f} GB")
    if size < exp: problems.append(f".pat is shorter than its header says ({size} < {exp} bytes): incomplete copy (Box still downloading?)")
    h = "".join(l for l in open(ang, errors="replace") if l.startswith("#"))
    hn, hr = int(re.search(r"NCOLS:\s*(\d+)", h).group(1)), int(re.search(r"NROWS:\s*(\d+)", h).group(1))
    if (hn, hr) != (nc, nr): problems.append(f".ANG grid {hn} x {hr} differs from the .pat grid {nc} x {nr}")
    V = rows_of(ang)
    if len(V) != hn * hr: problems.append(f".ANG has {len(V)} rows, expected {hn*hr}")
    nan = ~np.isfinite(V[:, :3]).all(1)
    notes.append(f".ANG: {nan.sum()} not-indexed rows (NaN), {((V[:, 6] <= 0) & ~nan).sum()} more with CI <= 0; indexed with CI > 0: {(~nan & (V[:, 6] > 0)).mean():.1%}")
    xs = re.search(r"x-star\s+([-\d.eE+]+)", h)
    if xs and abs(float(xs.group(1))) < 1e-9: notes.append("the .ANG pattern-centre header is 0 (not stored): the PC comes from --pc (the hand-edited report)")
    nm, lat = re.search(r"MaterialName\s+(.*)", h), re.search(r"LatticeConstants\s+(.*)", h)
    notes.append(f".ANG phase '{nm.group(1).strip() if nm else '?'}', lattice {lat.group(1).strip() if lat else '?'}: replaced by '{a.phase}' / material '{a.material}' for OpenXY")
    if a.pc is None: problems.append("--pc X Y L (Pathfinder mm, from the hand-edited report) is required")
    elif np.allclose(a.pc, DEFAULT_PC, atol=1e-6): problems.append("PC equals the instrument default 15.8 / 7.9 / 16.0: that means 'not entered' and is never correct")
    else:
        x, y, z = brady(*a.pc); notes.append(f"PC X,Y,L = {a.pc} -> x*,y*,z* = {x:.5f}, {y:.5f}, {z:.5f}")
        if not (0.3 < x < 0.7 and 0.4 < y < 0.9 and 0.6 < z < 1.3): problems.append("PC outside the plausible range (check units / typing)")
    if not os.path.exists(RMAT): problems.append(f"calibration file missing: {RMAT}")
    return dict(pat=pat, ang=ang, nc=nc, nr=nr, pw=pw, ph=ph, mperpix=PHOSPHOR_UM / ph), problems, notes


def sh(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode: sys.exit("FAILED: " + " ".join(cmd) + "\n" + r.stdout[-1500:] + r.stderr[-1500:])
    return r.stdout


def prepare(a, info, w):
    raw = w + "_raw"; os.makedirs(os.path.join(w, "Corrected_EBSPs"), exist_ok=True)
    mk = [PY, "-I", "tools/make_pat_subset.py", "--pat", info["pat"], "--ang", info["ang"], "--out", raw, "--stem", a.name, "--phase", a.phase, "--ref-margin", str(a.ref_margin)]
    mk += ["--n", str(a.block)] if a.block else ["--full"]
    print(" ", sh(mk).strip().splitlines()[-1])
    print(" ", sh([PY, "-I", "tools/correct_ang.py", "--in", f"{raw}/{a.name}.ang", "--grain", f"{raw}/{a.name}.txt", "--rmat", RMAT, "--az", str(a.az), "--out", w, "--stem", a.name]).strip())
    for f in glob.glob(os.path.join(raw, "*.tif")): os.replace(f, os.path.join(w, "Corrected_EBSPs", os.path.basename(f)))
    for f in glob.glob(os.path.join(raw, "*_index.txt")): os.replace(f, os.path.join(w, os.path.basename(f)))
    open(os.path.join(w, "first_image.txt"), "w").write(os.path.basename(sorted(glob.glob(os.path.join(w, "Corrected_EBSPs", "*_r0c0.tif")))[0]))


def run(a, info, w):
    x, y, z = brady(*a.pc)
    env = dict(os.environ, OPENXY_HEADLESS_MATERIAL=a.material, OPENXY_NPROC=str(a.nproc), OPENXY_GND="1" if a.gnd else "0", OPENXY_CAMELEV=str(a.elev),
               OPENXY_MPERPIX=f"{info['mperpix']:.4f}", OPENXY_PC=f"{x:.6f} {y:.6f} {z:.6f}", OPENXY_FIRSTIMAGE=open(os.path.join(w, "first_image.txt")).read().strip())
    os.makedirs("logs", exist_ok=True)
    subprocess.run(["powershell", "-File", "tools/run_case.ps1", "-Log", "hrebsd_" + a.name, "-Code", CODE, "-Data", w, "-Stem", a.name, "-TimeoutSec", "1800" if a.block else "10800"], env=env, capture_output=True)
    log = open(os.path.join("logs", "hrebsd_" + a.name + ".log"), errors="replace").read()
    ok = "HREBSDMain returned" in log and "TIMEOUT" not in log      # run_case.ps1 leaves the MATLAB exit code blank
    print("  OpenXY " + ("finished" if ok else "FAILED"))
    if not ok: sys.exit("OpenXY did not finish cleanly, see logs\\hrebsd_%s.log" % a.name)


def post(a, info, w):
    import scipy.io as sio
    from scipy.linalg import polar
    from scipy.spatial.transform import Rotation as SR
    mats = sorted(glob.glob(os.path.join(w, "out", "AnalysisParams_*.mat")), key=os.path.getmtime)
    if not mats: sys.exit("no OpenXY result in %s\\out" % w)
    S = sio.loadmat(mats[-1], squeeze_me=True, struct_as_record=False)["Settings"]
    F = np.transpose(S.data.F, (2, 0, 1)); n = len(F)
    nc = a.block or info["nc"]; nr = n // nc
    grain = np.loadtxt(os.path.join(w, a.name + ".txt"), usecols=8, dtype=int)
    meta = np.loadtxt(os.path.join(w, a.name + "_index.txt"))               # original index, instrument CI, original column
    ci_work = np.loadtxt(os.path.join(w, a.name + ".txt"), usecols=6)
    near = (ci_work == 0) & (meta[:, 1] > 0)                                  # prepare zeroes CI within --ref-margin px of a grain boundary
    U = np.array([polar(f)[1] - np.eye(3) for f in F]); dev = U - np.trace(U, axis1=1, axis2=2)[:, None, None] / 3 * np.eye(3)
    rot = np.array([np.degrees(np.linalg.norm(SR.from_matrix(polar(f)[0]).as_rotvec())) for f in F])
    finite = np.isfinite(F).all((1, 2)); mag = np.abs(dev).reshape(n, 9).max(1)
    col0, notidx = meta[:, 2] == 0, meta[:, 1] <= 0
    usable = finite & ~col0 & ~notidx & (mag <= 5e-3)
    with open(os.path.join(w, "results.csv"), "w", newline="") as f:
        wr = csv.writer(f)
        wr.writerow(["index", "row", "col", "original_index", "grain", "instrument_CI", "col0", "not_indexed", "near_boundary", "usable", "rot_deg", "e11", "e12", "e13", "e22", "e23", "e33"])
        for i in range(n):
            wr.writerow([i, i // nc, i % nc, int(meta[i, 0]), grain[i], meta[i, 1], int(col0[i]), int(notidx[i]), int(near[i]), int(usable[i]), f"{rot[i]:.5f}"]
                        + [f"{dev[i][p, q]:.6e}" for p, q in ((0, 0), (0, 1), (0, 2), (1, 1), (1, 2), (2, 2))])
    with open(os.path.join(w, "summary.md"), "w") as f:
        f.write(f"# {a.name}\n\n- points {n} ({nr} x {nc}); OpenXY result `{os.path.basename(mats[-1])}`\n")
        f.write(f"- usable (finite, not column 0, instrument CI > 0, max |dev| <= 5e-3): {usable.sum()} ({usable.mean():.1%})\n")
        f.write(f"- flagged: column 0 {col0.sum()}, not indexed {notidx.sum()}, within {a.ref_margin} px of a grain boundary {near.sum()} (still in the table), non-finite {(~finite).sum()}\n")
        if usable.any(): f.write(f"- usable points: median max-component |dev| {np.median(mag[usable])*1e3:.2f}e-3, median rotation vs reference {np.median(rot[usable]):.3f} deg\n")
        f.write("- strain is deviatoric, OpenXY frame, relative to each grain's reference pattern; not interpreted\n")
    xy = np.array([l.split()[3:5] for l in open(os.path.join(w, a.name + ".ang")) if not l.startswith("#") and l.strip()], float)
    rv = np.array([SR.from_matrix(polar(f)[0]).as_rotvec() for f in F]) * 180 / np.pi
    grid_ = lambda v: np.asarray(v).reshape((nr, nc) + np.shape(v)[1:])
    save = dict(nrows=nr, ncols=nc, row=grid_(np.arange(n) // nc), col=grid_(np.arange(n) % nc), x_um=grid_(xy[:, 0]), y_um=grid_(xy[:, 1]), original_index=grid_(meta[:, 0].astype(int)),
                grain=grid_(grain), instrument_CI=grid_(meta[:, 1]), usable=grid_(usable), col0=grid_(col0), not_indexed=grid_(notidx), near_boundary=grid_(near), finite=grid_(finite),
                rot_deg=grid_(rot), rotx=grid_(rv[:, 0]), roty=grid_(rv[:, 1]), rotz=grid_(rv[:, 2]), strain_norm=grid_(np.linalg.norm(dev.reshape(n, 9), axis=1)), F=grid_(F), dev_strain=grid_(dev))
    for k, (p_, q_) in dict(e11=(0, 0), e12=(0, 1), e13=(0, 2), e22=(1, 1), e23=(1, 2), e33=(2, 2)).items(): save[k] = grid_(dev[:, p_, q_])
    ad = sio.loadmat(mats[-1], squeeze_me=True, struct_as_record=False).get("alpha_data")       # present when run with --gnd
    if ad is not None and hasattr(ad, "alpha_total5"):
        for k in ("alpha_total3", "alpha_total5", "alpha_total9"): save["gnd" + k[-1]] = grid_(np.asarray(getattr(ad, k), float))   # OpenXY GND estimates from 3 / 5 / 9 Nye components
        save["nye_alpha"] = grid_(np.transpose(np.asarray(ad.alpha_filt, float), (2, 0, 1))); save["gnd_filtered"] = grid_(np.asarray(ad.filteredInds, bool))
        save["gnd_misorientation"] = grid_(np.asarray(ad.misang, float))
        from gnd_binned import gnd_binned                                                   # k x k binned GND from the rotation field (noise floor falls steeply with k)
        tmp = dict(save)
        for k in (2, 3, 4): save["gnd_bin%d" % k] = gnd_binned(tmp, k)
    np.savez_compressed(os.path.join(w, "maps.npz"), **save)
    from replot_maps import replot
    replot(os.path.join(w, "maps.npz"))
    print(open(os.path.join(w, "summary.md")).read())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("map"); ap.add_argument("--pc", nargs=3, type=float); ap.add_argument("--name")
    ap.add_argument("--block", type=int, default=0, help="N: only an N x N block in the largest grain (quick test)")
    ap.add_argument("--phase", default="Platinum"); ap.add_argument("--material", default="platinum")
    ap.add_argument("--az", type=float, default=90.0, help="extra azimuth (deg) about the sample normal; empirical, see NOTES.md section 12")
    ap.add_argument("--elev", type=float, default=11.04); ap.add_argument("--ref-margin", type=int, default=4); ap.add_argument("--nproc", type=int, default=3); ap.add_argument("--gnd", action="store_true")
    ap.add_argument("--stage", default="all", choices=["check", "prepare", "run", "post", "all"])
    a = ap.parse_args()
    if not a.name: a.name = re.sub(r"\W+", "", os.path.basename(os.path.normpath(a.map)).split(".")[0]) + ("b" if a.block else "f")
    w = os.path.join("work", "hrebsd_" + a.name)
    info, problems, notes = check(a)
    print("== check"); [print("  " + n) for n in notes]
    if problems: print("PROBLEMS:"); [print("  - " + p) for p in problems]; sys.exit(1)
    if a.stage == "check": sys.exit(0)
    if a.stage in ("prepare", "all"): print("== prepare"); prepare(a, info, w)
    if a.stage in ("run", "all"): print("== run"); run(a, info, w)
    if a.stage in ("post", "all"): print("== post"); post(a, info, w)
