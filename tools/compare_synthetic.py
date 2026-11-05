"""Compare OpenXY output F (AnalysisParams_*.mat) with the imposed truth (truth.csv).

Both are reduced to  F = R U  (polar decomposition):  strain = U - I (reported as its deviatoric part,
because HR-EBSD cannot measure the dilatation - OpenXY fixes it with the sigma33=0 condition),
rotation = axis-angle vector of R (rad).  Frames: imposed F is in the crystal frame, which is the frame
OpenXY's F is expected in (confirmed by the e11 cases landing on F11).

usage: python tools/compare_synthetic.py --mat X.mat --truth truth.csv --label name [--md out.md] [--csv out.csv]
"""
import argparse
import csv

import numpy as np
import scipy.io as sio
from scipy.linalg import polar
from scipy.spatial.transform import Rotation


def decompose(F):
    R, U = polar(F)           # F = R U
    eps = U - np.eye(3)
    dev = eps - np.trace(eps) / 3 * np.eye(3)
    rv = Rotation.from_matrix(R).as_rotvec()
    return dev, rv, np.trace(eps)


def load(mat, truth):
    S = sio.loadmat(mat, squeeze_me=True, struct_as_record=False)["Settings"]
    Frec = np.transpose(S.data.F, (2, 0, 1))
    rows = list(csv.DictReader(open(truth)))
    Ftrue = np.array([[float(r[f"F{i}{j}"]) for i in (1, 2, 3) for j in (1, 2, 3)] for r in rows]).reshape(-1, 3, 3)
    return S, Frec, Ftrue, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mat", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--label", default="run")
    ap.add_argument("--md")
    ap.add_argument("--csv")
    ap.add_argument("--shift-truth", help="CSV from tools/shift_truth.m (idx + F column-major) giving the exact expected F for rigid shifts")
    a = ap.parse_args()
    S, Frec, Ftrue, rows = load(a.mat, a.truth)
    have_shift_truth = False
    if a.shift_truth:
        st = np.loadtxt(a.shift_truth, delimiter=",")
        for r in np.atleast_2d(st):
            Ftrue[int(r[0])] = r[1:].reshape(3, 3, order="F")
        have_shift_truth = True
    n = len(rows)
    out = []
    P = out.append
    P(f"## {a.label}\n")
    P(f"OpenXY reference index (1-based): {sorted(set(np.atleast_1d(S.RefInd).tolist()))}; PC used (x*,y*,z*) = "
      f"({np.mean(S.XStar):.5f}, {np.mean(S.YStar):.5f}, {np.mean(S.ZStar):.5f}); scan {S.Nx}x{S.Ny}; PixelSize {S.PixelSize}\n")

    D = []
    for i in range(n):
        dt, rt, _ = decompose(Ftrue[i])
        dr, rr, trr = decompose(Frec[i])
        D.append((dr - dt, rr - rt, dt, dr, rt, rr, trr))
    dE = np.array([d[0] for d in D])          # (n,3,3) dev-strain error
    dW = np.array([d[1] for d in D])          # (n,3) rotation error
    grp = np.array([r["group"] for r in rows])

    def maxabs(x):
        return float(np.max(np.abs(x))) if x.size else float("nan")

    P("| group | n | max abs dev-strain error | RMS dev-strain error | max abs rotation error (rad) | imposed range |")
    P("|---|---|---|---|---|---|")
    for g in ("ref", "strain", "rot", "mixed", "noise", "shift"):
        m = grp == g
        if not m.any():
            continue
        if g == "shift" and not have_shift_truth:
            P(f"| {g} | {m.sum()} | n/a (no analytic truth) | n/a | n/a | rigid shift 0.1-5 px |")
            continue
        tr = np.array([[np.abs(decompose(Ftrue[i])[0]).max(), np.abs(decompose(Ftrue[i])[1]).max()] for i in np.where(m)[0]])
        rng = f"strain<= {tr[:,0].max():.1e}, rot<= {tr[:,1].max():.1e} rad" if g != "noise" else "F=I, noise 1-20% of range"
        if g == "shift":
            rng = "rigid shift 0.1-5 px -> expected F from OpenXY forward model (" + rng + ")"
        P(f"| {g} | {m.sum()} | {maxabs(dE[m]):.2e} | {np.sqrt(np.mean(dE[m]**2)):.2e} | {maxabs(dW[m]):.2e} | {rng} |")

    # per-component gain for strain cases
    P("\n### Strain cases: recovered vs imposed (deviatoric component), per imposed component\n")
    P("| imposed | amplitude | imposed dev | recovered dev (same comp) | error | largest other component |")
    P("|---|---|---|---|---|---|")
    idx = {"e11": (0, 0), "e22": (1, 1), "e33": (2, 2), "e12": (0, 1), "e13": (0, 2), "e23": (1, 2)}
    for i in np.where(grp == "strain")[0]:
        nm = rows[i]["name"]
        k = nm.split("=")[0]
        p, q = idx[k]
        dt, dr = D[i][2], D[i][3]
        other = np.abs(dr - dt).copy()
        other[p, q] = 0
        other[q, p] = 0
        P(f"| {k} | {nm.split('=')[1]} | {dt[p,q]:+.2e} | {dr[p,q]:+.2e} | {dr[p,q]-dt[p,q]:+.1e} | {other.max():.1e} |")

    P("\n### Rotation cases: recovered vs imposed rotation vector (rad)\n")
    P("| case | imposed rotvec | recovered rotvec | max abs error | max dev-strain leak |")
    P("|---|---|---|---|---|")
    for i in np.where(grp == "rot")[0]:
        P(f"| {rows[i]['name']} | {np.array2string(D[i][4], precision=5, suppress_small=True)} | "
          f"{np.array2string(D[i][5], precision=5, suppress_small=True)} | {np.abs(D[i][1]).max():.1e} | {np.abs(D[i][3]).max():.1e} |")

    P("\n### Rigid pattern shifts (no analytic truth; first-order expectation: rotation ~ shift / (z* x P) rad, strain ~ 0)\n")
    zpx = float(np.mean(S.ZStar) * S.PixelSize)
    P(f"z* x P = {zpx:.1f} px\n")
    if have_shift_truth:
        P("| shift (dx,dy) px | expected rotvec (rad) | recovered rotvec (rad) | expected max dev-strain | recovered max dev-strain |")
        P("|---|---|---|---|---|")
        for i in np.where(grp == "shift")[0]:
            P(f"| ({rows[i]['dx']},{rows[i]['dy']}) | {np.array2string(D[i][4], precision=5, suppress_small=True)} | "
              f"{np.array2string(D[i][5], precision=5, suppress_small=True)} | {np.abs(D[i][2]).max():.1e} | {np.abs(D[i][3]).max():.1e} |")
        P("\n(Expected values are exact under OpenXY's own forward model, `tools/shift_truth.m`.)")
    else:
        P("| shift (dx,dy) px | expected |rot| ~ | recovered rotvec (rad) | recovered |rot| | max dev-strain |")
        P("|---|---|---|---|---|")
    for i in (np.where(grp == "shift")[0] if not have_shift_truth else []):
        dx, dy = float(rows[i]["dx"]), float(rows[i]["dy"])
        rv = D[i][5]
        P(f"| ({dx:g},{dy:g}) | {np.hypot(dx,dy)/zpx:.2e} | {np.array2string(rv, precision=5, suppress_small=True)} | "
          f"{np.linalg.norm(rv):.2e} | {np.abs(D[i][3]).max():.1e} |")

    P("\n### Noise floor (F = I, Gaussian noise as fraction of full range; reference pattern is noise-free)\n")
    P("| noise | n | mean max-abs dev-strain | max max-abs dev-strain | mean |rot| (rad) |")
    P("|---|---|---|---|---|")
    for lev in sorted({float(rows[i]["noise"]) for i in np.where(grp == "noise")[0]}):
        ii = [i for i in np.where(grp == "noise")[0] if float(rows[i]["noise"]) == lev]
        ms = [np.abs(D[i][3]).max() for i in ii]
        P(f"| {lev:.2f} | {len(ii)} | {np.mean(ms):.1e} | {np.max(ms):.1e} | {np.mean([np.linalg.norm(D[i][5]) for i in ii]):.1e} |")

    txt = "\n".join(out)
    print(txt)
    if a.md:
        open(a.md, "w", encoding="utf8").write(txt + "\n")
    if a.csv:
        with open(a.csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["idx", "group", "name"] + [f"dev_true_{p}{q}" for p in (1, 2, 3) for q in (1, 2, 3)]
                       + [f"dev_rec_{p}{q}" for p in (1, 2, 3) for q in (1, 2, 3)] + ["rot_true_x", "rot_true_y", "rot_true_z", "rot_rec_x", "rot_rec_y", "rot_rec_z"])
            for i in range(n):
                w.writerow([i, rows[i]["group"], rows[i]["name"], *D[i][2].ravel(), *D[i][3].ravel(), *D[i][4], *D[i][5]])


if __name__ == "__main__":
    main()
