"""Apply a specimen-frame correction  R' = Rz(az) @ Rcorr @ R  (R = g^T = Bunge active matrix; same as MTEX 'Rcorr * ori') to the Euler columns
of an .ang and its grain .txt. Rcorr from Deba's EBSD_Rcorr_calibration .mat.  usage: correct_ang.py --in X.ang --grain X.txt --rmat M.mat --az 90 --out dir --stem Name"""
import argparse, os, numpy as np, scipy.io as sio
from scipy.spatial.transform import Rotation as SR
ap = argparse.ArgumentParser()
ap.add_argument("--in", dest="ang", required=True); ap.add_argument("--grain"); ap.add_argument("--rmat", required=True)
ap.add_argument("--az", type=float, default=0.0); ap.add_argument("--out", required=True); ap.add_argument("--stem", required=True)
a = ap.parse_args()
M = np.array(sio.loadmat(a.rmat, squeeze_me=True, struct_as_record=False)["calib"].Rmat, float)
C = SR.from_euler("z", a.az, degrees=True).as_matrix() @ M
os.makedirs(a.out, exist_ok=True)
def run(src, dst):
    out = []
    for L in open(src).read().splitlines():
        if L.startswith("#") or not L.strip(): out.append(L); continue
        v = L.split(); e0 = np.array(v[:3], float)
        if not np.all(np.isfinite(e0)): out.append(L); continue          # not-indexed rows (NaN) stay untouched
        R = SR.from_euler("ZXZ", e0).as_matrix()
        e = SR.from_matrix(C @ R).as_euler("ZXZ"); e[0] %= 2 * np.pi; e[2] %= 2 * np.pi
        v[:3] = [f"{x:.6f}" for x in e]; out.append("  ".join(v))
    open(dst, "w").write("\n".join(out) + "\n")
run(a.ang, os.path.join(a.out, a.stem + ".ang"))
if a.grain: run(a.grain, os.path.join(a.out, a.stem + ".txt"))
print("wrote", a.out, "correction angle %.3f deg" % np.degrees(np.linalg.norm(SR.from_matrix(C).as_rotvec())))
