"""Compare two hrebsd.py runs on the points usable in both: difference of deviatoric strain (norm) and rotation.
usage: compare_maps.py runA/maps.npz runB/maps.npz"""
import sys, numpy as np
A, B = np.load(sys.argv[1]), np.load(sys.argv[2]); u = A["usable"].astype(bool) & B["usable"].astype(bool)
d = np.linalg.norm((A["dev_strain"] - B["dev_strain"])[u].reshape(-1, 9), axis=1) * 1e3
r = np.abs(A["rot_deg"] - B["rot_deg"])[u]
print("points usable in both: %d; |dev(A) - dev(B)| (1e-3) median %.3f p90 %.3f max %.3f; |rot A - rot B| (deg) median %.4f max %.4f" % (u.sum(), np.median(d), np.percentile(d, 90), d.max(), np.median(r), r.max()))
for g in np.unique(A["grain"][u]):
    m = u & (A["grain"] == g)
    if m.sum() < 30: continue
    dg = np.linalg.norm((A["dev_strain"] - B["dev_strain"])[m].reshape(-1, 9), axis=1) * 1e3
    print("  grain %d: n %d, median diff %.3fe-3, max %.3fe-3" % (g, m.sum(), np.median(dg), dg.max()))
