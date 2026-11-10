"""What mixed patterns do to OpenXY. Reads a run made with make_mix_map.py (+ hrebsd.py --stage run, --stage post; truth = zero).
--levels runs: per grain and mixing fraction f, the apparent rotation and deviatoric strain norm (median, p90) and the systematic rotation
(norm of the mean rotation vector), against the f = 0 band of the same run.
--measured runs made with --gnd: OpenXY GND (gnd5) against distance into the grain (d_edge from gb_mixing.py), next to optional other runs.
usage: mix_bias.py --run work\hrebsd_Pt6mix --out DIR [--mixing work\gbmix_Pt6\mixing.npz --compare work\hrebsd_Pt6gnd work\hrebsd_Pt6null]
Error numbers only: the truth is zero, nothing here is a measured strain."""
import argparse, json, os, numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--run", required=True); ap.add_argument("--out", required=True); ap.add_argument("--mixing"); ap.add_argument("--compare", nargs="*", default=[])
ap.add_argument("--src", help="source work dir of the mixture map (only for runs whose mix_truth.npz predates the 'replaced' mask)")
a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
m = np.load(os.path.join(a.run, "maps.npz")); t = np.load(os.path.join(a.run, "mix_truth.npz"))
f, grains, levels = t["f"], [int(g) for g in t["grains"]], [float(x) for x in t["levels"]]
if "replaced" in t.files: replaced = t["replaced"]                       # points whose pattern make_mix_map.py rebuilt (others keep real patterns)
elif a.src: replaced = np.load(os.path.join(a.src, "maps.npz"))["usable"] & np.isin(m["grain"], grains)
else: replaced = np.ones_like(f, bool); print("warning: no 'replaced' mask; untouched real points are included")
grain, ok = m["grain"], m["finite"] & ~m["col0"] & ~m["not_indexed"] & replaced
rv = np.stack([m["rotx"], m["roty"], m["rotz"]], -1); res = dict(run=a.run)
if levels:
    tab = []
    for g in grains:
        for lv in levels:
            s = ok & (grain == g) & np.isclose(f, lv)
            if not s.any(): continue
            tab.append(dict(grain=g, f=lv, n=int(s.sum()), rot_median_deg=float(np.median(m["rot_deg"][s])), rot_p90_deg=float(np.percentile(m["rot_deg"][s], 90)),
                            rot_systematic_deg=float(np.linalg.norm(rv[s].mean(0))), strain_norm_median=float(np.median(m["strain_norm"][s])),
                            strain_norm_p90=float(np.percentile(m["strain_norm"][s], 90))))
    res["levels"] = tab
    print(f"{'grain':>5} {'f':>5} {'n':>4} {'rot med':>9} {'rot p90':>9} {'rot syst':>9} {'strain med':>11} {'strain p90':>11}   (deg, deg, deg, norm)")
    for x in tab: print(f"{x['grain']:5d} {x['f']:5.2f} {x['n']:4d} {x['rot_median_deg']:9.4f} {x['rot_p90_deg']:9.4f} {x['rot_systematic_deg']:9.4f} {x['strain_norm_median']:11.2e} {x['strain_norm_p90']:11.2e}")
    first = {}
    for g in grains:                                     # smallest f whose median error exceeds twice the f = 0 value of the same grain
        z = [x for x in tab if x["grain"] == g and x["f"] == 0]
        if not z: continue
        z = z[0]
        for x in sorted((x for x in tab if x["grain"] == g and x["f"] > 0), key=lambda x: x["f"]):
            if x["rot_median_deg"] > 2 * z["rot_median_deg"] or x["strain_norm_median"] > 2 * z["strain_norm_median"]: first[g] = x["f"]; break
    res["first_f_above_2x_null"] = first; print("smallest f with median error > 2x the unmixed band:", first)
if a.mixing:
    mx = np.load(a.mixing); de = mx["d_edge"]
    runs = [a.run] + list(a.compare); prof = {}
    for rdir in runs:
        q = np.load(os.path.join(rdir, "maps.npz"))
        if "gnd5" not in q.files: continue
        gnd, filt = q["gnd5"], q["gnd_filtered"] if "gnd_filtered" in q.files else np.zeros_like(q["gnd5"], bool)
        okq = q["finite"] & ~q["col0"] & ~q["not_indexed"] & ~filt & np.isin(q["grain"], grains) & (gnd > 0)
        prof[os.path.basename(rdir)] = [dict(d_edge=k, n=int((okq & (np.round(de) == k)).sum()),
                                             median=float(np.median(gnd[okq & (np.round(de) == k)])) if (okq & (np.round(de) == k)).any() else None)
                                        for k in range(1, 7)] + [dict(d_edge="interior (>8)", n=int((okq & (de > 8)).sum()),
                                                                      median=float(np.median(gnd[okq & (de > 8)])) if (okq & (de > 8)).any() else None)]
    res["gnd_vs_d_edge"] = prof
    for k, v in prof.items(): print(f"GND5 median by d_edge, {k}:", [(x["d_edge"], x["n"], None if x["median"] is None else f"{x['median']:.2e}") for x in v])
    err = []                                             # this run's own error (truth zero) by distance, and against the planted fraction
    for k in list(range(1, 7)) + ["interior (>8)"]:
        s = ok & np.isin(grain, grains) & ((de > 8) if isinstance(k, str) else (np.round(de) == k))
        if s.any(): err.append(dict(d_edge=k, n=int(s.sum()), f_median=float(np.median(f[s])), rot_median_deg=float(np.median(m["rot_deg"][s])),
                                    rot_p90_deg=float(np.percentile(m["rot_deg"][s], 90)), strain_norm_median=float(np.median(m["strain_norm"][s])),
                                    strain_norm_p90=float(np.percentile(m["strain_norm"][s], 90))))
    res["error_vs_d_edge"] = err
    print("error by d_edge (d_edge, n, planted f median, rot median / p90 deg, strain norm median / p90):")
    for x in err: print(f"  {x['d_edge']!s:>14} {x['n']:4d} {x['f_median']:.3f}  {x['rot_median_deg']:.4f} / {x['rot_p90_deg']:.4f}  {x['strain_norm_median']:.2e} / {x['strain_norm_p90']:.2e}")
json.dump(res, open(os.path.join(a.out, "mix_bias.json"), "w"), indent=1)
