"""Plot the maps saved by hrebsd.py (maps.npz) - everything is stored, so any figure can be redone without rerunning OpenXY.
usage: replot_maps.py work\hrebsd_NAME\maps.npz [--out DIR] [--clim 2e-3] [--rot-clim 0.5] [--show-unusable] [--cmap RdBu_r]
Writes strain_components.png (e11 e12 e13 e22 e23 e33, deviatoric, OpenXY frame), rotation_components.png (rotation vector, deg),
overview.png (strain norm, rotation angle, grain id, usable mask). Unusable points are grey unless --show-unusable.
Load in your own scripts:  d = np.load('maps.npz'); d['e11'] is (nrows, ncols); keys: see d.files."""
import argparse, os, numpy as np

def replot(npz, out=None, clim=None, rot_clim=None, show_unusable=False, cmap="RdBu_r"):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    d = np.load(npz); out = out or os.path.dirname(npz); u = d["usable"].astype(bool)
    ext = [d["x_um"].min(), d["x_um"].max(), d["y_um"].max(), d["y_um"].min()] if "x_um" in d.files else None
    mk = (lambda v: v) if show_unusable else (lambda v: np.where(u, v, np.nan))
    def grid(keys, titles, scale, lim, fname, cm, sup):
        fig, ax = plt.subplots(2 if len(keys) > 3 else 1, 3, figsize=(13, 7 if len(keys) > 3 else 3.8), squeeze=False)
        for A, k, t in zip(ax.ravel(), keys, titles):
            v = mk(d[k] * scale); L = lim or np.nanpercentile(np.abs(v), 98) or 1
            c = plt.get_cmap(cm).copy(); c.set_bad("0.8"); im = A.imshow(v, cmap=c, vmin=-L, vmax=L, extent=ext, interpolation="nearest")
            A.set_title(t, fontsize=10); plt.colorbar(im, ax=A, fraction=0.04)
        fig.suptitle(sup, fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(out, fname), dpi=130); plt.close(fig)
    grid(["e11", "e12", "e13", "e22", "e23", "e33"], ["e11", "e12", "e13", "e22", "e23", "e33"], 1e3, clim and clim * 1e3, "strain_components.png", cmap,
         "deviatoric strain (x1e-3), OpenXY frame, relative to each grain's reference; grey = flagged")
    grid(["rotx", "roty", "rotz"], ["rotation x (deg)", "rotation y (deg)", "rotation z (deg)"], 1.0, rot_clim, "rotation_components.png", cmap, "rotation vector vs reference (deg)")
    if "gnd5" in d.files:                                    # GND maps (hrebsd.py --gnd): OpenXY estimate from 3 / 5 / 9 Nye components, log10
        ok = u & ~d["gnd_filtered"].astype(bool) & (d["gnd5"] > 0); fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        for A, k in zip(ax, ("gnd3", "gnd5", "gnd9")):
            v = np.where(ok if not show_unusable else (d[k] > 0), np.log10(np.where(d[k] > 0, d[k], 1)), np.nan); lo, hi = np.nanpercentile(v, [2, 98]) if np.isfinite(v).any() else (0, 1)
            c = plt.get_cmap("viridis").copy(); c.set_bad("0.8"); plt.colorbar(A.imshow(v, cmap=c, vmin=lo, vmax=hi, extent=ext, interpolation="nearest"), ax=A, fraction=0.04); A.set_title("log10 " + k, fontsize=10)
        fig.suptitle("GND estimate (OpenXY, units as OpenXY); grey = flagged or filtered (grain/misorientation)", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(out, "gnd.png"), dpi=130); plt.close(fig)
    if "gnd_bin3" in d.files:                                # binned GND (tools/gnd_binned.py), own 3-component formula, log10 1/m^2
        fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
        for A, k in zip(ax, ("gnd_bin2", "gnd_bin3", "gnd_bin4")):
            v = np.log10(np.where(d[k] > 0, d[k], np.nan)); lo, hi = np.nanpercentile(v, [2, 98]) if np.isfinite(v).any() else (0, 1)
            c = plt.get_cmap("viridis").copy(); c.set_bad("0.8"); plt.colorbar(A.imshow(v, cmap=c, vmin=lo, vmax=hi, extent=ext, interpolation="nearest"), ax=A, fraction=0.04); A.set_title("log10 %s (%dx%d bins)" % (k, int(k[-1]), int(k[-1])), fontsize=10)
        fig.suptitle("binned GND (grey = unusable / grain edge); noise floor falls ~1/k^2, see NOTES 13.3", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(out, "gnd_binned.png"), dpi=130); plt.close(fig)
    fig, ax = plt.subplots(1, 4, figsize=(16, 3.8))
    for A, v, t, cm in ((ax[0], mk(d["strain_norm"] * 1e3), "|deviatoric strain| (1e-3)", "viridis"), (ax[1], mk(np.minimum(d["rot_deg"], 2)), "rotation angle (deg, clipped 2)", "magma"),
                        (ax[2], d["grain"].astype(float), "grain id", "tab20"), (ax[3], d["usable"].astype(float), "usable (1) / flagged (0)", "gray")):
        plt.colorbar(A.imshow(v, cmap=cm, extent=ext, interpolation="nearest"), ax=A, fraction=0.04); A.set_title(t, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(out, "overview.png"), dpi=130); plt.close(fig)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("npz"); ap.add_argument("--out"); ap.add_argument("--clim", type=float); ap.add_argument("--rot-clim", type=float)
    ap.add_argument("--show-unusable", action="store_true"); ap.add_argument("--cmap", default="RdBu_r"); a = ap.parse_args()
    replot(a.npz, a.out, a.clim, a.rot_clim, a.show_unusable, a.cmap); print("plots written next to", a.npz)
