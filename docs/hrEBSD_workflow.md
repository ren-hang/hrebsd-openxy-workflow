# hrEBSD with OpenXY on Pathfinder maps: step-by-step workflow

Scope: a scan from the Thermo Fisher Pathfinder (backscatter geometry, 70° sample tilt, 20 kV) stored as `<name>.pat` + `<name>.ANG`, analysed with OpenXY (MATLAB) to get per-point deformation gradient, deviatoric strain, rotation and (optionally) dislocation-density maps. Everything below is run from the repository root `C:\Users\Hang\Projects\hrEBSD` and reproduces the maps in `work\`. The full record of how each choice was reached is `NOTES.md`; this file only says what to do.

**Rules that apply to every step.** Box is read-only: copy the `*.MAP.XTAL` folder to `inputs\` first and work only on the copy. Strain values are produced, not interpreted, by this workflow; see section 7 for what the numbers do and do not mean.

---

## 1. One-time setup (already done on this PC)

| Item | Where | Notes |
|---|---|---|
| MATLAB R2024a | `C:\Program Files\MATLAB\R2024a` | `--nproc 3` uses `parfor`, so the Parallel Computing Toolbox is needed. |
| OpenXY source | `work\OpenXY_BYU\Code` | BYU-MicrostructureOfMaterials/OpenXY `master` cloned 2026-05-07 (+ quaternions submodule) **plus two local patches** in `tools\openxy_patches\` (patched file next to its `.orig`). Re-apply them if OpenXY is ever re-cloned. |
| Python env | `work\kpenv` | `python -m venv work\kpenv`, then `work\kpenv\Scripts\python -m pip install -r docs\requirements.txt`. Pinned versions in that file. |
| Orientation calibration | `inputs\Deba_calibration\EBSD_Rcorr_calibration_20260505.mat` | Deba's 11.04° camera-tilt correction `Rcorr`. |

The two OpenXY patches (both only matter for headless runs):
- `ValidatePhase.m`: if the environment variable `OPENXY_HEADLESS_MATERIAL` is set, unknown phase names (the grain file's `NotIndexed`, the .ANG's `Generic FCC`) are mapped to that material instead of opening a dialog.
- `HREBSDMain.m`: `ppm = [];` so the `parfor` body works when `DisplayGUI` is false.

---

## 2. What you need for each map

1. The map folder from the instrument (`<name>.MAP.XTAL`) copied to `inputs\`. It must contain `<name>.pat` and `<name>.ANG`. Check that the `.pat` has finished downloading from Box; the tool refuses a truncated file.
2. **The pattern centre X, Y, L in Pathfinder millimetres for that session.** It is typed by hand into the acquisition report (the .ANG header stores 0). The instrument default **15.8 / 7.9 / 16.0 means "not entered" and is never correct.** Confirmed values so far:

   | Map | X | Y | L |
   |---|---|---|---|
   | Pt(2), 25 Apr 2026 | 16.22 | 5.88 | 20.35 |
   | ORR(2), 18 Mar 2026 | 16.80 | 6.41 | 20.50 |
   | Pt(6), 9 May 2026 | 16.31 | 6.86 | 19.09 |
   | Pt(7), 9 May 2026 | 15.84 (typed) | 6.96 (typed) | 18.70 (typed): **does not fit the patterns**, pattern fit gives 17.60 / 6.58 / 19.18. Recheck the report. |

3. The material is cubic (the scripts assume an FCC metal, Pt by default). For another material pass `--material <OpenXY materials-folder name>` and `--phase <name>`.

Nothing else needs to be typed: pattern size, scan size, step, pixel size (`mperpix`), indexing, phase and the camera geometry are read or derived (section 4).

---

## 2b. Data organisation and file naming

Folder layout (the tool creates everything under `work\` and `logs\`; only `inputs\` is filled by hand). None of these folders is in git; back the raw maps up on Box.

```
hrEBSD\                          repository (tools\, docs\, NOTES.md, run_gui.bat)
  inputs\
    <date>_<sample>\<name>.MAP.XTAL\   COPY of the instrument folder: <name>.pat + <name>.ANG (never work inside Box)
    Deba_calibration\EBSD_Rcorr_calibration_<date>.mat
  work\
    OpenXY_BYU\Code\             patched OpenXY
    kpenv\                       Python environment
    hrebsd_<RUN>\                one folder per run (below)
  logs\hrebsd_<RUN>.log          MATLAB log of each run
```

One run folder `work\hrebsd_<RUN>\` (written by `prepare`/`run`/`post`):

```
<RUN>.ang                  orientation-corrected scan file (10 columns), x,y re-zeroed
<RUN>.txt                  grain file (11 columns), same row order as the .ang
<RUN>_index.txt            original scan index, instrument CI, original column of each row
first_image.txt            name of the pattern OpenXY opens first
Corrected_EBSPs\<RUN>_r<row>c<col>.tif     one 16-bit tif per scan point
out\AnalysisParams_hrebsd_<RUN>.mat, out\Corr_*.ang   OpenXY results
maps.npz, results.csv, summary.md, *.png   results of `post`
```

**What OpenXY requires, and what the tool already produces (so nothing has to be renamed by hand):**
- **Pattern files:** `<Prefix>_r<row>c<col>.tif` in a folder called `Corrected_EBSPs` next to the `.ang`. `r` = scan row (y), `c` = scan column (x), both **0-based**, e.g. `ORR2_r0c0.tif`, `ORR2_r41c63.tif`. OpenXY is pointed at the first pattern (`..._r0c0.tif`) and builds every other file name from the same prefix and the row/column of each .ang entry, so every point of the scan must exist, with no gaps. The prefix is the run name. 16-bit tif, one scan point each.
- **.ang:** same stem as the run (`<RUN>.ang`), TSL header with `NCOLS`, `NROWS`, `XSTEP`, `YSTEP`, columns phi1 PHI phi2 (radians), x, y, IQ, CI, fit, … in **row-major order, x fastest**, i.e. point i = row·ncols + col. The tool shifts x,y so the first point is (0,0); it sets CI/IQ/Fit near boundaries (reference trap) and writes Euler angles corrected by Rcorr and the azimuth.
- **Grain file:** `<RUN>.txt`, same stem, 11 columns (phi1 PHI phi2 x y IQ CI Fit grain-id, 1, phase name), same row order as the .ang. The phase name must be a material in OpenXY's `Materials` folder; the tool writes `Platinum` (`--phase`) and maps it to `platinum` (`--material`).
- Not-indexed points: their orientation is NaN in the instrument file; they stay in the files (OpenXY needs the complete grid) with grain 0 and CI 0, and are flagged `not_indexed` in the results.

**Compared with the older manual route** (ThermoPAT extraction to tifs named `X=…;Y=….tif`, `FileRenaming.m`, the ImageJ background macro, `analysisScript.m`, then the OpenXY GUI): `hrebsd.py` reads the patterns straight from the `.pat` file, names them as above, corrects and re-zeroes the .ang and writes the grain file itself. It does **not** apply the ImageJ background correction; patterns go to OpenXY as stored and OpenXY applies its own 9/90 pixel filter. If a different pre-processing is ever wanted, put the processed tifs in `Corrected_EBSPs\` under the same names before the `run` stage (`--stage run`).

## 3. The one command

```
work\kpenv\Scripts\python.exe -I tools\hrebsd.py "inputs\<folder>\<name>.MAP.XTAL" --pc X Y L --name NAME [--gnd]
```

It runs four stages (`--stage check|prepare|run|post|all`, default `all`) and writes everything to `work\hrebsd_NAME\`. Typical times: 28 s for a 10×10 block (`--block 10`), 2–3 min for a 2688-point map with 928×580 or 372×232 patterns, about 6 min with `--gnd` on the 928×580 maps. The scan size does not matter beyond that.

Useful options: `--block N` (quick test on an N×N block inside the largest grain), `--gnd`, `--ref-margin 4` (default), `--nproc 3` (default), `--az 90`, `--elev 11.04` (both defaults, section 5).

Example, ORR(2), whole map, with GND (verified end to end, 1 m 43 s):

```
work\kpenv\Scripts\python.exe -I tools\hrebsd.py "inputs\ORR_20260318\ORR(2).MAP.XTAL" --pc 16.80 6.41 20.50 --name ORR2 --gnd
```

---

### GUI
Double-click `run_gui.bat` (or `work\kpenv\Scripts\python.exe tools\hrebsd_gui.py`). Pick the map folder, enter the run name and the PC X, Y, L, optionally tick GND or give a block size for a quick test, then **Check only** or **Run everything**. The window shows the same output as the command line; **Open results folder** opens `work\hrebsd_<name>`. It only builds and runs the command above.

## 4. What each stage does and why

### Stage 1: `check` (stops with exit code 1 if anything below is wrong)
- finds `.pat` and `.ANG`; compares the `.pat` size with its header (truncated copy);
- compares the grid in the .pat header with the .ANG header and row count;
- counts not-indexed points (NaN orientations) and points with confidence index ≤ 0;
- reports the .ANG phase (always "Generic FCC", a = 4.0 Å, which is not the real material; it is replaced for OpenXY);
- rejects a missing PC, the default PC, or a PC outside the plausible range; prints the OpenXY (TSL) x\*, y\*, z\* it will use. Conversion (Brady App. E): x\* = 0.5 + 0.5·(X − 15.8)/20, y\* = 1 − Y/20, z\* = L/20.

### Stage 2: `prepare` (`make_pat_subset.py`, then `correct_ang.py`)
- each pattern is read from the `.pat` (20-byte header of five uint32; pattern k is uint16 little-endian at byte 25 + k·(w·h·2 + 8)) and written as a 16-bit tif named `<NAME>_r<row>c<col>.tif`, which is what OpenXY expects. OpenXY crops the central h × h square of each (580 × 580 or 232 × 232), and the PC refers to that square;
- grains are labelled by 5° flood fill on cubic misorientation (not-indexed points get grain 0) and written to the 11-column grain file OpenXY reads, with the phase column set to the OpenXY material name;
- **reference selection fix**: OpenXY chooses each grain's reference by voting on max CI, min Fit and max IQ. Points within `--ref-margin` pixels (default 4) of a grain boundary get IQ = CI = 0, Fit = 99 in the working files, so a boundary pattern is never picked as the reference;
- orientations are corrected: R′ = Rz(90°) · Rcorr(11.04°) · R. The 11.04° is Deba's calibration of the instrument's camera tilt; the 90° azimuth was found empirically (Hough indexing agrees with the .ANG at 90° and not at 0°) and holds for every 2026 map tested. The matching `CameraElevation = 11.04°` is set in OpenXY. Orientation and elevation must both be corrected;
- `first_image.txt` names the pattern OpenXY opens first; `<NAME>_index.txt` records for every point its original scan index, instrument CI and original column, so results can always be tied back to the instrument file.

### Stage 3: `run` (`run_case.ps1` → `run_openxy_headless.m`)
Starts MATLAB in batch mode. Settings: 20 kV, 70° sample tilt, camera elevation 11.04°, pixel size from the pattern height (`mperpix` = 93542.4 µm / pattern height, i.e. 161.28 for 580 px, 80.64 for 1160 px, 403.2 for 232 px), PC from `--pc`, material from `--material`, one reference pattern per grain chosen with `GrainRefImageType = 'IQ > Fit > CI'`, default image filter 9/90 px, three workers; GND off unless `--gnd` (then `CalcDerivatives` with the `Full` method). The MATLAB log is `logs\hrebsd_NAME.log`. The run is a success only if the log contains `HREBSDMain returned` and no `TIMEOUT`; the tool checks that (MATLAB's own exit code is not reliable here).

### Stage 4: `post`
Reads OpenXY's `out\AnalysisParams_*.mat`, polar-decomposes each F = R·U, and writes:

| File | Content |
|---|---|
| `maps.npz` | **everything, on the scan grid**, for replotting (keys below) |
| `results.csv` | one row per point: index, row, col, original index, grain, instrument CI, flags, usable, rotation, e11 … e33 |
| `summary.md` | counts and the flags |
| `overview.png`, `strain_components.png`, `rotation_components.png` | maps |
| `gnd.png`, `gnd_binned.png` | with `--gnd` |

`maps.npz` keys: `nrows ncols row col x_um y_um original_index grain instrument_CI`; flags `usable col0 not_indexed near_boundary finite`; `rot_deg rotx roty rotz strain_norm`; `F` and `dev_strain` (rows × cols × 3 × 3) and `e11 e12 e13 e22 e23 e33` (deviatoric strain components); with `--gnd`: `gnd3 gnd5 gnd9 nye_alpha gnd_filtered gnd_misorientation gnd_bin2 gnd_bin3 gnd_bin4`.

**Flags.** A point is `usable` if it is finite, not scan column 0, indexed (instrument CI > 0) and max |deviatoric component| ≤ 5×10⁻³. Column 0 (the first point of every scan line) is dropped because those patterns are anomalous in all 2026 maps (cause unknown; dropping is accepted). `near_boundary` points are kept in the table but flagged. Unusable points are grey in the plots.

**Replotting any time**, with no rerun:
```
work\kpenv\Scripts\python.exe -I -c "import sys; sys.path.insert(0,'tools'); from replot_maps import replot; replot(r'work\hrebsd_NAME\maps.npz', clim=2e-3)"
```
Arguments: `clim` (strain colour limit), `rot_clim`, `show_unusable`, `cmap`.

---

## 5. Settings that are decided, and why

| Setting | Value | Why |
|---|---|---|
| Orientation correction | Rz(90°)·Rcorr(11.04°) | instrument tilt calibration + empirical azimuth; without it OpenXY sees the wrong crystal frame (≈1.8×10⁻³ error on the ORR block) |
| CameraElevation | 11.04° | same physical effect as Rcorr |
| Image filter | OpenXY default 9/90 | recovers ~95% of an imposed 5×10⁻³ strain vs ~98% with no filter, with much better robustness; keep default |
| Reference | OpenXY per-grain, kept ≥ 4 px from boundaries | boundary references give garbage |
| Column 0 | dropped | anomalous patterns, cause unknown |
| Remap second pass | **off** | at 0.01–0.05° rotations it adds noise (NOTES 13.5); `remap_full.py --min-rot 0.5` exists for strongly rotated data |
| GND | binned maps for interpretation | native-step GND is at the noise floor (section 6) |

---

## 6. Noise floors (what the maps can resolve)

Measured on rendered strain-free maps (every pattern = one real pattern + measured pixel noise; truth = 0) run through the same pipeline:

| Quantity | Floor (median) |
|---|---|
| Deviatoric strain norm | ≈ 1×10⁻⁴ |
| Rotation vs reference | ≈ 0.005° |
| OpenXY GND (`gnd5`, 0.27 µm step) | ≈ 5–6×10¹² m⁻² |
| Binned GND, 3×3 / 4×4 (0.27 µm step) | ≈ 4–5×10¹¹ / 2×10¹¹ m⁻² |

GND units are m⁻² (OpenXY source: step in metres, Burgers vector from the material file, α = curl β / b). The floor scales as 1/step, and falls ≈ 1/k² when k×k points are binned. At a 0.27 µm step the native-resolution GND map is mostly noise; use `gnd_bin3`/`gnd_bin4`. Real values from the Pt(6)/Pt(7) maps exceed the binned floor by 3–15×.

Whole-map results through this tool: Pt(6) 88.1% usable, median strain norm 0.31×10⁻³; Pt(7) 93.0%, 0.41×10⁻³; ORR(2) 95.3%, 0.57×10⁻³ (median max-component). OpenXY recovers an imposed pure strain with a ratio of 1.001 on exactly rendered synthetic patterns, and under-recovers rotation by 3–7%.

---

## 7. Cautions: what the numbers mean

- OpenXY gives the deformation gradient F = R·U **relative to a reference pattern inside each grain**, not absolute strain. Strain is only comparable within a grain, and any residual strain already in the reference is invisible.
- Only the deviatoric part of the strain is measured. The dilatation printed by OpenXY comes from the assumption σ33 = 0.
- Rotation relative to the reference matters for accuracy: a pattern-centre error leaks into strain in proportion to that rotation (≈ rotation [rad] × δPC / z\*). At the 0.01–0.05° rotations in these maps the effect is below 0.1×10⁻³ (Pt(7) typed vs fitted PC: median 0.04×10⁻³).
- Points within a few pixels of a grain boundary (overlapping patterns) are unreliable and are flagged; the tool does not try to separate the two patterns.
- The 90° azimuth is empirical. If a new session gives a Hough or kikuchipy check that disagrees with the .ANG, rerun `tools\pat_az_scan.py` (section 8) before trusting any map.

---

## 8. Optional checks (use when something looks off)

| Task | Command / tool |
|---|---|
| Does the PC in the report fit the patterns? | `tools\make_pt_master.py` (once; Pt master pattern with `ebsdsim`, saved as `work\pt_master\Pt_20kV_hw250.npz`), then `tools\map_pc_fit.py --pat P.pat --ang orig.ANG --corr work\hrebsd_NAME\NAME.ang --grain work\hrebsd_NAME\NAME.txt --pc0 x* y* z*` (fits the PC to patterns spread over the map's grains with kikuchipy; prints the Pathfinder-unit equivalent). Fit and typed PC agreeing within ~0.02 in x\*, y\*, z\* is fine. |
| Is the azimuth right for a new session? | `tools\pat_az_scan.py --pat P.pat --ang A.ang --rmat inputs\Deba_calibration\EBSD_Rcorr_calibration_20260505.mat --pc X Y L` (indexes a few patterns with PyEBSDIndex at several azimuths; the right one gives a median misorientation of ~1–2° to the .ANG). |
| How different are two runs? | `tools\compare_maps.py runA\maps.npz runB\maps.npz` |
| What is the noise floor of this map? | `tools\make_null_map.py --src work\hrebsd_NAME --stem NAME --out work\hrebsd_NAMEnull --name NAMEnull`, then run `hrebsd.py` on that name (`--stage run`, then `--stage post`) and compare to the real run. |
| GND at coarser effective step | already in `maps.npz` (`gnd_bin2/3/4`); or `tools\gnd_binned.py maps.npz --k 1 2 3 4` for medians per grain. |
| Large rotations / deformed samples | `tools\remap_full.py` (second pass, F = R1·F2). |
| Does OpenXY itself recover a known deformation? | `tools\gen_synthetic.py` renders a 10×10 map with known F from a kikuchipy master pattern; `tools\compare_synthetic.py` compares OpenXY's output with the truth. |

---

## 9. When it fails

| Symptom | Cause | Fix |
|---|---|---|
| `check` says PC is the default | report not filled in | get the real X, Y, L |
| `check` says `.pat` is shorter than its header | incomplete copy | wait for Box, copy again |
| MATLAB error "Material detection failed … notindexed" | unpatched `ValidatePhase.m` | apply `tools\openxy_patches\ValidatePhase.m` |
| MATLAB error "Unrecognized function or variable 'ppm'" | unpatched `HREBSDMain.m` | apply `tools\openxy_patches\HREBSDMain.m` |
| Run takes far longer than ~3 min per 2700 points | MATLAB serial or parpool failing | check `--nproc`, and `logs\hrebsd_NAME.log` |
| Huge strain in a grain with a large rotation or garbage along a boundary | bad reference | rerun with a larger `--ref-margin` |
| Strain looks odd everywhere | wrong PC, azimuth or elevation | section 8 checks |

Removed development scripts (sweeps, one-off warp tests, kinematic PC fit, the Hough-index scan variants, block-only remap) are archived in `work\removed_tools_2026-11-05.zip`; `NOTES.md` still describes what they showed.
