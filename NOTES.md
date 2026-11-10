# OpenXY / hrEBSD on this PC — working notes

Status: steps 1-4 done and step 6 (synthetic validation, §9) done (2026-11-03). No strain values from the lab data have been interpreted. Box is read-only: nothing was written into it.

## 0. Layout

| Path | What |
|---|---|
| `Examples\` | Original files you supplied (protocol .docx, Brady PDF, lab .m/.ijm scripts). Untouched. |
| `inputs\` | Verbatim copies from Box (see §1). Untouched; treat as read-only. |
| `work\OpenXY\Code` | Local OpenXY copy used for runs (lineage `Deba_OpenXY_05_01_2024\OpenXY-master\OpenXY-master`, the one the protocol names). Two patched files, each with a `.orig` next to it (§5). |
| `work\OpenXY_DebaOpenXYfolder_Code` | Unpatched copy of the other lineage (`Deba_OpenXY_05_01_2024\OpenXY\Code`), only used to reproduce failure F2. |
| `work\pt10x10\` | 10x10 practice subset: `Pt10x10.ang`, `Pt10x10.txt` (grain file), `Corrected_EBSPs\` (100 tifs), `bkg1.tif`, `out\` (OpenXY output). |
| `tools\` | `preflight.py`, `make_subset.py`, `run_openxy_headless.m`, `run_case.ps1`, `headless_shims\` (dialog-to-error shims). |
| `logs\` | `preflight_*.log`, `run0N_*.log`, `matlab_ver.log`. |

## 1. What was copied from Box (`C:\Users\Hang\Box\Ren Lab Shared`)

| Local (`inputs\`) | Box source | Size | Note |
|---|---|---|---|
| `HREBSD_Strain_files_06_02_2024\OpenXY-master`, `OpenXY_updated_05_07_2024`, `Sample_Data`, `analysisScript.m`, `FileRenaming.m`, `Kikuchi_Correction_Method2.ijm.txt`, `IMG-20240524-WA0011.jpg` | `Results and Data\Data\Deba\HREBSD_Strain_files_06_02_2024\` | ~189 MB | |
| `...\Practice_05_02_2024_Strain-analysis\` (4 small files + two 100-pattern subsets) | same, `Practice_05_02_2024_Strain-analysis` | 4.6 MB + 2 x 411 MB | Full folder is **132.7 GB** (3 copies of 10,752 patterns x ~4.1 MB: `Corrected_EBSPs`, `EBSPs - Copy`, `Renames-EBSPs`). Only patterns r0-9,c0-9 of `Renames-EBSPs` (raw) and `Corrected_EBSPs` were copied. |
| `Deba_OpenXY_05_01_2024\` (`OpenXY`, `OpenXY-master`, `Files`) | `1. Group Members\G7. Deba\OpenXY_05_01_2024\` | ~243 MB | One file failed on first copy (Box not yet hydrated), recopied OK. |
| `Deba_EBSD\Open_XY_Protocol.pdf`, `Open_XY_data.pptx` | `...\G7. Deba\EBSD\` and `...\EBSD\2026\05-09-2025\Pt_map1_ORR\` | 1.5 MB | pptx not opened yet (earlier OpenXY output to compare against). |
| `Alumni_Munoz-Castaneda_OpenXY_REFERENCE\` | `1. Group Members\Z99. Alumni\U1. Jesus Munoz-Castaneda\3. Protocols and code\OpenXY` | ~17 MB | Reference only. |

Missing vs. your list: the `HREBSD_Strain_files_06_02_2024` folder has **no `Code` or `Documentation` subfolder** and no `Run_OpenXY.m` / `Batch_Analyze_Run_OpenXY.m` at its top level (they live inside the OpenXY trees, `...\Code\`). The practice folder is `Practice_05_02_2024_Strain-analysis`, not `practice_05_02_2024`. I also copied `Sample_Data`, `analysisScript.m`, `FileRenaming.m` and the `.jpg`/`.ijm.txt` that sit beside the OpenXY folders, which you did not list.

## 2. The pipeline (protocol + Brady Appendix E)

1. **Acquire** (ThermoFisher, EER PC): 20 kV, 26 nA, 70 deg tilt, 1x1 binning, zero gain. Save background tif; tick "Save Patterns to Disc"; record pattern centre (X, Y, L) from Pathfinder (View > Attributes > EBSD); convert .ebsd -> .ang with EBSDFileConverter; upload `.ang` + `.pat`.
2. **Extract** patterns: ThermoPATextractV2, tick "Force 8-bit EBSP to 16-bit Word" -> folder `EBSPs\` with names containing `X=..;Y=...tif`.
3. **Rename** with `FileRenaming.m` -> `Platinum_r<Y>c<X>.tif` (r = row = y, c = column = x), into a new folder.
4. **Correct + background-subtract** in Fiji with the `.ijm`: bkg to 16-bit + auto B&C (manual, saved); macro auto-contrasts each pattern, subtracts bkg (32-bit), converts to 16-bit, auto-contrasts, saves into `Corrected_EBSPs`.
5. **Grain file**: `analysisScript_share_3_modified_chatgpt.m` (MTEX) loads the .ang, `calcGrains` at 5 deg, calls `export_GrainFile` -> `Platinum.txt` (11 columns; must have the same number of rows as patterns; x steps first, then y).
6. **OpenXY** (`MainGUI`): scan file = original `.ang` (lowercase ext); grain file; material; first image = first corrected pattern; output file. Settings (Brady E): Square scan, Material "Scan File", 3 processors; Microscope 20 kV / 70 deg tilt / sample az 0 / camera elev 0 / camera az 0 / 80.64 um per pixel; Advanced: strain on, Real-Grain ref, Std Dev 2, Misorientation tol 5 deg, IQ>Fit>CI, Grain ID = Grain File, GND on (Full cross-correlation), skip 0, IQ cutoff 0; ROI 25 %, 48 ROIs, Annular, ROI filter 2/50/1/1, image filter standard 9/90/0/0; PC from Pathfinder via E.1-E.3 (below), optional strain-minimisation PC refinement, check with Test Geometry, Run.

PC conversion (Brady E.1-E.3; the PDF extraction of the equations is garbled, so **verify against the dissertation page 110 image**):
`x* = 0.5 + 0.5*((X_TF - 15.8)/20)`, `y* = 1 - Y_TF/20`, `z* = L_TF/20`.
For the protocol PDF's practice values X=16.31, Y=6.86, L=19.09 -> (0.51275, 0.65700, 0.95450). The `.docx` omits these practice values (only the older PDF has them).

## 3. Hardcoded paths and settings (lab PC used `C:\Users\renlab\...`)

| File | Line | Hardcoded | Local handling |
|---|---|---|---|
| `FileRenaming.m` | 8, 9 | `C:\Users\renlab\Desktop\Deba\05-09-2026_ORR_HREBSD\Pt_map1_ORR\EBSPs` / `...\Renamed_EBSPs` | not needed for the practice subset (already renamed) |
| `FileRenaming.m` | 18 | specimen name `'Platinum_r%sc%s.tif'` (protocol says "change in the for loop"); assumes source names contain `X=...;` and `Y=....tif`; `movefile` **moves**, doesn't copy | |
| `analysisScript_...chatgpt.m` | 14, 15 | `C:\Users\renlab\Desktop\Deba\05-09-2026_ORR_HREBSD\Pt_map1_ORR` and `'\Pt(corrected).ANG'` — **uppercase ANG** (protocol says lowercase) | |
| same | 8 | `crystalSymmetry('cubic',[3 3 3],'Generic FCC')` — lattice constant 3 (Pt is 3.92); phase 1 = notIndexed | |
| same | 391-395 | `phaseId==2 -> 'Platinum'` else `'NotIndexed'` (hardcoded) | source of the `NotIndexed` phase name, see F1 |
| same | 377, 400, many `print(...)` | output names `N205DF_Area10_Compiled_OpenXY_...` (from a Nb dataset), `Platinum.txt`; written to the current folder, `export_ang` would overwrite earlier files | |
| same | 72, 107, 122 | `min_g_size = 8` while `export_GrainFile` header says "Minimum Grain Size: 5" | |
| `export_GrainFile.m` | 62-76 | writes a `#` header, yet `PtGrain3.txt` has **no header** (protocol: "header isn't included in the current code") — the provided file is not from the current exporter, or the header was stripped | |
| `Kikuchi_Correction_Method_updated_10_02_2024.ijm` | 4, 5, 6, 11 | `C:/Users/renlab/Desktop/Deba/05-09-2026_ORR_HREBSD/Pt_map1_ORR/Bkg_HR-EBSD.tif` (twice) and two directory strings passed to `getDirectory()` (only used as dialog titles) | protocol says change "line 12"; the bkg is actually opened at lines 4 **and** 11 |
| `Batch_Analyze_Run_OpenXY.m` (OpenXY) | 12, 19, 58 | `load Batch_Analyze_Settings`, material `tantalum`, `*.ctf`, GUI folder picker (`Batch_Analyze_uipickfiles`), first image `0_0.*` | not used; headless driver replaces it |
| Protocol docx | 5, 24, 40 | `C:\Users\renlab\Box\...\EmilyBrady_Dissertation_20211105`, `C:\Users\renlab\Desktop\Deba\Codes_HREBSD_strain\Codes`, `C:\Users\renlab\Desktop\OpenXY-master\Files 05-02-2024\Deba Raw File\OpenXY-master\OpenXY-master\Code` | the last one matches the `OpenXY-master\OpenXY-master` lineage (used here) |

Headless driver paths are all parameters (`tools\run_openxy_headless.m`); nothing is hardcoded to `renlab`.

## 4. Protocol vs. code inconsistencies

1. Protocol/Brady cite `Kikuchi_Correction_Method2.ijm` and `RenameTifFiles.m`; the files are `Kikuchi_Correction_Method_updated_10_02_2024.ijm` and `FileRenaming.m`. Bkg line is "12" in text but lines 4/11 in code.
2. Protocol says lowercase `.ang`; `analysisScript` reads `Pt(corrected).ANG`.
3. OpenXY looks for `<angname>.txt` as the grain file (`GetOIMGrainVals.m`); the script writes `Platinum.txt`, so a manual-select dialog appears (blocks headless runs).
4. The phase name `NotIndexed` written by the script is not in OpenXY's materials list -> "Material detection failed" dialog (protocol says pick Platinum by hand). In code this is `ValidatePhase.m`.
5. `GetOIMGrainVals.m`/`ReadGrainFile.m` read a **11-column** grain file; row count must equal the ang row count (`ImportScanInfo` errors otherwise). Verified OK for the practice data (10,752 = 128 x 84).
6. **Two OpenXY lineages behave differently on the same ang.** `OpenXY-master\...\ReadAngFile.m` derives `NumColsOdd` from unique x; `Deba_OpenXY_05_01_2024\OpenXY\Code\ReadAngFile.m` requires `# NCOLS_ODD:` in the header. `Pt_plasticstrain.ang` has only `# NCOLS:` -> second lineage fails (F2). The MTEX-exported `..._ebsdG_20210817.ang` has `NCOLS_ODD/EVEN`. Four OpenXY trees exist (MainGUI.m sizes 40,959 / 40,971 / 40,971 / 40,973 bytes); which one is canonical is not stated.
7. Settings: protocol "Camera elevation 0 (will need to change)" vs code default 10 deg; "80.64 um/pixel" gives PhosphorSize = 1160 px x 80.64 = **93.5 mm**, vs the code's own default `30000/max(x,y)` = 16.2 um/px (30 mm). Needs checking against the detector.
8. **Non-square patterns**: patterns are 1856x1160; `ReadEBSDImage` crops to a square and `HREBSDMain.m` line 12 has a commented warning "change PC if the image is cropped". Nothing in the protocol says whether the Pathfinder PC (E.1-E.3) refers to the full or the cropped frame. `Batch_Analyze_Run_OpenXY.m` converts tiff-tag PCs with `x = (x_tag - (1-VH)/2)/VH`, `z = z_tag/VH`, which E.1-E.3 does not obviously reproduce. **PC consistency check: unresolved.**
9. `bkg1.tif` is 8-bit; protocol says convert to 16-bit before subtraction. Raw extracted patterns are uint16 with max 4095 (12-bit range), corrected ones span 0-65535.
10. The ang header `x-star/y-star/z-star` is 0 in both practice ang files (and copied unchanged into the OpenXY output ang); the PC only enters through the GUI/driver.
11. Protocol .docx numbering skips step 1 and merges steps; the .pdf version is the same text plus the practice PC.
12. `Batch_Analyze_Run_OpenXY.m` defaults (tantalum, Radial, 20 ROIs, `DoStrain=0`) differ from the protocol settings.

## 5. Preflight (`tools\preflight.py`) results

Run: `python tools\preflight.py --ang work\pt10x10\Pt10x10.ang --grain work\pt10x10\Pt10x10.txt --images work\pt10x10\Corrected_EBSPs --bkg work\pt10x10\bkg1.tif --pc-thermo 16.31 6.86 19.09 --expect-prefix Platinum` -> **0 FAIL, 7 WARN, 20 PASS** (`logs\preflight_subset.log`).
Passing: lowercase `.ang`; rows = 10x10; x-then-y order in both files; grain file 11 columns, same x,y as ang; 100 names `Platinum_r{0-9}c{0-9}.tif`, r0c0 present; all 1856x1160 uint16.
Warnings: no grain header; grain file name != ang name (fixed in the subset by naming both `Pt10x10`); PC 0 in header; non-square images; bkg 8-bit vs 16-bit patterns; `NCOLS` only (lineage-dependent, F2).
On the original full-size practice files (`logs\preflight_originals.log`) the same checks pass for the full 128x84 grid; the only FAIL there is "pattern count" because only 100 patterns were copied.
Not checked mechanically: whether E.1-E.3 PC matches the code's convention (see §4 item 8).

## 6. MATLAB environment

MATLAB R2024a Update 6 (24.1.0.2689473), Windows x64; `matlab -batch` works. MTEX 6.1.1 at `C:\mtex-6.1.1` (auto-initialised by startup). Toolboxes present: Image Processing, Parallel Computing, Statistics and Machine Learning, Signal Processing, Computer Vision, Curve Fitting, Deep Learning, Symbolic Math, PDE (full list: `logs\matlab_ver.log`). Brady used R2019b / MTEX 5.6.0. Parallel pool ('Processes' profile, 3 workers) starts fine.

## 7. Failure log and fixes (each in LOCAL copies only)

Harness (not an OpenXY change): `tools\headless_shims\*.m` replace `questdlg/errordlg/warndlg/msgbox/inputdlg/uigetfile/uiputfile/uigetdir` with errors so a GUI prompt becomes a logged failure instead of hanging MATLAB (the very first attempt hung silently on such a dialog and had to be killed). `tools\run_openxy_headless.m` reproduces `MainGUI` setup (ImportScanInfo, pattern provider, PixelSize/PhosphorSize, PC, `HREBSDMain`) with the Brady settings. `tools\run_case.ps1` runs it with a timeout and captures stdout/stderr.

| # | Run / log | Symptom | Cause | Fix |
|---|---|---|---|---|
| F1 | `run01_repro` | `Dialog errordlg() blocked ... Auto material detection failed. "notindexed" not found in list of known materials` (ValidatePhase <- GetOIMGrainVals <- ReadAngFile) | grain file phase column contains `NotIndexed` (written by analysisScript lines 391-395); not in `Materials\` | **C1**: `work\OpenXY\Code\ValidatePhase.m` — if env var `OPENXY_HEADLESS_MATERIAL` is set, unknown phases are mapped to that material (here `platinum`, with a warning), mirroring the manual "select Platinum" step. Original kept as `ValidatePhase.m.orig`. |
| F2 | `run03_debaOpenXY_unpatched` | `Unrecognized field name "NumColsOdd"` (ReadAngFile line 66), other lineage | `Pt_plasticstrain.ang` header has `NCOLS:` but no `NCOLS_ODD:` | Not patched: the `OpenXY-master` lineage handles it. Alternative would be adding `# NCOLS_ODD:`/`NCOLS_EVEN:` to the ang header. |
| F3 | `run04_parallel3` | `Unrecognized function or variable 'ppm'` in `HREBSDMain` with 3 workers | `parfor` body references `ppm`, only defined when `Settings.DisplayGUI` is true (MainGUI always sets it true; headless sets it false) | **C2**: `work\OpenXY\Code\HREBSDMain.m` line 61 `ppm = [];` (original kept as `.orig`). |

Result after C1 (+C2 for 3 workers): `run02` (serial) and `run05_parallel3_patched` (3 workers) both complete on the 10x10 subset, write `out\Corr_<name>.ang` (100 data rows) and `AnalysisParams_<name>.mat`, and the two `.ang` outputs are **byte-identical** (same SHA-256), so the pipeline is deterministic across worker counts. Values were not inspected.

**The original Feb 2026 Pt failure is still unconfirmed.** It was never recorded; F1 (every Pt grain file written by the current script has `NotIndexed` unless all points index) and F2 (a Pt ang with `NCOLS` only, run on the Deba/OpenXY lineage) are the two I could reproduce from the files in Box, but I cannot tell which, if either, was hit. Ask Deba which OpenXY folder and which ang/grain files were used.

## 8. Open items before step 6

- Resolve the PC convention (full vs cropped frame; E.1-E.3 vs the code) and the 80.64 um/px vs 1160 px phosphor size. §9.3 shows an x/y ambiguity of ~1 % is harmless at small deformations, a z* unit error is not.
- Confirm E.1-E.3 against page 110 of Brady (PDF text is garbled).
- Decide which OpenXY lineage is canonical; `OpenXY_updated_05_07_2024` was not exercised.
- `Open_XY_data.pptx` not yet opened.
- Fiji correction step (step 4 of the pipeline) was not re-run; the supplied `Corrected_EBSPs` were used as-is.
- Only `Pt_plasticstrain.ang`'s first 10x10 block was run; no test of a different scan size, hexagonal grid, or an `.ANG` extension.

## 9. Step 6 - validation with kikuchipy-simulated patterns of known strain / shift

Your step-6 message was cut off after "plus a ..."; this covers only the kikuchipy known-strain/known-shift part. Anything after "plus a" is **not done**.

### 9.1 Method (all in `tools\`, data in `work\syn_A` and `work\syn_C`, reports in `logs\val0*.md`)
- Environment: `work\kpenv` (Python 3.11 venv: kikuchipy 0.13.1, orix 0.15, numpy 2.4, scipy, tifffile). Master pattern: kikuchipy's shipped **Ni, 20 kV** master pattern (FCC; 401x401 Lambert). Pt was not available, so the material in OpenXY is `nickel` (from the grain-file phase), not platinum.
- `gen_synthetic.py` renders a 10x10 "scan" of one orientation (Bunge 20, 35, 10 deg), 1160x1160 uint16 (the size OpenXY sees after cropping), PC (x*,y*,z*) = (0.51275, 0.657, 0.9545), 70 deg tilt, camera elevation 0, written as `Ni_r{Y}c{X}.tif` + `.ang` + grain file. Pixel -> ray uses **OpenXY's own geometry** (`CalcF.m` / `Theoretical_Pixel_Shift.m`, pixel (col,row) 1-based = view coords), intensity from the kikuchipy sampler. A deformation F (crystal frame) is imposed as I'(ray) = I_ref(F^-1 ray), i.e. the same convention as OpenXY's forward model. Point 0 is the reference (highest IQ/CI; OpenXY picked index 1 as reference in every run).
- 100 points: 1 reference, 24 pure strains (e11..e23 at +-1e-3, +-5e-3), 12 rotations (+-0.1, +-0.5 deg about x,y,z), 16 rigid sub-pixel to 5 px shifts (Fourier shift), 12 random mixed (strain ~1.5e-3 + rotation <0.3 deg), 35 noise cases (F=I, Gaussian noise 1-20 % of range). `syn_C` replaces the noise cases with larger misorientations (+-1, 2, 4 deg and 17 random mixes up to 3 deg).
- Exact truth for shifts: a rigid shift is exactly representable by an F under OpenXY's forward model; `shift_truth.m` solves for it at the 48 annular ROI centres (residual 1e-13 px).
- Runs: `run_openxy_headless.m`, the patched local OpenXY (§7), 3 workers, Real/grain reference, 25 % ROI, 48 annular ROIs, strain on, GND off (to keep the check on F). Comparison (`compare_synthetic.py`): polar decomposition F = R U; strain reported as the **deviatoric** part of U-I (OpenXY fixes the dilatation with sigma33=0, so it is not measurable); rotation as rotation vector.

### 9.2 Results with the true pattern centre (`logs\val01_A_truePC.md`)
| group | max abs dev-strain error | max abs rotation error (rad) |
|---|---|---|
| reference pattern | 0 | 0 |
| pure strain, up to 5e-3 | 4.9e-4 (RMS 1.7e-4) | 5.5e-4 |
| rotation <= 0.5 deg | 6.2e-4 (strain leakage) | 6.1e-4 |
| mixed | 4.4e-4 | 5.2e-4 |
| rigid shift 0.1-5 px | 2.4e-4 | 2.7e-4 |
| noise 1 / 3 / 5 / 10 / 20 % | mean max-abs strain 3e-5 / 1e-4 / 2e-4 / 4e-4 / 1.2e-3 | rotation of the same order |

- **Frames/conventions agree**: an imposed crystal-frame F comes back in the crystal frame with the same sign (rotations keep their sign, so F is not returned inverted or transposed); sign-antisymmetric for +-shifts; the shift cases match the OpenXY-forward-model prediction to <= 2.7e-4 (rad / strain), i.e. cross-correlation + fit reproduce a known shift.
- **Accuracy**: on noise-free patterns the error is ~1e-4 to 5e-4, with a small systematic under-recovery (~7-15 %: e.g. e11 = 1e-3 recovered as 0.86e-3 in dev terms, e12 = 5e-3 as 4.63e-3) and 1-4e-4 cross-talk into other components. Noise at 5 % of the intensity range already gives ~2e-4.
- **Rotation leaks into strain, linearly**: with the exact PC, a pure rotation of 0.5 / 1 / 2 / 4 deg about the crystal x axis produces spurious dev-strain of ~6e-4 / 1.4e-3 / 2.8e-3 / 6e-3 and the rotation is under-recovered by ~8 % (`logs\val03_C0_bigrot_table.md`). Likely cause: `Real` reference runs `CalcFShift` once, with no remap/iteration (the iteration loop in `GetDefGradientTensor.m` belongs to the Simulated branch). I did not isolate this further.

### 9.3 Pattern-centre sensitivity (`logs\val02_pc_sensitivity.md` small deformations, `logs\val03_pc_sensitivity_bigrot.md` up to 4 deg)
Same patterns, OpenXY given a wrong PC:
- dx* = +0.010 (11.6 px), dy* = +0.010, dx* = +0.002, half-pixel offset in x and y (4e-4), and the "full-frame fraction used on the cropped frame" case dx* = -0.0085: changes of the errors are **below the 5e-4 floor** for deformations up to 0.5 deg / 5e-3, and <= ~10 % of the baseline rotation error for up to 4 deg.
- dz* = +0.020 (22 px): also within the floor.
- z* scaled by 1160/1856 (a units mix-up between full-frame and cropped width): **dev-strain error rises to 2.9e-3** for the <= 0.5 deg rotations and 1.5e-3 for mixed cases. So an x/y PC ambiguity of about 1 % (the crop question in §4 item 8) is harmless at these deformation levels, but a z* (detector-distance) unit error is not.
- Caveat: this measures OpenXY's response to a PC error, not which PC is physically right for the lab detector; E.1-E.3 and the 80.64 um/pixel phosphor size remain unverified (§8).

### 9.4 What this does and does not show
- It validates OpenXY's cross-correlation, F fit, frame handling and the headless pipeline against ground truth **under OpenXY's own geometry model** and a smooth, low-resolution Ni master pattern (no detector blur beyond the 401-px master, no experimental noise, no pattern preprocessing, no real background subtraction, a single orientation, single grain).
- It does not check that OpenXY's geometry (phosphor frame, tilt, camera elevation, 1-based pixel convention) matches the physical microscope, nor GND outputs, nor the PC-minimisation tools, nor `OpenXY_updated_05_07_2024`, nor non-square frames (the 1856x1160 crop was reasoned about, not simulated). I did not cross-check my renderer against kikuchipy's own `EBSDDetector` projection, so an independent confirmation of OpenXY's frame definitions is still open.

### 9.5 Reproduce
```
work\kpenv\Scripts\python.exe tools\gen_synthetic.py --out work\syn_A            # --profile bigrot --out work\syn_C
matlab -batch "shift_truth('<repo>\work\OpenXY\Code','<repo>\work\syn_A\truth.csv','<repo>\work\syn_A\shift_truth.csv',1160,[0.51275 0.657 0.9545],[20 35 10],70,0)"
set OPENXY_HEADLESS_MATERIAL=platinum, OPENXY_NPROC=3, OPENXY_FIRSTIMAGE=Ni_r0c0.tif, OPENXY_GND=0, OPENXY_PC="0.51275 0.657 0.9545"
powershell -File tools\run_case.ps1 -Log val01_A_truePC -Code work\OpenXY\Code -Data work\syn_A -Stem Syn10x10
work\kpenv\Scripts\python.exe tools\compare_synthetic.py --mat work\syn_A\out\AnalysisParams_val01_A_truePC.mat --truth work\syn_A\truth.csv --shift-truth work\syn_A\shift_truth.csv --label A
```
Each OpenXY run on 100 patterns of 1160x1160 takes ~35 s with 3 workers.

## 10. Camera tilt (~10-11 deg) and uncorrected orientations (added after step 6)

Context: Debasree reported (email 2026-10-01) that the instrument has a ~11 deg error, a camera (detector) tilt as in Britton et al. 2016 frames figure; the Pt .ang files do NOT contain the correction (confirmed by the user).

### 10.1 Sign / convention check (independent of my renderer)
- OpenXY `Qps` uses `alpha = pi/2 - SampleTilt + CameraElevation`. kikuchipy `EBSDDetector(sample_tilt=70, tilt=t).sample_to_detector` times OpenXY `Qps` equals the **same constant matrix diag(-1,1,-1) for t = -10,0,5,10,11** (checked numerically). So kikuchipy detector `tilt` == OpenXY `CameraElevation`, same sign (+ = same direction), differing only by a fixed frame flip.
- **OpenXY's own default is CameraElevation = 10 deg** (`GetHROIMDefaultSettings.m:65`); my headless runner defaulted to 0 (`run_openxy_headless.m`, env `OPENXY_CAMELEV`). All earlier Pt-subset runs therefore used 0 - they were plumbing tests only. Real runs must set the elevation explicitly.

### 10.2 Sensitivity (synthetic Ni, rendered with true elevation 10 deg, true PC; max dev-strain err / max rot err, rad) - `logs/val04_camera_elevation.md`
| analysis elevation | strain cases | rot cases | mixed |
|---|---|---|---|
| 10 (true) | 3.8e-4 / 4.0e-4 | 6.6e-4 / 7.1e-4 | 3.9e-4 / 4.5e-4 |
| 11 | 3.3e-4 / 4.0e-4 | 6.5e-4 / 6.1e-4 | 3.9e-4 / 4.4e-4 |
| 15 | 4.4e-4 / 4.1e-4 | 6.4e-4 / 1.0e-3 | 4.0e-4 / 4.2e-4 |
| 5 | 9.7e-4 / 3.8e-4 | 6.5e-4 / 1.2e-3 | 5.1e-4 / 4.8e-4 |
| 0 (tilt ignored) | 1.6e-3 / 3.7e-4 | 6.4e-4 / 1.7e-3 | 7.2e-4 / 6.4e-4 |
(The "shift" column in that file is meaningless here: no shift truth was computed for elevation 10.)

### 10.3 Uncorrected orientations (g' = g * Rx(+-10 deg) in .ang and grain file, elevation 10) - `logs/val05_orientation_error.md`
- strain 1.4e-3 (+10) / 1.8e-3 (-10) vs 3.8e-4; rot-case rotation error 1.6e-3 / 1.5e-3 vs 7.1e-4; mixed 7.1e-4 / 7.5e-4 vs 3.9e-4.
- Caveat: the way the 10-deg error enters real data (rotation about which axis, in which frame) is my assumption; Debasree's correction code decides this.

### 10.4 Conclusion / limits
- Ignoring a 10 deg tilt, or feeding 10 deg-wrong orientations, costs about 1-2e-3 in strain at imposed strains <=5e-3 and rotations <=0.5 deg (errors ~4x the correct-geometry floor). Being within ~1 deg is harmless; 5 deg off is acceptable for strain, not for rotation.
- Real-data usage: set CameraElevation to the calibrated tilt AND use orientations from the same geometry (or apply Debasree's correction to the .ang, not both). Scripts: `tools/rotate_ang.py` (makes the wrong-orientation variants), `tools/summarize_pc_runs.py`.
- Still using Ni master pattern, OpenXY's own renderer geometry, one orientation, no experimental noise.

## 11. Fresh BYU OpenXY, independent geometry check, PC refinement, preprocessing sweep

### 11.1 Fresh BYU OpenXY (`work/OpenXY_BYU/Code`, clone of github.com/BYU-MicrostructureOfMaterials/OpenXY master, 2026-05-07, + quaternions submodule)
- Unpatched, on the Pt 10x10 subset: still fails with **F1** (unknown phase "notindexed" -> errordlg). `ReadAngFile` derives NumColsOdd from unique x, so **F2 does not occur** in BYU master (F2 is specific to the Deba `OpenXY\Code` tree).
- `ValidatePhase.m` and `HREBSDMain.m` are byte-identical to the Deba-derived copy, so patches C1 and C2 apply unchanged (originals kept as `.orig`).
- With C1+C2 the F array on the Pt subset is **identical (max difference 0)** to run05. Core numerical files (Theoretical_Pixel_Shift, CalcF, CalcFShift, GetDefGradientTensor, CropSquare, GetROIs, ReadEBSDImage, HREBSDPrep, euler2gmat, custimfilt, ...) are identical to the Deba-derived copy ignoring line endings. **Decision: BYU master + C1 + C2 is the working copy.**
- Kacher fork (`Bethany` branch, 30 commits since the common ancestor f65ab08, last code change 2023): adds EMsoft dynamic simulated references (needs separate EMsoft install), XASGO/amoeba, `SwitchF`, `ResolveFandDelta`, ang reader deriving the grid from unique x/y. Not used.

### 11.2 Independent geometry check (`tools/kp_geometry_check.py`)
- Clean reference pattern rendered with OpenXY's geometry vs kikuchipy projection: **kikuchipy 'tsl' convention with the same (x*,y*,z*), Bunge Euler as-is, sample tilt 70, detector tilt = OpenXY CameraElevation (10) gives NCC 0.998.** Equivalent Bruker form (x*, 1-y*, z*) also 0.998. Flipped x (1-x*) gives 0.50; inverted rotation gives <=0.09. So OpenXY's PC convention == kikuchipy/EDAX-TSL convention and my renderer's frames (incl. elevation sign) are confirmed by an independent code.
- Still open: mapping Pathfinder (X,Y,L) -> (x*,y*,z*) (Brady E.1-E.3) is NOT tested by this.

### 11.3 PC refinement (`tools/kp_pc_refine.py`, NCC + Nelder-Mead vs master pattern; optional orientation correction)
- Synthetic (noisy, background-corrected, 5 patterns, one orientation): start (0.55,0.62,0.90) -> (0.5132,0.6557,0.9515); truth (0.51275,0.657,0.9545): errors 5e-4, 1.3e-3, 3e-3. So kikuchipy-based refinement can recover the PC to ~0.003 provided the orientation and geometry are right.
- **Real Pt patterns: no match** (NCC 0.02-0.04, also with orientation refinement and with a +-16 deg orientation / flip grid, best 0.13). Not usable yet: the Ni master pattern stands in for Pt, and the .ang orientations are not in a frame I can reproduce.
- Independent Hough indexing (PyEBSDIndex 0.3.10.1, `tools/pei_index.py`, vendor EDAX, PC from Brady conversion, elevation 10): on synthetic patterns it reproduces the true orientation to <1 deg (convention verified). On the real Pt patterns r0c0/r0c4/r2c7 it indexes consistently (fit 0.6 deg) but the orientation differs from the .ang by a constant **~42 deg** (axis ~(-0.84,-0.54,-0.09) in the sample frame), not ~11 deg; image flips (ud, lr, rot180, T) do not reduce it to ~11. Patterns r5c3, r8c8 give other values (3-4 bands matched, unreliable). Interpretation open: the .ang (Thermo-exported) frame convention, the camera-tilt correction, or an incorrect PC could each contribute. **Needs Debasree's correction code and the raw .ang / calibrated PC.**

### 11.4 Preprocessing sweep (synthetic Ni, elev 10, true PC; degraded patterns: smooth background + vignetting, band contrast 15 % of background, 1 px blur, Gaussian noise sigma = 0.3 x band amplitude; `tools/degrade_synthetic.py`, `tools/sweep_pre.ps1`, `tools/sweep_table.py`; table `logs/val06_preprocessing_n30.md`)
- **Bug found in my own runner and fixed:** with the pattern-provider route (`Settings.patterns`) the image filter is `pats.filter` (default lower 9 / upper 90 px), **not `Settings.ImageFilter`** (legacy path only). `Settings.patterns` is a value object, so the filter must be changed before it is stored. The first two sweeps were invalid (filters had no effect) and were discarded. All earlier runs used the default 9/90, i.e. Brady's setting, so they are unaffected.
- Background handling (none / divide / subtract / divide+contrast stretch) at filter 9/90: **no difference** (scatter 1.5e-3, bias 1.6e-3 in all four). The 9-px lower radius of the image filter acts as a high-pass that removes the slowly varying background anyway; contrast stretching only rescales and cannot matter to a normalised cross-correlation.
- Filter radii (strain scatter RMS / max error): none 3.5e-3 / 1.7e-2; 9-90 (default) 1.5e-3 / 6e-3; 20-90 2.4e-3; 15-120 2.7e-3; 9-150 3.3e-3; **5-60 9.3e-4 / 3.5e-3; 9-45 1.0e-3 / 3.5e-3**. A tighter upper radius is better here (it removes noise).
- ROI size 15-50 % and number of ROIs 24/48/96: no systematic effect at this noise level (scatter 1.5e-3 throughout).
- Common bias 1-2e-3 in every run comes from the noisy reference pattern (r0c0), as expected for a single-pattern 'Real' reference; a cleaner or averaged reference would remove it.
- **Limits:** the synthetic master pattern is smooth (401-px Lambert), has no high-frequency detail beyond that, so the optimal upper radius for real patterns (which contain sharper bands and different noise) cannot be read off from this; the filter optimum must be re-checked on real data, ideally from the noise scatter in a strain-free region. The noise model is Gaussian, white, single orientation.

## 12. Deba's calibration (`correct_EBSD_hr.m`, `EBSD_Rcorr_calibration_20260505.mat`; copies in `inputs/Deba_calibration/`)

### 12.1 What it is
- `correct_EBSD_hr.m` (MTEX): `loadEBSD('Au111.ANG')`, loads `calib`, `Rcorr = rotation.byMatrix(calib.Rmat)`, applies **`ori_corr = Rcorr * ori` (left multiplication, specimen frame) to indexed points only**, then grain reconstruction/plots (GND section commented out). Needs MTEX only for this; the numbers are a plain matrix.
- `calib` (made 2026-05-05 in MATLAB R2025a from **Si111.ANG and Si100.ANG wafer scans**): `Rmat` = rotation of **11.04 deg about the specimen X axis** (axis (1.000, 0.012, 0)); per-wafer 10.51 deg (111) and 11.57 deg (100); orientation error before 10.5 / 11.6 deg, after 0.55 deg. Its own warning: corrects the common wafer-normal tilt; **does not determine the azimuth about the wafer normal** (notch/flat not tracked).
- Reproduced through MTEX 6.1.1 (`tools/apply_rcorr_mtex.m`; MTEX's `export` to .ang fails on this file, so Euler angles are dumped instead): the applied rotation equals `Rmat` exactly (0.000 deg); my numpy version (`tools/correct_ang.py`, az=0) reproduces MTEX's corrected Euler angles to 6 decimals.

### 12.2 Test on the Pt practice patterns (PyEBSDIndex orientations from the real patterns vs the .ang; `tools/pei_vs_ang.py`, `tools/pei_pc_scan.py`)
- Original Pt .ang vs Hough: 37.7 deg (camera elev 0) / 42.8 deg (elev 11).
- After Deba's Rcorr alone: 27.9 deg (elev 0) / 34.2 deg (elev 11): moves by ~10 deg as expected, but a large residual remains with axis ~(-0.7..-0.85, -0.4..-0.7, ~-0.2) in the sample frame.
- After **Rz(90 deg) * Rcorr** and Hough at elev 11: **2.55 deg**; with a PC scan around Brady's value the residual falls to **0.9 deg** at PC ~(0.503-0.513, 0.677, 0.915) (scan edge; PC and orientation are correlated for a single orientation, so this does not pin the PC). So the Pt sample sits **90 deg about the surface normal** from the Si wafers (or the pattern frame differs by 90 deg); this azimuth is not in Deba's calibration and was found empirically here. It must be confirmed (ask Deba: sample mounting/azimuth for Pt vs the Si wafers).
- With that corrected orientation (`work/pt10x10_corr/`, built by `tools/correct_ang.py --az 90`; NaN/NotIndexed rows left untouched) the earlier failing real-data PC refinement now matches: NCC 0.02 -> 0.24 (Ni master pattern as stand-in for Pt, Fiji-corrected patterns; synthetic reference value 0.48). Refined PC (orientation fixed): (0.4974, 0.6894, 0.9594); with a free small orientation correction (1.8 deg): (0.4928, 0.6551, 0.9693). Brady E.1-E.3 gave (0.51275, 0.657, 0.9545) for X=16.31, Y=6.86, L=19.09: differences 0.015-0.03 in x*, <=0.03 in y*, <=0.015 in z*, i.e. the Brady conversion is plausible for this detector. By §9.3 these x/y differences are harmless; z* differences of ~1.5 % should be too, but were not simulated.

### 12.3 OpenXY on the real Pt 10x10 subset (BYU master + C1 + C2, elev 11.04, default filter; plumbing / sensitivity only, no strain interpretation)
| comparison | max abs difference in deviatoric strain / rotation |
|---|---|
| Brady PC vs kikuchipy-refined PC (same corrected ang) | 1.0e-4 / 6.3e-5 rad |
| corrected vs uncorrected ang (same PC, elev 11.04) | 2.9e-3 / 2.7e-3 rad |
So at this scan the PC choice is immaterial, but the orientation/tilt correction is not.

### 12.4 How to use for real data (proposed, pending confirmation)
1. Correct orientations: `Rz(az) * Rcorr` left-multiplied on the TSL Euler matrix (R = g^T) - `tools/correct_ang.py --az <azimuth>`; Rcorr from the .mat; **azimuth per sample must be established** (here 90 deg for Pt).
2. Set `CameraElevation` = 11.04 deg (the tilt in Rcorr) in OpenXY at the same time. Orientation and elevation must describe the same geometry (see §10).
3. Pattern centre: Brady's conversion is consistent with kikuchipy refinement to ~0.03; refine per sample with `tools/kp_pc_refine.py` when a master pattern for the material exists (Ni stand-in used here).
Caveats: one pattern set, one orientation, Ni master pattern for Pt, azimuth chosen by best match (2.55 -> 0.9 deg residual), Hough indexing accuracy ~0.5-1 deg.

### 12.5 Is the 90 deg azimuth a mounting effect or an instrument convention? (`tools/pei_az_scan.py`)
- User: the Pt sample is a polycrystalline rectangle mounted with one edge roughly parallel to the holder; azimuth of the Si wafers was not tracked.
- Reasoning: indexing is referenced to the detector-fixed sample frame (RD/TD/ND of the scan), so mounting azimuth should not change the Euler angles of a given pattern; a **fixed 90 deg** between the Pathfinder .ang frame and the TSL/PyEBSDIndex frame would be an instrument/export convention, identical for all samples, and invisible to Deba's wafer-normal check (rotation about ND leaves the wafer normal unchanged). This is a hypothesis; the Pt data alone cannot separate it from a mounting effect.
- Evidence: patterns from both Pt grains in the subset (Euler ~(250,27,22) and ~(257,20,23) deg, ~8 deg apart) give best azimuth 89-90 deg with residual 2.2-2.5 deg; without the azimuth the residual is 27-34 deg. r5c3 is a NotIndexed point in the .ang (NaN) and was skipped. The two orientations are close, so this is a weak independence test.
- To settle: one more dataset from the same instrument with raw patterns and a very different orientation (e.g. Si111/Si100/Au111 patterns with their .ang), run through `pei_az_scan.py`. If the same 90 deg appears, it is an instrument convention and can be applied to all data; if it varies with the sample, it must be fitted per dataset (as a QC step: Hough-index 3-5 patterns, compare to the corrected .ang, expect <~3 deg).

### 12.6 Second and third Pt maps with raw patterns: the 90 deg azimuth reproduces (2026-11-04; `tools/pat_az_scan.py`)
- Box Drive is mounted locally (`C:\Users\Hang\Box\...`): read directly, copy to `inputs\`; no HTTP needed. Pt PC values are in the acquisition reports (`*.docx`, "Pattern Center (mm)"), not in the .ang (x/y/z-star = 0 there).
- **.pat layout (reverse-engineered, verified by frame-to-frame correlation):** 20-byte header = 5 x uint32 (1, ncols, nrows, pw, ph); pattern k = uint16 LE, pw*ph pixels at byte 25 + k*(pw*ph*2 + 8). Pt(2) and Si: 372x232 (5x binned); Pt_Hole: 184x116 (10x binned). Not full-resolution, so these are fine for orientation tests, not for HR-EBSD strain.
- Maps used (all copied to `inputs\`; not checked for PCs beyond the reports):
  - `Pt_HREBSD_20260427\Pt(2).MAP.XTAL` (G7. Deba\EBSD\2026\04-27-2026-HREBSD, 25 Apr 2026, 128x84, 0.47 um step). Report PC: X=16.22 Y=5.88 L=20.35 (a second report of the same session lists the default 15.8/7.9/16.0; Pt(1) of that session: 16.10/6.13/20.5). Practice-folder PC given by user: 16.17/5.76/20.43.
  - `Pt_Hole1_20260506\Pt_Hole(1)` (05-06-2026_Pt_repeat_calibration, 128x84, 2.52 um step, report PC is the uncalibrated default 15.8/7.9/16.0).
  - `Si111_20260504\111(2)` (64x42; Si is not indexable by PyEBSDIndex with the FCC/BCC reflector lists: 0-3 of 8 patterns index; the 3 that did give no azimuth consistency -> inconclusive, not used).
- Result (12 random Pt(2) patterns; 10 random Pt_Hole(1) patterns; orientation from PyEBSDIndex vs .ang corrected by Rz(az)*Rcorr, cubic symmetry):
  - Pt(2), PC 16.22/5.88/20.35: 8/12 patterns give 0.9-1.2 deg at az=90 (best az 90-91); 0 deg azimuth gives 41-54 deg. The 4 other patterns have weak Hough fits (nmatch 3-5) and scatter.
  - Pt_Hole(1), same PC: 10/10 patterns at az=90 give 1.0-2.3 deg (best common az 91, median 1.5 deg); az=0 gives 17-48 deg. With the default PC 15.8/7.9/16.0 the Hough indexing is poor (median 33 deg) -> **a calibrated PC is needed even for the orientation QC**.
- Interpretation: the Rz(90)*Rcorr correction holds for three Pt maps from different dates (Apr 2024 practice map, Apr 2026, May 2026), different sample regions, 100s of orientations. This supports "fixed Pathfinder-export convention" over "mounting azimuth", but all maps were mounted/acquired by the same group on the same instrument, so a differently mounted sample (e.g. a Si wafer indexed with a diamond-cubic reflector list, or another user's map) is still the decisive test.
- **Confirmed by the student (2026-11-04):** Pt(2) (04-27-2026-HREBSD, 25 Apr 2026) was acquired with Pathfinder PC X=16.22, Y=5.88, L=20.35, i.e. the 21:06 report value (the 21:08 report's default 15.8/7.9/16.0 was not used). Hough QC on 30 random patterns: median 15.8 deg at az=90 (report PC) vs 17.0 (practice-map PC) vs 37.1 (default); the first two are not separable by Hough. Pattern-centre sources: the session .docx reports; AcqProps.tnp holds only the default.
- **Azimuth check with the confirmed PC (log `logs/val07_pt2_azimuth.log`; 80 random Pt(2) patterns, nbands 11):** 78 indexed; 43 land within 3 deg at az=90 (median of those 1.0 deg, best azimuth 90-92 deg); 9 of the 10 patterns with nmatch>=5 are within 3 deg. No pattern is better than 22 deg without the azimuth. The rest are weak Hough fits (nmatch 3-4) on the binned 372x232 patterns. Whole-sample median is 1.2 deg at az=90. Conclusion unchanged: Rz(90)*Rcorr holds for Pt(2).
- Report-to-map assignment (first image of each .docx): 19.34 report = Pt(1) (PC 16.10/6.13/20.5); 21.06 and 21.08 reports are both Pt(2); the 21.08 one shows the default PC, the 21.06 one the calibrated PC confirmed by the student.

### 12.7 Candidate maps for a real HR-EBSD test (2026-11-04; headers read from Box Drive, nothing copied)
- `G7. Deba\EBSD\2026\07_17_2026\Pt(1..6)`: 256x170 maps, 1.64 um step, patterns 184x116 (10x binned), 1.77 GB each; all reports list the default PC 15.8/7.9/16.0 and WD changes between maps (9.9-11.1 mm). **Not suitable** (low resolution, uncalibrated PC).
- `G7. Deba\EBSD\2026\05-09-2025\` (report dated 2026-05-09; has Bkg_HR-EBSD.tif): Pt(6) and Pt(7) are 64x42 maps with 928x580 patterns (2x binned, 2.9 GB), step 0.273 / 0.328 um. The single report lists maps 256x170 (default PC), 128x84 (16.31/6.86/19.09), 64x42 at mag 12000 (no PC listed) and 64x42 at mag 10000 (16.72/6.54/18.89). Which entry belongs to Pt(6)/Pt(7) is not yet matched; PC to be confirmed by the student.
- Full-resolution (1856x1160) maps are the 46 GB .pat files (e.g. the 2024 practice map), too large to copy casually.
- Box Drive hydrates a whole .pat on first read (1-3 GB takes ~minutes).

### 12.8 ORR(2), 18 Mar 2026 (`G7. Deba\EBSD\2026\March_26\03_18_2026_ORR_HREBSD\ORR(2).MAP.XTAL`, copy in `inputs\ORR_20260318`)
- User info (student): the report .docx files are written **by hand**; the instrument always gives the default PC X=15.8 Y=7.9 L=16.0, which is never correct. A default value in a report therefore means "not entered", not a measured PC. Pt(2) of 25 Apr 2026 has two reports only because one report did not load all images. This folder has no report at all (PC unknown).
- ORR(2): 64x42, step 0.504 um, patterns 372x232 (5x binned), Generic FCC in the .ang. Hough QC (30 random patterns, nbands 11, `tools/pat_az_scan.py`): az=90 median 1.1 deg for PC 16.22/5.88/20.35 (April Pt(2) PC), 1.4 deg for 16.72/6.54/18.89, 2.1 deg for 16.31/6.86/19.09, but 33 deg for the default PC. az=0 gives 44-47 deg for all PCs. The Rz(90)*Rcorr correction therefore also holds for a **March 2026** map (before Deba's calibration date), another sample/session. Hough cannot separate the three calibrated PCs well (median-of-sum 8.6 / 15.0 / 13.6 deg, the April PC slightly best).
- Patterns are binned (372x232), as for Pt(2): good for orientation checks, limited for HR-EBSD precision.

### 12.9 OpenXY on a 10x10 block of ORR(2) (run11; plumbing / scatter only, strain values not interpreted)
- Tools: `tools/make_pat_subset.py` (cuts an N x N block from .pat + .ANG, patterns as uint16 tifs, 11-column grain file with grain ids from 5 deg flood fill; here rows 15-24, cols 12-21, entirely inside the largest grain, 1 grain), then `tools/correct_ang.py --az 90` (Rz(90)*Rcorr, 90.53 deg total), then `run_case.ps1` (BYU master + C1 + C2).
- Settings: PC (TSL) 0.5105 0.706 1.0175 (= Pathfinder 16.22/5.88/20.35, the **Pt(2) April value used as a stand-in; ORR(2)'s own PC is unknown**), CameraElevation 11.04, tilt 70, phosphor 403.2 um per binned pixel (5 x 80.64; 232 px = 93.5 mm as for 1160 px), default filter 9/90, GND off, `platinum` material, 3 workers.
- Result: 100/100 finite F. Scatter of the deviatoric strain over the block: std 0.1-0.3e-3 per component, max 1.5e-3, rotation std 0.01-0.02 deg. This is a noise/consistency figure (same grain, grain reference), not a strain measurement; the PC is not verified for this map.
- Not yet done: PC for ORR(2) from the student; comparison with uncorrected ang (as in 12.3); a full-resolution map.
- **Confirmed PCs (student, 2026-11-04):** ORR(2) (18 Mar 2026) X=16.80 Y=6.41 L=20.50; Pt(6) (05-09-2025 folder, 9 May 2026) X=16.31 Y=6.86 L=19.09; Pt(7) X=15.84 Y=6.96 L=18.70. (So the 16.72/6.54/18.89 and 16.31/6.86/19.09 entries in that report belonged to other maps / Pt(6); the Pt(7) value was not in the report.)
- ORR(2) with its own PC (TSL 0.525 0.6795 1.025): Hough QC 30 patterns: az=90 median 1.5 deg, az=0 44.7 deg. OpenXY 10x10 block (`tools/compare_runs_real.py`, logs `run11/12/13`):

| run | finite | mean per-component std (1e-3) | max abs (1e-3) | max abs diff vs run12 (1e-3) | mean abs diff vs run12 (1e-3) |
|---|---|---|---|---|---|
| run12 corrected ang, own PC | 100 | 0.23 | 1.46 | 0 | 0 |
| run11 corrected ang, April Pt(2) PC | 100 | 0.23 | 1.46 | 0.07 | 0.01 |
| run13 uncorrected ang, own PC | 100 | 0.24 | 1.58 | 1.81 | 0.34 |

  Within one grain (grain reference) the PC change of ~0.015 in x*/z* moves the strain by <=0.07e-3, whereas dropping the orientation correction moves it by up to 1.8e-3, consistent with the sensitivity measured on synthetic data (10.4). Values not interpreted.

### 12.10 Full ORR(2) map (run14) and a reference-selection trap (PC = confirmed 16.80/6.41/20.50, corrected ang, GND off; strain values not interpreted)
- Folder rename: the 9 May session folder is now `G7. Deba\EBSD\2026\05-09-2026` (the Box search index still says `05-09-2025`); Pt(6)/Pt(7) are copied to `inputs\Pt_HREBSD_20260509`.
- `make_pat_subset.py --full` builds the whole 64x42 map (13 grains by the 5 deg flood fill, two big ones split by a vertical boundary at column ~36; grain 1: 1485 pts, grain 2: 1153).
- Result: 2688/2688 finite, but **grain 2 is garbage** (mean F-I has F21 = -0.178, median deviatoric component 7.7e-2), while grain 1 gives 0.6e-3 (median). The ang says every grain-2 point is within 0.2 deg of OpenXY's reference point, so a 10 deg pattern rotation is not real.
- Cause: OpenXY's reference selection `IQ > Fit > CI` chose the reference of both grains at column 35, i.e. right at the grain boundary (grain 2: row 26, col 35), where the pattern is a superposition of two crystals. A 10x10 block in the grain-2 interior (rows 15-24, cols 50-59, run15; reference chosen inside the block) gives normal results (std 0.21e-3, max 0.96e-3).
- **Rule for real maps:** check where the reference lies (`Settings.RefInd`, row/col) and keep it a few pixels (>= ~3 steps) from any grain boundary; otherwise override it (or erode the grain mask before the reference selection) and compare NCC of the reference to the grain points (here 0.85 for the bad grain vs 0.95 for the good one).

### 12.11 Pt(6) and Pt(7) (9 May 2026, 928x580 patterns = 2x binned; copies in `inputs\Pt_HREBSD_20260509`; PCs confirmed by the student)
- Pt(6) PC 16.31/6.86/19.09 (TSL 0.51275 0.657 0.9545); Pt(7) PC 15.84/6.96/18.70 (TSL 0.501 0.652 0.935). Both maps 64x42; 44 / 32 grains by the 5 deg flood fill.
- Hough QC (30 random patterns, 11 bands): Pt(6) az=90 median 0.9 deg (az=0: 31.0 deg); Pt(7) az=90 median 3.7 deg (az=0: 35.8 deg, best common az 92). Correction holds on a third session (9 May 2026, after Deba's calibration), with the confirmed PCs.
- OpenXY 10x10 blocks inside the largest grain (Pt(6): rows 20-29, cols 12-21; Pt(7): rows 18-27, cols 42-51), corrected ang (Rz(90)*Rcorr), elevation 11.04, 161.28 um per binned pixel (580 px square = 93.5 mm), default 9/90 filter, GND off, reference inside the block (runs 16, 17): 100/100 finite; per-component scatter 0.05e-3 (max 0.23e-3) and 0.07e-3 (max 0.41e-3), i.e. ~3-4x lower than the 372x232 ORR(2) block (0.23e-3, 10.4 px vs 2.5 um of resolution differ by 2.5x in pattern size). Noise-floor figures only; values not interpreted.

### 12.12 Filters on real Pt(6)/Pt(7), warped-pattern recovery test, corrected full maps (strain values not interpreted)
- **Correction to 12.10:** `Settings.RefInd` is **1-based**. The ORR(2) grain-2 reference was therefore (26, 34), the max-CI point and the first column of grain 2 at the boundary; its NCC with the grain points is 0.15 (not 0.85). Conclusion unchanged, the reference sat in the boundary mix.
- **Reference selection** (`GetRefImageInds.m`): votes on max CI, min Fit and max IQ in each grain and takes the best trade-off. Making only IQ unattractive near boundaries (first attempt) did not work. `make_pat_subset.py --ref-margin 4` now sets IQ=0, CI=0, Fit=99 within 4 px of a boundary (ang and grain file); the references are then interior (IQ 790-800, CI 0.86-0.92) for the large grains.
- **Scatter in the 10x10 blocks vs filter** (`sweep_real_filters.ps1`, `compare_runs_real.py`; mean per-component std, 1e-3): Pt(6): 9/90 0.05, 5/60 0.06, 9/45 0.07, 5/45 0.07, 15/120 0.04, none 0.03; Pt(7): 0.07, 0.08, 0.08, 0.08, 0.06, none 0.04. Opposite to the synthetic sweep (where no filter was worst), so scatter alone is ambiguous (a strong common background can suppress real shifts).
- **Recovery test** (`make_strained_block.py`, `sweep_warp_filters.ps1`, `warp_recovery.py`): every non-reference pattern of each block is warped by a known deviatoric F (3e-3 x [[1,.7,0],[.7,-1,0],[0,0,0]], norm 5.18e-3, ray map around the PC, bicubic); recovered change of |dev| relative to the unwarped run, same filter (median over 99 points): Pt(6): 9/90 0.95, 5/60 0.94, 9/45 0.93, 5/45 0.94, 15/120 0.97, none 0.98; Pt(7): 0.95, 0.94, 0.93, 0.94, 0.96, 0.98 (IQR of the ratio <= 0.01). So on these real patterns no filter suppresses the signal, all recover 93-98 %, the unfiltered / wide-band filters slightly most, and also give the lowest scatter. Differences are <= 5 % of 5e-3 (0.25e-3); the 9/90 default is acceptable. Caveats: one strain tensor, one block per map, same PC used to warp and to analyse (a PC error is not tested here), norm-only comparison (the sample-frame rotation hides the components).
- **Full maps with the interior reference** (runs 26, 27): all 2688 finite. Pt(6): per-grain median |dev| 0.26e-3 (grain 1, 1220 pts, now fixed; was ~320e-3), 0.50e-3 (grain 3); Pt(7): 0.60e-3, 0.34e-3. Fraction of points above 5e-3: 2.8 % (Pt(6)), 3.5 % (Pt(7)), consisting of
  - **scan column 0**: 41 of 42 rows in both maps. The col-0 patterns differ from col 1 (NCC 0.97 vs 0.999 between other neighbours, mean intensity +13 counts) - a first-point-of-line artifact (beam settling / drift?), not OpenXY. Exclude column 0 (or ask the instrument people why).
  - 35 (Pt(6)) and 54 (Pt(7)) other points, mostly in noisy sub-grains with CI 0 / IQ 0 and along a diagonal line in Pt(7) (grain-boundary or scratch).

### 12.13 Column-0 artefact, full-map statistics without it, and analysis-PC error test (Pt(6)/Pt(7) blocks; strain values not interpreted)
- **Column 0** (first point of each scan line) is anomalous in every Pt/ORR map from March-May 2026 (Pt(6), Pt(7), Pt(2), ORR(2), Pt_Hole(1)): NCC(col0, col1) = 0.88-0.98 vs 0.985-0.999 between other neighbours; on 9 May mean intensity +14 counts (0.8 %) vs column 1; **no translation** (phase correlation < 0.1 px); OpenXY strain vs an interior reference exceeds 5e-3 in 41/42 rows. Last columns (62, 63) are normal. Si(111) 4 May shows no difference (near-identical patterns, uninformative). Questions for Deba: (1) is the first pattern of each line acquired differently (beam settling / flyback delay, auto-gain, different exposure)? (2) is there a settling / line-delay setting in the map set-up? (3) is the background correction applied the same way to column 0? (4) treat column 0 as invalid, or add a dummy column?
- **Full maps excluding column 0** (`tools/full_map_stats.py`, reference >= 4 px from boundaries): Pt(6): median 0.31e-3, 1.3 % above 5e-3, usable (col > 0, instrument CI > 0, |dev| <= 5e-3) 2368/2688 = 88 %; bad points are in the tiny unindexed grains 2 and 11. Pt(7): median 0.41e-3, 2.0 % above 5e-3, usable 2499/2688 = 93 %; all 54 bad points are in grain 4, indexed by the instrument and lying within ~4 px of the grain boundary (pattern mixture). Instrument CI taken from the original .ang (the working copy has CI = 0 near boundaries on purpose).
- **Analysis-PC error, no rotation** (`sweep_pc_warp.ps1`, `pc_error_table.py`; known strain 5.18e-3 norm): moving x*, y* by +-0.01 (~6 px) or z* by +-0.02 changes the strain of the unwarped block by <= 0.01e-3 and the recovered fraction stays 0.95 on both maps. Without relative rotation the PC hardly matters.
- **With a 1 deg rotation imposed on top of the strain** (`sweep_pc_rot.ps1`, Pt(6)): extra strain from the PC error (mean over non-reference points): x* +-0.01: 0.24e-3; y* +-0.01: 0.12e-3; z* +-0.02: 0.48e-3 per degree of rotation, in line with the first-order estimate (rotation [rad] x dPC/z* ~ 0.18e-3 per 0.01). So PC accuracy matters in proportion to the rotation between pattern and reference; a PC good to ~0.003 gives <= ~0.1e-3 per degree.
- **Caveat found:** in this warp test OpenXY recovers only ~0.95 of the imposed strain (0.80 when 1 deg of rotation is added) and 0.945 deg of the 1 deg rotation. Either my warp model (ray map around the PC without the camera elevation / bicubic smoothing) differs from OpenXY's geometry by ~5 %, or OpenXY's recovered F is ~5 % low for these 580 px patterns; the synthetic kikuchipy tests (independent geometry) are the better calibration and should be compared before quoting absolute accuracy.

**Decision (user, 2026-11-04):** column 0 of every scan line is dropped from all full-map results (the acquisition cause is unknown; with maps this large the loss is acceptable). Deba question is optional, not blocking. 	ools\full_map_stats.py already excludes it; apply the same mask to any later map or figure.

### 12.14 Recovery deficit vs synthetic truth (`tools/syn_recovery_ratio.py`; medians of recovered/imposed norms, exact PC, OpenXY-geometry synthetic Ni patterns `work\syn_A`, `work\syn_C`)
- Pure strain up to 5e-3 (24 cases): **1.001**. Rotation <= 0.5 deg: rotation 0.959. Rotation 1-4 deg: rotation 0.967. Mixed small: strain 1.024, rotation 0.935. Large mixed (up to 3 deg): rotation 0.926, strain norm 1.66 (rotation leaks into strain, §9.2).
- So OpenXY itself shows **no strain deficit** on exact-geometry patterns, but under-recovers rotation by ~3-7 %. The real-pattern rotation recovery (0.945) agrees with that.
- The 0.95 strain recovery on real patterns therefore most likely comes from my warp (bicubic smoothing, no camera elevation in the ray map) or the real-pattern content, not from OpenXY. The 0.80 with 1 deg rotation mixes the rotation under-recovery and rotation-to-strain leakage. Do not quote an absolute accuracy from the warp test; use the synthetic figures (~1e-4 to 5e-4 absolute, §9.2). Not resolved by a dedicated test (e.g. warp the synthetic patterns with my warp code).

### 12.15 My warp applied to synthetic patterns (`tools\make_strained_block.py`, `run_syn_warp.ps1`, `syn_warp_ratio.py`; work\syn_warp, syn_warpR; Ni, 1160 px, elevation 0, true PC, all 100 patterns = the reference pattern before warping, 3e-3 strain, norm 5.18e-3)
- Strain only: recovered/imposed **0.962**, rotation 0.04 deg. With 1 deg rotation: strain **0.847**, rotation 0.995 deg (imposed 1.00).
- Real Pt(6) (580 px): 0.95 / 0.80 / 0.945 deg. The strain deficit is therefore reproduced on synthetic patterns, so it comes from my warp model (not from real-pattern content, filters or elevation, which was 0 here), while OpenXY's own forward-model patterns give 1.001 (§12.14). The 1 deg rotation case also loses ~15 %, i.e. warp + rotation interact.
- Not isolated yet: the difference between my ray-map warp (bicubic, F applied to the ray around the PC) and `gen_synthetic.py`'s rendering (kikuchipy sampler, OpenXY geometry) - interpolation or a frame/sign detail in the shear/rotation terms. The warp-test recovery fractions are not a property of OpenXY; the real-pattern strain figures are unaffected by it.

### 12.16 Cause of the warp/rotation discrepancy: y* convention in my warp (`make_strained_block.py` usage, `warp_recovery_y.py`, `run_pt6_warpY.ps1`)
- `gen_synthetic.py` (OpenXY geometry): ray y = (1-y*)P - row, i.e. y* is measured **from the bottom**; my warp put the PC at y*·P from the top (mirrored). **Pass `--pc x (1-y*) z` to `make_strained_block.py`** (e.g. 0.51275 0.343 0.9545 for y* = 0.657).
- With that fixed: synthetic strain only 0.954 (was 0.962), strain + 1 deg rotation **0.996** (was 0.847), rotation 0.996 deg. Real Pt(6): strain only 0.950, strain + 1 deg **0.999** (was 0.80), rotation 0.91 deg. So the earlier 0.80 was my mirrored PC; the earlier PC-error-with-rotation table (§12.13, x*/y*/z* extra strain) was computed with the mirrored warp and should be redone before it is quoted (the no-rotation PC results are unaffected).
- Still unexplained: strain-only recovery 0.95 on both synthetic and real, whereas OpenXY-rendered crystal-frame patterns give 1.001 (§12.14). Not the dilatation closure (adding d*I cannot change the deviator; checked). Candidates not yet tested: my mixed e11/e22/e12 detector-frame tensor (out-of-plane shear components in the crystal frame), bicubic resampling, the x-mirroring of my ray map.

### 12.17 Remapping pass around OpenXY (`tools\remap_pass.py`, `remap_eval.py`, `run_remap_syn.ps1`; synthetic syn_C, true PC, elevation 0)
- Method (black box, no MATLAB change): pass 1 = normal OpenXY run -> F = R1 U per pattern; un-rotate each test pattern by R1 in the crystal frame (I_new(r) = I_test(R1 r) = I_ref(U^-1 r)); pass 2 = OpenXY on the un-rotated patterns -> F2; final F = R1 . F2. OpenXY's own 'Remapping' method (`rotateImage.m`, GUI-selectable) exists but was never used here; its code calls `getPattern` without Settings (broken as shipped) and loops over every pixel, so I did not use it.
- Error vs truth (max / median dev-strain error 1e-3; max rotation error mrad):

| group | n | pass 1 | remap | rot pass 1 | rot remap |
|---|---|---|---|---|---|
| pure strain | 24 | 0.90 / 0.43 | 0.84 / 0.43 | 0.62 | 0.60 |
| rotation <= 0.5 deg | 12 | 1.22 / 0.42 | 0.38 / 0.19 | 0.91 | 0.25 |
| small mixed | 12 | 0.80 / 0.40 | 0.50 / 0.35 | 0.57 | 0.33 |
| rotation 1-4 deg | 18 | 11.75 / 3.17 | 1.89 / 0.40 | 7.31 | 1.37 |
| large mixed (<= 3 deg) | 17 | 8.26 / 2.96 | 1.54 / 0.49 | 6.29 | 1.10 |
- So one remap pass removes most of the rotation-to-strain leakage and the rotation under-recovery. Not yet tested on real maps, with a second pass, with elevation 11.04, or with noisy patterns.

### 12.18 Remap pass on the real Pt(6)/Pt(7) 10x10 blocks (`run_remap_real.ps1`, `remap_real_eval.py`; no truth, strain values not interpreted)
- Rotation of every point relative to the reference is tiny (median 0.005-0.007 deg, max 0.013-0.016 deg), so the remap pass changes the deviatoric strain by only 0.06-0.10e-3 (median), <= 0.25e-3 (max); the within-block scatter rises slightly (0.050 -> 0.065e-3 Pt(6), 0.066 -> 0.080e-3 Pt(7), extra resampling noise). Remapping is therefore not needed inside a grain with a grain reference; it is for points with larger rotation (sub-grains, near boundaries, large-rotation maps). Not yet run on a full map.

### 12.19 Pattern-centre refinement with a kinematic Pt band model (`tools\kin_pc_refine.py`; no EMsoft and no Pt master pattern available here)
- Model: exact FCC Pt plane positions (a = 3.924 A, 20 kV), Gaussian band profile, geometry of `gen_synthetic.build_rays_crystal`; only band positions are trusted. Orientation fixed (corrected ang, Rz(90) Rcorr). NCC is low in absolute terms (max ~0.14-0.16, the Ni master gave 0.24), because the band shape is crude.
- Pt(6), 2 clean patterns from the 2 large grains (margin 6): refined PC (0.5014, 0.6501, 0.9511) vs Brady/confirmed (0.51275, 0.657, 0.9545): differences (-0.011, -0.007, -0.003). With a third pattern from a CI = 0 grain the answer moves by -0.010 in z* (that point should be excluded).
- Free camera elevation (start offset 0.02/-0.02/0.03): elevation -> 10.41 deg (from 11.04) while y* moves +0.011: elevation and PC are strongly coupled; with 2 patterns of 2 orientations it is **not identifiable**. Same NCC (0.0974) as the fixed-elevation fit.
- Pt(7), 2 patterns: NCC at the confirmed PC is ~0 (-0.003); refined PC differs by (+0.044, +0.019, +0.021) and reaches NCC 0.138 (elevation free: 10.40 deg). Conflicts with the confirmed Pt(7) values by more than the Pt(6) result does; either the corrected orientation for these two points, the typed PC for Pt(7), or the crude model is off. Unresolved; only two patterns.
- Take-away: Brady/typed PCs are good to ~0.01 for Pt(6); do not use these fits to replace them. Tilt as a free variable needs more grains (map with >= 6 well-indexed grains) and a better band model.

### 12.20 Real Pt master pattern (ebsdsim) and PC checks with kikuchipy (`tools\make_pt_master.py`, `kp_pc_refine.py --master`)
- `pip install ebsdsim` (0.2.4) into `work\kpenv`; WebGPU works on this PC (Quadro RTX 4000, Vulkan). **Cell units are Angstrom** (first attempt with 0.3924 was wrong). `make_pt_master.py 250` makes a dynamical Pt (Fm-3m, a = 3.924 A, 20 kV, 501x501) master pattern in ~2-3 s: `work\pt_master\Pt_20kV_hw250.npz`; kikuchipy 0.13.1 reads it with `kp.load(path)` (upper hemisphere only in the file; get_patterns worked).
- NCC against 3 real Pt(6) block patterns (corrected ang, typed PC): Pt master 0.342 start / 0.484 refined; Ni stand-in 0.217 / 0.278. The Pt master is clearly the better model.
- **Pt(6)**: refined PC (0.5009, 0.6519, 0.9644) vs typed (0.51275, 0.657, 0.9545): (-0.012, -0.005, +0.010); agrees with the kinematic fit (§12.19) and the Ni fit (-0.013, -0.005, +0.005). In Pathfinder units (inverse Brady): X 15.84, Y 6.96, L 19.29 (typed 16.31, 6.86, 19.09).
- **Pt(7)**: NCC at the typed PC is -0.05 (no match); refined PC (0.5448, 0.6714, 0.9601), NCC 0.49, i.e. (+0.044, +0.019, +0.025) from typed (0.501, 0.652, 0.935); with a free orientation correction (-1.1 deg about x) (0.5248, 0.6667, 0.9616). Pathfinder units X 17.59, Y 6.57, L 19.20 (typed 15.84, 6.96, 18.70). Same answer as the kinematic fit, so the Pt(7) typed values (or its orientation) do not fit the data. Note: Pt(6)'s refined X,Y (15.84, 6.96) equal Pt(7)'s typed X,Y - possibly a transcription mix-up; not established.

## 13. One-command workflow (`tools\hrebsd.py`)
`python tools\hrebsd.py "<MAP.XTAL folder>" --pc X Y L [--block 10] [--name N] [--phase Platinum --material platinum --az 90 --elev 11.04 --ref-margin 4 --gnd]`, stages `check | prepare | run | post | all`; run from the repo root with the local copy of the map folder.
- **check** stops (exit 1) on: no .pat/.ANG, truncated .pat (Box not fully downloaded), .pat/.ANG grid mismatch, wrong row count, PC missing / equal to the instrument default 15.8/7.9/16.0 / implausible; reports pattern size, not-indexed (NaN) and CI <= 0 counts, the wrong phase in the .ANG header ('Generic FCC', a = 4.0) that is replaced by the Platinum grain file, and that the .ANG PC header is 0.
- **prepare**: exports patterns, grain file (grains by cubic misorientation < 5 deg), zeroes CI/IQ/Fit within `--ref-margin` px of a boundary so the reference is chosen inside the grain, applies Rz(az) Rcorr to the orientations (NaN rows untouched), writes `<stem>_index.txt` (original index, instrument CI, column).
- **run**: OpenXY headless (patched BYU copy), mperpix = 93542.4 / pattern height, PC by Brady E.1-E.3. **post**: `results.csv` (per point flags col0 / not_indexed / near_boundary / usable, rotation, deviatoric strain), `summary.md`, `maps.png`.
- Acceptance test: `Pt(6)` 10x10 block through the tool reproduces the earlier run16 result exactly (max difference 0.00e-3).
- Not yet tested through the tool: full maps, other materials/phases, ORR, Si; `--az 90` is still the empirical azimuth (NOTES 12.6); GND post-processing is not included; `run_case.ps1` leaves the MATLAB exit code blank, so success is read from the log.
- **Saved maps** (added to `hrebsd.py post`): `maps.npz` holds every per-point quantity as (nrows, ncols) arrays (e11 e12 e13 e22 e23 e33 deviatoric strain, `dev_strain` (r,c,3,3), full `F` (r,c,3,3), rotation vector rotx/roty/rotz in deg, rot_deg, strain_norm, grain, instrument_CI, usable/col0/not_indexed/near_boundary/finite flags, x_um/y_um, original_index) plus `results.csv`. `python tools\replot_maps.py work\hrebsd_NAME\maps.npz [--clim 2e-3] [--rot-clim 0.5] [--show-unusable] [--cmap ...]` redraws `strain_components.png`, `rotation_components.png`, `overview.png` without rerunning OpenXY. GND maps are not produced (run with `--gnd` is available but not post-processed).
- Full Pt(6), Pt(7) maps through the tool: `tools\run_full_via_tool.ps1` (output `work\hrebsd_Pt6full`, `hrebsd_Pt7full`), started 2026-11-04.

### 13.1 Full maps through the tool, PC check, GND (2026-11-05; strain/GND values not interpreted)
- `hrebsd.py` full Pt(6) / Pt(7): identical statistics to the earlier hand-driven runs (Pt(6) 88.1 % usable, median 0.31e-3; Pt(7) 93.0 %, 0.41e-3). A full map takes ~3 min.
- **Pt(7) fitted vs typed PC** (X 17.60, Y 6.58, L 19.18 vs 15.84, 6.96, 18.70; `compare_maps.py`): on the 2499 points usable in both, |dev difference| median 0.042e-3, p90 0.10e-3, max 0.64e-3; rotation difference median 0.0004 deg. A 0.04 PC error is harmless here because rotations relative to the reference are only ~0.03 deg (consistent with the first-order rule, NOTES 12.13). It would matter for larger rotations.
- **GND** (`--gnd`, OpenXY 'Full' method = `DislocationDensityCalculate`, MisoTol 5 deg, points with different grain/high misorientation filtered, 298 of 2688 on Pt(6)): saved as `gnd3/gnd5/gnd9` (estimates from 3/5/9 Nye components), `nye_alpha` (r,c,3,3), `gnd_filtered`, `gnd_misorientation`; `gnd.png` (log10). Pt(6) usable points: median ~1.1-1.3e13 (p10 4-6e12, p90 2.7-3.3e13) in OpenXY units (presumably 1/m^2; units not verified). Map is pixel-to-pixel noisy (~0.5 decade): likely close to the noise floor for a 0.27 um step; a noise floor was not measured (needs GND of the synthetic strain-free block or a repeat-pattern test) - do not interpret absolute GND values yet.

### 13.2 GND noise floor (Pt(6), `work\hrebsd_Pt6gnd\maps.npz`)
- Per-point rotation noise (MAD of residual to a 3x3 median, in-grain): 0.003-0.004 deg per component (grain 1: 0.0028/0.0031/0.0016; grain 3: 0.0041/0.0037/0.0029).
- White-noise GND from that noise: ~sqrt(2)*sigma/(b*step) = 1-1.6e12 per component (b = 0.2775 nm, step = 0.273 um); OpenXY sums several components with prefactors, so its scalar is larger.
- Shuffle test (random permutation of rotations inside each grain, same finite-difference metric): real field gives 2.3e12 (grain 1) and 3.9e12 (grain 3); shuffled gives 3.2e12 and 8.0e12. So real GND is only 1.4-2x below pure noise.
- Lag-1 correlation of log10(OpenXY gnd5) between neighbours: 0.12 (grain 1), 0.19 (grain 3). Noise-dominated.
- Conclusion: OpenXY GND of ~0.9-1.4e13 in Pt(6) IS the noise floor at this step. Do not interpret. Floor scales as 1/step, so binning 2x2 / using every 2nd point halves it at the cost of resolution; only a field clearly above ~1.5e13 (e.g. near boundaries, indents) is resolvable.
- Not done: floor from a rendered strain-free synthetic map through OpenXY itself (this estimate uses the real rotation noise instead); units of OpenXY's value still unverified (assumed 1/m^2).

### 13.3 GND floor through OpenXY itself, and binned GND
- `tools/make_null_map.py`: copies the Pt(6) work dir and replaces every usable pattern in grains 1 and 3 by ONE pattern of that grain plus independent Gaussian noise (sigma 15.5 counts, measured between real neighbours; signal std ~690). True strain, rotation and GND = 0. Run through `hrebsd.py --gnd` as `work\hrebsd_Pt6null` (5.7 min).
- Null result (median, unfiltered usable points): strain norm 1.2e-4, rotation 0.005 deg, OpenXY gnd5 6.2e12 (grain 1) / 6.4e12 (grain 3), p10-p90 2.8e12-1.3e13. Real Pt(6): gnd5 9.4e12 / 1.4e13, strain norm 4.9e-4 / 9.5e-4, rotation 0.011 / 0.029 deg. Real/null gnd5 = 1.5 and 2.2. So the floor measured through OpenXY is ~6e12 and real GND is only 1.5-2x above it at the native step. (The noise sigma includes a little real rotation, so the floor is slightly conservative; the null has no pattern-shape variation either.)
- `tools/gnd_binned.py`: GND from the saved rotation field, k x k binned (block-mean rotation, forward differences between blocks in the same grain, alpha13, alpha23, alpha33, x3 / b). k = 1 is ~0.63x OpenXY gnd3 (different stencil) but the real/null ratio matches OpenXY's (1.6, 2.1), so it is a usable surrogate.
- Median GND (1/m^2, my formula) real vs null and ratio, grain 1 / grain 3:
  k=1 6.1e12 / 3.9e12 (1.6) ; 9.8e12 / 4.7e12 (2.1)
  k=2 2.0e12 / 9.6e11 (2.1) ; 4.0e12 / 1.2e12 (3.3)
  k=3 1.4e12 / 4.7e11 (3.0) ; 2.6e12 / 5.4e11 (4.9)
  k=4 9.9e11 / 2.3e11 (4.4) ; 2.1e12 / 2.6e11 (8.0)
  Null floor falls ~1/k^2 (noise averaging x 1/k from the longer baseline); real values flatten (~1e12, ~2e12), i.e. real signal emerges above the floor from k=3. Spatial resolution at k=3 is 0.82 um.
- `hrebsd.py post --gnd` now also saves `gnd_bin2/3/4` in maps.npz and writes `gnd_binned.png`. Units 1/m^2 assumed (b = 0.2775 nm).
- Caveat: the Nye tensor uses only the 3 components obtainable from a 2-D scan (lower bound); no interpretation of values here.

### 13.4 GND units (OpenXY source) and Pt(7) null test
- `DislocationDensityCalculate.m`: stepsize_orig = (XData(3)-XData(2))/1e6 (m); Burgers from the material file (platinum 2.7719e-10 m); beta = (F - I)/step; alpha = curl(beta)/b, so alpha and alpha_total* are in 1/m^2. alpha_total3 = 30/10 * sum|alpha_i3|, total5 = 30/14 * (5 comps), total9 = 30/20 * (9 comps). The script's own noise cutoff is log10(0.006 deg / (step*b)) = 1.4e12 per component at our step: consistent with the floor measured here.
- Pt(7) (fitted PC 17.60/6.58/19.18, grains 1 and 4), real vs strain-free null through OpenXY: gnd5 median 1.05e13 / 1.00e13 vs 5.1e12 / 5.3e12 (ratio 2.1, 1.9). Floor in Pt(7) ~5e12, Pt(6) ~6e12: same, so the floor is a property of the noise/step, not of the map.
- Binned (my 3-component formula) real/null (ratio): grain 1: k1 9.5e12/3.5e12 (2.7), k2 4.4e12/8.4e11 (5.2), k3 3.3e12/3.9e11 (8.5), k4 3.3e12/2.1e11 (15.5); grain 4: k1 6.6e12/3.3e12 (2.0), k2 2.5e12/7.8e11 (3.3), k3 1.6e12/3.6e11 (4.5), k4 1.3e12/2.1e11 (6.1). Real values plateau at k>=3 (3e12 and 1.3-1.6e12), floor keeps falling.
- Pt(7) grain 1 also has larger rotations vs its reference (median 0.049 deg vs 0.014 deg in grain 4), consistent with more real lattice curvature; not interpreted.

### 13.5 Full-map remap on real data (tools\remap_full.py)
- `remap_full.py` takes a finished hrebsd.py run, un-rotates every non-reference pattern by R1 of its pass-1 F (rays in the crystal frame of its grain's reference, per-grain Euler from the corrected .ang), reruns OpenXY as `<NAME>_remap` and writes `work\hrebsd_<NAME>_remap\maps.npz` with F = R1 . F2 (pass-1 F and strain kept as `F_pass1`, `dev_strain_pass1`). Pass 2 takes ~5 min.
- Pt(6) full map (remap every point, min-rot 0): median |final - pass 1| deviatoric strain 0.27e-3 (p90 0.97e-3, max 3.4e-3). Median strain norm pass 1 -> final: grain 1 0.50 -> 0.61, grain 3 0.95 -> 1.37 (x1e-3). Pass-2 residual rotation 0.006 deg (was 0.011), so the un-rotation direction is right.
- Control on the strain-free null map (truth = 0): strain norm 0.12 -> 0.18 (grain 1) and 0.12 -> 0.23 (grain 3) x1e-3, rotation 0.005 -> 0.008 deg. So at ~0.01-0.05 deg rotations the remap pass ADDS noise (it combines two noisy fits) and cannot help; the Pt(6) increase above is the same effect.
- Real maps have almost no rotation to remove: Pt(6) p99 0.06 deg, max 0.28 deg; Pt(7) p99 0.10 deg, max 0.34 deg; none above 0.5 deg. Leakage matters only for large rotations (synthetic: 11.75 -> 1.89e-3 at large rotation, NOTES 12.17).
- Decision: remap is OFF for these maps. `remap_full.py` has `--min-rot` (default 0.5 deg): only points rotated at least that much are remapped, everything else keeps pass 1. Use it on deformed samples / reference-far points, not on these.

## 14. Cleanup (2026-11-05)
- ORR(2) run end to end through `hrebsd.py` (`work\hrebsd_ORR2`, PC 16.80/6.41/20.50, `--gnd`): check/prepare/run/post in 1 m 43 s; 95.3% usable, 13 grains (two large ones), median max-component strain 0.57e-3, median rotation 0.041 deg. A 10x10 block (`ORR2blk`, 28 s) after the cleanup also runs.
- 40 development scripts removed from `tools\` and archived in `work\removed_tools_2026-11-05.zip` (sweeps, one-off warp/recovery tests, pei_* Hough variants, kinematic PC fit, block-only remap, old subset/preflight tools). The sections above that cite them describe what they showed; they are no longer runnable from `tools\`. Kept: hrebsd.py, make_pat_subset.py, correct_ang.py, replot_maps.py, run_case.ps1, run_openxy_headless.m (+ headless_shims), pat_az_scan.py, gnd_binned.py, compare_maps.py, make_null_map.py, remap_full.py, map_pc_fit.py, make_pt_master.py, gen_synthetic.py, compare_synthetic.py.
- OpenXY patches C1/C2 copied to `tools\openxy_patches\`; pinned Python packages in `docs\requirements.txt`.
- The how-to document is `docs\hrEBSD_workflow.md`.

## 15. General EBSD lessons from the Zn-on-HOPG project (2026-11-10)
That project is kept separate: code, findings and figures are in the private repo `ren-hang/hopg-zn-orientation`; data and work are in
`C:\Users\Hang\Projects\HOPG_Zn`. Its data are Elias's SECCM Zn grids on HOPG: Thermo Pathfinder maps of 16 and 18 Apr 2025, 256x170 at
0.546 um nominal step with 184x116 patterns, and Nov 2024 / Feb 2025 maps with 372x232 patterns. This section keeps what applies to EBSD in
general. Tools named `HOPG:...` live in that repo. Only `tools\register_grid.py` was added here.

### 15.1 Thermo .pat, beyond 12.6
- Record 0 is dark in every map seen; leave it out of backgrounds and means.
- The 8 bytes after each pattern are a per-record counter, not pixel data.
- The file holds one record more than ncols x nrows.
- A 256x170 map of 184x116 patterns is 1.86 GB. With `np.memmap` at `25 + k*(pw*ph*2 + 8)`, the mean of every pattern takes ~1 min.
- That per-pattern brightness map shows surface features (deposits, particles, indents) as dots, so it is the map to register against
  (15.6).

### 15.2 PC and detector tilt without indexing (`HOPG:tools/hopg_twist_map.py --fit-pc`)
- Graphite patterns (10x binned, weak) could not be Hough-indexed with PyEBSDIndex.
- With the orientation family constrained (c along ND, only phi1 free), the PC and the detector tilt were fitted instead. The fit is
  Nelder-Mead on NCC against ebsdsim simulations, using 8 spread points.
- The fit reproduced within a session (Nov 2024: X 16.0-16.3, Y 6.7-6.8, L 20.3-20.4 mm, tilt 11.0-11.3 deg). Across sessions it moved
  (Feb 2025: X 16.8-17.0, Y 5.3-5.9; Apr 2025: 16.0-16.8 / 5.2-5.5 / 20.5-20.6, tilt 11.1-11.7).
- The acquisition reports showed the default 15.8/7.9/16.0 or inconsistent values (see 12.8).
- In kikuchipy, sample tilt and detector tilt enter the projection together, as alpha = 90 - sample_tilt + detector_tilt. A fit therefore
  fixes alpha (~31 deg here), not the two tilts separately. This is the same coupling seen in 12.19, and it matters for 15.4.

### 15.3 ebsdsim for other phases, and an average-Z bug
- ebsdsim 0.2.4 makes hexagonal master patterns as easily as cubic ones (graphite P6_3/mmc with 2 sites, Zn P6_3/mmc), each in a few
  seconds on the GPU.
- Each pattern is a Monte Carlo electron run at the chosen `voltage_kv` and specimen tilt `sigma_deg` (default 70), followed by a
  Bloch-wave solve.
- `mp.bin_voltages_kv` and `mp.bin_weights` give the backscattered-electron energy distribution. At 20 kV, the share in the top 1 keV bin is
  Pt 52% and Zn 31%; light phases contribute many low-energy (deeper, more scattered) electrons. The graphite value (4%) comes from the
  wrong Z described next.
- **Bug:** for structures with more than one atomic site, ebsdsim 0.2.4 stores the average Z and A divided by the number of sites
  (graphite 3.0 / 6.0 instead of 6 / 12.01; the density is right). Single-site phases (Pt, Zn, and in principle Ni, Si) are correct. The
  Monte Carlo step presumably uses these values, so the energy weighting of multi-site phases is off; band geometry is not affected.
- Before using a new phase, check `mp.metadata['cell']['average_atomic_number']`. Not yet reported upstream.

### 15.4 The map's physical step is not the nominal one (measured with a fiducial grid; important for GND)
- Both Apr 2025 sessions were registered to the SECCM deposit grid: 5 um pitch, axes 89.3-90.2 deg apart and scales isotropic within 2%
  in a 0 deg-tilt SEM image.

| session | pixel along columns | pixel along rows | angle between map axes |
|---|---|---|---|
| 16 Apr | 0.553 um (1.01x nominal) | 0.451 um (0.83x) | 86.2 deg |
| 18 Apr | 0.547 um (1.00x) | 0.464 um (0.85x) | 86.4 deg |

- So the rows, which run along the tilt direction, cover 15-17% less sample than nominal, and the axes are ~3.7 deg from orthogonal.
- If the scan's tilt correction assumes 70 deg, that row scale corresponds to a surface tilt of ~66 deg. The PC/elevation fits cannot tell
  (they see only alpha, 15.2), so a sample at ~66 deg with a ~7 deg camera would match the same patterns as 70 + 11 deg. **Not
  established**; it would also bear on what Deba's 11.04 deg Rcorr represents physically.
- Consequences for hrEBSD if the same holds for the 2026 sessions:
  1. OpenXY's GND uses one nominal step for both directions (13.4: `XData(3)-XData(2)`), so terms with y-derivatives would be
     underestimated by ~15%.
  2. Map shapes, and grain sizes along y, are stretched ~1.2x.
  3. Overlaying EBSD maps on SEM images needs a full affine, not a scale factor.
- **To check:** map one fiducial grid (TEM grid, FIB marks, an SECCM or indent array) in a typical hrEBSD session and run
  `tools\register_grid.py`; it reports both pixel sizes and the axis angle. If rows come out at 0.83x again, consider a step correction in
  GND, and measure the true stage tilt.
- **Update (15.13):** checked with Deba's Pt hole maps and SEM images. The row scale is not a fixed instrument property: it varies between
  sessions (0.77-1.19x), and Jan 2026 is undistorted.

### 15.5 Mounting and SEM metadata
- The same sample was mounted 90 deg (16 Apr) and 180 deg (18 Apr) rotated in EBSD relative to the SEM images. Both are proper rotations,
  never a mirror.
- An FEI zoom series of 76 images is unreliable in its metadata:
  - StageX/Y jumped mid-series (stage re-home), by ~430 um;
  - ScanRotation/StageR changed in groups while the image orientation stayed the same (measured 0.0-0.6 deg; compensated);
  - pixel size varied between images of one series (6.3-7.8 nm).
- So: use per-image metadata, place things in um from the image centre, and verify geometry against features, not against stage
  coordinates.

### 15.6 Registration recipe (`tools\register_grid.py`)
1. **Dots.** Take local maxima of the high-passed map above 4 robust SDs (the per-pattern brightness map of 15.1, or any image).
2. **Dot lattice.** Fit the map's own dot lattice from nearest-neighbour vectors; this needs no prior on scale or rotation.
3. **Orientation search.** Index the fiducials on their own lattice. Score all 8 rotations/mirrors x all integer shifts on
   deposit/no-deposit agreement (agree - disagree - unexplained dots).
   - Missing fiducials, an extra point or a deliberately asymmetric marker make the answer unique. A periodic grid alone does not: whole-pitch
     shifts tie, and fitting only on dot positions picked the wrong shift and orientation here.
   - The margin to the runner-up is reported.
4. **Refinement.** Refine the affine by least squares on the matched dots.
5. **Validation.** Every absent fiducial must land on a blank spot.

Results:
- 16 Apr: rot90, 153/154 present fiducials matched, RMS 0.63 px, score 163 vs 133.
- 18 Apr: rot180, 178/190, RMS 0.79 px, 155 vs 143.
- 0 absent fiducials with a dot in either map.
- Same numbers as the project's own code (`HOPG:tools/ebsd_grid_coloc.py`). The SEM side (`HOPG:tools/seccm_sem_coloc.py`) uses the same
  idea.

Uses here: EBSD <-> SEM overlays, indent arrays, and the step check of 15.4.

### 15.7 An independent frame check from crystal shape (bears on the 90 deg azimuth, 12.5-12.6)
- Hexagonal Zn plates were imaged in SEM (edge directions) and indexed in EBSD (kikuchipy phi1). After registering the frames (15.6):
  phi1 - SEM angle - mount rotation = **27.6 deg (16 Apr) and 25.9 deg (18 Apr)** (mod 60). The mounts differed by 93 deg.
- A constant across two mounts means the relation between the scan/kikuchipy frame and the physical sample is fixed, not mount-dependent.
  That is the "instrument convention" reading of the 90 deg azimuth in 12.5.
- (Mod 60, the constant cannot separate a hexagon-vertex convention from a 90 deg frame rotation.)
- For cubic Pt the same test is possible with any orientation-defined feature visible in SEM: facets, etch pits, slip or twin traces, or a
  cleaved Si edge along <110>. Register SEM to EBSD with fiducials, then compare. This would test Rz(90) independently of the .ANG export.

### 15.8 Mixed patterns: unmixing with measured neighbour references (near boundaries, inclusions, thin overlayers)
- **Method** (`HOPG:tools/deposit_deconv.py`):
  - flat-field and high-pass the patterns (both linear);
  - take the mean of pure patterns on a ring 3-5 px from the site as the reference for the host phase;
  - fit each pattern as alpha*reference + beta*simulated(second phase), and combine the beta-weighted residuals into an extracted pattern
    of the second phase.
- **Measured references beat simulated ones.** Against an independent truth (SEM shape of 9-11 clean plates), the agreement rose from
  R = 0.90 with a simulated host to 0.99.
- **Controls are essential.**
  - Detection thresholds come from the same procedure run on host-only control pixels (99th percentile).
  - Template crosstalk exists: host-only pixels give "second-phase" angles within 5 deg of the host angle 22-24% of the time, against 17% for
    random angles.
  - Sensitivity comes from injecting a synthetic second component into real patterns. An equal-weight component 5-30 deg away is found in
    82-98% of cases; a 25% minor one in <= 29%.
- **For hrEBSD.** Near-boundary points are masked by a pixel margin (`--ref-margin`, `near_boundary` flag) and never separated (workflow
  section 7). Unmixing each near-boundary pattern against interior patterns of the two adjacent grains would give a **measured mixing
  fraction per point**. That makes a principled mask (HR-EBSD cross-correlation is biased by even a few % of a second pattern) and may
  recover points that the margin discards. Not implemented here.
- A related recent method is iterative simulated-pattern subtraction (Cios, Winkelmann et al., arXiv:2601.14155, Jan 2026).

### 15.9 Probe footprint vs step
- Zn plates of 0.08-0.2 um^2, mapped at 0.55 x 0.45 um pixels (20 kV, 70 deg, on graphite), were detected more steeply with size than point
  sampling predicts: less often than the hit probability below ~0.11 um^2, more often above ~0.14 um^2.
- So each pattern integrates roughly 0.1-0.2 um around the beam position for this light-element case; a high-Z phase like Pt is smaller.
- For hrEBSD: a 4 px `--ref-margin` at 0.27 um (1.1 um) is conservative against that scale. A measured mixing fraction (15.8) would show how
  far mixing actually reaches.

### 15.10 Statistics when correlating per-point quantities
- Acquisition-order drift is common (SECCM: open-circuit potential rho = 0.7 with hop order). Detrend against scan order (running median)
  before correlating.
- Use FDR (Benjamini-Hochberg) over the family of tests.
- Use permutation nulls for angle statistics (6-fold or 4-fold circular R).
- Use partial correlations when a confounder (feature size) drives both quantities.
- Report injection-test sensitivity alongside any "not detected".

### 15.11 Data access
- Elias's data: Box Drive `C:\Users\Hang\Box\Ren Lab Shared\Results and Data\Data\Elias\Zn Deposition on HOPG`. `robocopy` from Box Drive
  copied ~2 GB in under a minute.
- The Box connector's download links expire after 60 s and are single-use, so it is only practical for a few files.
- Box stays read-only.

### 15.12 Methods worth knowing (literature check 2026-11-10; details in the HOPG findings and the chat record)
- **Depth.** The EBSD information depth grows with kV and shrinks with Z and with tilt. Example: at 15 kV, 70 deg, >80% of contributing
  electrons scattered within ~40 nm (material-dependent; a textbook figure).
  - Orientation was still measurable through ~116 nm of amorphous Si at 30 kV (Fraunhofer wedge study).
  - A kV series (with Monte Carlo, e.g. ebsdsim's MC or CASINO, to interpret it) gives coarse depth weighting, not tomography. True depth
    needs FIB serial sectioning or cross-section TKD.
- **Resolution.** Low kV (5-7.5 kV) shrinks the interaction volume (~10 nm reported for steel at 7.5 kV). Energy-filtered EBSD with
  direct electron detectors (event counting, filtering on the primary-energy window) sharpens patterns and is reported to improve
  resolution (Della Ventura et al., Ultramicroscopy 2025).
- **Overlapping patterns.**
  - Dictionary indexing and spherical indexing (EMSphInx) are robust to noise and can handle overlaps near boundaries.
  - Constrained NMF (Chauniyal et al.) and iterative subtraction (Cios et al. 2026) separate overlaps.
- **HR-EBSD.**
  - Integrated/global DIC (Vermeij & Hoefnagels; Ernould et al.) and simulation-based references (Alkorta et al.) avoid the
    reference-pattern problem.
  - Global-DIC refinement of PC and sample/detector angles (Griesbach & Kochmann, arXiv:2604.25869, 2026) targets exactly the
    PC/tilt coupling of 12.19 and 15.2.
- **Denoising.** Autoencoders and diffusion models restore noisy patterns. All of them can hallucinate, so validate against controls
  before trusting orientations.

### 15.13 The map geometry varies per session (Deba's Pt hole maps, 2026-11-10)
No fiducial grid has been mapped in an hrEBSD session yet, so three methods were used on existing data. The data are 45 `.ANG` maps of
Deba's Pt foil holes (Box `G7. Deba\EBSD`, 2024-2026, copies in `inputs\Deba_hole_maps`) and 97 top-view SEM images of the holes (Box
`G7. Deba\SEM + EDX`, copies in `inputs\Deba_hole_SEM`). The SEM images are all 0 deg tilt with tilt correction off, so their pixels are
square on the sample.
1. **Grain-boundary registration** (`tools\register_sem_ebsd.py`).
   - The EBSD boundary map (neighbour misorientation, cubic, plus the unindexed hole edge) is registered to the SEM edge map with a full
     affine: rotation and row-scale search, masked NCC for the shift, then Powell on all six parameters.
   - The fit is repeated on map quadrants for an uncertainty, and a result counts as "reliable" only if the peak beats other rotations by
     > 0.1 and the quadrants agree.
   - **Self-test**: a known distortion (rows 0.83x, axes 86.3 deg, rotation 37 deg) is recovered to 0.001 / 0.02 deg.
2. **Round features** (`tools\round_feature_check.py`): from the second moments of a round hole in IQ, the row/column pixel ratio and the
   axis angle follow.
   - The holes are round in SEM: major/minor 1.00-1.03 for 7 holes (Sep/Nov 2024), 1.04 (Feb 2025 hole 1), 1.07 (Jan 2025 hole 5); one
     outlier at 1.16 (Nov 2024 hole 5).
   - Segmentation fails where a dark grain touches the hole.
3. **Contamination dots** left by an earlier scan, measured in SEM (`tools\scan_dot_lattice.py`, experimental): the dot lattice gives the
   true step along both scan axes.

**Jan 2026: three methods agree, no distortion.**

| method | rows/columns | axis angle |
|---|---|---|
| registration, map 1 vs SEM 13 Jan (NCC 0.52, reliable) | 1.025 +- 0.034 | 88.6 +- 1.7 deg |
| dot lattice, SEM 15 Jan | 1.020 | 90.0 deg |
| holes, 2 maps | 0.98 | |

The registration's pixel sizes were columns 1.238 um and rows 1.269 um, against a nominal 1.260. The dot lattice was 9.91 x 9.71 um, with
its axes within 0.4-1.7 deg of the EBSD scan axes after the -90 deg stage rotation, so it was left by an EBSD scan on that mount.

**Per-session row scale from round holes** (median per map, then per session), with the surface tilt that would produce it:

| session | maps | rows/columns | surface tilt implied |
|---|---|---|---|
| 2026-05-06 repeat calibration | 5 | 0.77 (0.73-0.80) | 63.5 deg |
| 2025-02-19 | 2 | 0.79 | 64.3 deg |
| 2025-07-02 | 3 | 0.80 | 64.8 deg |
| 2025-01-22 | 2 | 0.82 | 65.3 deg |
| 2024-11-01, 2024-11-22, 2025-01-10, 2025-01-16, 2025-05-06, 2026-01 | 2-6 each | 0.95-1.04 | 69-71 deg |
| 2025-01-08 | 6 | 1.07 | 71.4 deg |
| 2024-12-06 | 2 | 1.13 | 72.4 deg |
| 2025-04-25 | 4 | 1.19 | 73.2 deg |

- The spread within a session is small (e.g. 0.73-0.80 over 5 maps of 2026-05-06), so the scale is a property of the session or mount.
- The 2025 Zn/HOPG sessions (15.4: 0.83 / 0.85) fit the same picture.
- **Registration on a compressed session.** Feb 2025, map `Pt_hole1(1)` (19 Feb) vs SEM `hole-1` (21 Feb, stage rotated 90 deg): rows/columns
  **0.846**, axes 89.6 deg, NCC 0.53. Boundaries and the hole rim overlay visually (`work\reg3_feb2025_h1m1\register_sem_ebsd.png`). The same
  map's round-hole estimate is 0.81; the hole is 1.04 out of round in SEM.
  - It is flagged not reliable only because the SEM field is narrower than the map, so no quadrant has the 60% coverage needed for an
    uncertainty.
  - So grain boundaries alone give 0.85 here, against 1.03 in Jan 2026: the session dependence is real, not a hole-shape artefact.
  - The Jan 2025 hole-5 pairs did not register (partial overlap). The coarse search needs the row scale as a search dimension: without
    it, compressed sessions gave no peak.

**Model.**
- **Row scale.** The scan's tilt correction assumes 70 deg. With the surface really at theta, the row step on the sample is
  nominal x cos 70 / cos theta. That is ~5% per degree near 70 deg, so a foil sitting 4-6 deg off (tape, bending, holder) gives
  0.8-0.85x or 1.1-1.2x.
- **Shear.** If the stage tilt axis is rotated by phi from the scan x-axis, the tilt correction shears the map axes by about
  sin^2(theta)/cos(theta) x phi = 2.6 phi at 70 deg. The ~3.7 deg shear of the 2025 maps corresponds to phi ~ 1.4 deg.
- **Pattern geometry.** The same surface-tilt deviation enters the pattern geometry: alpha = 90 - sample tilt + camera elevation (15.2).
  - For Apr 2025, the PC/tilt fits (alpha 31.4 / 31.1 deg) with the measured surface tilt (65.5 / 66.3 deg) give a camera elevation of
    6.9 / 7.4 deg; the same within 0.5 deg.
  - That would mean the ~11 deg found with the sample assumed at 70 deg is camera elevation plus mount tilt. **Not established**; it needs
    the pilot below.
- **For hrEBSD**, if the model holds:
  1. GND needs the session's real row step.
  2. Deba's Rcorr (11.04 deg, calibrated 2026-05-05, the day before the most compressed session) may be partly mount-specific; applied
     to a differently mounted sample, the orientation correction could be off by the mount difference (several deg).
  3. The sample frame of the strain/rotation results is tilted by the same amount.

**Next pilot: per-session geometry calibration from data that already exist.** For sessions that have both an EBSD map with raw patterns
and a top-view SEM image of the same area:
1. Measure the row scale and shear (register_sem_ebsd or round features), giving theta and phi.
2. Fit the PC and alpha on many grains (`map_pc_fit.py --free-elev`).
3. Compute the camera elevation = alpha - 90 + theta.

If it comes out constant across mounts, each session's geometry and Rcorr can be built from that constant plus the measured theta. Then
rerun Pt(6)/(7) with session-specific corrections and compare strain and rotation.
