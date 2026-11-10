# hrEBSD with OpenXY on Thermo Pathfinder maps

One-command workflow (check, prepare, run OpenXY headless, post-process) for `.pat` + `.ANG` maps, with a small GUI.

- How to use it, data organisation and file naming: `docs/hrEBSD_workflow.md`
- Why each setting is what it is (validation, noise floors, failures): `NOTES.md`
- Code: `tools/` (`hrebsd.py` is the entry point, `hrebsd_gui.py` the GUI, `openxy_patches/` two local fixes to OpenXY)
- Python packages: `docs/requirements.txt`. MATLAB R2024a and a local OpenXY copy (BYU master + the two patches) are required.
- General EBSD lessons from other projects: `NOTES.md` section 15. It covers .pat details, PC fitting without indexing, the measured map step, registration (`tools/register_grid.py`), unmixing of mixed patterns, and methods to consider. The Zn-on-HOPG project itself is kept separate in `ren-hang/hopg-zn-orientation`.

Lab data (`.pat`, `.ANG`, patterns) and result folders are not in the repository.
