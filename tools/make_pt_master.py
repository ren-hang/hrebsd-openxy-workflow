"""Dynamical Pt (Fm-3m, a = 3.924 Angstrom; ebsdsim cell units are Angstrom) master pattern at 20 kV with ebsdsim (WebGPU). usage: make_pt_master.py [halfw=250]"""
import sys, time, ebsdsim as es
halfw = int(sys.argv[1]) if len(sys.argv) > 1 else 250
pt = es.Material(cell=es.Cell(a=3.924, b=3.924, c=3.924, space_group="Fm-3m"), atoms=[es.Atom("Pt", x=0.0, y=0.0, z=0.0)], name="Pt")
t = time.time(); mp = es.master_pattern(pt, voltage_kv=20.0, halfw=halfw); print("done in %.1f s" % (time.time() - t), mp.pattern.shape)
mp.save("work/pt_master/Pt_20kV_hw%d.npz" % halfw)
