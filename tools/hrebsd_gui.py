"""Small GUI around tools\\hrebsd.py (no new logic: it builds the same command line and shows its output).
Start from the repository root:  work\\kpenv\\Scripts\\python.exe tools\\hrebsd_gui.py   (or double-click run_gui.bat)"""
import os, queue, subprocess, sys, threading, tkinter as tk
from tkinter import filedialog, messagebox, ttk

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PY = os.path.join(ROOT, "work", "kpenv", "Scripts", "python.exe")


class App(tk.Tk):
    def __init__(self):
        super().__init__(); self.title("hrEBSD / OpenXY workflow"); self.geometry("820x640"); self.q = queue.Queue(); self.proc = None
        f = ttk.Frame(self, padding=8); f.pack(fill="x")
        self.v = {k: tk.StringVar(value=d) for k, d in dict(map="", name="", X="", Y="", L="", block="", material="platinum", phase="Platinum").items()}
        self.gnd = tk.BooleanVar(value=False)
        ttk.Label(f, text="Map folder (*.MAP.XTAL copy)").grid(row=0, column=0, sticky="w"); ttk.Entry(f, textvariable=self.v["map"], width=70).grid(row=0, column=1, columnspan=5, sticky="we")
        ttk.Button(f, text="Browse", command=self.browse).grid(row=0, column=6, padx=4)
        ttk.Label(f, text="Run name").grid(row=1, column=0, sticky="w"); ttk.Entry(f, textvariable=self.v["name"], width=16).grid(row=1, column=1, sticky="w")
        ttk.Label(f, text="PC (Pathfinder mm)  X").grid(row=2, column=0, sticky="w")
        for i, k in enumerate("XYL"):
            ttk.Label(f, text=k if k != "X" else "").grid(row=2, column=1 + 2 * i, sticky="e"); ttk.Entry(f, textvariable=self.v[k], width=8).grid(row=2, column=2 + 2 * i, sticky="w")
        ttk.Label(f, text="Material / phase").grid(row=3, column=0, sticky="w"); ttk.Entry(f, textvariable=self.v["material"], width=16).grid(row=3, column=1, sticky="w"); ttk.Entry(f, textvariable=self.v["phase"], width=16).grid(row=3, column=2, columnspan=2, sticky="w")
        ttk.Label(f, text="Quick test: N x N block (blank = whole map)").grid(row=4, column=0, sticky="w"); ttk.Entry(f, textvariable=self.v["block"], width=6).grid(row=4, column=1, sticky="w")
        ttk.Checkbutton(f, text="Dislocation density (GND, +3 min)", variable=self.gnd).grid(row=4, column=2, columnspan=4, sticky="w")
        b = ttk.Frame(self, padding=(8, 0)); b.pack(fill="x")
        ttk.Button(b, text="1. Check only", command=lambda: self.go("check")).pack(side="left", padx=2)
        ttk.Button(b, text="2. Run everything", command=lambda: self.go("all")).pack(side="left", padx=2)
        ttk.Button(b, text="Re-make plots only (post)", command=lambda: self.go("post")).pack(side="left", padx=2)
        ttk.Button(b, text="Open results folder", command=self.open_results).pack(side="left", padx=2)
        ttk.Button(b, text="Stop", command=self.stop).pack(side="right", padx=2)
        self.log = tk.Text(self, wrap="word", height=26, font=("Consolas", 9)); self.log.pack(fill="both", expand=True, padx=8, pady=8)
        self.after(100, self.pump)

    def browse(self):
        d = filedialog.askdirectory(title="Select the *.MAP.XTAL folder", initialdir=os.path.join(ROOT, "inputs"))
        if d:
            self.v["map"].set(d)
            if not self.v["name"].get(): self.v["name"].set("".join(c for c in os.path.basename(d).split(".")[0] if c.isalnum()))

    def cmd(self, stage):
        v = {k: x.get().strip() for k, x in self.v.items()}
        if not v["map"] or not v["name"]: raise ValueError("Choose the map folder and give the run a name.")
        try: pc = [float(v[k]) for k in "XYL"]
        except ValueError: raise ValueError("Enter the pattern centre X, Y, L (mm) from the acquisition report.")
        c = [PY, "-I", os.path.join("tools", "hrebsd.py"), v["map"], "--pc", *map(str, pc), "--name", v["name"], "--material", v["material"], "--phase", v["phase"], "--stage", stage]
        if v["block"]: c += ["--block", v["block"]]
        if self.gnd.get(): c.append("--gnd")
        return c

    def go(self, stage):
        if self.proc and self.proc.poll() is None: return messagebox.showinfo("Busy", "A run is already in progress.")
        try: c = self.cmd(stage)
        except ValueError as e: return messagebox.showerror("Missing input", str(e))
        self.log.delete("1.0", "end"); self.say("> " + " ".join(c) + "\n")
        threading.Thread(target=self.work, args=(c,), daemon=True).start()

    def work(self, c):
        self.proc = subprocess.Popen(c, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
        for line in self.proc.stdout: self.q.put(line)
        self.q.put("\n--- finished (exit code %s) ---\n" % self.proc.wait())

    def say(self, t): self.log.insert("end", t); self.log.see("end")

    def pump(self):
        try:
            while True: self.say(self.q.get_nowait())
        except queue.Empty: pass
        self.after(100, self.pump)

    def stop(self):
        if self.proc and self.proc.poll() is None:
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(self.proc.pid)], capture_output=True); self.say("\n--- stopped ---\n")

    def open_results(self):
        d = os.path.join(ROOT, "work", "hrebsd_" + self.v["name"].get().strip())
        os.startfile(d) if os.path.isdir(d) else messagebox.showinfo("Not yet", "No results folder yet: " + d)


if __name__ == "__main__":
    App().mainloop()
