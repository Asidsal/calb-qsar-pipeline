"""Benchmark or run CAL-028 in 1 ns segments on one CUDA GPU."""
import argparse
import csv
import hashlib
import json
import math
import time
from pathlib import Path
from zipfile import ZipFile

import openmm as mm
from openmm import app, unit as u

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--ns", type=int, choices=(1, 50, 100), default=1,
                    help="total ns from the equilibrated state (default 1)")
args = parser.parse_args()

archive = Path("CAL028_CORRECTED_AMBER.zip")
initial_state = Path("state.xml")
with ZipFile(archive) as bundle:
    topology = bundle.read("cal028.prmtop")
if hashlib.sha256(archive.read_bytes()).hexdigest() != "25507e5910b207ea7a28be67ec34f18f530e6dfff4dfab85d0dc0c2d3948f06b":
    raise ValueError("Use the verified CAL028_CORRECTED_AMBER.zip from the pilot.")
Path("cal028.prmtop").write_bytes(topology)
top = app.AmberPrmtopFile("cal028.prmtop")
system = top.createSystem(nonbondedMethod=app.PME,
                          nonbondedCutoff=1*u.nanometer, constraints=app.HBonds)
barostat = mm.MonteCarloBarostat(1*u.bar, 300*u.kelvin, 25)
barostat.setRandomNumberSeed(20261124)
system.addForce(barostat)
integrator = mm.LangevinMiddleIntegrator(300*u.kelvin, 1/u.picosecond, 0.002*u.picoseconds)
integrator.setRandomNumberSeed(20261123)
sim = app.Simulation(top.topology, system, integrator,
                     mm.Platform.getPlatformByName("CUDA"), {"Precision": "mixed"})

atoms = list(top.topology.atoms())
fe, = (a for a in atoms if a.name == "FE" and a.residue.name == "HEM")
neighbors = {b.name: b.index for a,b in top.topology.bonds() if a == fe}
neighbors.update({a.name: a.index for a,b in top.topology.bonds() if b == fe})
if set(neighbors) != {"SG", "NA", "NB", "NC", "ND"} or len(atoms) != 75277:
    raise ValueError("Unexpected heme connectivity or atom count.")

def check_state():
    state = sim.context.getState(getPositions=True, getEnergy=True)
    pos = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
    e = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
    distances = {name: math.dist(pos[fe.index], pos[idx]) for name,idx in neighbors.items()}
    if not math.isfinite(e) or not all(map(math.isfinite, distances.values())):
        raise ValueError("Nonfinite energy or heme distance.")
    if not (2.1 <= distances["SG"] <= 3.2 and
            all(1.7 <= distances[n] <= 2.5 for n in ("NA", "NB", "NC", "ND"))):
        raise ValueError(f"Heme coordination outside the pilot screening range: {distances}")
    print("Fe distances (Å):", {k:round(v,3) for k,v in distances.items()}, flush=True)

out = Path("CAL028_5090_results")
out.mkdir(exist_ok=True)
metadata = {
    "openmm": mm.__version__,
    "topology_sha256": hashlib.sha256(topology).hexdigest(),
    "state_sha256": hashlib.sha256(initial_state.read_bytes()).hexdigest(),
    "temperature_K": 300, "pressure_bar": 1, "dt_fs": 2,
    "precision": "mixed", "steps_per_segment": 500000,
    "report_every_steps": 5000,
}
manifest = out/"protocol.json"
if manifest.exists() and json.loads(manifest.read_text()) != metadata:
    raise ValueError("The existing results use different inputs or a different OpenMM version.")
manifest.write_text(json.dumps(metadata, indent=2))

done = sorted(out.glob("ns_[0-9][0-9][0-9]"))
expected = [f"ns_{i:03}" for i in range(1,len(done)+1)]
if [p.name for p in done] != expected:
    raise ValueError("Completed segments must form a contiguous sequence.")
if done:
    sim.loadCheckpoint(str(done[-1]/"restart.chk"))
    print(f"Resumed from {done[-1]}", flush=True)
else:
    sim.loadState(str(initial_state))
    state_time = sim.context.getTime().value_in_unit(u.picosecond)
    if sim.currentStep != 550000 or not math.isclose(state_time,1100,abs_tol=0.001):
        raise ValueError("Use the 03_npt/state.xml equilibration endpoint (step 550000).")
    sim.currentStep = 0
    sim.context.setTime(0)
    check_state()
    print("Running a 10 ps warm-up before the timed segment.", flush=True)
    sim.step(5000)
    sim.currentStep = 0
    sim.context.setTime(0)

for i in range(len(done)+1,args.ns+1):
    directory = out/f"ns_{i:03}"
    work = directory.with_suffix(".partial")
    work.mkdir(exist_ok=True)
    sim.reporters = [
        app.DCDReporter(str(work/"trajectory.dcd"),5000,enforcePeriodicBox=False),
        app.StateDataReporter(str(work/"state.csv"),5000,step=True,time=True,
                              potentialEnergy=True,temperature=True,density=True),
    ]
    ticks = 500000
    start = time.monotonic()
    sim.step(ticks)
    check_state()
    elapsed = time.monotonic()-start
    sim.reporters.clear()
    sim.saveCheckpoint(str(work/"restart.chk"))
    sim.saveState(str(work/"state.xml"))
    work.rename(directory)
    actual_ns = ticks * 0.002 / 1000
    throughput = actual_ns*86400/elapsed
    with (out/"timing.csv").open("a") as handle:
        csv.writer(handle).writerow((i, actual_ns, elapsed, throughput))
    print(f"Saved {i} ns | {throughput:.1f} ns/day | "
          f"1 ns equivalent: {86400/throughput:.1f} s | "
          f"100 ns estimate: {2400/throughput:.2f} h", flush=True)
