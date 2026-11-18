import fcntl
import json
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit as u

workspace = Path("/workspace")
inputs = workspace / "CAL-015/N14_frame_corrected"
output = workspace / "CAL-015/production_50ns_1fs"
output.mkdir(exist_ok=True)

with (output / "queue.lock").open("a") as queue_lock, open("/tmp/cal028-production.lock", "a") as gpu_lock:
    fcntl.flock(queue_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    for name in ("system.prmtop", "system.rst7", "minimized_state.xml", "minimization.json"):
        assert (inputs / name).stat().st_size > 0, name
    assert not (output / "npt_state.xml").exists(), "NPT state already saved"
    print("CAL-015 staged. Waiting for CAL-028 to release the GPU lock.", flush=True)
    fcntl.flock(gpu_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    complete = any(line.startswith("Saved 50 ns |") for line in
                   (workspace / "CAL028_50ns.log").read_text().splitlines())
    assert complete and (workspace / "CAL028_5090_results/ns_050/restart.chk").stat().st_size > 0, "CAL-028 released the lock without a verified 50 ns completion"

    top = app.AmberPrmtopFile(str(inputs / "system.prmtop"))
    crd = app.AmberInpcrdFile(str(inputs / "system.rst7"))
    start = mm.XmlSerializer.deserialize((workspace / "CAL-015/npt_retry_1fs/npt_state.xml").read_text())
    system = top.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1*u.nanometer, constraints=app.HBonds)
    atoms = list(top.topology.atoms())
    ligand = [a.index for a in atoms if a.residue.name == "LIG"]
    heavy = [i for i in ligand if atoms[i].element.atomic_number > 1]
    ca = [a.index for a in atoms if a.name == "CA" and a.residue.name != "LIG"]
    assert len(heavy) == 32 and len(ca) == 484
    nb = next(f for f in system.getForces() if isinstance(f, mm.NonbondedForce))
    charges = np.array([nb.getParticleParameters(i)[0].value_in_unit(u.elementary_charge) for i in range(len(atoms))])
    assert abs(charges.sum()) < 0.002 and abs(charges[ligand].sum()-1) < 0.002
    fe, = [a for a in atoms if a.name == "FE" and a.residue.name == "HEM"]
    neighbors = {b if a == fe else a for a,b in top.topology.bonds() if fe in (a,b)}
    assert {a.name for a in neighbors} == {"SG","NA","NB","NC","ND"}
    barostat = mm.MonteCarloBarostat(1*u.bar, 300*u.kelvin, 25)
    barostat.setRandomNumberSeed(15039)
    system.addForce(barostat)
    integ = mm.LangevinMiddleIntegrator(300*u.kelvin, 1/u.picosecond, 0.001*u.picoseconds)
    integ.setRandomNumberSeed(15038)
    sim = app.Simulation(top.topology, system, integ, mm.Platform.getPlatformByName("CUDA"), {"Precision":"mixed"})
    assert len(start.getPositions()) == system.getNumParticles() == len(crd.positions)
    sim.context.setPeriodicBoxVectors(*start.getPeriodicBoxVectors())
    sim.context.setPositions(start.getPositions())
    sim.context.applyConstraints(1e-5)
    sim.context.setVelocities(start.getVelocities())
    reference = start.getPositions(asNumpy=True).value_in_unit(u.angstrom)

    def inspect():
        state = sim.context.getState(getPositions=True, getEnergy=True)
        x = state.getPositions(asNumpy=True).value_in_unit(u.angstrom)
        e = state.getPotentialEnergy().value_in_unit(u.kilojoule_per_mole)
        assert np.isfinite(e) and np.isfinite(x).all()
        distances = {a.name: float(np.linalg.norm(x[a.index]-x[fe.index])) for a in neighbors}
        assert 2.1 <= distances["SG"] <= 3.2 and all(1.7 <= distances[n] <= 2.5 for n in ("NA","NB","NC","ND")), distances
        center = x[ca].mean(0)
        ref_center = reference[ca].mean(0)
        left, _, right = np.linalg.svd((x[ca]-center).T @ (reference[ca]-ref_center))
        handed = np.diag([1, 1, np.linalg.det(left @ right)])
        aligned = (x[heavy]-center) @ (left @ handed @ right) + ref_center
        rmsd = float(np.sqrt(np.mean(np.sum((aligned-reference[heavy])**2, axis=1))))
        return {"step":sim.currentStep, "potential_kJ_mol":e, "Fe_distances_A":distances, "ligand_RMSD_A_from_minimized":rmsd}


    done = sorted(output.glob("ns_[0-9][0-9][0-9]"))
    assert [p.name for p in done] == [f"ns_{i:03d}" for i in range(1,len(done)+1)]
    if done:
        sim.loadCheckpoint(str(done[-1]/"restart.chk"))
        assert sim.currentStep == 1000000*len(done)
    for i in range(len(done)+1,51):
        dest = output/f"ns_{i:03d}"
        partial = dest.with_suffix(".partial")
        partial.mkdir(exist_ok=True)
        sim.reporters = [
            app.DCDReporter(str(partial/"trajectory.dcd"), 10000, enforcePeriodicBox=False),
            app.StateDataReporter(str(partial/"state.csv"), 10000, step=True, time=True,
                                  potentialEnergy=True, temperature=True, density=True),
        ]
        sim.step(1000000)
        sim.reporters.clear()
        report = inspect()
        sim.saveCheckpoint(str(partial/"restart.chk"))
        sim.saveState(str(partial/"state.xml"))
        (partial/"screen.json").write_text(json.dumps(report,indent=2))
        partial.rename(dest)
        print(f"Saved {i} ns | {json.dumps(report)}",flush=True)
