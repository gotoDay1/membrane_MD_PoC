"""
POC: POPC二重膜だけを OpenMM で組んで回す (CPUラップトップ想定)

環境:  pip install openmm   (または conda install -c conda-forge openmm)
実行:  python popc_bilayer_poc.py            # 既定 0.1 ns
       python popc_bilayer_poc.py --ns 1.0
力場:  Amber Lipid17 (amber14-all.xml に同梱) + TIP3P, HMR 4 fs
出力:  bilayer_init.pdb, bilayer_final.pdb, traj.dcd, log.csv(APL・膜厚の時系列)
判定:  APL(1脂質あたり面積) ≈ 0.62–0.68 nm², P–P膜厚 ≈ 3.7–4.0 nm (POPC, 303 K)
"""
import argparse
import sys
import time

import numpy as np
from openmm import Vec3, app, unit, LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform

ap = argparse.ArgumentParser()
ap.add_argument("--ns", type=float, default=0.1, help="本計算の長さ [ns]")
ap.add_argument("--pad", type=float, default=2.5, help="膜パッチの半幅の目安 [nm] (2.5で約5 nm角)")
args = ap.parse_args()

TEMP = 303 * unit.kelvin
DT = 4 * unit.femtoseconds
REPORT_PS = 10

# ---- 1. 系の構築 ----
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
# addMembrane は溶質を必須とするため、水1分子を「芯」として原点に置く
# (膜の中心に入る分子なので、addMembraneによって除去されるか水として残るだけで影響は無視できる)
seed = app.Topology()
chain = seed.addChain()
res = seed.addResidue("HOH", chain)
o = seed.addAtom("O", app.element.oxygen, res)
h1 = seed.addAtom("H1", app.element.hydrogen, res)
h2 = seed.addAtom("H2", app.element.hydrogen, res)
seed.addBond(o, h1)
seed.addBond(o, h2)
seed_pos = [Vec3(0, 0, 0), Vec3(0.09572, 0, 0), Vec3(-0.024, 0.0927, 0)] * unit.nanometer
modeller = app.Modeller(seed, seed_pos)
t0 = time.time()
modeller.addMembrane(ff, lipidType="POPC",
                     minimumPadding=args.pad * unit.nanometer,
                     positiveIon="K+", negativeIon="Cl-",
                     ionicStrength=0.15 * unit.molar)
app.PDBFile.writeFile(modeller.topology, modeller.positions, open("bilayer_init.pdb", "w"))

atoms = list(modeller.topology.atoms())
p_idx = [a.index for a in atoms if a.element is not None and a.element.symbol == "P"]
n_lipid = len(p_idx)
print(f"原子数: {len(atoms)}  POPC: {n_lipid}  構築: {time.time() - t0:.1f} s")

# ---- 2. System ----
system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds, hydrogenMass=4 * unit.amu)
system.addForce(MonteCarloMembraneBarostat(
    1 * unit.bar, 0 * unit.bar * unit.nanometer, TEMP,
    MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
integrator = LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)

platform = max((Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())),
               key=lambda p: p.getSpeed())
print("Platform:", platform.getName())
sim = app.Simulation(modeller.topology, system, integrator, platform)
sim.context.setPositions(modeller.positions)

print("エネルギー最小化...")
sim.minimizeEnergy(maxIterations=2000)
sim.context.setVelocitiesToTemperature(TEMP)


def bilayer_metrics(state):
    """APL [nm²] と P–P膜厚 [nm] を返す"""
    box = state.getPeriodicBoxVectors(asNumpy=True).value_in_unit(unit.nanometer)
    pos = state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    area = box[0][0] * box[1][1]
    z = pos[p_idx, 2]
    zc = z.mean()
    thick = z[z > zc].mean() - z[z <= zc].mean()
    return area / (n_lipid / 2), thick


# ---- 3. 本計算 ----
steps_per_report = int(REPORT_PS * unit.picosecond / DT)
n_reports = max(1, int(args.ns * 1000 / REPORT_PS))
sim.reporters.append(app.DCDReporter("traj.dcd", steps_per_report))
sim.reporters.append(app.StateDataReporter(sys.stdout, steps_per_report * 5,
                     time=True, potentialEnergy=True, temperature=True,
                     speed=True, remainingTime=True,
                     totalSteps=steps_per_report * n_reports))

with open("log.csv", "w") as f:
    f.write("time_ps,apl_nm2,thickness_nm\n")
    for i in range(n_reports):
        sim.step(steps_per_report)
        apl, th = bilayer_metrics(sim.context.getState(getPositions=True))
        f.write(f"{(i + 1) * REPORT_PS},{apl:.4f},{th:.3f}\n")
        f.flush()

state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
app.PDBFile.writeFile(modeller.topology, state.getPositions(), open("bilayer_final.pdb", "w"))
apl, th = bilayer_metrics(state)
print(f"最終: APL = {apl:.3f} nm²  P–P膜厚 = {th:.2f} nm")
