"""
POPC二重膜 MD(OpenMM, CPUラップトップ想定)

  平衡化: python run_md.py equil --ns 2      # 構築 → 最小化 → NPT (約126脂質, 約3.2万原子)
  本計算: python run_md.py prod  --ns 8      # equil の続きから DCD を出力
  再開:   python run_md.py prod  --ns 8 --resume   # prod.chk から続行(中断対策)

力場: Amber Lipid17 (amber14-all.xml) + TIP3P, 150 mM KCl, 303 K, 1 bar
      HMR(H質量 4 amu)+ 4 fs, PME 1.0 nm, MonteCarloMembraneBarostat(xy等方・z独立)
"""
import argparse
import os
import sys
import time

from openmm import (Vec3, app, unit, LangevinMiddleIntegrator,
                    MonteCarloMembraneBarostat, Platform, XmlSerializer)

TEMP = 303 * unit.kelvin
DT = 4 * unit.femtoseconds
FRAME_PS = 10          # DCD・ログの出力間隔
CHK_PS = 100           # チェックポイント間隔

ap = argparse.ArgumentParser()
ap.add_argument("stage", choices=["equil", "prod"])
ap.add_argument("--ns", type=float, default=None, help="計算長 [ns](既定: equil 2, prod 8)")
ap.add_argument("--pad", type=float, default=1.0,
                help="膜パッチの大きさの目安 [nm]。addMembraneの最小パッチは約126脂質(2.5以下では変わらない)")
ap.add_argument("--threads", type=int, default=None, help="CPUスレッド数(既定: 全コア)")
ap.add_argument("--resume", action="store_true", help="prod.chk から再開")
args = ap.parse_args()
ns = args.ns if args.ns is not None else (2.0 if args.stage == "equil" else 8.0)

ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")


def build_system():
    # addMembrane は溶質が必須のため、水1分子を芯として原点に置く
    seed = app.Topology()
    res = seed.addResidue("HOH", seed.addChain())
    o = seed.addAtom("O", app.element.oxygen, res)
    h1 = seed.addAtom("H1", app.element.hydrogen, res)
    h2 = seed.addAtom("H2", app.element.hydrogen, res)
    seed.addBond(o, h1)
    seed.addBond(o, h2)
    pos = [Vec3(0, 0, 0), Vec3(0.09572, 0, 0), Vec3(-0.024, 0.0927, 0)] * unit.nanometer
    m = app.Modeller(seed, pos)
    m.addMembrane(ff, lipidType="POPC", minimumPadding=args.pad * unit.nanometer,
                  positiveIon="K+", negativeIon="Cl-", ionicStrength=0.15 * unit.molar)
    return m.topology, m.positions


def make_simulation(topology):
    system = ff.createSystem(topology, nonbondedMethod=app.PME,
                             nonbondedCutoff=1.0 * unit.nanometer,
                             constraints=app.HBonds, hydrogenMass=4 * unit.amu)
    system.addForce(MonteCarloMembraneBarostat(
        1 * unit.bar, 0 * unit.bar * unit.nanometer, TEMP,
        MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree))
    integrator = LangevinMiddleIntegrator(TEMP, 1 / unit.picosecond, DT)
    platform = max((Platform.getPlatform(i) for i in range(Platform.getNumPlatforms())),
                   key=lambda p: p.getSpeed())
    props = {}
    if args.threads and platform.getName() == "CPU":
        props["Threads"] = str(args.threads)
    print("Platform:", platform.getName(), props)
    return app.Simulation(topology, system, integrator, platform, props)


def run(sim, prefix, nsteps, append=False):
    frame = int(FRAME_PS * unit.picosecond / DT)
    sim.reporters.append(app.DCDReporter(f"{prefix}.dcd", frame, append=append))
    sim.reporters.append(app.StateDataReporter(
        f"{prefix}_log.csv", frame, append=append, step=True, time=True,
        potentialEnergy=True, temperature=True, volume=True, speed=True))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, frame * 10, time=True, temperature=True, speed=True,
        remainingTime=True, totalSteps=sim.currentStep + nsteps))
    sim.reporters.append(app.CheckpointReporter(f"{prefix}.chk",
                         int(CHK_PS * unit.picosecond / DT)))
    sim.step(nsteps)
    sim.saveCheckpoint(f"{prefix}.chk")


t0 = time.time()
if args.stage == "equil":
    top, pos = build_system()
    n_lip = sum(1 for r in top.residues() if r.name == "POP")
    print(f"原子数: {top.getNumAtoms()}  POPC: {n_lip}")
    app.PDBFile.writeFile(top, pos, open("init.pdb", "w"))
    sim = make_simulation(top)
    sim.context.setPositions(pos)
    print("エネルギー最小化...")
    sim.minimizeEnergy(maxIterations=5000)
    sim.context.setVelocitiesToTemperature(TEMP)
    run(sim, "equil", int(ns * 1e6 * unit.femtoseconds / DT))
    state = sim.context.getState(getPositions=True, getVelocities=True, enforcePeriodicBox=True)
    open("equil_state.xml", "w").write(XmlSerializer.serialize(state))
    app.PDBFile.writeFile(top, state.getPositions(), open("equil.pdb", "w"))

else:
    top = app.PDBFile("equil.pdb").topology
    sim = make_simulation(top)
    if args.resume and os.path.exists("prod.chk"):
        sim.loadCheckpoint("prod.chk")
        done = sim.currentStep
        print(f"prod.chk から再開: {done * DT.value_in_unit(unit.picosecond) / 1000:.2f} ns 済み")
    else:
        sim.loadState("equil_state.xml")
        sim.currentStep = 0
        done = 0
    total = int(ns * 1e6 * unit.femtoseconds / DT)
    if total - done > 0:
        run(sim, "prod", total - done, append=args.resume and done > 0)
    state = sim.context.getState(getPositions=True, enforcePeriodicBox=True)
    app.PDBFile.writeFile(top, state.getPositions(), open("prod_final.pdb", "w"))

print(f"完了: {(time.time() - t0) / 60:.1f} 分")
