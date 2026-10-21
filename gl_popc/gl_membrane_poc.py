"""
POC: グリチルリチン酸(GL) + POPC二重膜 を OpenMM で回せるか確認する

環境 (conda-forge):
  conda create -n glmd -c conda-forge python=3.11 openmm openmmforcefields \
      openff-toolkit rdkit ambertools requests
構成:
  膜 = Amber Lipid17 (amber14-all.xml に同梱) / 水 = TIP3P / GL = GAFF2 + AM1-BCC
  GL は膜上の水相に置き、拘束なしで膜へ近づくか・吸着するかを見る
出力:
  system.pdb, traj.dcd, gl_z.csv (GL重心z − 膜中心z [nm])
"""
import sys
import requests
import numpy as np
from rdkit import Chem
from openff.toolkit import Molecule
from openmmforcefields.generators import GAFFTemplateGenerator
from openmm import app, unit, LangevinMiddleIntegrator, MonteCarloMembraneBarostat, Platform
from openmm import Vec3

# ---- 設定 ----
CHARGE_STATE = -3        # -3: 3つのCOOHすべて脱プロトン化(pH5〜7想定) / 0: 中性
MEMBRANE_OFFSET = -4.0   # 膜中心z [nm]。GLを原点に置き、膜はその下に配置
SIM_NS = 0.2             # 本計算の長さ [ns] (CPU向けに短縮。GPUなら1〜10に)
TEMP = 303 * unit.kelvin
DT = 4 * unit.femtoseconds  # 水素質量再分配(HMR)前提の4 fs
REPORT_PS = 10

# ---- 1. GLの構造 (PubChemから立体化学つきSMILESを取得) ----
url = ("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
       "glycyrrhizic%20acid/property/IsomericSMILES/TXT")
smiles = requests.get(url, timeout=30).text.strip().splitlines()[0]
rdmol = Chem.MolFromSmiles(smiles)

if CHARGE_STATE == -3:  # カルボン酸 -C(=O)OH -> -C(=O)O-
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    for match in rdmol.GetSubstructMatches(patt):
        o = rdmol.GetAtomWithIdx(match[2])
        o.SetFormalCharge(-1)
        o.SetNoImplicit(True)
        o.SetNumExplicitHs(0)
    Chem.SanitizeMol(rdmol)
print("GL SMILES:", Chem.MolToSmiles(rdmol), "charge:", Chem.GetFormalCharge(rdmol))

mol = Molecule.from_rdkit(Chem.AddHs(rdmol), allow_undefined_stereo=True)
mol.name = "GLY"
mol.generate_conformers(n_conformers=1)
print("AM1-BCC charges (数分かかる)...")
mol.assign_partial_charges("am1bcc")

# ---- 2. 力場 ----
ff = app.ForceField("amber14-all.xml", "amber14/tip3p.xml")
gaff = GAFFTemplateGenerator(molecules=[mol], forcefield="gaff-2.11")
ff.registerTemplateGenerator(gaff.generator)

# ---- 3. 系の構築: GLを原点、膜をその下に ----
top = mol.to_topology().to_openmm()
xyz = mol.conformers[0].m_as("nanometer")
xyz = xyz - xyz.mean(axis=0)
modeller = app.Modeller(top, [Vec3(*p) for p in xyz] * unit.nanometer)
print("膜を構築中...")
modeller.addMembrane(ff, lipidType="POPC",
                     membraneCenterZ=MEMBRANE_OFFSET * unit.nanometer,
                     minimumPadding=1.0 * unit.nanometer,  # CPU向けに膜パッチを小さく
                     positiveIon="K+", negativeIon="Cl-",
                     ionicStrength=0.15 * unit.molar)
app.PDBFile.writeFile(modeller.topology, modeller.positions, open("system.pdb", "w"))
print("原子数:", modeller.topology.getNumAtoms())

# ---- 4. System / Simulation ----
system = ff.createSystem(modeller.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=1.0 * unit.nanometer,
                         constraints=app.HBonds,
                         hydrogenMass=4 * unit.amu)  # HMRで4 fs刻みを可能に
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
sim.minimizeEnergy()
sim.context.setVelocitiesToTemperature(TEMP)

# ---- 5. 解析用のインデックス ----
atoms = list(modeller.topology.atoms())
gl_idx = [a.index for a in atoms if a.residue.chain.index == 0]  # GLは最初のchain
p_idx = [a.index for a in atoms if a.residue.name in ("PC", "POPC") and a.element and a.element.symbol == "P"]
if not p_idx:  # Lipid17の残基名差に備える
    p_idx = [a.index for a in atoms if a.element and a.element.symbol == "P"]

steps_per_report = int(REPORT_PS * unit.picosecond / DT)
n_reports = int(SIM_NS * 1000 / REPORT_PS)
sim.reporters.append(app.DCDReporter("traj.dcd", steps_per_report))
sim.reporters.append(app.StateDataReporter(sys.stdout, steps_per_report * 10,
                     step=True, time=True, potentialEnergy=True, temperature=True,
                     volume=True, speed=True))

with open("gl_z.csv", "w") as f:
    f.write("time_ps,dz_nm,min_dist_P_nm\n")
    for i in range(n_reports):
        sim.step(steps_per_report)
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
        gl_com = pos[gl_idx].mean(axis=0)
        p = pos[p_idx]
        mem_z = p[:, 2].mean()
        dmin = np.min(np.linalg.norm(pos[gl_idx][:, None, :] - p[None, :, :], axis=2))
        f.write(f"{(i + 1) * REPORT_PS},{gl_com[2] - mem_z:.3f},{dmin:.3f}\n")
        f.flush()
print("done")
