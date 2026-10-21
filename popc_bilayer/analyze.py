"""
POPC二重膜の物性解析

  python analyze.py                       # equil.pdb + prod.dcd を解析
  python analyze.py --traj equil.dcd --skip 0.5   # 平衡化の様子を見る
  --skip [ns]: 先頭を捨てる長さ

出力:
  summary.txt        平均値 ± 標準誤差(ブロック平均)
  timeseries.csv/png APL・P–P膜厚の時系列
  density.csv/png    成分別質量密度プロファイル(膜中心基準)
  scd.csv/png        アシル鎖の秩序パラメータ -S_CD(sn-1 パルミトイル / sn-2 オレオイル)
"""
import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--top", default="equil.pdb")
ap.add_argument("--traj", default="prod.dcd")
ap.add_argument("--skip", type=float, default=0.0, help="先頭を捨てる長さ [ns]")
ap.add_argument("--frame-ps", type=float, default=10.0, help="DCDのフレーム間隔 [ps]")
ap.add_argument("--dz", type=float, default=0.1, help="密度プロファイルのビン幅 [nm]")
args = ap.parse_args()

u = mda.Universe(args.top, args.traj)
lip = u.select_atoms("resname POP")
P = lip.select_atoms("name P")
n_lip = P.n_atoms
start = int(args.skip * 1000 / args.frame_ps)
frames = range(start, len(u.trajectory))
if len(frames) < 2:
    raise SystemExit(f"フレーム不足: 全 {len(u.trajectory)} フレーム, skip後 {len(frames)}")
t_ns = np.array(frames) * args.frame_ps / 1000
print(f"POPC {n_lip}分子, 解析フレーム {len(frames)} ({t_ns[0]:.2f}–{t_ns[-1]:.2f} ns)")

# ---- 成分の定義 (CHARMM式の原子名) ----
groups = {
    "water": u.select_atoms("resname HOH"),
    "choline": lip.select_atoms("name N C11 C12 C13 C14 C15 or name H11* H12* H13* H14* H15*"),
    "phosphate": lip.select_atoms("name P O11 O12 O13 O14"),
    "glycerol+ester": lip.select_atoms("name C1 HA HB C2 HS C3 HX HY O21 C21 O22 O31 C31 O32"),
}
used = sum(groups.values(), u.atoms[[]])
groups["acyl chains"] = lip - used
groups["ions"] = u.select_atoms("resname K CL")

# ---- S_CD 用の C–H ペア (先頭フレームで距離から決める) ----
def chain_carbons(prefix):  # sn-2: C22..C218, sn-1: C32..C316
    names = sorted({a.name for a in lip if re.fullmatch(rf"C{prefix}\d+", a.name)},
                   key=lambda s: int(s[2:]))
    return [n for n in names if int(n[2:]) >= 2]

u.trajectory[0]
pairs = {}
for chain, prefix in (("sn-1", "3"), ("sn-2", "2")):
    for cname in chain_carbons(prefix):
        ci, hi = [], []
        for res in lip.residues:
            c = res.atoms.select_atoms(f"name {cname}")
            h = res.atoms.select_atoms("element H or name H*")
            d = np.linalg.norm(h.positions - c.positions[0], axis=1)
            for idx in h.indices[d < 1.2]:
                ci.append(c.indices[0])
                hi.append(idx)
        pairs[(chain, int(cname[2:]))] = (np.array(ci), np.array(hi))

# ---- ループ ----
apl, thick = [], []
edges = np.arange(-40, 40.001, args.dz * 10)  # Å
centers = 0.5 * (edges[1:] + edges[:-1])
dens = {k: np.zeros(len(centers)) for k in groups}
cos2 = {k: [] for k in pairs}

for ts in u.trajectory[start:]:
    Lx, Ly, Lz = ts.dimensions[:3]
    # 膜中心: P原子zの周期境界を考慮した平均
    th = P.positions[:, 2] / Lz * 2 * np.pi
    zc = (np.arctan2(np.sin(th).mean(), np.cos(th).mean()) / (2 * np.pi) * Lz) % Lz

    def zrel(z):
        return (z - zc + Lz / 2) % Lz - Lz / 2

    zp = zrel(P.positions[:, 2])
    apl.append(Lx * Ly / (n_lip / 2) / 100)                      # nm²
    thick.append((zp[zp > 0].mean() - zp[zp <= 0].mean()) / 10)  # nm

    vol = Lx * Ly * (edges[1] - edges[0])                        # Å³
    for k, g in groups.items():
        h, _ = np.histogram(zrel(g.positions[:, 2]), bins=edges, weights=g.masses)
        dens[k] += h / vol * 1.66054                             # amu/Å³ → g/cm³

    pos = u.atoms.positions
    for k, (ci, hi) in pairs.items():
        v = pos[hi] - pos[ci]
        cos2[k].append(np.mean((v[:, 2] / np.linalg.norm(v, axis=1)) ** 2))

nf = len(apl)
for k in dens:
    dens[k] /= nf
apl, thick = np.array(apl), np.array(thick)


def block_se(x, nblock=5):
    if len(x) < 2 * nblock:
        return float("nan")  # フレーム不足
    b = np.array_split(x, nblock)
    m = np.array([bb.mean() for bb in b])
    return m.std(ddof=1) / np.sqrt(nblock)


# ---- 出力 ----
pk = centers[np.argmax(np.where(centers > 0, dens["phosphate"], 0))] - \
     centers[np.argmax(np.where(centers < 0, dens["phosphate"], 0))]
scd = {k: -(1.5 * np.mean(v) - 0.5) for k, v in cos2.items()}

lines = [
    f"frames: {nf}  ({t_ns[0]:.2f}–{t_ns[-1]:.2f} ns), POPC {n_lip}",
    f"APL            = {apl.mean():.4f} ± {block_se(apl):.4f} nm²   (実験目安 POPC 303 K: 0.643)",
    f"P–P膜厚        = {thick.mean():.3f} ± {block_se(thick):.3f} nm",
    f"D_PP(密度ピーク) = {pk / 10:.2f} nm",
    "-S_CD sn-1: " + " ".join(f"C{n}:{s:.3f}" for (c, n), s in scd.items() if c == "sn-1"),
    "-S_CD sn-2: " + " ".join(f"C{n}:{s:.3f}" for (c, n), s in scd.items() if c == "sn-2"),
]
open("summary.txt", "w").write("\n".join(lines) + "\n")
print("\n".join(lines))

np.savetxt("timeseries.csv", np.c_[t_ns, apl, thick], delimiter=",",
           header="time_ns,apl_nm2,thickness_nm", comments="")
np.savetxt("density.csv", np.c_[centers / 10, np.array([dens[k] for k in dens]).T],
           delimiter=",", header="z_nm," + ",".join(dens), comments="")
with open("scd.csv", "w") as f:
    f.write("chain,carbon,minus_scd\n")
    for (c, n), s in scd.items():
        f.write(f"{c},{n},{s:.4f}\n")

fig, ax = plt.subplots(2, 1, figsize=(6, 5), sharex=True)
ax[0].plot(t_ns, apl); ax[0].set_ylabel("APL [nm²]")
ax[0].axhline(0.643, ls="--", c="gray", lw=0.8)
ax[1].plot(t_ns, thick); ax[1].set_ylabel("P–P thickness [nm]"); ax[1].set_xlabel("time [ns]")
fig.tight_layout(); fig.savefig("timeseries.png", dpi=150)

fig, ax = plt.subplots(figsize=(6, 4))
for k, d in dens.items():
    ax.plot(centers / 10, d, label=k)
ax.set_xlabel("z from bilayer center [nm]"); ax.set_ylabel("mass density [g/cm³]")
ax.legend(fontsize=8); fig.tight_layout(); fig.savefig("density.png", dpi=150)

fig, ax = plt.subplots(figsize=(6, 4))
for chain in ("sn-1", "sn-2"):
    xs = [n for (c, n) in scd if c == chain]
    ax.plot(xs, [scd[(chain, n)] for n in xs], "o-", label=chain)
ax.set_xlabel("carbon number"); ax.set_ylabel("-S_CD"); ax.legend()
fig.tight_layout(); fig.savefig("scd.png", dpi=150)
