# popc_bilayer

POPC二重膜だけをOpenMMで組んで回すPOC。CPUのラップトップで回すことを想定している。

```bash
python popc_bilayer_poc.py            # 既定 0.1 ns
python popc_bilayer_poc.py --ns 1.0 --pad 2.5
```

- 力場:Amber Lipid17(`amber14-all.xml` に同梱)+ TIP3P、150 mM KCl
- 水素質量再分配(HMR)を使い、時間刻みは4 fs
- 温度と圧力:303 K、`MonteCarloMembraneBarostat`(xyは等方的、zは独立)
- `addMembrane` は溶質を必須とするため、水1分子を芯として置いている

## 出力
- `bilayer_init.pdb`, `bilayer_final.pdb`, `traj.dcd`
- `log.csv`:1脂質あたりの面積(APL)とP–P膜厚の時系列

## 妥当性チェック(POPC, 303 K)
- APL ≈ 0.62–0.68 nm²
- P–P膜厚 ≈ 3.7–4.0 nm

## 実行記録
- 2コアのクラウドCPUで確認:系の構築に約20秒。`--pad 2.5` で POPC 126分子・32,191原子。
