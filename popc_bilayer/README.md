# popc_bilayer

POPC二重膜の全原子MDと物性評価。4コア程度のCPUラップトップで回すことを想定している。

## 手順

```bash
pip install openmm MDAnalysis matplotlib

python run_md.py equil --ns 2        # 構築 → 最小化 → NPT平衡化
python run_md.py prod  --ns 8        # 本計算 (10 ps毎にDCD出力)
python run_md.py prod  --ns 8 --resume   # 中断したら prod.chk から再開
python analyze.py --skip 1           # 先頭 1 ns を捨てて解析
```

- 系はPOPC 126分子・約3.2万原子。`addMembrane` の最小パッチがこの大きさなので、`--pad` で小さくはできない(大きくはできる)
- `--threads N` でCPUスレッド数を指定できる
- 再開時、最後のチェックポイント以降に書かれたDCDフレームは重複する。チェックポイントは100 ps毎に保存する

## 設定
- 力場:Amber Lipid17(`amber14-all.xml` に同梱)+ TIP3P、150 mM KCl
- 303 K、1 bar、`MonteCarloMembraneBarostat`(xyは等方的、zは独立)
- 水素質量再分配(HMR)を使い、時間刻みは4 fs。PMEのカットオフは1.0 nm
- `addMembrane` は溶質を必須とするため、水1分子を芯として原点に置いている(初期状態では膜中心に水が1分子だけ入る)

## 解析出力(`analyze.py`)
| ファイル | 内容 |
|---|---|
| `summary.txt` | APL・P–P膜厚の平均とブロック平均による標準誤差、D_PP、−S_CD |
| `timeseries.csv/png` | APLと膜厚の時系列(平衡に達したかの判断に使う) |
| `density.csv/png` | 成分別の質量密度プロファイル(水、コリン、リン酸、グリセロールとエステル、アシル鎖、イオン) |
| `scd.csv/png` | アシル鎖の秩序パラメータ −S_CD(sn-1はパルミトイル、sn-2はオレオイル) |

## 妥当性の目安(POPC, 303 K)
- APL:実験値は約0.64 nm²(Kučerka et al. 2011)。Lipid17はこれに近い値を再現するとされる
- P–P膜厚:約3.7–4.0 nm
- −S_CD:sn-1のC3–C8付近に約0.2のプラトーがあり、鎖の末端に向かって減少する。sn-2はC9–C10(二重結合)でくぼむ

## 実行記録
| 環境 | 速度 | メモ |
|---|---|---|
| クラウド 2コア CPU | 2.76 ns/day | 最小化に約13分。4コアのラップトップなら4〜5 ns/日程度の見込み → 平衡化2 ns + 本計算8 ns で約2〜3日 |

動作テスト(平衡化 40 ps → 本計算 40 ps → 再開して 50 ps):APL = 0.663 nm²、P–P膜厚 = 3.66 nm、D_PP = 3.70 nm、sn-1の−S_CDのプラトー ≈ 0.20、sn-2はC10(二重結合)でくぼむ。ごく短い計算だが、値はいずれも妥当な範囲に入っている。
