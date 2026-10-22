# membrane_MD_PoC

OpenMMで脂質二重膜のMDを試すリポジトリ。
最終目標:物質が酵母表面近傍でどのように安定化されているかを見る。

## 構成

計算系ごとに1ディレクトリを作り、スクリプトとREADMEをまとめる。新しいMDを追加するときも同じ形で増やす。

| ディレクトリ | 内容 | 状態 |
|---|---|---|
| [`popc_bilayer/`](popc_bilayer/) | POPC単成分膜(Amber Lipid17 + TIP3P)のMDと物性評価(APL・膜厚・密度プロファイル・S_CD) | 平衡化・本計算・再開・解析の動作確認済み |
| [`gl_popc/`](gl_popc/) | GL 1分子 + POPC膜(GAFF2 + AM1-BCC) | 未実行 |

## 環境

```bash
pip install openmm MDAnalysis matplotlib   # popc_bilayer
conda install -c conda-forge openmm openmmforcefields openff-toolkit rdkit ambertools   # gl_popc
```

スクリプトはすべてCPU計算を既定にしている(GPUは基本的に使わない前提)。

計算の出力(`*.dcd`, `*.pdb`, `log.csv` など)はコミットしない(`.gitignore` 参照)。
