# gl_popc

グリチルリチン酸(GL)1分子をPOPC膜の上の水相に置き、拘束なしで膜との相互作用を見るPOC。**まだ実行していない。**

- GL:PubChemから立体化学つきSMILESを取得し、GAFF2 + AM1-BCCで力場パラメータを付与(openmmforcefieldsを使用)
- 電荷:既定は −3(3つのCOOHをすべて脱プロトン化)。`CHARGE_STATE` で変更できる
- 出力 `gl_z.csv`:GL重心のz座標と膜中心の差、GLとリン原子の最短距離

膜部分は `../popc_bilayer` と同じ設定。`addMembrane` に溶質を渡す必要があるので、GL自体が芯の役割を兼ねる。
