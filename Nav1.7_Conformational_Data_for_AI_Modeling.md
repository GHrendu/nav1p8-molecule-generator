# Nav1.7 构象结构数据手册 (AI建模专用)
## 数据来源: Wu et al., Nature Communications 2023 (DOI: 10.1038/s41467-023-38942-3)

---

## 1. 基础结构参数

### 1.1 整体架构
| 参数 | 值 |
|------|-----|
| 复合物 | Human Nav1.7-β1-β2 |
| 分辨率范围 | 2.6 – 3.2 Å |
| α亚基残基数 | ~1500–2200 |
| 跨膜拓扑 | 4 × 6-TM (S1–S6) 伪对称 |
| VSD | S1–S4 (电压感受域) |
| PD | S5–S6 (孔道域) |
| ECL | 多对二硫键稳定的糖基化胞外环 |

### 1.2 结构存档编号
| 状态 | PDB ID | EMDB ID |
|------|--------|---------|
| Apo | 7W9K | — |
| + Bupivacaine (BPV) | 8I5B | EMD-35193 |
| + Lacosamide (LCM) | 8S9B | EMD-40238 |
| + Carbamazepine (CBZ) | 8S9C | EMD-40239 |
| + Vinpocetine (VPC) | 8I5X | EMD-35197 |
| + Hardwickiic acid (HDA) | 8J4F | EMD-35975 |
| + Vixotrigine (VXT) | 8I5Y | EMD-35198 |
| + PF-05089771 | 8I5G | EMD-35194 |

---

## 2. 孔道域(PD)构象变化核心数据

### 2.1 细胞内门控(Intracellular Gate)构象

#### Apo状态门控轮廓 (Oval Contour)
构成残基（按重复域分类）：
- **Repeat I (S6I)**: Leu398
- **Repeat II (S6II)**: Leu960, Phe963
- **Repeat III (S6III)**: Ile1453
- **Repeat IV (S6IV)**: Val1752, Tyr1755

#### 药物结合后收缩态 (BPV/LCM/CBZ)
关键构象变化：
1. **S6IV α→π 螺旋转变** (α-to-π helical transition)
   - 位置：S6IV中段
   - 效应：导致门控显著收缩
2. **S6II/S6III轻微内向移动**
3. **门控残基重排**：
   - Leu398 (S6I) → 向中心移动
   - Leu964 (S6II) → 替代Leu960成为门控接触残基
   - Ile1457 (S6III) → 向中心移动
   - Ile1756 (S6IV) → 替代Val1752成为门控接触残基

#### 孔道半径量化 (HOLE计算)
| 参数 | Apo | BPV/LCM/CBZ结合态 |
|------|-----|-------------------|
| 收缩位点半径 | ~2 Å | **~1 Å** |
| 半径变化 | 基准 | 缩短 >1 Å |
| 门控区域 | 较宽 | 显著收窄 |

---

## 3. 药物结合位点原子级构象数据

### 3.1 Site BIG (Beneath Intracellular Gate)
**通用特征**: 位于S6四螺旋束胞质侧，紧邻门控残基

#### 3.1.1 Bupivacaine (BPV) 结合构象
```
结合模式: 芳香环+哌啶环共平面，丁基尾指向胞质
手性: R型和S型均可结合（丁基尾取向相反）

疏水接触残基:
├── S6II: Leu964, Leu967, Leu968, Phe971
├── S6III: Ile1457
└── S6IV: Ile1756, Leu1760

关键几何:
├── 3,5-二甲基苯基环 → 嵌入疏水核心
├── 哌啶环 → 与苯环近似共面
└── 丁基尾 → 指向胞质出口方向
```

#### 3.1.2 Lacosamide-1 (LCM-1) 结合构象
```
结合模式: 与BPV/CBZ共享腔隙，额外氢键稳定

疏水接触: 同BPV核心区域
氢键网络:
├── 供体: LCM-1 中央酰胺 N-H
├── 受体: Glu406 侧链羧基 (S6I)
└── 键型: N–H···O (常规氢键)
```

#### 3.1.3 Carbamazepine (CBZ) 结合构象
```
结合模式: 三环结构嵌入疏水核心

疏水接触: 同Site BIG核心
氢键网络:
├── 供体: Asn1461 侧链酰胺 N-H (S6III)
├── 受体: CBZ 羰基氧 (C=O)
└── 键型: N–H···O=C
```

### 3.2 Lacosamide-2 (LCM-2) — 中央腔/选择性过滤器下方
**位置**: 中央腔顶部，偏离中心轴，偏向Repeats I & IV
**构象**: Y形分子构象

```
疏水口袋 (苄基环):
├── Repeat IV: Ser1697, Ile1744, Phe1748, Val1751, Val1752
└── 深度: 苯环嵌入疏水口袋

氢键锚点网络:
├── Anchor 1 (双氢键):
│   ├── LCM-2 酰胺 N-H ··· O=C Thr1696 (主链, P1IV)
│   └── LCM-2 酰胺 C=O ··· H-N Thr1696 (主链, P1IV)
│   └── 位置: Repeat IV P1 loop
├── Anchor 2:
│   ├── LCM-2 乙酰氨基 C=O ··· H-N? Lys1406 (DEKA motif, S6III)
│   └── 备注: DEKA选择性过滤器基序组成部分
└── Anchor 3:
    ├── LCM-2 酰胺 N-H ··· Oε Gln360 (S6I)
    └── 稳定残基: Ser390 侧链羟基 (氢键稳定Gln360)
```

### 3.3 III-IV Fenestration (Site F3)
**VPC与HDA共享位点，但结合深度不同**

```
共享疏水接触面:
├── Repeat III:
│   ├── P1III: Thr1404
│   ├── S5III: Trp1332
│   └── S6III: Thr1448, Leu1449, Phe1452
└── Repeat IV:
    ├── S6IV: Ser1697, Ile1744, Phe1748

构象差异:
├── VPC: 分子延伸更深，部分进入中央腔
└── HDA: 结合较浅，与VPC仅部分重叠
```

### 3.4 IV-I Fenestration (Site F4) — Vixotrigine
**构象特征**: 线性分子贯穿窗孔，从中央腔延伸至脂双层

```
分子走向:
├── 酰氨基部分 → 位于SF下方 (中央腔侧)
├── 三芳环骨架 → 穿过窗孔通道
└── 苯基环末端 → 接触脂双层 (胞外侧)

固定残基:
├── Repeat IV (S6IV/P1IV):
│   ├── Phe1692 (P1IV)
│   ├── Thr1695 (P1IV)
│   ├── Thr1696 (P1IV)
│   ├── Val1751 (S6IV)
│   └── Tyr1755 (S6IV) ← 关键锚点
└── Repeat I (S6I):
    ├── Ile386
    ├── Phe387
    └── Phe391

关键氢键 (稳定锚点):
├── 供体: Tyr1755 酚羟基 (S6IV)
├── 受体: VXT 吡咯烷酰胺 C=O 和 N-H
└── 效应: 形成双氢键锚定

状态依赖性结构基础:
├── S6IV α-helix → IV-I窗孔开放 (允许VXT进入)
├── S6IV π-helix → IV-I窗孔闭合 (阻断VXT结合)
└── 结论: VXT结合严格依赖S6IV构象状态
```

### 3.5 VSDIV 细胞外腔 (Site V4EC) — PF-05089771
```
结合位置: VSDIV S1-S4 细胞外半区形成的腔隙
验证结果: 野生型Nav1.7与NavAb嵌合体结构高度一致
功能意义: 亚型选择性抑制剂靶点
```

---

## 4. 功能状态相关的构象开关

### 4.1 S6IV α/π 螺旋转变 (核心构象开关)
| 特征 | α-螺旋态 | π-螺旋态 |
|------|----------|----------|
| 门控状态 | 较开放 | 收缩/关闭 |
| 细胞内门半径 | ~2 Å | ~1 Å |
| IV-I窗孔 | **开放** (VXT可进入) | **闭合** (VXT被排除) |
| Site BIG可及性 | 允许药物进入 | 药物被锁定 |
| 相关药物 | VXT结合态 | BPV/CBZ/LCM结合态 |

### 4.2 窗孔状态开关
| 窗孔 | 开放条件 | 闭合条件 | 可结合药物 |
|------|----------|----------|-----------|
| III-IV (F3) | 激活/失活态 | 静息态 | VPC, HDA |
| IV-I (F4) | S6IV α-helix | S6IV π-helix | VXT |

### 4.3 中央腔体积变化
- 中央腔体积随功能状态波动
- 激活/失活态：腔体较大，可容纳Site C药物
- 不同状态下S6段轴向旋转导致腔内残基重定位

---

## 5. 3D结合位点坐标化命名系统

### 5.1 孔道域(PD)位点层级
```
PD垂直层级 (从胞外到胞内):
├── Site E (Extracellular loops dome)
│   └── ECL形成的穹顶，μ-芋螺毒素KIIIA
├── Site S (Selectivity Filter)
│   └── SF外口，TTX/STX结合
├── Site C (Cavity) — 可细分:
│   ├── UC (Upper Cavity) — LCM-2 (UC14, 偏I/IV)
│   ├── CC (Central Cavity) — 奎尼丁等
│   └── LC (Lower Cavity)
├── Site G (Gate) — 细胞内门控
├── Site I (Inactivation motif) — CBD结合
└── Site BIG (Beneath Intracellular Gate)
    └── BPV, LCM-1, CBZ

PD侧向位点 (Fenestrations):
├── Site F1 (I-II fenestration) — Bulleyaconitine A
├── Site F3 (III-IV fenestration) — VPC, HDA
└── Site F4 (IV-I fenestration) — VXT
```

### 5.2 VSD位点命名
```
Vn[位置][环境]
├── V4EC: VSDIV Extracellular Cavity — PF-05089771
├── V4EM: VSDIV Extracellular Membrane — HWTX-IV, LqhIII
├── V2EP: VSDII Extracellular Pore — Dc1a
└── [预留] VnIC/M/P: 若发现胞内侧配体
```

---

## 6. AI建模关键结构特征摘要

### 6.1 必须精确建模的构象元素
1. **S6IV α→π转变**: 这是多个位点的核心开关，需精确建模中段二面角变化
2. **门控残基侧链旋转**: Leu964, Ile1457, Ile1756的χ1/χ2角在药物结合后显著变化
3. **氢键网络方向性**: 
   - LCM-1与Glu406 (S6I)
   - LCM-2与Thr1696主链 (P1IV)
   - VXT与Tyr1755 (S6IV)
4. **窗孔几何**: III-IV和IV-I窗孔在不同状态下的开口尺寸
5. **中央腔非对称性**: LCM-2偏向Repeats I/IV的偏移结合

### 6.2 多构象系综建议
为准确预测药物结合，建议至少准备以下构象系综：
- **开放态 (α-S6IV)**: 用于窗孔药物对接 (VXT, VPC, HDA)
- **收缩态 (π-S6IV)**: 用于Site BIG药物对接 (BPV, LCM, CBZ)
- **失活态**: 中央腔扩大，Site C药物可及
- **静息态**: 窗孔关闭，排除状态依赖性药物

### 6.3 保守性警示 (AI训练注意)
- S6段在Nav亚型间高度保守 → 结合Site BIG/C/G/F的药物难以获得亚型选择性
- VSD和ECL区域保守性较低 → 选择性药物设计的结构基础
- 状态依赖性药物的评分函数需考虑构象自由能代价

---

## 7. 关键残基速查表 (按重复域)

### Repeat I
| 残基 | 位置 | 功能/相互作用 |
|------|------|--------------|
| Gln360 | S6I | LCM-2氢键受体 |
| Ser390 | S6I | 稳定Gln360 |
| Leu398 | S6I | 门控残基 |
| Ile386 | S6I | VXT疏水接触 |
| Phe387 | S6I | VXT疏水接触 |
| Phe391 | S6I | VXT疏水接触 |

### Repeat II
| 残基 | 位置 | 功能/相互作用 |
|------|------|--------------|
| Leu960 | S6II | Apo门控 |
| Leu964 | S6II | 收缩态门控 + BPV接触 |
| Leu967 | S6II | BPV疏水接触 |
| Leu968 | S6II | BPV疏水接触 |
| Phe971 | S6II | BPV疏水接触 |

### Repeat III
| 残基 | 位置 | 功能/相互作用 |
|------|------|--------------|
| Ile1453 | S6III | Apo门控 |
| Ile1457 | S6III | 收缩态门控 + BPV接触 |
| Thr1404 | P1III | VPC/HDA接触 |
| Trp1332 | S5III | VPC/HDA接触 |
| Thr1448 | S6III | VPC/HDA接触 |
| Leu1449 | S6III | VPC/HDA接触 |
| Phe1452 | S6III | VPC/HDA接触 |
| Lys1406 | S6III | DEKA motif, LCM-2接触 |
| Asn1461 | S6III | CBZ氢键供体 |

### Repeat IV
| 残基 | 位置 | 功能/相互作用 |
|------|------|--------------|
| Val1752 | S6IV | Apo门控 |
| Tyr1755 | S6IV | Apo门控 + VXT氢键锚点 |
| Ile1756 | S6IV | 收缩态门控 + BPV接触 |
| Leu1760 | S6IV | BPV疏水接触 |
| Ser1697 | S6IV/P1IV | LCM-2/VPC/HDA接触 |
| Ile1744 | S6IV | LCM-2/VPC/HDA接触 |
| Phe1748 | S6IV | LCM-2/VPC/HDA接触 |
| Val1751 | S6IV | LCM-2/VXT接触 |
| Thr1696 | P1IV | LCM-2双氢键锚点 |
| Phe1692 | P1IV | VXT接触 |
| Thr1695 | P1IV | VXT接触 |

---

*文档生成时间: 2026-09-07*
*用途: AI辅助药物设计与分子动力学模拟*
