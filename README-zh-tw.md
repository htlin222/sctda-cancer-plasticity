# scTDA 癌症可塑性

> 拓撲資料分析揭示 EGFR 突變肺癌藥物耐受性背後的環狀細胞狀態可塑性

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## 這是什麼？

本專案將持續同調（persistent homology）和 Mapper 演算法應用於經藥物處理的 EGFR 突變非小細胞肺癌（NSCLC）細胞株的單細胞 RNA 定序資料。我們假設藥物耐受性是透過**環狀**（loop-like）的狀態轉換運作——而非線性分支樹狀結構——這些環狀結構可被拓撲資料分析偵測，但標準軌跡推斷方法無法發現。

## 主要發現

🚧 分析進行中——結果將在此更新。

## 快速開始

```bash
# 複製儲存庫
git clone https://github.com/<user>/sctda-cancer-plasticity
cd sctda-cancer-plasticity

# 建立環境
mamba env create -f envs/environment.yml
conda activate sctda

# 下載資料並執行完整流程
make pipeline
```

## 專案結構

```
├── CLAUDE.md                    # Claude Code 指令
├── Makefile                     # 流程指令
├── plan/                        # 研究設計與決策
│   ├── 00-overview.md           # 完整專案計畫
│   ├── 01-todo.md               # 任務清單
│   ├── 02-decisions.md          # 架構決策
│   ├── 03-analysis-protocol.md  # 逐步分析方案
│   └── 04-figure-plan.md        # 論文圖表計畫
├── data/                        # 資料（不納入版控）
├── notebooks/                   # Jupyter 筆記本
├── src/sctda_plasticity/        # Python 套件
│   ├── data.py                  # 載入與前處理
│   ├── topology.py              # PH、Mapper、置換檢定
│   └── visualize.py             # 出版品質圖表
├── scripts/                     # 流程腳本
├── results/                     # 輸出結果（可重新產生）
└── tests/                       # 單元測試
```

## 資料集

| 資料集 | GEO | 說明 | 角色 |
|--------|-----|------|------|
| PC9 + erlotinib | [GSE134839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134839) | 時間序列 Drop-seq | 發現組 |
| PC9 + osimertinib | [GSE150949](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE150949) | 第三代 TKI | 驗證組 |

## 方法摘要

1. **標準 scRNA-seq**：QC → 正規化 → HVG → PCA → UMAP → 分群 → RNA velocity
2. **持續同調**：在 PCA 空間上計算每個處理時間點的 H₀、H₁
3. **統計檢定**：置換檢定（n=500）驗證 H₁ 顯著性
4. **細胞週期消融**：移除 CC 基因 → 重跑 PH → 確認環非細胞週期所致
5. **基因歸因**：Mapper 圖上的基因連結性 → 通路富集分析
6. **跨資料集驗證**：在獨立 TKI 資料集上重複分析

## 引用

論文準備中。

## 授權

MIT © 2026
