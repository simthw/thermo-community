# Analysis Module

[![Julia](https://img.shields.io/badge/Julia-1.11+-9558B2?logo=julia)](https://julialang.org/)

Analysis 是一个专门用于 ANMBR（厌氧膜生物反应器）微生物生态系统组装模拟的 Julia 模块。该模块基于 Assembly 项目构建，通过内嵌的 `AssemblyCook` 子模块封装了完整的微生物群落动力学模型，专注于研究厌氧消化系统中的微生物组装过程、热力学约束和代谢途径选择。

本文档整合了模型介绍、复现步骤与当前审计状态。所有命令均从工作区根目录 `D:\demo\um3097_test` 的 PowerShell 运行。

---

## 目录

1. [项目概述](#1-项目概述)
2. [项目结构](#2-项目结构)
3. [核心逻辑](#3-核心逻辑)
4. [当前数据集与状态](#4-当前数据集与状态)
5. [复现与执行](#5-复现与执行)
6. [派生分析与出图](#6-派生分析与出图)
7. [从 candidate 到 canonical 的门槛](#7-从-candidate-到-canonical-的门槛)
8. [已知限制](#8-已知限制)
9. [API 使用](#9-api-使用)

---

## 1. 项目概述

Analysis 模块的核心目标是模拟和分析微生物群落的组装动力学，特别关注：

1. **热力学约束**：通过热力学抑制因子 θ 实现反应的可逆性控制
2. **蛋白组分配**：动态调节代谢酶、核糖体和管家蛋白的比例
3. **能量耦合**：ATP 产量与生长率的耦合机制
4. **种间相互作用**：竞争、互养、促进和污染等相互作用类型
5. **稳健性分析**：不同参数配置下的群落稳定性

### 核心特性

- **内嵌 Assembly 代码**：通过 `AssemblyCook` 子模块封装 Assembly 核心功能，无需外部依赖
- **ANMBR 专用配置**：针对厌氧膜生物反应器优化的参数集（60 菌株、15 代谢物、22 反应）
- **自动热力学计算**：从代谢物标准生成自由能自动计算反应 ΔG0
- **批处理模拟**：支持多组独立模拟的自动化运行和结果归档（schema-3 provenance 元数据）
- **完整分析工具**：热力学分析、相互作用量化、稳健性测试、符号 Jacobian 计算

### 研究设计

**Paired-seed counterfactual**：每个 seed 生成一个固定的组内功能实现（反应分配 + 动力学参数 + 蛋白组分配），同一 seed 在 5 档 SO4 浓度（10/30/100/300/1000 mg/L）下重复模拟，SO4 是唯一受控变化的实验变量；不同 seed 是不同组内功能实现。SO4 的作用通过两条耦合通道实现：(a) 电子流的热力学可行性；(b) 共享底物的生长竞争；哪条占主导取决于具体组装出的菌株组合。

---

## 2. 项目结构

```
Analysis/
├── README.md                 # 本文档
├── src/
│   ├── Analysis.jl           # 主模块：配置管理和批处理流程
│   ├── constants.jl          # 物理常数（Rgas、ΔGATP、NA）及辅助函数（Keq、Q）
│   ├── simulation_registry.jl # run-set 注册表：数据集解析与校验
│   ├── assembly_cook/        # 内嵌 Assembly 核心代码
│   │   ├── AssemblyCook.jl   # 子模块入口，导出所有核心函数
│   │   ├── fullParameters.jl # 数据结构定义（Reaction、MicrobeP、FullParameters）
│   │   ├── make_full.jl      # 参数生成函数（initialise_anmbr 等）
│   │   ├── make_full_robust.jl    # 稳健性测试参数生成器
│   │   ├── simulate_full.jl  # 动力学模拟核心（ODE 求解）
│   │   ├── thermodynamics.jl # 热力学计算函数
│   │   ├── interactions.jl   # 种间相互作用分析
│   │   ├── robustness.jl     # 稳健性测试
│   │   └── analytic.jl       # 解析分析工具（SymPy 符号 Jacobian）
│   ├── data/
│   │   ├── reactions.dat           # 反应数据库（22 个厌氧消化反应）
│   │   └── metabolite_energies.csv # 代谢物标准生成自由能（ΔGf0）
│   ├── figures/              # 图脚本（见 §6）
│   └── scripts/
│       ├── run_anmbr.jl            # 串行批处理脚本
│       ├── run_anmbr_parallel.jl   # 并行批处理脚本
│       ├── run_ablation_theta.jl   # θ-off 消融运行脚本
│       └── run_ablation_adm1.jl    # deterministic five-guild baseline
├── test/
│   └── runtests.jl           # 单元测试（守恒/RHS 纯度/registry 校验等）
├── simulations/              # 原始模拟输出（含注册表 0_rule/）
│   └── 0_rule/runs.toml      # run-set 注册表（数据集的规范来源）
├── results/                  # 派生 CSV 与图表（见下）
└── .ai/                      # 项目记忆（写作约束与解释边界）
```

### results/ 布局

```
results/
├── <run_set>/            # 本次运行的全部结果（图 + 表格 CSV），按 run-set 名（与 runs.toml 一致）
└── z_archive/            # 历史（旧模型/旧数据）结果，保留但不参与当前叙述
```

### 命名规则

- **原始目录**：`simulations/yyyy-mm-ddTHHMMSSsssZ/`（θ-on）、`simulations/ablation_theta/...`（θ-off）、`simulations/ablation_adm1/...`（baseline）。原始模拟目录一经生成不重命名。
- **Run-set 名称**：`<purpose>_<YYYY_MM_DD>`（如 `full_ensemble_2026_08_01`），在 `runs.toml` 注册。
- **派生输出**：`results/<run_set>/<figure>/...`（与 `runs.toml` 的 `outputs` 路径一致）；消融结果在 `results/<run_set>/ablation_theta/`、`results/<run_set>/ablation_adm1/`；每个 artifact 目录有 `source.toml`（含生成脚本相对路径与 SHA-256）。
- **归档**：不再使用的历史结果移入 `results/z_archive/`，不删除、不参与当前分析。

---

## 3. 核心逻辑

### 3.1 数据模型

#### Reaction（反应）

```julia
struct Reaction
  ID::Int64                    # 反应标识符
  Rct::Vector{Int64}           # 底物代谢物索引
  Prd::Vector{Int64}           # 产物代谢物索引
  stoichR::Vector{Float64}     # 底物化学计量系数
  stoichP::Vector{Float64}     # 产物化学计量系数
  ΔG0::Float64                 # 标准自由能变化 (J/mol)
end
```

#### MicrobeP（微生物参数）

```julia
struct MicrobeP
  MC::Int64                    # 细胞质量（氨基酸数）
  γm::Float64                  # 最大延伸率（氨基酸/分钟）
  ρ::Float64                   # ATP/合成步骤
  Kγ::Float64                  # 能量阈值
  Pb::Float64                  # 结合核糖体比例
  d::Float64                   # 死亡率
  ϕH::Float64                  # 管家蛋白比例
  KΩ::Float64                  # 蛋白组饱和常数
  fd::Float64                  # 蛋白组重排所需倍增次数
  R::Int64                     # 反应数量
  Reacs::Vector{Int64}         # 反应索引列表
  η::Vector{Float64}           # ATP 产量（每摩尔反应）
  kc::Vector{Float64}          # 催化速率常数
  KS::Vector{Float64}          # 半饱和常数
  kr::Vector{Float64}          # 可逆性因子
  n::Vector{Int64}             # 蛋白类型氨基酸数 [r,p,h]
  ϕP::Vector{Float64}          # 代谢蛋白比例（sum ≈ 1.0）
end
```

#### FullParameters（系统参数）

```julia
struct FullParameters
  N::Int64                     # 菌株数量
  M::Int64                     # 代谢物数量
  O::Int64                     # 反应数量
  T::Float64                   # 温度 (K)
  κ::Vector{Float64}           # 外部资源供给率
  δ::Vector{Float64}           # 资源衰减率
  reacs::Vector{Reaction}      # 反应列表
  mics::Vector{MicrobeP}       # 菌株列表
end
```

### 3.2 热力学模型

#### 热力学抑制因子 θ

```julia
θ(S, P, T, η, r) = clamp(Q / Keq, 0.0, 1.0)
```

- `Q = ∏(P^νP) / ∏(S^νS)`：反应商
- `Keq = exp((-ΔG0 - η·ΔGATP) / (R·T))`：平衡常数（含 ATP 修正）
- θ ∈ [0, 1]，θ = 1 表示完全抑制（反应处于/越过平衡），θ = 0 表示无热力学抑制
- 当底物 S ≤ 0 时 θ = 1（无底物，完全抑制）；当产物 P ≤ 0 时 θ = 0（无热力学抑制）
- 缓冲代谢物（`BUFFERED_METABOLITES = Set([5, 8])`，即 H₂O 和 H⁺）浓度恒定

#### 反应速率 q

```julia
q = kc·E·S·(1-θ) / (KS + S·(1+kr·θ))
```

- 热力学抑制通过 `(1-θ)` 项实现；可逆性通过 `kr·θ` 项实现
- `S` 取底物向量中的限制性底物（`minimum(S ./ stoichR)`）
- `thermodynamic_feedback=false` 时 θ ≡ 0，退化为纯 Monod 速率

#### 酶量计算 Eα

```julia
Eα(ϕR, mic, i) = MC·(1-ϕR-ϕH)·ϕP[i] / n[2]
```

#### 生长与蛋白组

```julia
γ(a, mic)   = γm·a / (Kγ + a)          # 延伸率
λ(a, ϕR, mic) = γ·ϕR·Pb / n[1]          # 生长率
ϕR(a, mic)  = (1-ϕH)·a / (KΩ + a)       # 核糖体比例
```

### 3.3 动力学方程

系统状态向量：`x = [populations(N), metabolites(M), energies(N), fractions(N)]`

```
d(pop_i)/dt   = (λ_i - d_i) · pop_i
d(energy_i)/dt = J_i - (MC·ρ + energy_i) · λ_i     # J_i = Σ(η_j · rate_ij)
d(ϕ_i)/dt     = (ϕR - ϕ_i) / τ,  τ = fd / λ
d(conc_j)/dt  = κ_j - δ_j·conc_j + Σ(rate_ij · stoich · pop_i / NA)
```

边界行为：

- 种群 `pop ≤ 1e-10` 或能量 `energy ≤ 0.33` 时该菌株的种群/能量/蛋白组动力学置零
- 代谢物与种群在非正边界只抑制外向负导数，正产生仍可从零池发生（RHS 不修改 solver 状态）
- 缓冲代谢物（H₂O、H⁺）浓度导数恒为 0

### 3.4 参数生成（initialise_anmbr）

1. **解析反应数据库**：从 `reactions.dat` 读取 22 个厌氧消化反应，从 `metabolite_energies.csv` 的 ΔGf0 自动计算 ΔG0
2. **菌株类型分配**：60 株按固定比例分配到 5 个 guild（`STRAIN_TYPE_RATIOS = [0.35, 0.20, 0.25, 0.10, 0.10]`）：

| Guild                         | 反应 ID | 菌株数 | 占比 |
| ----------------------------- | ------- | -----: | ---: |
| acidifier（酸化发酵）         | 1-6     |     21 |  35% |
| hydrogen_producer（产氢）     | 14-20   |     12 |  20% |
| methanogen（产甲烷）          | 21-22   |     15 |  25% |
| sulfate_reducer（硫酸盐还原） | 7-11    |      6 |  10% |
| acetogen（产乙酸）            | 12-13   |      6 |  10% |

3. **随机化**：每株从所属 guild 随机抽 1-5 个反应；`kc/KS/kr` 乘以 `2^Normal(0,1)`；`η` 按热力学可行范围抽样（放能反应 η ∈ [1/3, ηh]；吸能产氢反应默认 η=0，与可靠 H₂ 消费反应形成放能耦合时才给保守正 η）；`ϕP` 随机后归一化
4. **H2 耦合**：产氢反应的静态 η 资格按 producer H2 产量与 consumer H2 化学计量之比缩放（`compute_coupled_dG0`），仅绑定可靠消费反应（rx21 氢营养产甲烷）；条件性消费（rx10 硫酸盐还原）交由运行时 θ 控制

### 3.5 代谢物定义

| 索引 | 化学式               | 名称           |
| ---- | -------------------- | -------------- |
| 1    | C₆H₁₂O₆(aq)          | 葡萄糖         |
| 2    | HCO₃⁻(aq)            | 碳酸氢根       |
| 3    | CH₃COO⁻(aq)          | 乙酸           |
| 4    | CO₂(g)               | 二氧化碳       |
| 5    | H₂O(l)               | 水（缓冲）     |
| 6    | SO₄²⁻(aq)            | 硫酸根         |
| 7    | HS⁻(aq)              | 硫化氢         |
| 8    | H⁺(aq)               | 氢离子（缓冲） |
| 9    | CH₃CH(OH)COO⁻(aq)    | 乳酸           |
| 10   | CH₃CH₂CH₂COO⁻(aq)    | 丁酸           |
| 11   | CH₃CH₂COO⁻(aq)       | 丙酸           |
| 12   | CH₃CH₂CH₂CH₂COO⁻(aq) | 戊酸           |
| 13   | H₂(g)                | 氢气           |
| 14   | CH₄(g)               | 甲烷           |
| 15   | CO₃²⁻(aq)            | 碳酸根         |

### 3.6 运行默认值

| 参数                 |                    默认值 | 说明                        |
| -------------------- | ------------------------: | --------------------------- |
| repeats              | 5（串行）/ 50（并行 CLI） | 独立模拟次数                |
| kc / KS / kr         |  10.0 / 1.375e-3 M / 10.0 | 动力学均值（KS = ¼×5.5e-3） |
| Tmax                 |                   2.0e6 s | 最大积分时间                |
| initial_population   |            1.0e5 cells/株 | 初始生物量                  |
| initial_energy       |                     1.0e5 | 初始 ATP 储备               |
| initial_fraction     |                     0.128 | 初始核糖体比例              |
| glucose_mM           |                       5.0 | 初始葡萄糖（mM）            |
| alkalinity_mgCaCO3_L |                    1200.0 | 碱度                        |
| sulfate_mg_L         |                      30.0 | 硫酸盐（mg/L）              |
| seed                 |                         1 | 基种子（repeats 内递增）    |
| δ（稀释/衰减）       |                6.0e-5 s⁻¹ | 所有代谢物                  |
| T                    |                     310 K | 温度                        |

进料换算：`κ[1]=3.3e-7 M/s`（葡萄糖）、`κ[2]=δ·alkalinity/100000`（HCO₃⁻）、`κ[6]=δ·sulfate/96000`（SO₄²⁻）。初始条件：葡萄糖 `glucose_mM/1000`、HCO₃⁻ `alkalinity/100000`、SO₄²⁻ `sulfate/96000`、CO₂ `1.5e-5 M`、H₂O `1.0`（活度≡1，buffered）、H⁺ `1e-7`、CO₃²⁻ `alkalinity/1000000`。

求解器：`Rodas5(autodiff=AutoFiniteDiff())`，`reltol=1e-6, abstol=1e-8, maxiters=1e7`。

---

## 4. 当前数据集与状态

### 4.1 注册数据集（simulations/0_rule/runs.toml）

| run set                     | 内容                                                   | 状态      |
| --------------------------- | ------------------------------------------------------ | --------- |
| `full_ensemble_2026_08_01`  | θ-on 100-seed × 5-SO4（10/30/100/300/1000 mg/L）       | candidate |
| `theta_ablation_2026_08_03` | θ-off 100-seed × 5-SO4（thermodynamic_feedback=false） | candidate |

数据为 schema-3 输出：config/JLD2 记录 `thermodynamic_feedback`、`model_revision`、`input_sha256`（源码/输入文件 SHA-256）、solver retcode/success、`reached_tmax`，并写 `COMPLETED.toml`。两组均通过 `validate_run_set` + `validate_paired_run_directories`（同 seed 跨条件的 reaction/strain 参数身份一致）。

### 4.2 当前 regime 分布（dissipation_rate.jl 分类）

分类判据（按 SO4 升序的终态 CH4 曲线）：`NO_CH4` = 五档全 < 1e-6；否则 `r = CH4@1000/CH4@10`（下限 1e-12）：`r ≥ 1.5` → SYNTROPHY_ACT，`r ≤ 0.5` → ACET_COMP，否则 NO_RESP。

| Regime        | θ-on（2026-08-01） | θ-off（2026-08-03） |
| ------------- | -----------------: | ------------------: |
| SYNTROPHY_ACT |                 22 |                   7 |
| ACET_COMP     |                 48 |                  19 |
| NO_RESP       |                 19 |                  68 |
| NO_CH4        |                 11 |                   6 |

θ-off 使 SYNTROPHY_ACT 从 22 降至 7（约 -68%），NO_RESP 上升；四类标签在 θ-off 下仍全部存在。θ 改变分类频率，但不使四类结构消失。seed 级归属用 `compare_syntrophy_seeds.jl` 核对。

### 4.3 验证状态

- `Pkg.test()` 通过 **51/51** 断言（模型不变量 28 + 元数据/solver 6 + registry 10 + schema-3 目录校验 7）。
- 两组 schema-3 数据全部通过 `validate_run_set` 与 `validate_paired_run_directories`（含 1000 个参数 JLD2 的 paired-community 身份核对）。
- `dissipation_rate.jl` 对两组数据均成功运行；CSV 存于 `results/<set>/dissipation_rate/`。

---

## 5. 复现与执行

### 5.1 预检

```powershell
Set-Location D:\demo\um3097_test
julia --project=Analysis -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

### 5.2 小规模 pilot

```powershell
julia --project=Analysis -p 2 Analysis/src/scripts/run_anmbr_parallel.jl --repeats 1 --seed 1 --sulfate 10 --tmax 2.0e6
julia --project=Analysis -p 2 Analysis/src/scripts/run_anmbr_parallel.jl --repeats 1 --seed 1 --sulfate 1000 --tmax 2.0e6
```

检查每条命令打印的 `Results directory`，确认 `config.toml` 的 `schema_version = 3`、`model_revision = "2026-07-22"`、存在 `input_sha256`；`summary.csv` 的 `solver_success=true`、`reached_tmax=true`；`COMPLETED.toml` 的 `all_solver_success=true`、`all_reached_tmax=true`。

### 5.3 Theta-on 100-seed × 5-SO4 sweep

相同 `--seed 1 --repeats 100` 保证每个 SO4 条件使用 seeds 1:100，可按 seed 配对。**不要用 `--no_save_trajectories`**（图脚本需要轨迹重建指标）。

```powershell
$so4Levels = 10, 30, 100, 300, 1000
foreach ($so4 in $so4Levels) {
  julia --project=Analysis -p 23 Analysis/src/scripts/run_anmbr_parallel.jl `
    --repeats 100 --seed 1 --sulfate $so4 --tmax 2.0e6
}
```

记录五个 `Results directory`，按 SO4 10/30/100/300/1000 排列。

### 5.4 Theta-off 100-seed × 5-SO4 sweep

```powershell
$env:ABLATION_THETA_SEEDS = "100"
$env:ABLATION_THETA_WORKERS = "23"
$env:ABLATION_THETA_SO4 = "10,30,100,300,1000"
julia --project=Analysis Analysis/src/scripts/run_ablation_theta.jl
Remove-Item Env:ABLATION_THETA_SEEDS, Env:ABLATION_THETA_WORKERS, Env:ABLATION_THETA_SO4
```

该脚本会打印五个有序目录，并将 `thermodynamic_feedback=false` 写入 config/JLD2。

### 5.5 Deterministic five-guild baseline

```powershell
julia --project=Analysis Analysis/src/scripts/run_ablation_adm1.jl
```

它不是完整 ADM1，也不是单因素 diversity ablation（一次改变菌株数、反应覆盖、动力学随机性、eta 和 proteome allocation），作为确定性参考曲线使用。

### 5.6 注册 run sets

将实际路径填入 `Analysis/simulations/0_rule/runs.toml`：

```toml
[sets.full_ensemble_YYYY_MM_DD]
description = "Theta-on 100-seed x 5-SO4 ensemble."
kind = "simulation_sweep"
status = "candidate"
thermodynamic_feedback = true
paired_seeds = true
varying_config_keys = ["sulfate_mg_L"]
dirs = [
  { label = "SO4_10",   so4_mg_L = 10.0,   path = "simulations/<timestamp-10>" },
  # ... 五个档位
]

[sets.full_ensemble_YYYY_MM_DD.outputs]
dissipation_rate_dir = "results/full_ensemble_YYYY_MM_DD/dissipation_rate"
dissipation_csv = "results/full_ensemble_YYYY_MM_DD/dissipation_rate/dissipation_rate.csv"
rf_attribution_dir = "results/figures/full_ensemble_YYYY_MM_DD/rf_attribution"
```

Theta-off set 增加 `control_set = "full_ensemble_YYYY_MM_DD"` 与 `thermodynamic_feedback = false`。新 set 先用 `candidate`，不要直接标 `canonical`。

### 5.7 Registry 校验

```powershell
julia --project=Analysis -e 'using Analysis; include("Analysis/src/simulation_registry.jl"); name="<set>"; validate_run_set(name); validate_paired_run_directories(run_set_paths(name))'
```

> **必须带 `using Analysis`**，否则 JLD2 走类型重建路径可能产生"community differs"误报。校验失败时不要生成图，也不要把 set 改为 `canonical`。

---

## 6. 派生分析与出图

所有图脚本**无默认 run set**：必须显式 `--set <name>`（或手动传目录），未指定时报错并列出注册表中可用集。`--list-sets` 可列出全部注册集。

### 6.1 核心指标

```powershell
julia --project=Analysis Analysis/src/figures/dissipation_rate.jl --set full_ensemble_2026_08_01
julia --project=Analysis Analysis/src/figures/dissipation_rate.jl --set theta_ablation_2026_08_03
```

### 6.2 θ-on/off 配对对比

```powershell
$env:ANALYSIS_SYNTH_ON_SET  = "full_ensemble_2026_08_01"
$env:ANALYSIS_SYNTH_OFF_SET = "theta_ablation_2026_08_03"
julia --project=Analysis Analysis/src/figures/compare_syntrophy_seeds.jl
Remove-Item Env:ANALYSIS_SYNTH_ON_SET, Env:ANALYSIS_SYNTH_OFF_SET
```

### 6.3 RF 特征归因

RF 归因分两条独立链路（共享 `rf_common.jl` 的 CV 机制，分别用结构/参数特征）：

```powershell
$env:ANALYSIS_FEATURE_RUN_SET = "full_ensemble_2026_08_01"
$env:ANALYSIS_RF_RUN_SET      = "full_ensemble_2026_08_01"

# 链路 1（结构特征）：extract_features.jl → features.csv → rf1_structural.jl
# （diagnose_features.jl / parametric_probe.jl 读 features.csv，紧跟 v1）
julia --project=Analysis Analysis/src/figures/extract_features.jl
julia --project=Analysis Analysis/src/figures/diagnose_features.jl
julia --project=Analysis Analysis/src/figures/parametric_probe.jl
julia --project=Analysis Analysis/src/figures/rf1_structural.jl

# 链路 2（参数特征）：extract_features_v2.jl → features_v2.csv → rf2_parametric.jl
# （v2 内部自带 param probe v2，无需单独跑 parametric_probe）
julia --project=Analysis Analysis/src/figures/extract_features_v2.jl
julia --project=Analysis Analysis/src/figures/rf2_parametric.jl

Remove-Item Env:ANALYSIS_FEATURE_RUN_SET, Env:ANALYSIS_RF_RUN_SET
```

RF 输出是一次分层 5-fold 结果；不要把未超过 fold-wise baseline 解释成不可约高维机制。

### 6.4 单 seed 图与轨迹诊断

```powershell
# 跨 SO4 机制对比（原 Fig1）：传两个 SO4 目录（低 vs 高），同 seed 画 θ/ΔG/通量
julia --project=Analysis Analysis/src/figures/cross_so4_thermodynamics.jl `
  Analysis/simulations/2026-08-01T075236234Z Analysis/simulations/2026-08-01T084218698Z
# 通量结构时序（原 Fig2）：跨条件 × 关键路径通量随时间
julia --project=Analysis Analysis/src/figures/flux_timeseries.jl --set full_ensemble_2026_08_01 --seed 1
# 效率-生长（原 Fig4）
julia --project=Analysis Analysis/src/figures/efficiency_growth.jl --set full_ensemble_2026_08_01 --condition SO4_10
# 轨迹诊断
julia --project=Analysis Analysis/src/figures/trajectory_diagnostics.jl --set full_ensemble_2026_08_01 --condition SO4_10 --seed 1
```


---

## 7. 从 candidate 到 canonical 的门槛

只有同时满足以下条件才把 `status` 改为 `canonical`：

1. Registry semantic validation 通过。
2. 所有 summary 行 `solver_success=true` 且 `reached_tmax=true`。
3. Paired config 只改变声明过的实验变量；`validate_paired_run_directories` 确认同 seed 跨 SO4（以及 theta 对照间）的 reaction/strain 参数身份一致。
4. 明确定义并通过 endpoint 稳态或收敛诊断。
5. 水、H⁺、气体标准态和 operational sigma 的解释边界已冻结。
6. 执行所用源码处于干净、已提交的 git revision，并把 commit 写入 registry notes。

在第 4-5 项解决前，推荐保持 `candidate`，即使程序运行成功。

---

## 8. 已知限制

1. **气相标准态**：CO2/H2/CH4 仍按液相浓度处理，无 fugacity、Henry 定律或气液传质（因 [H2]^4 出现在氢营养反应中，此项影响最大）。
2. **熵产生口径**：σ = Σ|J·ΔG|/T 是项目定义的 operational gross Gibbs-throughput / Clausius 指标，尚未建立含 ATP 储能、进出料和环境边界的完整反应器熵平衡。不要将其等同于完整 reactor entropy production，也不要静默改用 Cook 式 `Fd = -(ΔG + η·ΔGATP)`。
3. **无稳态诊断**：固定积分到 `Tmax` 不等于稳态；summary 无末端导数范数、质量平衡残差或 perturb-and-return 测试。
4. **确定性 baseline 非单因素**：同时改变菌株数、反应覆盖、动力学随机性、eta 和 proteome allocation，不能作为 diversity-only 因果估计。
5. **θ-off 保留四类标签**：数据只支持"θ 改变分类频率"，不支持"θ 对激活必要"或"四类结构依赖 θ"。
6. **RF 不超 baseline**：当前低阶特征与一次分层 5-fold 未检测到稳定预测信号，不能证明机制不可约；仍需 repeated/nested CV。
8. **git commit 未自动记录**：schema-3 只记录源码 SHA-256；正式重跑应在干净已提交 revision 上执行并登记 commit。

---

## 9. API 使用

### 基本用法

```julia
using Analysis

config = Analysis.anmbrRunConfig(repeats = 3, kc = 10.0, KS = 5.5e-3/4, kr = 10.0,
                                 Tmax = 2.0e6, save_trajectories = true)
run_dir = Analysis.run_anmbr_analysis(config)
println("Results stored in: $run_dir")
```

### 自定义参数与模拟

```julia
using Analysis

ps = Analysis.initialise_anmbr(10.0, 5.5e-3/4, 10.0; κ_sulfate=1.875e-8)
N, M = ps.N, ps.M
pop0 = fill(1e5, N)
conc0 = zeros(M)
conc0[1] = 5e-3         # 葡萄糖 5 mM
conc0[5] = 1.0          # H₂O 活度 1（buffered）
energy0 = fill(1e5, N)
phi0 = fill(0.128, N)

C, T, metadata = Analysis.full_simulate(ps, 2.0e6, pop0, conc0, energy0, phi0;
                                        return_metadata=true)
final_pop = C[end, 1:N]
survivors = count(final_pop .> 1e-9)
println("Survivors: $survivors/$(ps.N)")
```

### 热力学分析

```julia
final_conc = C[end, ps.N+1:ps.N+ps.M]
Analysis.analyze_thermodynamic_feasibility(ps, final_conc)
σ_total = Analysis.total_entropy_production(C, ps, size(C,1))
println("Operational σ: $σ_total J/(K·s)")
```

### 相互作用分析

```julia
result = Analysis.quantify_ints(ps, C[end, :])   # 计算昂贵，仅小规模
stats  = Analysis.analyse_ints(result, ps)
```

### 稳健性测试

```julia
result = Analysis.robustness_test(60, 15, 22, 1, 5, 10.0, 5.5e-3/4, 10.0, true, 5e6, 0, 2.0e6)
results = Analysis.batch_robustness_test(60, 15, 22, 1, 5, 10.0, 5.5e-3/4, 10.0, true, 5e6; Tmax=2.0e6)
```
