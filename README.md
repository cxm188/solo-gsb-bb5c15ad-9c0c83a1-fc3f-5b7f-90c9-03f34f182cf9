# PCR 引物互补性筛选台

在网页上粘贴引物序列 JSON、设置盐浓度与退火温度，后端计算反向互补的局部
配对，分别检测 **发卡结构（Hairpin）** 与 **引物二聚体（含自二聚体）**，
按简化 Tm 阈值筛出高风险位点，返回 **可定位的配对区间**（1-based，5′→3′）
与 **建议淘汰的引物/组合**。任务历史持久化到 SQLite，结果可一键导出 JSON。

技术栈：**FastAPI + SQLite + 原生 HTML/CSS/JS（零构建）**，Docker / Compose
一键启动。

---

## 一、快速启动

### 方式 A：Docker Compose（推荐）

```bash
docker compose up -d --build
```

启动后打开：

```
http://localhost:8000
```

停止：

```bash
docker compose down          # 停止并移除容器（数据保留在 volume）
docker compose down -v       # 同时删除任务历史数据
```

### 方式 B：直接用 Docker

```bash
docker build -t pcr-primer-screen:latest .
docker run -d --name pcr-primer-screen \
  -p 8000:8000 \
  -v pcr-primer-screen-data:/data \
  pcr-primer-screen:latest
```

### 方式 C：本地 Python 运行（开发用）

需要 Python 3.11+：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 端口与配置

| 项目 | 默认值 | 说明 |
| --- | --- | --- |
| 服务端口 | `8000` | Compose 中通过 `ports: "8000:8000"` 映射，改宿主机端口写成 `"18000:8000"` |
| `PCR_DB_PATH` | `/data/jobs.db`（容器内） | SQLite 数据库文件路径，挂载到 named volume `pcr-primer-screen-data` 持久化 |
| 盐浓度 `[Na⁺]` | `0.05` M | 页面/API 可调，范围 0.001–1.5 M |
| 退火温度 / Tm 阈值 | `25.0` °C | 页面/API 可调，范围 0–95 °C |
| 最小互补 bp 数 | `4` | 粗筛门槛，范围 3–30 |

本地直接运行时可用环境变量覆盖数据库位置：

```bash
PCR_DB_PATH=./data/jobs.db uvicorn app.main:app --reload
```

健康检查：`GET http://localhost:8000/healthz` → `{"status":"ok"}`

---

## 二、使用方法

1. 打开页面，点击 **「载入示例」**（或自行粘贴引物 JSON）；
2. 设置 **盐浓度** 与 **退火温度（Tm 阈值）**；
3. 点击 **「开始筛选」**；
4. 查看概览卡片、建议淘汰清单、发卡/二聚体位点（含坐标与三行配对展示）；
5. 点击 **「导出 JSON」** 下载完整结构化结果；历史任务可在页面底部重新载入或删除。
6. （可选）在结果页底部的 **「条件敏感性分析」** 中输入盐浓度列表与 Tm 阈值
   列表（各 ≤8 项、组合 ≤64 格），点击「运行敏感性分析」：后端复用当前任务的
   引物逐格重算风险，页面以 **热力表** 展示每格高风险数量，★ 标出有事件
   **首次达到 critical** 的条件格，蓝框为基准条件（原任务参数）；点击任意格子
   可查看该格的淘汰引物与相对基准的风险变化，「下载 JSON」可导出完整网格结果。

### 输入 JSON 格式

数组形式（推荐）：

```json
[
  {"id": "Primer_F", "sequence": "ATGCGTACGTTAGCGGATCCGGATCCGGAT"},
  {"id": "Primer_R", "sequence": "TACGATCGATTACGGGATCCGGATCCGGAT"}
]
```

或扁平对象形式（页面会自动归一化）：

```json
{ "Primer_F": "ATGCGTACGTTAGCGGATCCGGATCCGGAT",
  "Primer_R": "TACGATCGATTACGGGATCCGGATCCGGAT" }
```

也可以直接提交完整请求体（含参数与备注，见 `/api/example` 的返回结构）：

```json
{
  "name": "IL-6 qPCR 初筛",
  "salt_conc": 0.05,
  "temperature": 25.0,
  "min_bp": 4,
  "primers": [ {"id": "F", "sequence": "ACGTACGTACGTACGT"} ]
}
```

序列只允许 A/C/G/T（大小写、空白均可，会自动规范化），长度 4–120 nt，
引物 id 需唯一（字母/数字/`._-`，1–64 字符）。

### 示例数据

点击「载入示例」或调用 `GET /api/example` 可取得内置示例；仓库中也放了一份
可直接提交的 [`examples/sample_request.json`](examples/sample_request.json)：

```bash
curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' \
  -d @examples/sample_request.json | python3 -m json.tool
```

示例包含：

| 引物 | 设计意图 |
| --- | --- |
| `Hairpin_F1` | 14 bp 茎 + `TTT` 环、3′ 臂贴 3′ 端 → critical 发卡 |
| `Primer_F2` / `Primer_R2` | 共享 16 nt 含回文（GGATCC）的 3′ 尾 → critical 的 3′ 锚定异源二聚体 + 自二聚体 |
| `Normal_F3` | 常规干净引物，作为对照 |

---

## 三、结果如何解读

每个位点都给出 **可定位区间**（1-based、闭区间，方向为引物自身 5′→3′）
和三行配对展示（顶行 5′→3′、中间互补标记、底行 3′→5′）。

### 发卡结构（Hairpin）

- `arm_5prime` / `arm_3prime`：形成茎区的两臂坐标；
- `loop_length`：两臂之间未配对的环长度（模型强制 ≥ 3 nt）；
- `at_edge`：茎区是否贴在 `3-prime` / `5-prime` 端，或位于链内部。
  **3′ 端发卡最危险**——聚合酶无法从被锁住的 3′ 端起延。

### 引物二聚体（Dimer）

- `primer_a` / `primer_b`：参与配对的两条引物（相同则为自二聚体）；
- `interval_a` / `interval_b`：两条引物各自的配对区间；
- `anchored_three_prime`：**双方配对区间都终止于各自 3′ 端**，即 3′ 锚定
  二聚体——聚合酶可直接沿自由 3′ 端延伸，是最差的组合；
- 自二聚体天然对称，结果已做去重，不会重复计数。

### 风险分级（简化 Tm 阈值模型）

对每个局部互补区间：

```
Tm_base = 4·(G+C) + 2·(A+T)                       # Wallace 经验式（基准 ~50 mM）
Tm      = Tm_base + 16.6·log10( f([Na⁺]) / f(0.05) )
          f(x) = x / (1 + 0.7x)                   # 相对 50 mM 基准的 SantaLucia 盐校正
```

| 级别 | 判定 |
| --- | --- |
| `critical` | 3′ 锚定（二聚体）或 3′ 端发卡，且 Tm ≥ 阈值 |
| `high` | Tm ≥ 阈值 |
| `medium` | Tm 在阈值以下 5 °C 以内且配对 ≥ 8 bp |
| `low` | 其余随机短互补，不进入报告 |

- 任引物关联到 `high` / `critical` 事件 → 进入 **建议淘汰** 清单并附原因；
- 每类事件最多保留 200 条（按风险、Tm 排序），被截断时会在响应 `capped`
  字段和页面上注明。

> 注：这是教学/初筛用途的 **简化 Tm 模型**，未使用完整最近邻热力学参数，
> 绝对 Tm 会偏保守；请以相对风险分级为主，关键引物请用专业工具复核。

### JSON 导出

页面「导出 JSON」下载的是完整响应（与 `POST /api/screen` 返回一致）：

```jsonc
{
  "job_id": "a02d91240a1b",
  "result": {
    "salt_conc": 0.05, "temperature": 25.0, "min_bp": 4,
    "primer_count": 4,
    "hairpins": [
      { "type": "hairpin", "primer_id": "Hairpin_F1",
        "arm_5prime": [1, 14], "arm_3prime": [18, 31],
        "loop_length": 3, "bp": 14, "gc_pairs": 7,
        "tm_base": 42.0, "tm": 42.0,
        "at_edge": "3-prime", "risk": "critical",
        "recommendation": "...",
        "alignment": { "top_5to3": "...", "pairs": "||||...", "bottom_3to5": "..." } }
    ],
    "dimers": [
      { "type": "dimer", "primer_a": "Primer_F2", "primer_b": "Primer_R2",
        "interval_a": [17, 30], "interval_b": [17, 30],
        "anchored_three_prime": true, "bp": 14, "gc_pairs": 8,
        "tm_base": 44.0, "tm": 44.0, "risk": "critical",
        "recommendation": "...", "alignment": { "..." : "..." } }
    ],
    "discard": ["Hairpin_F1", "Primer_F2", "Primer_R2"],
    "summary": { "...": "每引物聚合统计与淘汰原因" },
    "events_total": 12,
    "capped": { "hairpins": 0, "dimers": 0 }
  }
}
```

---

## 四、HTTP API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/api/screen` | 提交引物与参数，执行筛选并落库 |
| `POST` | `/api/jobs/{job_id}/sensitivity` | 条件敏感性分析：复用该任务的引物，在「盐浓度 × Tm 阈值」网格上逐格重算风险 |
| `GET` | `/api/example` | 获取内置示例请求体 |
| `GET` | `/api/jobs` | 历史任务列表（最近 50 条） |
| `GET` | `/api/jobs/{job_id}` | 某任务的请求 + 完整结果 |
| `DELETE` | `/api/jobs/{job_id}` | 删除历史任务 |
| `GET` | `/healthz` | 健康检查 |
| `GET` | `/docs` | FastAPI 自动生成的 Swagger 文档 |

命令行示例：

```bash
curl -s http://localhost:8000/api/example | curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' -d @- | python3 -m json.tool
```

### 条件敏感性分析

对某个历史任务（`/api/screen` 返回的 `job_id`）执行：

```bash
JOB_ID=$(curl -s -X POST http://localhost:8000/api/screen \
  -H 'Content-Type: application/json' -d @examples/sample_request.json \
  | python3 -c "import sys,json;print(json.load(sys.stdin)['job_id'])")

curl -s -X POST "http://localhost:8000/api/jobs/${JOB_ID}/sensitivity" \
  -H 'Content-Type: application/json' \
  -d '{"salt_concs": [0.03, 0.05, 0.1], "tm_thresholds": [20.0, 25.0, 30.0]}' \
  | python3 -m json.tool
```

请求体：

| 字段 | 约束 |
| --- | --- |
| `salt_concs` | 1–8 项，每项 0.001–1.5 M，不得重复 |
| `tm_thresholds` | 1–8 项，每项 0–95 °C，不得重复 |
| 组合上限 | 两轴乘积 ≤ 64 格 |

要点：

- **复用同一批引物**：配对区间只与序列和 `min_bp` 有关，后端先扫描一次事件
  骨架，再逐格重算 Tm（盐校正）与风险等级（阈值比较），网格计算不重复扫描；
- **事件对齐**：每个事件以「类型 + 引物 + 配对区间」为键（如
  `hairpin|Hairpin_F1|1-14|18-31`），`events[].cell_risks` 与 `cells` 同序，
  可跨格追踪同一事件的风险变化；
- **基准条件**：原任务的 `salt_conc` / `temperature`，每格 `delta` 均相对它计算
  （高风险数变化、新增/解除淘汰、事件升级/降级/新出现/消失）；
- **首次 critical**：`events[].first_critical_cell` 给出该事件首次达到 critical
  的格子索引（取最温和的临界条件：阈值最高者优先，并列时盐浓度最低者），
  每格 `first_critical_events` 为对应计数，前端以 ★ 标注；
- 进入对齐矩阵的事件上限为 500 个（按最严苛格风险排序），超出部分计入
  `events_capped`；在所有格均为 low 的背景事件不进入矩阵。

响应节选：

```jsonc
{
  "job_id": "a02d91240a1b",
  "base": { "salt_conc": 0.05, "temperature": 25.0, "min_bp": 4 },
  "salt_concs": [0.03, 0.05, 0.1],
  "tm_thresholds": [20.0, 25.0, 30.0],
  "events_considered": 33,
  "events_capped": 0,
  "cells": [
    {
      "salt_conc": 0.05, "tm_threshold": 25.0, "is_base": true,
      "high_risk_hairpins": 1, "high_risk_dimers": 8, "high_risk_total": 9,
      "critical_events": 4,
      "discard": ["Hairpin_F1", "Primer_F2", "Primer_R2"],
      "first_critical_events": 0,
      "delta": { "high_risk_total": 0, "discard_added": [], "discard_removed": [],
                 "events_new": 0, "events_gone": 0,
                 "events_upgraded": 0, "events_downgraded": 0 }
    }
  ],
  "events": [
    {
      "key": "dimer|Primer_F2|Primer_F2|17-30|17-30",
      "type": "dimer",
      "primers": ["Primer_F2", "Primer_F2"],
      "locus": { "primer_a": "Primer_F2", "interval_a": [17, 30], "...": "..." },
      "bp": 14, "tm_base": 44.0,
      "base_risk": "critical",
      "cell_risks": ["critical", "critical", "..."],
      "first_critical_cell": 2
    }
  ]
}
```

错误返回（均为 4xx + 明确 `detail`）：

| 情形 | 状态码 |
| --- | --- |
| 任务不存在 | `404` |
| 列表为空 / 超过 8 项 / 数值越界 / 含重复值 / 非数值 | `422` |
| 组合超过 64 格 | `422` |

---

## 五、算法说明（简版）

两条核酸链 Watson-Crick 配对时方向相反。固定链 A 按 5′→3′、链 B 在复合
物中按 3′→5′ 摆放，则 A 上位置 `i` 与 B 上位置 `k` 的互补碱基满足
`i + k = const`，即落在同一 **反对角线** 上：

1. 枚举所有反对角线（允许两链带悬挂端相对滑动）；
2. 沿对角线逐格判断互补，切出所有极大连续互补段（无 gap，结果可解释）；
3. 发卡检测：链与自身配对，丢弃两臂重叠或环长 < 3 nt 的平凡自配对；
4. 二聚体检测：自身（自二聚体，不要求环长）与两两异源配对，对称结果去重；
5. 对每个互补段计算 Wallace + 盐校正 Tm，按上表分级并聚合成淘汰建议。

纯标准库实现，无 numpy/生物信息学重依赖；30 条约 20–28 nt 的引物全量
（含所有两两组合）筛选约 0.4 秒。

## 项目结构

```
.
├── app/
│   ├── main.py          # FastAPI 路由 + 静态托管
│   ├── screen.py        # 核心：反对角线扫描、发卡/二聚体、Tm、分级、事件骨架
│   ├── sensitivity.py   # 条件敏感性分析：网格重算、事件对齐、风险变化
│   ├── models.py        # Pydantic 请求模型与输入校验
│   ├── db.py            # SQLite 任务历史
│   └── static/          # 原生前端 index.html / styles.css / app.js
├── data/                # 本地运行时的 SQLite 目录（已 gitignore）
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```
