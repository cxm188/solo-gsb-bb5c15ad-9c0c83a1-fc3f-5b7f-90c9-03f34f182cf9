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
│   ├── screen.py        # 核心：反对角线扫描、发卡/二聚体、Tm、分级
│   ├── models.py        # Pydantic 请求模型与输入校验
│   ├── db.py            # SQLite 任务历史
│   └── static/          # 原生前端 index.html / styles.css / app.js
├── data/                # 本地运行时的 SQLite 目录（已 gitignore）
├── Dockerfile
├── docker-compose.yml
└── requirements.txt
```
