"""条件敏感性分析：同一批引物 × 「盐浓度 × Tm 阈值」网格，逐格重算风险。

复用 :func:`app.screen.collect_raw_events` 扫描出的事件骨架（与盐/温度无关），
在每个条件格上重算 Tm 与风险等级，并按「事件类型 + 引物 + 配对区间」的对齐键
跨格比较，输出每格的高风险数量、淘汰引物与相对基准条件的风险变化。

单调性说明（预筛事件集的依据）：
  - Tm 关于盐浓度单调不减（盐校正项随 [Na⁺] 递增）；
  - 风险等级关于 Tm 单调不减、关于阈值温度单调不增。
因此「最高盐浓度 × 最低阈值」是最严苛格：该格为 low 的事件在所有格都为 low，
可以安全地排除出对齐矩阵。
"""

from __future__ import annotations

from .screen import (
    DEFAULT_MIN_BP,
    _RISK_ORDER,
    collect_raw_events,
    evaluate_raw_event,
)

# 单轴最多取值数（与 API 校验一致，这里仅作文档化常量）。
MAX_AXIS_VALUES = 8
# 网格最大格数。
MAX_CELLS = MAX_AXIS_VALUES * MAX_AXIS_VALUES
# 进入对齐矩阵的事件上限（按最严苛格的风险与 Tm 排序保留）。
EVENT_CAP = 500

# 进入报告（非背景）的风险等级。
_REPORTED = ("medium", "high", "critical")
# 触发淘汰建议的风险等级。
_HIGH = ("high", "critical")


def _discard_set(event_risks: list[tuple[list[str], str]]) -> set[str]:
    """由 [(涉及引物, 风险)] 聚合出应淘汰的引物集合（worst risk ≥ high）。"""

    worst: dict[str, str] = {}
    for pids, risk in event_risks:
        for pid in pids:
            if _RISK_ORDER[risk] > _RISK_ORDER[worst.get(pid, "low")]:
                worst[pid] = risk
    return {pid for pid, r in worst.items() if r in _HIGH}


def _aggregate(kept, risks):
    """统计一组风险值（与 kept 事件同序）下的格子指标。"""

    high_hp = high_dm = critical = 0
    for ev, risk in zip(kept, risks):
        if risk in _HIGH:
            if ev.type == "hairpin":
                high_hp += 1
            else:
                high_dm += 1
        if risk == "critical":
            critical += 1
    discard = _discard_set([(ev.primer_ids, r) for ev, r in zip(kept, risks)])
    return {
        "high_risk_hairpins": high_hp,
        "high_risk_dimers": high_dm,
        "high_risk_total": high_hp + high_dm,
        "critical_events": critical,
        "discard": discard,
    }


def _delta(cell, base, cell_risks, base_risks) -> dict:
    """相对基准条件的风险变化。"""

    new = gone = upgraded = downgraded = 0
    for ev_type_risk, base_risk in zip(cell_risks, base_risks):
        cell_reported = ev_type_risk in _REPORTED
        base_reported = base_risk in _REPORTED
        if cell_reported and not base_reported:
            new += 1
        elif base_reported and not cell_reported:
            gone += 1
        elif cell_reported and base_reported:
            if _RISK_ORDER[ev_type_risk] > _RISK_ORDER[base_risk]:
                upgraded += 1
            elif _RISK_ORDER[ev_type_risk] < _RISK_ORDER[base_risk]:
                downgraded += 1
    return {
        "high_risk_total": cell["high_risk_total"] - base["high_risk_total"],
        "discard_added": sorted(cell["discard"] - base["discard"]),
        "discard_removed": sorted(base["discard"] - cell["discard"]),
        "events_new": new,
        "events_gone": gone,
        "events_upgraded": upgraded,
        "events_downgraded": downgraded,
    }


def run_sensitivity(
    primers: list[tuple[str, str]],
    *,
    salt_concs: list[float],
    tm_thresholds: list[float],
    min_bp: int = DEFAULT_MIN_BP,
    base_salt: float,
    base_temp: float,
) -> dict:
    """执行条件敏感性分析。

    参数
    ----
    primers: ``[(id, 5'->3' 序列), ...]``（与 /api/screen 相同，来自历史任务）。
    salt_concs / tm_thresholds: 网格两轴（各 1-8 项，组合 ≤ 64，由 API 层校验）。
    min_bp: 沿用原任务的最小互补 bp 数。
    base_salt / base_temp: 基准条件（原任务的盐浓度与 Tm 阈值），
        每格的风险变化均相对该基准计算。
    """

    salts = sorted(salt_concs)
    temps = sorted(tm_thresholds)
    cells = [(s, t) for s in salts for t in temps]

    # 1) 扫描事件骨架，并用「最严苛格」（最高盐 × 最低阈值，含基准条件）预筛：
    #    在该格仍为 low 的事件由单调性可知在所有格均为 low，不进入对齐矩阵。
    raws = collect_raw_events(primers, min_bp=min_bp)
    harshest_salt = max(list(salts) + [base_salt])
    lowest_temp = min(list(temps) + [base_temp])
    scored = []
    for ev in raws:
        tm, risk = evaluate_raw_event(
            ev, salt_conc=harshest_salt, temperature=lowest_temp
        )
        if risk != "low":
            scored.append((_RISK_ORDER[risk], tm, ev.key, ev))
    scored.sort(key=lambda x: (-x[0], -x[1], x[2]))
    events_capped = max(0, len(scored) - EVENT_CAP)
    kept = [ev for _, _, _, ev in scored[:EVENT_CAP]]

    # 2) 逐格重算每个事件的风险；基准条件同样重算（保证同一套事件集与算法）。
    cell_risks: list[list[str]] = []  # cell_risks[cell_idx][event_idx]
    for salt, temp in cells:
        cell_risks.append(
            [
                evaluate_raw_event(ev, salt_conc=salt, temperature=temp)[1]
                for ev in kept
            ]
        )
    base_risks = [
        evaluate_raw_event(ev, salt_conc=base_salt, temperature=base_temp)[1]
        for ev in kept
    ]
    base_stats = _aggregate(kept, base_risks)

    # 3) 每个事件的「首次达到 critical」条件：在所有为 critical 的格子中取
    #    最温和的一格——阈值温度最高者优先，并列时盐浓度最低者（临界条件）。
    first_critical: list[int | None] = []
    for ev_idx in range(len(kept)):
        crit_cells = [
            i for i, risks in enumerate(cell_risks) if risks[ev_idx] == "critical"
        ]
        if crit_cells:
            first_critical.append(
                min(crit_cells, key=lambda i: (-cells[i][1], cells[i][0]))
            )
        else:
            first_critical.append(None)

    first_crit_count = [0] * len(cells)
    for idx in first_critical:
        if idx is not None:
            first_crit_count[idx] += 1

    # 4) 汇总每格指标与相对基准的变化。
    cell_out = []
    for i, (salt, temp) in enumerate(cells):
        stats = _aggregate(kept, cell_risks[i])
        cell_out.append(
            {
                "salt_conc": salt,
                "tm_threshold": temp,
                "is_base": abs(salt - base_salt) < 1e-9
                and abs(temp - base_temp) < 1e-9,
                "high_risk_hairpins": stats["high_risk_hairpins"],
                "high_risk_dimers": stats["high_risk_dimers"],
                "high_risk_total": stats["high_risk_total"],
                "critical_events": stats["critical_events"],
                "discard": sorted(stats["discard"]),
                "first_critical_events": first_crit_count[i],
                "delta": _delta(stats, base_stats, cell_risks[i], base_risks),
            }
        )

    events_out = [
        {
            "key": ev.key,
            "type": ev.type,
            "primers": ev.primer_ids,
            "locus": ev.locus,
            "bp": ev.bp,
            "matches": ev.matches,
            "tm_base": round(ev.tm_base, 2),
            "base_risk": base_risks[i],
            "cell_risks": [risks[i] for risks in cell_risks],
            "first_critical_cell": first_critical[i],
        }
        for i, ev in enumerate(kept)
    ]

    return {
        "base": {
            "salt_conc": base_salt,
            "temperature": base_temp,
            "min_bp": min_bp,
        },
        "salt_concs": salts,
        "tm_thresholds": temps,
        "primer_count": len(primers),
        "events_considered": len(kept),
        "events_capped": events_capped,
        "cells": cell_out,
        "events": events_out,
    }
