"""请求/响应数据模型（Pydantic v2）。"""

from __future__ import annotations

import math
import re
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator

from .screen import VALID_BASES, clean_sequence

ID_PATTERN = re.compile(r"^[A-Za-z0-9_.\-]{1,64}$")


class PrimerIn(BaseModel):
    id: str = Field(..., description="引物名称，1-64 个字母/数字/._- 字符")
    sequence: str = Field(..., description="5'->3' 引物序列，仅允许 A/C/G/T")

    @field_validator("id")
    @classmethod
    def _check_id(cls, v: str) -> str:
        v = v.strip()
        if not ID_PATTERN.match(v):
            raise ValueError("引物 id 只能包含字母、数字、点、下划线、连字符（1-64 字符）")
        return v

    @field_validator("sequence")
    @classmethod
    def _check_sequence(cls, v: str) -> str:
        seq = clean_sequence(v)
        if not 4 <= len(seq) <= 120:
            raise ValueError("引物长度需在 4-120 nt 之间")
        bad = set(seq) - VALID_BASES
        if bad:
            raise ValueError(f"引物序列含非法碱基：{''.join(sorted(bad))}（仅允许 A/C/G/T）")
        return seq


class ScreenRequest(BaseModel):
    primers: list[PrimerIn] = Field(..., min_length=1, max_length=200)
    salt_conc: float = Field(0.05, ge=0.001, le=1.5, description="一价盐浓度（M Na+ 当量）")
    temperature: float = Field(25.0, ge=0.0, le=95.0, description="退火温度 / Tm 筛选阈值（°C）")
    min_bp: int = Field(4, ge=3, le=30, description="最小连续互补 bp 数（粗筛门槛）")
    name: str | None = Field(None, max_length=120, description="本次筛选任务的备注名")


class SensitivityRequest(BaseModel):
    """条件敏感性分析请求：盐浓度 × Tm 阈值网格（各 ≤8 项，组合 ≤64 格）。"""

    salt_concs: list[float] = Field(
        ..., min_length=1, max_length=8, description="盐浓度列表（M Na+ 当量），1-8 项"
    )
    tm_thresholds: list[float] = Field(
        ..., min_length=1, max_length=8, description="Tm 阈值列表（°C），1-8 项"
    )

    @field_validator("salt_concs")
    @classmethod
    def _check_salts(cls, v: list[float]) -> list[float]:
        for x in v:
            if not math.isfinite(x) or not 0.001 <= x <= 1.5:
                raise ValueError(f"盐浓度 {x} 非法：需在 0.001-1.5 M 之间")
        if len(set(v)) != len(v):
            raise ValueError("盐浓度列表存在重复值，请去重")
        return v

    @field_validator("tm_thresholds")
    @classmethod
    def _check_tms(cls, v: list[float]) -> list[float]:
        for x in v:
            if not math.isfinite(x) or not 0.0 <= x <= 95.0:
                raise ValueError(f"Tm 阈值 {x} 非法：需在 0-95 °C 之间")
        if len(set(v)) != len(v):
            raise ValueError("Tm 阈值列表存在重复值，请去重")
        return v

    @model_validator(mode="after")
    def _check_cell_count(self) -> "SensitivityRequest":
        cells = len(self.salt_concs) * len(self.tm_thresholds)
        if cells > 64:
            raise ValueError(
                f"条件组合共 {cells} 格，超过上限 64 格，请减少盐浓度或 Tm 阈值个数"
            )
        return self


class JobSummary(BaseModel):
    job_id: str
    name: str | None
    created_at: str
    salt_conc: float
    temperature: float
    min_bp: int
    primer_count: int
    events_total: int
    high_risk_hairpins: int
    high_risk_dimers: int
    discard: list[str]


class ScreenResponse(BaseModel):
    job_id: str
    result: dict[str, Any]


class ErrorResponse(BaseModel):
    detail: str
