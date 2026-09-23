"""请求/响应数据模型（Pydantic v2）。"""

from __future__ import annotations

import re
from typing import Any, Literal

from pydantic import BaseModel, Field, field_validator

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
