/* PCR 引物互补性筛选台 —— 原生 JS，无依赖、无构建。 */
"use strict";

const $ = (sel) => document.querySelector(sel);

const els = {
  jobName: $("#jobName"),
  salt: $("#salt"),
  temp: $("#temp"),
  minbp: $("#minbp"),
  json: $("#jsonInput"),
  msg: $("#inputMsg"),
  run: $("#btnRun"),
  example: $("#btnLoadExample"),
  format: $("#btnFormat"),
  exportBtn: $("#btnExport"),
  refreshJobs: $("#btnRefreshJobs"),
  empty: $("#emptyResult"),
  body: $("#resultBody"),
  overview: $("#overview"),
  discardBox: $("#discardBox"),
  hairpins: $("#hairpinList"),
  dimers: $("#dimerList"),
  jobList: $("#jobList"),
};

let lastResponse = null; // { job_id, result }

/* ---------------- 工具 ---------------- */

function setMsg(text, kind = "") {
  els.msg.textContent = text;
  els.msg.className = "msg" + (kind ? " " + kind : "");
}

const RISK_CN = { critical: "极高", high: "高", medium: "中", low: "低" };

function badge(text, cls) {
  return `<span class="badge ${cls}">${text}</span>`;
}

function esc(s) {
  return String(s).replace(/[&<>"]/g, (c) =>
    ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
}

/* 接受数组形式或扁平 id->seq 形式的 JSON。 */
function normalizePrimers(parsed) {
  if (Array.isArray(parsed)) {
    return parsed.map((p) => {
      if (!p || typeof p.id !== "string" || typeof p.sequence !== "string") {
        throw new Error("数组中每项需为 {\"id\": ..., \"sequence\": ...}");
      }
      return { id: p.id, sequence: p.sequence };
    });
  }
  if (parsed && typeof parsed === "object") {
    return Object.entries(parsed).map(([id, sequence]) => {
      if (typeof sequence !== "string") throw new Error(`引物 ${id} 的序列必须是字符串`);
      return { id, sequence };
    });
  }
  throw new Error("JSON 需为引物数组或 id->sequence 对象");
}

/* ---------------- 渲染 ---------------- */

function alignmentBlock(al) {
  return `<div class="align">` +
    `<span class="dir">5'→3'</span> ${esc(al.top_5to3)}\n` +
    `             ${esc(al.pairs)}\n` +
    `<span class="dir">3'←5'</span> ${esc(al.bottom_3to5)}` +
    `</div>`;
}

function renderHairpin(e) {
  const edgeTag =
    e.at_edge === "3-prime" ? badge("3′ 端锚定", "tag")
    : e.at_edge === "5-prime" ? badge("5′ 端", "tag")
    : badge("链内部", "tag");
  return `<div class="event ${e.risk}">
    <div class="event-title">
      ${badge(RISK_CN[e.risk] + "风险", e.risk)}
      <span class="who">${esc(e.primer_id)}</span>
      ${edgeTag}
    </div>
    <div class="event-meta">
      <span>5′ 臂 <b>${e.arm_5prime[0]}–${e.arm_5prime[1]}</b></span>
      <span>3′ 臂 <b>${e.arm_3prime[0]}–${e.arm_3prime[1]}</b></span>
      <span>环长 <b>${e.loop_length}</b> nt</span>
      <span>茎区 <b>${e.bp}</b> bp</span>
      <span>GC <b>${e.gc_pairs}</b></span>
      <span>Tm <b>${e.tm}°C</b> <span style="opacity:.7">(无盐 ${e.tm_base}°C)</span></span>
    </div>
    ${alignmentBlock(e.alignment)}
    <div class="advice">💡 ${esc(e.recommendation)}</div>
  </div>`;
}

function renderDimer(e) {
  const aTag = e.a_is_three_prime ? badge("A 侧 3′", "tag") : "";
  const bTag = e.b_is_three_prime ? badge("B 侧 3′", "tag") : "";
  const selfTag = e.same_primer ? badge("自二聚体", "tag") : badge("异源二聚体", "tag");
  const anchorTag = e.anchored_three_prime ? badge("3′ 锚定 · 可延伸", "critical") : "";
  return `<div class="event ${e.risk}">
    <div class="event-title">
      ${badge(RISK_CN[e.risk] + "风险", e.risk)}
      <span class="who">${esc(e.primer_a)}</span>
      <span style="color:var(--muted)">×</span>
      <span class="who">${esc(e.primer_b)}</span>
      ${selfTag}${anchorTag}${aTag}${bTag}
    </div>
    <div class="event-meta">
      <span>A 区间 <b>${e.interval_a[0]}–${e.interval_a[1]}</b></span>
      <span>B 区间 <b>${e.interval_b[0]}–${e.interval_b[1]}</b></span>
      <span>配对 <b>${e.bp}</b> bp</span>
      <span>GC <b>${e.gc_pairs}</b></span>
      <span>Tm <b>${e.tm}°C</b> <span style="opacity:.7">(无盐 ${e.tm_base}°C)</span></span>
    </div>
    ${alignmentBlock(e.alignment)}
    <div class="advice">💡 ${esc(e.recommendation)}</div>
  </div>`;
}

function renderResult(data) {
  lastResponse = data;
  const r = data.result;
  els.empty.hidden = true;
  els.body.hidden = false;
  els.exportBtn.disabled = false;

  const cards = [
    [r.primer_count, "引物数量"],
    [r.summary.hairpin_events, "发卡位点"],
    [r.summary.dimer_events, "二聚体位点"],
    [r.summary.high_risk_hairpins, "高风险发卡"],
    [r.summary.high_risk_dimers, "高风险二聚体"],
    [r.discard.length, "建议淘汰"],
  ];
  els.overview.innerHTML = cards
    .map(([num, lbl]) => `<div class="card"><div class="num">${num}</div><div class="lbl">${lbl}</div></div>`)
    .join("");

  if (r.discard.length) {
    const perPrimer = r.summary.primers.filter((p) => p.discard);
    els.discardBox.className = "discard-box";
    els.discardBox.innerHTML = `
      <h4>🚫 建议淘汰的引物（${r.discard.length}）</h4>
      <div class="ids">${r.discard.map(esc).join("，")}</div>
      <ul>${perPrimer
        .map((p) => `<li><b>${esc(p.primer_id)}</b>：${p.reasons.map(esc).join("；")}</li>`)
        .join("")}</ul>`;
  } else {
    els.discardBox.className = "discard-box safe";
    els.discardBox.innerHTML = `<h4>✅ 未检出高风险组合</h4>
      <div>所有引物在当前盐浓度（${r.salt_conc} M）与 Tm 阈值（${r.temperature}°C）下均可保留。</div>`;
  }

  els.hairpins.innerHTML = r.hairpins.length
    ? r.hairpins.map(renderHairpin).join("")
    : `<div class="empty" style="padding:10px 0">未检出超过阈值的发卡结构。</div>`;
  els.dimers.innerHTML = r.dimers.length
    ? r.dimers.map(renderDimer).join("")
    : `<div class="empty" style="padding:10px 0">未检出超过阈值的引物二聚体。</div>`;

  if (r.capped && (r.capped.hairpins || r.capped.dimers)) {
    els.dimers.insertAdjacentHTML(
      "beforebegin",
      `<p class="hint">注：另有 ${r.capped.hairpins} 个发卡 / ${r.capped.dimers} 个二聚体事件因数量上限被截断。</p>`
    );
  }
}

/* ---------------- 动作 ---------------- */

async function runScreen() {
  let parsed;
  try {
    parsed = JSON.parse(els.json.value);
  } catch (e) {
    setMsg("JSON 解析失败：" + e.message, "err");
    return;
  }
  let primers;
  try {
    primers = normalizePrimers(parsed);
  } catch (e) {
    setMsg(e.message, "err");
    return;
  }
  const payload = {
    name: els.jobName.value.trim() || null,
    salt_conc: parseFloat(els.salt.value),
    temperature: parseFloat(els.temp.value),
    min_bp: parseInt(els.minbp.value, 10),
    primers,
  };
  els.run.disabled = true;
  setMsg("筛选中…");
  try {
    const resp = await fetch("/api/screen", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const body = await resp.json();
    if (!resp.ok) {
      setMsg("请求失败（" + resp.status + "）：" + (body.detail || ""), "err");
      return;
    }
    setMsg(`完成 ✓ 任务 ID ${body.job_id}（已存入 SQLite）`, "ok");
    renderResult(body);
    loadJobs();
  } catch (e) {
    setMsg("网络错误：" + e.message, "err");
  } finally {
    els.run.disabled = false;
  }
}

function exportJson() {
  if (!lastResponse) return;
  const blob = new Blob([JSON.stringify(lastResponse, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `primer-screen-${lastResponse.job_id}.json`;
  a.click();
  URL.revokeObjectURL(a.href);
}

async function loadExample() {
  const ex = await (await fetch("/api/example")).json();
  els.jobName.value = ex.name || "";
  els.salt.value = ex.salt_conc;
  els.temp.value = ex.temperature;
  els.minbp.value = ex.min_bp;
  els.json.value = JSON.stringify(ex.primers, null, 2);
  setMsg("已载入示例数据", "ok");
}

function formatJson() {
  try {
    els.json.value = JSON.stringify(JSON.parse(els.json.value), null, 2);
    setMsg("JSON 已格式化", "ok");
  } catch (e) {
    setMsg("JSON 解析失败：" + e.message, "err");
  }
}

/* ---------------- 历史任务 ---------------- */

async function loadJobs() {
  const body = await (await fetch("/api/jobs")).json();
  if (!body.jobs.length) {
    els.jobList.innerHTML = `<div class="empty" style="padding:10px 0">暂无历史任务。</div>`;
    return;
  }
  els.jobList.innerHTML = body.jobs
    .map((j) => {
      const risk = (j.high_risk_hairpins + j.high_risk_dimers > 0)
        ? `<span class="badge high">高风险 ${j.high_risk_hairpins + j.high_risk_dimers}</span>`
        : `<span class="badge low">通过</span>`;
      return `<div class="job-row">
        <span class="jid" data-id="${j.job_id}" title="点击查看完整结果">${j.job_id}</span>
        ${risk}
        <span class="jname">${esc(j.name || "（未命名任务）")}</span>
        <span class="jmeta">
          <span>${j.created_at.replace("T", " ").replace("+00:00", " UTC")}</span>
          <span>${j.primer_count} 条引物</span>
          <span>[Na⁺] ${j.salt_conc} M</span>
          <span>Tm 阈值 ${j.temperature}°C</span>
          <span>事件 ${j.events_total}</span>
          <span>淘汰 ${j.discard.length}</span>
        </span>
        <button class="btn tiny danger" data-del="${j.job_id}" type="button">删除</button>
      </div>`;
    })
    .join("");

  els.jobList.querySelectorAll(".jid").forEach((node) =>
    node.addEventListener("click", () => openJob(node.dataset.id)));
  els.jobList.querySelectorAll("[data-del]").forEach((node) =>
    node.addEventListener("click", () => removeJob(node.dataset.del)));
}

async function openJob(jobId) {
  const job = await (await fetch(`/api/jobs/${jobId}`)).json();
  const req = job.request;
  els.jobName.value = job.name || "";
  els.salt.value = req.salt_conc;
  els.temp.value = req.temperature;
  els.minbp.value = req.min_bp;
  els.json.value = JSON.stringify(req.primers, null, 2);
  renderResult({ job_id: job.job_id, result: job.result });
  setMsg(`已载入历史任务 ${job.job_id}`, "ok");
  window.scrollTo({ top: 0, behavior: "smooth" });
}

async function removeJob(jobId) {
  if (!confirm(`确定删除任务 ${jobId}？`)) return;
  const resp = await fetch(`/api/jobs/${jobId}`, { method: "DELETE" });
  if (resp.ok) {
    setMsg(`已删除 ${jobId}`, "ok");
    loadJobs();
  }
}

/* ---------------- 绑定 ---------------- */

els.run.addEventListener("click", runScreen);
els.exportBtn.addEventListener("click", exportJson);
els.example.addEventListener("click", loadExample);
els.format.addEventListener("click", formatJson);
els.refreshJobs.addEventListener("click", loadJobs);

loadJobs();
