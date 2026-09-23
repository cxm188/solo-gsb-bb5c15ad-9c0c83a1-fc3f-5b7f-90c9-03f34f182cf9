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
  sensSalts: $("#sensSalts"),
  sensTms: $("#sensTms"),
  sensRun: $("#btnSens"),
  sensExport: $("#btnSensExport"),
  sensMsg: $("#sensMsg"),
  sensResult: $("#sensResult"),
};

let lastResponse = null; // { job_id, result }
let currentJobId = null; // 当前结果所属任务（敏感性分析复用其引物）
let lastSens = null; // 最近一次敏感性分析响应

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
  currentJobId = data.job_id;
  const r = data.result;
  els.empty.hidden = true;
  els.body.hidden = false;
  els.exportBtn.disabled = false;
  resetSens(r);

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

/* ---------------- 条件敏感性分析 ---------------- */

const RISK_SHORT = { critical: "极", high: "高", medium: "中", low: "·" };

/* 新结果渲染时重置敏感性区块，并按基准参数预填默认网格（不覆盖用户已填内容）。 */
function resetSens(r) {
  lastSens = null;
  els.sensExport.disabled = true;
  els.sensResult.innerHTML = "";
  els.sensMsg.textContent = "";
  els.sensMsg.className = "msg";
  if (!els.sensSalts.value.trim()) {
    const around = [r.salt_conc * 0.6, r.salt_conc, r.salt_conc * 1.5]
      .map((x) => Math.min(1.5, Math.max(0.001, +x.toFixed(3))));
    els.sensSalts.value = [...new Set(around)].join(", ");
  }
  if (!els.sensTms.value.trim()) {
    const around = [r.temperature - 5, r.temperature, r.temperature + 5]
      .map((x) => Math.min(95, Math.max(0, +x.toFixed(1))));
    els.sensTms.value = [...new Set(around)].join(", ");
  }
}

function setSensMsg(text, kind = "") {
  els.sensMsg.textContent = text;
  els.sensMsg.className = "msg" + (kind ? " " + kind : "");
}

function parseNumList(text, label) {
  const parts = text.split(/[\s,，、;；]+/).filter(Boolean);
  if (!parts.length) throw new Error(`请填写${label}`);
  const nums = parts.map((p) => {
    const x = Number(p);
    if (!Number.isFinite(x)) throw new Error(`${label}含非法数值：${p}`);
    return x;
  });
  if (nums.length > 8) throw new Error(`${label}最多 8 项（当前 ${nums.length} 项）`);
  if (new Set(nums).size !== nums.length) throw new Error(`${label}存在重复值，请去重`);
  return nums;
}

async function runSensitivity() {
  if (!currentJobId) {
    setSensMsg("请先执行筛选或载入历史任务", "err");
    return;
  }
  let salts, tms;
  try {
    salts = parseNumList(els.sensSalts.value, "盐浓度列表");
    tms = parseNumList(els.sensTms.value, "Tm 阈值列表");
    if (salts.length * tms.length > 64) {
      throw new Error(`组合共 ${salts.length * tms.length} 格，超过上限 64 格`);
    }
  } catch (e) {
    setSensMsg(e.message, "err");
    return;
  }
  els.sensRun.disabled = true;
  setSensMsg("分析中…");
  try {
    const resp = await fetch(`/api/jobs/${currentJobId}/sensitivity`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ salt_concs: salts, tm_thresholds: tms }),
    });
    const body = await resp.json();
    if (!resp.ok) {
      const detail = typeof body.detail === "string"
        ? body.detail
        : (body.detail || []).map((d) => d.msg).join("；");
      setSensMsg(`请求失败（${resp.status}）：${detail || "参数非法"}`, "err");
      return;
    }
    lastSens = body;
    els.sensExport.disabled = false;
    setSensMsg(`完成 ✓ ${body.cells.length} 格 × ${body.events_considered} 个对齐事件`, "ok");
    renderSensitivity(body);
  } catch (e) {
    setSensMsg("网络错误：" + e.message, "err");
  } finally {
    els.sensRun.disabled = false;
  }
}

/* 高风险数 → 热力色（绿→黄→红）。 */
function heatColor(ratio) {
  const hue = 140 - 140 * Math.min(1, ratio);
  return `hsl(${hue}, 65%, 32%)`;
}

function locusText(ev) {
  const L = ev.locus;
  if (ev.type === "hairpin") {
    return `${L.primer_id} 5′臂${L.arm_5prime[0]}–${L.arm_5prime[1]} × 3′臂${L.arm_3prime[0]}–${L.arm_3prime[1]}`;
  }
  return `${L.primer_a} ${L.interval_a[0]}–${L.interval_a[1]} × ${L.primer_b} ${L.interval_b[0]}–${L.interval_b[1]}`;
}

function renderSensitivity(data) {
  const nT = data.tm_thresholds.length;
  const maxHigh = Math.max(1, ...data.cells.map((c) => c.high_risk_total));

  let html = `<div class="sens-summary">
    基准条件：[Na⁺] <b>${data.base.salt_conc} M</b> · Tm 阈值 <b>${data.base.temperature}°C</b> · min_bp <b>${data.base.min_bp}</b>；
    共 <b>${data.cells.length}</b> 格，对齐事件 <b>${data.events_considered}</b> 个` +
    (data.events_capped ? `（另有 ${data.events_capped} 个事件超出上限未纳入）` : "") +
    `。点击格子查看淘汰引物与风险变化。</div>`;

  html += `<div class="heat-wrap"><table class="heat"><thead><tr><th>[Na⁺] ＼ Tm</th>` +
    data.tm_thresholds.map((t) => `<th>${t}°C</th>`).join("") + `</tr></thead><tbody>`;
  data.salt_concs.forEach((s, ri) => {
    html += `<tr><th>${s} M</th>`;
    data.tm_thresholds.forEach((t, ci) => {
      const idx = ri * nT + ci;
      const c = data.cells[idx];
      const star = c.first_critical_events
        ? `<span class="star" title="${c.first_critical_events} 个事件在此格首次达到 critical">★${c.first_critical_events}</span>`
        : "";
      html += `<td class="hcell${c.is_base ? " base" : ""}" data-cell="${idx}"
        style="background:${heatColor(c.high_risk_total / maxHigh)}"
        title="高风险 ${c.high_risk_total}（critical ${c.critical_events}）· 淘汰 ${c.discard.length}">
        <div class="hnum">${c.high_risk_total}</div>
        <div class="hsub">极${c.critical_events} 淘${c.discard.length}</div>${star}</td>`;
    });
    html += `</tr>`;
  });
  html += `</tbody></table></div><div id="sensCellDetail" class="sens-detail"></div>`;

  // 事件对齐矩阵（折叠，最多展示前 150 行）
  const MAX_ROWS = 150;
  const rows = data.events.slice(0, MAX_ROWS);
  html += `<details class="sens-matrix"><summary>事件对齐矩阵（${data.events.length} 个事件 × ${data.cells.length} 格，按最严苛条件排序` +
    (data.events.length > MAX_ROWS ? `，仅显示前 ${MAX_ROWS} 行` : "") + `）</summary>`;
  html += `<div class="heat-wrap"><table class="heat matrix"><thead><tr><th>事件（类型 · 引物 · 配对区间）</th>` +
    data.cells.map((c) => `<th>${c.salt_conc}M<br>${c.tm_threshold}°C</th>`).join("") +
    `<th>基准</th></tr></thead><tbody>`;
  rows.forEach((ev) => {
    const icon = ev.type === "hairpin" ? "🔻" : "🔗";
    html += `<tr><td class="mlabel" title="${esc(ev.key)}">${icon} ${esc(locusText(ev))}</td>`;
    ev.cell_risks.forEach((risk, i) => {
      const star = ev.first_critical_cell === i ? `<span class="star in-cell">★</span>` : "";
      html += `<td class="rcell ${risk}">${RISK_SHORT[risk]}${star}</td>`;
    });
    html += `<td class="rcell ${ev.base_risk} basecol">${RISK_SHORT[ev.base_risk]}</td></tr>`;
  });
  html += `</tbody></table></div></details>`;

  els.sensResult.innerHTML = html;
  els.sensResult.querySelectorAll(".hcell").forEach((td) =>
    td.addEventListener("click", () => showCellDetail(data, +td.dataset.cell)));
}

function showCellDetail(data, idx) {
  const c = data.cells[idx];
  const d = c.delta;
  const sign = (n) => (n > 0 ? `+${n}` : `${n}`);
  const parts = [];
  if (d.events_upgraded) parts.push(`升级 ${d.events_upgraded}`);
  if (d.events_downgraded) parts.push(`降级 ${d.events_downgraded}`);
  if (d.events_new) parts.push(`新出现 ${d.events_new}`);
  if (d.events_gone) parts.push(`消失 ${d.events_gone}`);
  const deltaEvents = parts.length ? parts.join(" · ") : "事件风险无变化";
  const discardHtml = c.discard.length
    ? `<span class="badge high">淘汰 ${c.discard.length}</span> ${c.discard.map(esc).join("，")}`
    : `<span class="badge low">无淘汰</span>`;
  const added = d.discard_added.length ? `新增淘汰：${d.discard_added.map(esc).join("，")}` : "";
  const removed = d.discard_removed.length ? `解除淘汰：${d.discard_removed.map(esc).join("，")}` : "";

  document.getElementById("sensCellDetail").innerHTML = `
    <h4>条件格：[Na⁺] ${c.salt_conc} M × Tm 阈值 ${c.tm_threshold}°C${c.is_base ? ' <span class="badge tag">基准条件</span>' : ""}</h4>
    <div class="sens-detail-grid">
      <span>高风险事件 <b>${c.high_risk_total}</b>（发卡 ${c.high_risk_hairpins} / 二聚体 ${c.high_risk_dimers}，critical ${c.critical_events}，相对基准 ${sign(d.high_risk_total)}）</span>
      <span>相对基准：${esc(deltaEvents)}</span>
      <span>${discardHtml}</span>
      ${added ? `<span>➕ ${added}</span>` : ""}
      ${removed ? `<span>➖ ${removed}</span>` : ""}
      ${c.first_critical_events ? `<span>★ ${c.first_critical_events} 个事件在本格首次达到 critical</span>` : ""}
    </div>`;
}

function exportSensJson() {
  if (!lastSens) return;
  const blob = new Blob([JSON.stringify(lastSens, null, 2)], { type: "application/json" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = `primer-sensitivity-${lastSens.job_id}.json`;
  a.click();
  URL.revokeObjectURL(a.href);
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
els.sensRun.addEventListener("click", runSensitivity);
els.sensExport.addEventListener("click", exportSensJson);

loadJobs();
