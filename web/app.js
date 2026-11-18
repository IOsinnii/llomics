// The page: tabs, the benchmark gallery, the build form and the payload inspector. All Python runs
// in worker.js; this file only sends {fn, args} to web/bridge.py and renders what comes back.
// Every value from a payload or an uploaded file is inserted as text, never as HTML.

// ------------------------------------------------------------------------------------------------
// Python in a worker
// ------------------------------------------------------------------------------------------------

const worker = new Worker("worker.js", { type: "module" });
const pending = new Map();
let nextId = 0;
let manifest = null;
const engine = document.getElementById("engine");

class PyError extends Error {}

function failed(text) {
  engine.textContent = "Python failed to load";
  engine.className = "engine failed";
  document.getElementById("gallery").replaceChildren(
    h("p", { class: "error" }, "Python could not be started in this browser. " + text),
  );
}

worker.onerror = (event) => {
  event.preventDefault();
  failed(event.message || "The worker script could not be loaded.");
};

worker.onmessage = (event) => {
  const m = event.data;
  if (m.type === "status") {
    engine.textContent = m.text + "…";
    return;
  }
  if (m.type === "ready") {
    manifest = m.manifest;
    onReady();
    return;
  }
  if (m.type === "failed") {
    failed(m.text);
    return;
  }
  const p = pending.get(m.id);
  if (!p) return;
  pending.delete(m.id);
  if (m.error) {
    p.reject(new Error(m.error));
    return;
  }
  const r = JSON.parse(m.out);
  if (r.ok) p.resolve(r.result);
  else p.reject(new PyError(r.error));
};

function py(fn, args = {}) {
  return new Promise((resolve, reject) => {
    const id = ++nextId;
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, fn, args });
  });
}

// ------------------------------------------------------------------------------------------------
// Small DOM helpers
// ------------------------------------------------------------------------------------------------

function h(tag, props = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(props || {})) {
    if (v === undefined || v === null || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of children.flat(Infinity)) {
    if (c === undefined || c === null || c === false) continue;
    el.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return el;
}

function fmt(x) {
  if (x === null || x === undefined) return "—";
  if (typeof x !== "number") return String(x);
  if (x === 0) return "0";
  const a = Math.abs(x);
  if (a < 1e-3 || a >= 1e6) return x.toExponential(2);
  if (Number.isInteger(x)) return x.toLocaleString("en");
  return String(Number(x.toPrecision(3)));
}

function kb(chars) {
  return chars < 1024 ? `${chars} characters` : `${(chars / 1024).toFixed(1)} KB`;
}

function download(name, text, type = "application/json") {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const a = h("a", { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function copyButton(label, getText) {
  const btn = h("button", { type: "button", class: "secondary small" }, label);
  btn.addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(getText());
      btn.textContent = "Copied";
    } catch {
      btn.textContent = "Copy failed";
    }
    setTimeout(() => (btn.textContent = label), 1500);
  });
  return btn;
}

const STATUS = {
  pass: ["✓", "pass"],
  warn: ["!", "warn"],
  fail: ["✕", "fail"],
  not_run: ["–", "not run"],
  not_applicable: ["○", "not applicable"],
};

function badge(status, label) {
  const [icon, text] = STATUS[status] || ["?", status];
  return h("span", { class: `status ${status}` }, h("span", { class: "icon", "aria-hidden": "true" }, icon), label ?? text);
}

function verdict(kind, text) {
  const icons = { good: "✓", bad: "✕", pending: "…", neutral: "→", warn: "!" };
  return h("span", { class: `status ${kind}` }, h("span", { class: "icon", "aria-hidden": "true" }, icons[kind]), text);
}

// What each gate reason means, for people reading the page.
const REASONS = {
  too_few_hits: "Fewer than 10 significant features: nothing to interpret.",
  technical_signature: "The hits are dominated by a technical panel (sex, stress, cell cycle, haemoglobin, mitochondrial or ribosomal genes).",
  signal_is_not_the_comparison: "The model's answer, agreed across runs, was that the signal is not the comparison asked about.",
};
const CHECKS = {
  replicates: "replicates per group",
  pvalue_distribution: "shape of the p-value histogram",
  batch_confounding: "group-by-batch confounding",
  multiple_testing: "multiple-testing scope",
  mean_variance: "mean-variance modelling",
  pseudoreplication: "cells used as replicates",
};

function explainReason(code) {
  if (REASONS[code]) return REASONS[code];
  const m = code.match(/^trust:([^:]+):(.+)$/);
  if (m) {
    const [, check, status] = m;
    const what = CHECKS[check] || check;
    if (status === "absent") return `The trust check on ${what} is missing from the payload.`;
    return `The trust check on ${what} is "${status.replace("_", " ")}", which the gate settings treat as blocking.`;
  }
  if (code.startsWith("no_answer")) return "Some model requests returned no usable answer.";
  if (code.startsWith("fact:")) return "A fact the gates need is missing from the payload.";
  return "";
}

function reasonList(reasons) {
  return h(
    "div",
    { class: "decision-reasons" },
    reasons.map((r) => h("p", { class: "explain" }, h("span", { class: "reason" }, r), " ", explainReason(r))),
  );
}

// ------------------------------------------------------------------------------------------------
// Tabs
// ------------------------------------------------------------------------------------------------

const TABS = ["problems", "build", "inspect"];

function showTab(name) {
  if (!TABS.includes(name)) name = "problems";
  for (const t of TABS) {
    document.getElementById("tab-" + t).setAttribute("aria-selected", String(t === name));
    document.getElementById("panel-" + t).hidden = t !== name;
  }
  if (location.hash !== "#" + name) history.replaceState(null, "", "#" + name);
}

for (const t of TABS) document.getElementById("tab-" + t).addEventListener("click", () => showTab(t));
document.querySelector(".tabs").addEventListener("keydown", (e) => {
  if (e.key !== "ArrowRight" && e.key !== "ArrowLeft") return;
  const current = TABS.findIndex((t) => document.getElementById("tab-" + t).getAttribute("aria-selected") === "true");
  const next = TABS[(current + (e.key === "ArrowRight" ? 1 : TABS.length - 1)) % TABS.length];
  showTab(next);
  document.getElementById("tab-" + next).focus();
});
showTab(location.hash.slice(1));
window.addEventListener("hashchange", () => showTab(location.hash.slice(1)));

// ------------------------------------------------------------------------------------------------
// Start-up
// ------------------------------------------------------------------------------------------------

async function onReady() {
  const info = await py("info");
  engine.textContent = `Python ${info.python} ready in your browser · llomics ${info.llomics}`;
  engine.className = "engine ready";
  document.getElementById("versions").textContent =
    `llomics ${info.llomics}, contract ${info.contract}, Python ${info.python} (Pyodide).`;
  fillExamples();
  updateBuildButton();
  renderGallery().catch((e) => showError(document.getElementById("gallery"), e));
}

function showError(container, error) {
  container.replaceChildren(h("p", { class: "error" }, error instanceof PyError ? error.message : "Error: " + error.message));
}

// ------------------------------------------------------------------------------------------------
// Planted problems
// ------------------------------------------------------------------------------------------------

const TITLES = {
  clean_interferon: "Clean: interferon response",
  clean_hypoxia: "Clean: hypoxia experiment",
  null: "No real difference",
  batch_aliased: "Batch aliased with group",
  no_batch_info: "No batch information",
  unmodelled_batch: "Batch effect not modelled",
  sex_imbalance: "Sex imbalance between groups",
  stress_handling: "Sample-handling stress",
  composition: "Cell-type composition shift",
  two_replicates: "Two replicates per group",
  technical_replicates: "Technical replicates counted as biological",
};
const LAYER = { trust: "a trust check", fact: "a fact gate", model: "the model", none: "nothing (known blind spot)" };

function gateVerdict(s) {
  const w = s.without_model;
  if (w === "caught") return verdict("good", "Sent to review by code, before any model");
  if (w.startsWith("ok so far")) return verdict("neutral", "Passes every gate: goes to the model");
  if (w.startsWith("pending")) return verdict("pending", "Passes every gate: only the model can catch this");
  if (w.startsWith("not caught")) return verdict("warn", "Not caught: a known blind spot of the payload");
  return verdict("bad", w);
}

async function renderGallery() {
  const { summary, scenarios } = await py("gallery");
  const tile = (big, total, what) =>
    h("div", { class: "tile" }, h("div", { class: "big" }, String(big), total !== null ? h("small", {}, ` of ${total}`) : null), h("div", { class: "what" }, what));
  document.getElementById("scoreboard").replaceChildren(
    tile(summary.caught_by_gates, summary.n_problems, "planted problems sent to review by code alone"),
    tile(summary.false_alarms, summary.n_clean, "false alarms on the clean analyses"),
    tile(summary.pending_model, null, "problem that only the model's judgment can catch"),
    tile(summary.blind_spots, null, "known blind spot, documented in the README"),
  );
  const note = summary.answers
    ? [`Model answers replayed from `, h("code", {}, `benchmark/answers/${summary.answers}`), `, gated again by this version's code: ${summary.with_model_correct} of ${summary.with_model_total} scenarios decided correctly.`]
    : ["Model answers for the benchmark have not been recorded yet, so this page shows the gates that need no model. Once ", h("code", {}, "python benchmark/score.py --ask"), " has been run and its answers committed, they are replayed here."];
  document.getElementById("scoreboard").after(h("p", { class: "answers-note" }, note));

  document.getElementById("gallery").replaceChildren(...scenarios.map(card));
}

function card(s) {
  const model = s.model;
  return h(
    "article",
    { class: "card" },
    h("div", { class: "head" }, h("h3", {}, TITLES[s.id] || s.id), h("span", { class: "expect" }, `correct decision: ${s.expected}`)),
    h("p", { class: "truth" }, s.truth),
    h(
      "div",
      { class: "genes" },
      `${fmt(s.n_hits)} hits (${fmt(s.n_up)} up, ${fmt(s.n_down)} down)`,
      s.top_up.length ? [". Top up: ", h("code", {}, s.top_up.slice(0, 5).join(", "))] : null,
    ),
    h(
      "div",
      { class: "checks" },
      s.trust.map((c) => badge(c.status, CHECKS[c.check] || c.check)),
    ),
    h(
      "div",
      { class: "decision" },
      h("span", { class: "label-sm" }, s.expected === "auto" ? "Gates (nothing to catch: it is clean)" : `Gates (should be caught by ${LAYER[s.caught_by ?? "none"]})`),
      gateVerdict(s),
      s.reasons.length ? reasonList(s.reasons) : null,
      model
        ? h(
            "div",
            {},
            h("span", { class: "label-sm" }, `Model (${model.model}, ${model.n_runs} runs)`),
            h(
              "p",
              { class: "explain" },
              verdict(model.verdict === "correct" ? "good" : "bad", `${model.status}: ${model.verdict}`),
              model.dominant_signal.top ? ` · dominant signal "${model.dominant_signal.top}" in ${Math.round(100 * model.dominant_signal.share)}% of runs` : "",
            ),
          )
        : null,
    ),
    h(
      "div",
      { class: "foot" },
      h("span", { class: "expect" }, `payload ${kb(s.payload_chars)} · model request ${kb(s.request_chars)}`),
      h("button", { type: "button", class: "secondary small", onclick: () => openInInspector(s.payload_path) }, "Inspect payload"),
    ),
  );
}

// ------------------------------------------------------------------------------------------------
// Inspector (shared by the build and inspect tabs)
// ------------------------------------------------------------------------------------------------

function renderInspection(container, r, opts = {}) {
  if (r.parse_error) {
    container.replaceChildren(h("div", { class: "panel" }, h("p", { class: "error" }, r.parse_error)));
    return;
  }
  const parts = [];
  const valid = r.errors.length === 0;
  const nUnits = r.units ? r.units.length : 0;
  parts.push(
    h(
      "div",
      { class: "panel" },
      h(
        "div",
        { class: "banner" },
        valid
          ? verdict("good", r.warnings.length ? `Valid against contract ${r.llomics}, with ${r.warnings.length} warning(s)` : `Valid against contract ${r.llomics}`)
          : verdict("bad", `${r.errors.length} contract error(s)`),
        h(
          "span",
          { class: "size" },
          "Payload ",
          h("strong", {}, kb(r.chars)),
          nUnits ? ` · ${nUnits} unit${nUnits > 1 ? "s" : ""}` : "",
          opts.inputChars ? [" · built from a ", h("strong", {}, kb(opts.inputChars)), ` table of ${fmt(opts.inputRows)} rows`] : "",
        ),
        opts.payloadText ? h("span", { class: "toolbar" },
          h("button", { type: "button", class: "small", onclick: () => download(opts.fileName || "payload.llomics.json", opts.payloadText) }, "Download payload"),
          copyButton("Copy JSON", () => opts.payloadText),
        ) : null,
      ),
      r.errors.length ? h("ul", { class: "messages error" }, r.errors.map((e) => h("li", {}, e))) : null,
      r.warnings.length ? h("ul", { class: "messages" }, r.warnings.map((e) => h("li", {}, e))) : null,
    ),
  );
  if (!valid) {
    container.replaceChildren(...parts);
    return;
  }
  if (opts.extra) parts.push(...opts.extra);

  const run = r.run;
  const ds = r.dataset;
  parts.push(
    h(
      "div",
      { class: "panel" },
      h("h3", {}, "Run card"),
      h(
        "div",
        { class: "grid2" },
        h(
          "dl",
          { class: "kv" },
          h("dt", {}, "Question"), h("dd", {}, run.question || "—"),
          h("dt", {}, "Method"), h("dd", {}, r.method),
          h("dt", {}, "Tool"), h("dd", {}, run.tool || "—"),
          h("dt", {}, "Analysis"), h("dd", {}, run.analysis_id || "—"),
          h("dt", {}, "Observations"), h("dd", {}, `${fmt(ds.n_observations)} ${ds.observation_type || ""}s`),
          h("dt", {}, "Features"), h("dd", {}, fmt(ds.n_features)),
        ),
        h(
          "div",
          {},
          h("span", { class: "label-sm" }, "Steps"),
          h("ul", { class: "plain" }, (run.steps || []).map((s) => h("li", {}, s))),
          h("span", { class: "label-sm" }, "Filters"),
          h("ul", { class: "plain" }, (run.filters || []).map((f) => h("li", {}, f.rule, f.why ? h("span", { class: "muted" }, ` (${f.why})`) : ""))),
        ),
      ),
      h("h4", {}, "Trust checks for the whole run"),
      trustTable(r.trust),
    ),
  );

  const unitPanel = h("div", { class: "panel" });
  const unitBody = h("div");
  const showUnit = (i) => unitBody.replaceChildren(...renderUnit(r.units[i], r));
  if (nUnits > 1) {
    const select = h(
      "select",
      { "aria-label": "Unit" },
      r.units.map((u, i) =>
        h("option", { value: String(i) }, `${u.unit.unit_id}${u.gate ? ` · ${u.gate.reasons.length ? "review" : "to the model"}` : ""}`),
      ),
    );
    select.addEventListener("change", () => showUnit(Number(select.value)));
    unitPanel.append(h("div", { class: "unit-picker" }, h("h3", { style: "margin:0" }, "Units"), select, h("span", { class: "muted" }, "one request per unit")));
  } else {
    unitPanel.append(h("h3", {}, "Unit"));
  }
  unitPanel.append(unitBody);
  if (nUnits) showUnit(0);
  parts.push(unitPanel);
  container.replaceChildren(...parts);
}

function trustTable(checks) {
  if (!checks || !checks.length) return h("p", { class: "muted" }, "No trust checks.");
  return h(
    "table",
    { class: "trust-table" },
    h("thead", {}, h("tr", {}, h("th", {}, "Check"), h("th", {}, "Status"), h("th", {}, "Note"))),
    h("tbody", {}, checks.map((c) => h("tr", {}, h("td", {}, c.check), h("td", {}, badge(c.status)), h("td", {}, c.note || "")))),
  );
}

const PREFERRED = ["feature", "symbol", "effect", "stat", "p_adj", "pct_in", "pct_out", "ave_expr"];

function featureTable(rows, caption) {
  if (!rows || !rows.length) return h("p", { class: "muted" }, `${caption}: none eligible.`);
  const keys = [...new Set(rows.flatMap((r) => Object.keys(r)))];
  const cols = [...PREFERRED.filter((k) => keys.includes(k)), ...keys.filter((k) => !PREFERRED.includes(k))];
  return h(
    "table",
    {},
    h("caption", { class: "label-sm", style: "text-align:left" }, caption),
    h("thead", {}, h("tr", {}, cols.map((c) => h("th", { class: c === "feature" || c === "symbol" ? "" : "num" }, c)))),
    h("tbody", {}, rows.map((r) => h("tr", {}, cols.map((c) => h("td", { class: typeof r[c] === "number" ? "num" : "" }, fmt(r[c])))))),
  );
}

function panelTable(panel) {
  const entries = Object.entries(panel || {});
  if (!entries.length) return null;
  return h(
    "table",
    {},
    h("caption", { class: "label-sm", style: "text-align:left" }, "Panel (reported whether or not it ranks; — means not measured)"),
    h("thead", {}, h("tr", {}, h("th", {}, "feature"), h("th", {}, "value"))),
    h(
      "tbody",
      {},
      entries.map(([name, v]) =>
        h(
          "tr",
          {},
          h("td", {}, name),
          h(
            "td",
            {},
            v === null ? "— not measured" : typeof v === "object" ? Object.entries(v).map(([k, x]) => `${k} ${fmt(x)}`).join(" · ") : fmt(v),
          ),
        ),
      ),
    ),
  );
}

function histogram(counts) {
  const tip = document.getElementById("tooltip");
  const max = Math.max(...counts, 1);
  const total = counts.reduce((a, b) => a + b, 0);
  const bars = counts.map((n, i) => {
    const label = `p ${(i / 10).toFixed(1)}–${((i + 1) / 10).toFixed(1)}: ${n.toLocaleString("en")} features`;
    const bar = h("div", { class: "bar", style: `height:${(100 * n) / max}%`, tabindex: "0", "aria-label": label }, h("span", { class: "hit" }));
    const show = () => {
      const box = bar.getBoundingClientRect();
      tip.textContent = label;
      tip.hidden = false;
      tip.style.left = `${Math.min(box.left, window.innerWidth - 200)}px`;
      tip.style.top = `${box.top - 30}px`;
    };
    const hide = () => (tip.hidden = true);
    bar.addEventListener("mouseenter", show);
    bar.addEventListener("focus", show);
    bar.addEventListener("mouseleave", hide);
    bar.addEventListener("blur", hide);
    return bar;
  });
  return h(
    "figure",
    { class: "hist", style: "margin:0" },
    h("div", { class: "bars", "aria-hidden": "false" }, bars),
    h("div", { class: "axis", "aria-hidden": "true" }, h("span", {}, "p = 0"), h("span", {}, "0.5"), h("span", {}, "1")),
    h("figcaption", { class: "caption" }, `Raw p-values of all ${total.toLocaleString("en")} tested features, in tenths. A flat floor with a spike near 0 is the healthy shape.`),
    h(
      "table",
      { class: "sr-only" },
      h("tr", {}, h("th", {}, "p-value bin"), h("th", {}, "features")),
      counts.map((n, i) => h("tr", {}, h("td", {}, `${(i / 10).toFixed(1)}–${((i + 1) / 10).toFixed(1)}`), h("td", {}, String(n)))),
    ),
  );
}

const CONTRAST_FACTS = [
  ["n_hits", "significant features"],
  ["n_up", "up"],
  ["n_down", "down"],
  ["n_near_miss", "near misses"],
  ["pi0", "estimated share of nulls (π₀)"],
  ["n_technical_flags", "technical flags"],
];

function renderUnit(entry, r) {
  const u = entry.unit;
  const out = [];
  out.push(
    h(
      "dl",
      { class: "kv" },
      h("dt", {}, "Unit"), h("dd", {}, `${u.unit_id} (${u.unit_type}) · ${kb(entry.chars)}`),
      u.definition ? [h("dt", {}, "Definition"), h("dd", {}, u.definition)] : null,
      u.size ? [h("dt", {}, "Size"), h("dd", {}, Object.entries(u.size).map(([k, v]) => `${k} ${typeof v === "object" ? JSON.stringify(v) : fmt(v)}`).join(" · "))] : null,
    ),
  );

  out.push(h("h4", {}, "Decision before asking a model"));
  if (entry.gate) {
    out.push(
      entry.gate.reasons.length
        ? h("div", {}, verdict("bad", "Goes to review: a gate failed"), reasonList(entry.gate.reasons))
        : verdict("neutral", "No trust or fact gate failed: this unit would be sent to the model"),
    );
  } else {
    out.push(h("p", { class: "explain" }, r.questions_note));
  }

  if (u.trust && u.trust.length) {
    out.push(h("h4", {}, "Trust checks for this unit"), trustTable(u.trust));
  }

  const facts = u.facts || {};
  const shown =
    u.unit_type === "contrast"
      ? CONTRAST_FACTS.filter(([k]) => k in facts)
      : Object.entries(facts).filter(([, v]) => typeof v === "number").map(([k]) => [k, k.replaceAll("_", " ")]);
  out.push(h("h4", {}, "Facts computed in code"));
  out.push(h("div", { class: "facts" }, shown.map(([k, label]) => h("div", { class: "fact" }, h("div", { class: "v" }, fmt(facts[k])), h("div", { class: "k" }, label)))));
  if (Array.isArray(facts.p_histogram)) out.push(histogram(facts.p_histogram));
  if (facts.technical_flags && facts.technical_flags.length) {
    out.push(h("p", { class: "explain" }, "Technical flags: ", facts.technical_flags.map((f) => h("span", { class: "reason" }, typeof f === "string" ? f : JSON.stringify(f)))));
  }

  out.push(h("h4", {}, "What the model would read"));
  const tf = u.top_features;
  if (Array.isArray(tf)) out.push(featureTable(tf, "Top features"));
  else if (tf) out.push(h("div", { class: "grid2" }, featureTable(tf.up, "Top up"), featureTable(tf.down, "Top down")));
  if (u.tail) {
    out.push(
      h("p", { class: "explain" }, `Showing ${fmt(u.tail.n_shown)} of ${u.tail.n_eligible === null ? "an unrecorded number of" : fmt(u.tail.n_eligible)} eligible features. Rule: ${u.tail.rule}.`),
    );
  }
  const pt = panelTable(u.panel);
  if (pt) out.push(pt);
  if (u.enrichment) {
    const sets = [...(u.enrichment.up || []).map((s) => ({ direction: "up", ...s })), ...(u.enrichment.down || []).map((s) => ({ direction: "down", ...s }))];
    out.push(h("h4", {}, `Enrichment (${fmt(u.enrichment.n_sets_tested)} sets tested)`), sets.length ? featureTable(sets, "Sets") : h("p", { class: "muted" }, "No set passed."));
  }

  out.push(h("details", {}, h("summary", {}, "Unit JSON"), h("pre", {}, JSON.stringify(u, null, 2))));
  if (entry.request) {
    let user = entry.request.user;
    try {
      user = JSON.stringify(JSON.parse(user), null, 2);
    } catch {}
    out.push(
      h(
        "details",
        {},
        h("summary", {}, `Request that would be sent to the model (${kb(entry.request.system.length + entry.request.user.length)}; nothing is sent from this page)`),
        h("span", { class: "label-sm" }, "System prompt"),
        h("pre", {}, entry.request.system),
        h("span", { class: "label-sm" }, "User message"),
        h("pre", {}, user),
        h("span", { class: "label-sm" }, "Answer schema"),
        h("pre", {}, JSON.stringify(entry.request.schema, null, 2)),
      ),
    );
  }
  return out;
}

// ------------------------------------------------------------------------------------------------
// Inspect tab
// ------------------------------------------------------------------------------------------------

const inspectResult = document.getElementById("inspect-result");
const inspectSelect = document.getElementById("inspect-example");

function fillExamples() {
  const payloads = manifest.files.filter((f) => f.endsWith(".llomics.json") || f === "examples/pbmc10k/payload.json");
  for (const f of payloads) {
    const label = f.startsWith("benchmark/")
      ? `Benchmark: ${TITLES[f.split("/").pop().replace(".llomics.json", "")] || f}`
      : "PBMC 10k: 19 clusters from a Scanpy pipeline";
    inspectSelect.append(h("option", { value: f }, label));
  }
}

async function inspectText(text, name) {
  inspectResult.replaceChildren(h("p", { class: "loading" }, "Checking…"));
  try {
    const r = await py("inspect", { text });
    renderInspection(inspectResult, r, { payloadText: text, fileName: name });
  } catch (e) {
    showError(inspectResult, e);
  }
}

async function openInInspector(path) {
  showTab("inspect");
  inspectSelect.value = path;
  const text = await (await fetch("repo/" + path)).text();
  await inspectText(text, path.split("/").pop());
  window.scrollTo({ top: document.querySelector(".tabs").offsetTop - 8, behavior: "smooth" });
}

inspectSelect.addEventListener("change", () => inspectSelect.value && openInInspector(inspectSelect.value));
document.getElementById("inspect-file").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (file) await inspectText(await file.text(), file.name);
});
document.getElementById("inspect-paste-toggle").addEventListener("click", () => {
  const box = document.getElementById("inspect-paste");
  box.hidden = !box.hidden;
});
document.getElementById("inspect-run").addEventListener("click", () =>
  inspectText(document.getElementById("inspect-text").value, "pasted.llomics.json"),
);

// ------------------------------------------------------------------------------------------------
// Build tab
// ------------------------------------------------------------------------------------------------

const form = document.getElementById("build-form");
const buildButton = document.getElementById("build-button");
const buildError = document.getElementById("build-error");
const buildResult = document.getElementById("build-result");
const tableStatus = document.getElementById("table-status");
let table = null; // {text, name}

function updateBuildButton() {
  buildButton.disabled = !(manifest && table);
}

function setTable(text, name, note) {
  table = { text, name };
  const rows = Math.max(0, text.trimEnd().split("\n").length - 1);
  tableStatus.replaceChildren(
    `Loaded ${name}: ${rows.toLocaleString("en")} rows, ${kb(text.length)}${note ? ` (${note})` : ""}. `,
    h("a", { href: "#build", onclick: (e) => { e.preventDefault(); download(name, text, "text/tab-separated-values"); } }, "Download it"),
  );
  updateBuildButton();
}

document.getElementById("table-file").addEventListener("change", async (e) => {
  const file = e.target.files[0];
  if (file) setTable(await file.text(), file.name);
});

document.getElementById("load-example").addEventListener("click", async () => {
  if (!manifest) return;
  const ex = await py("example");
  for (const [k, v] of Object.entries(ex.form)) {
    if (form.elements[k]) form.elements[k].value = v;
  }
  document.getElementById("table-file").value = "";
  setTable(ex.text, ex.filename, "simulated: about 50 interferon-stimulated genes up, nothing else changed");
});

form.addEventListener("submit", async (e) => {
  e.preventDefault();
  if (!table) return;
  buildError.textContent = "";
  buildButton.disabled = true;
  buildResult.replaceChildren(h("p", { class: "loading" }, "Building…"));
  const fields = Object.fromEntries(new FormData(form).entries());
  const gmtFile = document.getElementById("gmt-file").files[0];
  try {
    const r = await py("build", {
      form: fields,
      table_text: table.text,
      filename: table.name,
      gmt_text: gmtFile ? await gmtFile.text() : "",
      gmt_name: gmtFile ? gmtFile.name : "",
    });
    const t = r.table;
    const mapping = h(
      "div",
      { class: "panel" },
      h("h3", {}, "How the table was read"),
      h(
        "p",
        { class: "mapping" },
        "Columns: ",
        Object.entries(t.mapping).map(([std, col], i) => [i ? ", " : "", h("code", {}, col || "(unnamed first column)"), " → ", std]),
        t.unmapped.length ? `. Not used: ${t.unmapped.join(", ")}.` : ".",
        t.tool_guess ? ` Looks like ${t.tool_guess} output.` : "",
      ),
      t.notes.length ? h("ul", { class: "messages" }, t.notes.map((n) => h("li", {}, n))) : null,
      h(
        "details",
        {},
        h("summary", {}, "The same call in Python, for your pipeline"),
        h("pre", {}, r.python),
        copyButton("Copy code", () => r.python),
      ),
    );
    const id = fields.analysis_id || "payload";
    renderInspection(buildResult, r.inspect, {
      payloadText: r.payload_text,
      fileName: `${id}.llomics.json`,
      inputChars: r.input_chars,
      inputRows: t.n_rows,
      extra: [mapping],
    });
  } catch (err) {
    buildResult.replaceChildren();
    buildError.textContent = err instanceof PyError ? err.message : "Error: " + err.message;
  } finally {
    updateBuildButton();
  }
});
