import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import test from "node:test";

async function render() {
  const workerUrl = new URL("../dist/server/index.js", import.meta.url);
  workerUrl.searchParams.set("test", `${process.pid}-${Date.now()}`);
  const { default: worker } = await import(workerUrl.href);

  return worker.fetch(
    new Request("http://localhost/", { headers: { accept: "text/html" } }),
    { ASSETS: { fetch: async () => new Response("Not found", { status: 404 }) } },
    { waitUntil() {}, passThroughOnException() {} },
  );
}

test("server-renders the ODE-driven organ heatmap workspace", async () => {
  const response = await render();
  assert.equal(response.status, 200);
  assert.match(response.headers.get("content-type") ?? "", /^text\/html\b/i);

  const html = await response.text();
  assert.match(html, /<title>SINEUP Delivery Atlas<\/title>/i);
  assert.match(html, /疾病驱动的 AAV 空间递送设计/);
  assert.match(html, /人体多区域/);
  assert.match(html, /小鼠器官级/);
  assert.match(html, /PBPK–空间 CNS–SINEUP ODE/);
  assert.match(html, /1328 个模型输出/);
  assert.match(html, /疾病设计空间/);
  assert.match(html, /正在载入人体多区域 ODE 数据/);
  assert.doesNotMatch(html, /codex-preview|Your site is taking shape/);
});

test("includes lazy human model shell and social metadata", async () => {
  const response = await render();
  const html = await response.text();

  assert.match(html, /正在载入人体多区域 ODE 数据/);
  assert.match(html, /EN/);
  assert.match(html, /property="og:image" content="http:\/\/localhost(?::3000)?\/og\.png"/);
});

test("exports mass-balanced human multiregion trajectories", async () => {
  const file = new URL("../public/data/human-spatial-results.json", import.meta.url);
  const payload = JSON.parse(await readFile(file, "utf8"));
  assert.equal(payload.administration_routes.length, 6);
  assert.deepEqual(payload.administration_routes.map((route) => route.route_id), [
    "iv", "intrathecal", "intramuscular", "intracisternal", "intracerebroventricular", "inhaled",
  ]);
  assert.ok(payload.administration_routes.every((route) => route.capsids.length === 8));
  assert.equal(payload.region_ids.length, 24);
  assert.equal(payload.time_h.length, 287);
  assert.equal(payload.max_time_days, 730);
  assert.equal(payload.body_weight_kg, 70);
  assert.equal(payload.dose_vg_per_kg, 4e13);
  assert.equal(payload.dose_vg, 2.8e15);
  assert.equal(payload.state_count, 301);
  assert.match(payload.equation_family, /Q-PS-Kp-J_res-J_deg/);
  assert.equal(payload.default_route_id, "iv");
  assert.equal(payload.effective_flow_scale, 0.05);
  const iv = payload.administration_routes.find((route) => route.route_id === "iv");
  const intrathecal = payload.administration_routes.find((route) => route.route_id === "intrathecal");
  const intramuscular = payload.administration_routes.find((route) => route.route_id === "intramuscular");
  const intracisternal = payload.administration_routes.find((route) => route.route_id === "intracisternal");
  const intracerebroventricular = payload.administration_routes.find((route) => route.route_id === "intracerebroventricular");
  const inhaled = payload.administration_routes.find((route) => route.route_id === "inhaled");
  assert.match(iv.description, /left arm vein/i);
  assert.match(intrathecal.description, /lumbar CSF/i);
  assert.match(intramuscular.description, /deltoid depot/i);
  assert.match(intracisternal.description, /cisterna magna/i);
  assert.match(intracerebroventricular.description, /ventricular CSF/i);
  assert.match(inhaled.description, /airway depot/i);
  const aav9 = iv.capsids.find((capsid) => capsid.capsid_id === "aav9");
  assert.ok(aav9);
  for (const route of payload.administration_routes) {
    assert.ok(route.capsids.every((capsid) => capsid.max_mass_balance_error < 1e-6));
  }
  assert.equal(aav9.regions.brain_frontal.protein_au.length, payload.time_h.length);
  assert.ok(Math.min(...aav9.regions.liver.isf_concentration_vg_ml) >= 0);
  assert.ok(Math.min(...aav9.regions.brain_frontal.protein_au) >= 0);
  const sixHourIndex = payload.time_h.reduce((best, value, index) =>
    Math.abs(value - 6) < Math.abs(payload.time_h[best] - 6) ? index : best, 0);
  const sixHourIsf = payload.region_ids.map((regionId) =>
    aav9.regions[regionId].isf_concentration_vg_ml[sixHourIndex]);
  assert.ok(Math.max(...sixHourIsf) / Math.min(...sixHourIsf.filter((value) => value > 0)) > 50);
  const cnsSixHour = ["brain_frontal", "brain_parietal", "brain_temporal", "brain_occipital", "brain_deep_gray", "brain_cerebellum", "brainstem_spinal"]
    .map((regionId) => aav9.regions[regionId].isf_concentration_vg_ml[sixHourIndex]);
  assert.ok(Math.max(...cnsSixHour) / Math.min(...cnsSixHour) > 3);
  const ivTmax = payload.region_ids.map((regionId) => aav9.regions[regionId].tmax_isf_h);
  assert.ok(new Set(ivTmax.map((value) => value.toFixed(2))).size >= 15);
  const itAav9 = intrathecal.capsids.find((capsid) => capsid.capsid_id === "aav9");
  const imAav9 = intramuscular.capsids.find((capsid) => capsid.capsid_id === "aav9");
  const icmAav9 = intracisternal.capsids.find((capsid) => capsid.capsid_id === "aav9");
  const icvAav9 = intracerebroventricular.capsids.find((capsid) => capsid.capsid_id === "aav9");
  const inhaledAav9 = inhaled.capsids.find((capsid) => capsid.capsid_id === "aav9");
  assert.ok(itAav9.regions.brainstem_spinal.peak_isf_concentration_vg_ml > aav9.regions.brainstem_spinal.peak_isf_concentration_vg_ml * 100);
  assert.ok(imAav9.regions.muscle_injected_arm.peak_isf_concentration_vg_ml > imAav9.regions.muscle_contralateral_arm.peak_isf_concentration_vg_ml * 100);
  assert.ok(icmAav9.regions.brain_cerebellum.peak_isf_concentration_vg_ml > icmAav9.regions.brain_deep_gray.peak_isf_concentration_vg_ml);
  assert.ok(icvAav9.regions.brain_deep_gray.peak_isf_concentration_vg_ml > icmAav9.regions.brain_deep_gray.peak_isf_concentration_vg_ml);
  assert.ok(inhaledAav9.regions.lung_right.peak_isf_concentration_vg_ml > aav9.regions.lung_right.peak_isf_concentration_vg_ml);
  const php = iv.capsids.find((capsid) => capsid.capsid_id === "php-eb");
  assert.match(php.human_translation_note, /LY6A-dependent CNS gain removed/);
});

test("caches GTEx and HPA expression evidence for the disease library", async () => {
  const file = new URL("../public/data/gene-expression.json", import.meta.url);
  const payload = JSON.parse(await readFile(file, "utf8"));
  assert.ok(payload.gene_count >= 40);
  assert.equal(payload.genes.SCN1A.top_modeled_organ, "CNS");
  assert.ok(payload.genes.SCN1A.organ_median_tpm.CNS > payload.genes.SCN1A.organ_median_tpm.Liver);
  assert.match(payload.genes.SCN1A.hpa.entry_url, /proteinatlas\.org/);
});
