import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const projectRoot = "/Users/zhuzixuan/Documents/Spatial PK";
const outputDir = path.join(projectRoot, "outputs/aav_safety_20260813");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    if (quoted) {
      if (ch === '"' && text[i + 1] === '"') {
        field += '"';
        i += 1;
      } else if (ch === '"') {
        quoted = false;
      } else {
        field += ch;
      }
    } else if (ch === '"') {
      quoted = true;
    } else if (ch === ',') {
      row.push(field);
      field = "";
    } else if (ch === '\n') {
      row.push(field.replace(/\r$/, ""));
      rows.push(row);
      row = [];
      field = "";
    } else {
      field += ch;
    }
  }
  if (field.length || row.length) {
    row.push(field.replace(/\r$/, ""));
    rows.push(row);
  }
  if (rows[0]?.[0]) rows[0][0] = rows[0][0].replace(/^\uFEFF/, "");
  return rows.filter((r) => r.some((v) => v !== ""));
}

function records(rows) {
  const [headers, ...data] = rows;
  return data.map((row) => Object.fromEntries(headers.map((h, i) => [h, row[i] ?? ""])));
}

function colName(index) {
  let n = index + 1;
  let out = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    out = String.fromCharCode(65 + rem) + out;
    n = Math.floor((n - 1) / 26);
  }
  return out;
}

const marginRows = records(parseCsv(await fs.readFile(path.join(projectRoot, "model/data/aav_safety_organ_margins.csv"), "utf8")));
const evidenceRows = records(parseCsv(await fs.readFile(path.join(projectRoot, "model/data/aav_safety_evidence.csv"), "utf8")));

const wb = Workbook.create();
wb.comments.setSelf({ displayName: "Codex" });
const summary = wb.worksheets.add("Summary");
const margins = wb.worksheets.add("Organ_Margins");
const evidence = wb.worksheets.add("Evidence");
const codebook = wb.worksheets.add("Codebook");

for (const sheet of [summary, margins, evidence, codebook]) sheet.showGridLines = false;

const navy = "#17324D";
const teal = "#0F766E";
const lightBlue = "#E8F1F8";
const lightGreen = "#E8F5E9";
const lightAmber = "#FFF4D6";
const lightRed = "#FDECEC";
const gray = "#5E6B75";
const border = "#C9D4DC";

summary.getRange("A1:H1").merge();
summary.getRange("A1").values = [["AAV–SINEUP 安全性证据与模型内暴露裕度（2026-08-13）"]];
summary.getRange("A1:H1").format = { fill: navy, font: { bold: true, color: "#FFFFFF", size: 16 }, rowHeight: 30 };
summary.getRange("A3:B8").values = [
  ["当前人体模型剂量", 4.0e13],
  ["70 kg 总剂量", 2.8e15],
  ["ZOLGENSMA IV 剂量", 1.1e14],
  ["小鼠心脏组织学信号剂量", 7.9e13],
  ["ITVISMA 鞘内总剂量", 1.2e14],
  ["KEBILIDI 脑内总剂量", 1.8e11],
];
summary.getRange("A3:A8").format = { fill: lightBlue, font: { bold: true, color: navy } };
summary.getRange("B3:B8").format.numberFormat = "0.00E+00";
summary.getRange("D3:E8").values = [
  ["IV / ZOLGENSMA 剂量比", null],
  ["IV / 小鼠信号剂量比", null],
  ["IV 器官 AUC 裕度", null],
  ["IT / ITVISMA 剂量比", null],
  ["IT 器官 AUC 裕度", null],
  ["IV AAV9 可提供背景支持器官", null],
];
summary.getRange("E3:E8").formulas = [
  ["='Summary'!$B$3/'Summary'!$B$5"],
  ["='Summary'!$B$3/'Summary'!$B$6"],
  ["='Summary'!$B$6/'Summary'!$B$3"],
  ["='Summary'!$B$4/'Summary'!$B$7"],
  ["='Summary'!$B$7/'Summary'!$B$4"],
  ["=5/8"],
];
summary.getRange("D3:D8").format = { fill: lightBlue, font: { bold: true, color: navy } };
summary.getRange("E3:E8").format.numberFormat = "0.0%";
summary.getRange("E5:E7").format.numberFormat = "0.00x";
summary.getRange("A10:H10").merge();
summary.getRange("A10").values = [["结论边界"]];
summary.getRange("A10:H10").format = { fill: teal, font: { bold: true, color: "#FFFFFF" } };
summary.getRange("A11:H16").merge(true);
summary.getRange("A11:H16").values = [
  ["IV AAV9：当前 4.0×10^13 vg/kg 在同一 PBPK 方程中，相对 ZOLGENSMA 上市剂量的器官 AUC 裕度约 2.75；相对标签中 7.9×10^13 vg/kg 小鼠心脏组织学信号剂量约 1.98。该结果仅支持‘低于选定参照’，不是人体 NOAEL。"],
  ["器官证据：8 个器官组中，心、肝、脑/DRG 背景、肾和脾可获得有限的同途径/同衣壳背景支持；肺、肌肉和 rest-of-body 无器官特异毒性阈值，应标为证据不足。"],
  ["鞘内 AAV9：当前 2.8×10^15 vg 是 ITVISMA 上市剂量的 23.3 倍；按当前模型计算，各器官 AUC 约为 ITVISMA 剂量重跑结果的 23.3–23.7 倍，因此当前 CSF 剂量不受支持。"],
  ["质粒/SINEUP：现有直接证据是小鼠纹状体 7.0×10^9 vg 的 AAV9-miniSINEUP-GDNF，使内源 GDNF 约升高 2 倍、持续至少 6 个月且未报告体重/摄食信号。不能外推为当前未定序列质粒的人体安全性。"],
  ["判定规则：margin ≥ 1 只表示模型暴露低于参照；任何绿色或黄色状态均不得解释为临床安全、无不良反应或可直接选定人体剂量。"],
  ["建议：IV 场景可作为研究假设继续；鞘内先降到 1.2×10^14 vg total 或更低重求解；DRG、免疫/补体和目标蛋白过表达窗口需要新增模型与实验。"],
];
summary.getRange("A11:H16").format = { wrapText: true, verticalAlignment: "top" };
summary.getRange("A11:H12").format.fill = lightGreen;
summary.getRange("A13:H13").format.fill = lightRed;
summary.getRange("A14:H16").format.fill = lightAmber;
summary.getRange("A3:E8").format.borders = { preset: "outside", style: "thin", color: border };
summary.getRange("A10:H16").format.borders = { preset: "outside", style: "thin", color: border };
summary.getRange("A1:H16").format.font.name = "Arial";
summary.getRange("A:A").format.columnWidth = 28;
summary.getRange("B:B").format.columnWidth = 18;
summary.getRange("C:C").format.columnWidth = 3;
summary.getRange("D:D").format.columnWidth = 28;
summary.getRange("E:E").format.columnWidth = 16;
summary.getRange("F:H").format.columnWidth = 12;
summary.getRange("11:16").format.rowHeight = 43;
summary.freezePanes.freezeRows(1);

const marginHeaders = [
  "route_id", "route", "capsid_id", "capsid", "organ", "current_dose_vg",
  "current_dose_vg_per_kg", "current_auc_isf_amount_vg_h", "reference_id",
  "reference_product", "reference_auc_isf_amount_vg_h",
  "exposure_margin_reference_over_current", "assessment", "evidence_grade",
  "comparator_scope", "evidence_note", "source_url",
];
const marginMatrix = [marginHeaders, ...marginRows.map((r) => marginHeaders.map((h) => {
  if (["current_dose_vg", "current_dose_vg_per_kg", "current_auc_isf_amount_vg_h", "reference_auc_isf_amount_vg_h", "exposure_margin_reference_over_current"].includes(h)) {
    return r[h] === "" ? null : Number(r[h]);
  }
  return r[h];
}))];
const marginLastRow = marginMatrix.length;
const marginLastCol = colName(marginHeaders.length - 1);
margins.getRange(`A1:${marginLastCol}${marginLastRow}`).values = marginMatrix;
margins.getRange(`A1:${marginLastCol}1`).format = { fill: navy, font: { bold: true, color: "#FFFFFF" }, wrapText: true, rowHeight: 36 };
margins.getRange(`F2:L${marginLastRow}`).format.numberFormat = "0.00E+00";
margins.getRange(`L2:L${marginLastRow}`).format.numberFormat = "0.00x";
margins.getRange(`M2:M${marginLastRow}`).conditionalFormats.add("containsText", { text: "exceeds", format: { fill: lightRed, font: { color: "#9A1B1B", bold: true } } });
margins.getRange(`M2:M${marginLastRow}`).conditionalFormats.add("containsText", { text: "below-signal", format: { fill: lightGreen, font: { color: "#166534" } } });
margins.getRange(`M2:M${marginLastRow}`).conditionalFormats.add("containsText", { text: "context-only", format: { fill: lightAmber, font: { color: "#854D0E" } } });
margins.getRange(`N2:N${marginLastRow}`).conditionalFormats.add("containsText", { text: "insufficient", format: { fill: lightAmber, font: { color: "#854D0E" } } });
margins.tables.add(`A1:${marginLastCol}${marginLastRow}`, true, "OrganMarginsTable").style = "TableStyleMedium2";
margins.freezePanes.freezeRows(1);
margins.freezePanes.freezeColumns(5);
margins.getRange("A:Q").format.font.name = "Arial";
margins.getRange("A:A").format.columnWidth = 14;
margins.getRange("B:B").format.columnWidth = 24;
margins.getRange("C:D").format.columnWidth = 13;
margins.getRange("E:E").format.columnWidth = 13;
margins.getRange("F:L").format.columnWidth = 18;
margins.getRange("M:N").format.columnWidth = 30;
margins.getRange("O:P").format.columnWidth = 40;
margins.getRange("Q:Q").format.columnWidth = 42;
margins.getRange(`M2:Q${marginLastRow}`).format.wrapText = true;

const evidenceHeaders = ["reference_id", "product_or_study", "capsid_or_modality", "route", "dose_value", "dose_unit", "evidence_kind", "model_use", "interpretation_limit", "source_url"];
const evidenceMatrix = [evidenceHeaders, ...evidenceRows.map((r) => evidenceHeaders.map((h) => h === "dose_value" ? Number(r[h]) : r[h]))];
const evidenceLastRow = evidenceMatrix.length;
evidence.getRange(`A1:J${evidenceLastRow}`).values = evidenceMatrix;
evidence.getRange("A1:J1").format = { fill: navy, font: { bold: true, color: "#FFFFFF" }, wrapText: true, rowHeight: 32 };
evidence.getRange(`E2:E${evidenceLastRow}`).format.numberFormat = "0.00E+00";
evidence.tables.add(`A1:J${evidenceLastRow}`, true, "SafetyEvidenceTable").style = "TableStyleMedium2";
evidence.freezePanes.freezeRows(1);
evidence.getRange("A:J").format.font.name = "Arial";
evidence.getRange("A:A").format.columnWidth = 32;
evidence.getRange("B:B").format.columnWidth = 40;
evidence.getRange("C:D").format.columnWidth = 24;
evidence.getRange("E:F").format.columnWidth = 20;
evidence.getRange("G:I").format.columnWidth = 44;
evidence.getRange("J:J").format.columnWidth = 48;
evidence.getRange(`A2:J${evidenceLastRow}`).format.wrapText = true;
evidence.getRange(`A2:J${evidenceLastRow}`).format.rowHeight = 54;

codebook.getRange("A1:D1").values = [["字段/状态", "中文解释", "允许的结论", "禁止的结论"]];
codebook.getRange("A2:D9").values = [
  ["exposure_margin_reference_over_current", "参照剂量重跑所得器官 ISF amount AUC / 当前器官 AUC", "比较模型内暴露", "不良反应概率或人体安全窗"],
  ["below-signal-or-clinical-context-not-proven-safe", "低于选定同途径/同衣壳参照", "可作为起始剂量背景", "已证明安全"],
  ["exceeds-marketed-product-exposure-context", "高于上市产品同途径剂量下的模型暴露", "应降剂量并重求解", "只是保守、可以忽略"],
  ["cross-capsid-context-only", "用 AAV9 参照比较非 AAV9 衣壳", "风险排序假设", "等效安全阈值"],
  ["closest-route-context-only", "ICM/ICV 暂用 IT 作为最近途径背景", "提示局部 CNS 风险", "途径等价"],
  ["insufficient-evidence", "无器官或产品特异阈值", "标记缺口并设计实验", "默认安全"],
  ["AUC", "组织间液中载体量对时间积分，单位 vg·h", "相对暴露指标", "组织病理、ALT 或事件概率"],
  ["DRG", "背根神经节；当前 PBPK 未单独建模", "独立神经毒理终点", "用 whole-brain AUC 代替"],
];
codebook.getRange("A1:D1").format = { fill: navy, font: { bold: true, color: "#FFFFFF" } };
codebook.getRange("A2:D9").format = { wrapText: true, rowHeight: 48 };
codebook.getRange("A:D").format.font.name = "Arial";
codebook.getRange("A:A").format.columnWidth = 44;
codebook.getRange("B:D").format.columnWidth = 40;
codebook.getRange("A1:D9").format.borders = { preset: "outside", style: "thin", color: border };
codebook.freezePanes.freezeRows(1);

summary.getRange("B3").conditionalFormats.add("cellIs", { operator: "lessThan", formula: "='Summary'!$B$6", format: { fill: lightGreen } });
summary.getRange("E6").conditionalFormats.add("cellIs", { operator: "greaterThan", formula: 1, format: { fill: lightRed, font: { bold: true, color: "#9A1B1B" } } });

await fs.mkdir(outputDir, { recursive: true });
for (const sheetName of ["Summary", "Organ_Margins", "Evidence", "Codebook"]) {
  const preview = await wb.render({ sheetName, autoCrop: "all", scale: 1, format: "png" });
  await fs.writeFile(path.join(outputDir, `${sheetName}.png`), new Uint8Array(await preview.arrayBuffer()));
}
const inspection = await wb.inspect({ kind: "table", range: "Summary!A1:H16", include: "values,formulas", tableMaxRows: 20, tableMaxCols: 10 });
await fs.writeFile(path.join(outputDir, "inspection.ndjson"), inspection.ndjson, "utf8");
const errors = await wb.inspect({ kind: "match", searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A", options: { useRegex: true, maxResults: 300 }, summary: "final formula error scan" });
await fs.writeFile(path.join(outputDir, "formula_errors.ndjson"), errors.ndjson, "utf8");
const xlsx = await SpreadsheetFile.exportXlsx(wb);
await xlsx.save(path.join(outputDir, "aav_safety_margin_review_2026.xlsx"));
console.log(path.join(outputDir, "aav_safety_margin_review_2026.xlsx"));
