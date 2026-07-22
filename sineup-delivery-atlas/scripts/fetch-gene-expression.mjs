import { readFile, writeFile } from "node:fs/promises";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const root = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const diseaseSource = await readFile(resolve(root, "app/disease-data.ts"), "utf8");
const genes = [...new Set([...diseaseSource.matchAll(/gene: "([A-Z0-9-]+)"/g)].map((match) => match[1]))].sort();

const referenceUrl = new URL("https://gtexportal.org/api/v2/reference/gene");
for (const gene of genes) referenceUrl.searchParams.append("geneId", gene);
referenceUrl.searchParams.set("gencodeVersion", "v39");
referenceUrl.searchParams.set("itemsPerPage", "1000");
const referenceResponse = await fetch(referenceUrl);
if (!referenceResponse.ok) throw new Error(`GTEx reference request failed: ${referenceResponse.status}`);
const referencePayload = await referenceResponse.json();
const referenceByGene = new Map(referencePayload.data.map((entry) => [entry.geneSymbol, entry]));

const expressionUrl = new URL("https://gtexportal.org/api/v2/expression/medianGeneExpression");
for (const gene of genes) {
  const reference = referenceByGene.get(gene);
  if (reference) expressionUrl.searchParams.append("gencodeId", reference.gencodeId);
}
expressionUrl.searchParams.set("datasetId", "gtex_v10");
expressionUrl.searchParams.set("itemsPerPage", "100000");
const expressionResponse = await fetch(expressionUrl);
if (!expressionResponse.ok) throw new Error(`GTEx expression request failed: ${expressionResponse.status}`);
const expressionPayload = await expressionResponse.json();

const organMatchers = {
  CNS: (id) => id.startsWith("Brain_"),
  Heart: (id) => id.startsWith("Heart_"),
  Liver: (id) => id === "Liver",
  Kidney: (id) => id.startsWith("Kidney_"),
  Muscle: (id) => id === "Muscle_Skeletal",
  Lung: (id) => id === "Lung",
};

function median(values) {
  if (!values.length) return null;
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

async function fetchHpa(reference) {
  const ensembl = reference.gencodeId.split(".")[0];
  const response = await fetch(`https://www.proteinatlas.org/${ensembl}.json`);
  if (!response.ok) return null;
  const payload = await response.json();
  return {
    tissue_specificity: payload["RNA tissue specificity"] ?? null,
    tissue_distribution: payload["RNA tissue distribution"] ?? null,
    tissue_specificity_score: payload["RNA tissue specificity score"] ?? null,
    protein_tissue_specificity: payload["Protein tissue specificity"] ?? null,
    entry_url: `https://www.proteinatlas.org/${ensembl}`,
  };
}

const rowsByGene = new Map(genes.map((gene) => [gene, []]));
for (const row of expressionPayload.data) rowsByGene.get(row.geneSymbol)?.push(row);

const outputGenes = {};
for (let offset = 0; offset < genes.length; offset += 6) {
  const batch = genes.slice(offset, offset + 6);
  const hpaRows = await Promise.all(batch.map((gene) => {
    const reference = referenceByGene.get(gene);
    return reference ? fetchHpa(reference) : null;
  }));
  batch.forEach((gene, index) => {
    const reference = referenceByGene.get(gene);
    const rows = rowsByGene.get(gene) ?? [];
    const values = rows.map((row) => row.median).filter((value) => Number.isFinite(value));
    const max = Math.max(...values, 0);
    const tau = values.length > 1 && max > 0
      ? values.reduce((sum, value) => sum + (1 - value / max), 0) / (values.length - 1)
      : null;
    const top = [...rows].sort((a, b) => b.median - a.median)[0] ?? null;
    const organMedianTpm = Object.fromEntries(Object.entries(organMatchers).map(([organ, matches]) => [
      organ,
      median(rows.filter((row) => matches(row.tissueSiteDetailId)).map((row) => row.median)),
    ]));
    const rankedOrgans = Object.entries(organMedianTpm)
      .filter(([, value]) => value !== null)
      .sort(([, a], [, b]) => b - a);
    outputGenes[gene] = {
      gencode_id: reference?.gencodeId ?? null,
      organ_median_tpm: organMedianTpm,
      top_modeled_organ: rankedOrgans[0]?.[0] ?? null,
      top_modeled_organ_tpm: rankedOrgans[0]?.[1] ?? null,
      top_gtex_tissue: top?.tissueSiteDetailId ?? null,
      top_gtex_tissue_tpm: top?.median ?? null,
      tissue_tau: tau,
      hpa: hpaRows[index],
    };
  });
}

const payload = {
  schema_version: "1.0",
  generated_at: new Date().toISOString(),
  gene_count: Object.keys(outputGenes).length,
  method: {
    gtex: "Median TPM from GTEx v10; modeled-organ values are medians across matching tissue subregions.",
    tau: "Yanai tissue-specificity index across all GTEx tissue sites; 0 is ubiquitous and 1 is highly tissue-restricted.",
    limitation: "Adult bulk-tissue RNA is a target-presence prior, not cell-type expression, developmental expression, disease causality, or therapeutic efficacy.",
  },
  sources: {
    gtex_api: "https://gtexportal.org/api/v2/docs",
    hpa_downloads: "https://www.proteinatlas.org/about/download",
  },
  genes: outputGenes,
};

await writeFile(resolve(root, "public/data/gene-expression.json"), `${JSON.stringify(payload, null, 2)}\n`);
console.log(`Wrote GTEx/HPA evidence for ${payload.gene_count} genes.`);
